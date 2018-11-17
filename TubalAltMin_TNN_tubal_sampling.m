clc; 
clear all; 
close all;


%% data loading
load('T_synthetic_tuabl_rank_2.mat');   %加载我们自己合成的人工合成数据;
% load('T_synthetic_tuabl_rank_3.mat');   %加载我们自己合成的人工合成数据;
% load('traces_100_100_1000.mat');      %加载真实地震数据
% load('05HBC3D_ALL_POST_MIG_200_25_601_T.mat');  %加载真实地震数据


 
tubalRank = LowTubalCDF(T, 1);
r = tubalRank;

%%  Tubal sampling   %%%
samplingRate = 0.5;
T1 = T;
szT = size(T1);
omega = repmat((rand(szT(1:2)) > samplingRate), [1, 1, szT(3)]);
T2 = T1;
T2(omega) = 0; 
omega = abs(1 - omega);
 


%% observations
[m,n,k] = size(T2);
T_f = fft(T2, [], 3);

T_omega = omega .* T2;  
T_omega_f = fft(T_omega,[],3);
omega_f = fft(omega, [], 3);

%% Alternating Minimization
%% X: m * r * k
%% Y: r * n * k
%% Given Y, do LS to get X
    Y = rand(r, n, k);
    Y_f = fft(Y, [], 3);

%% do the transpose for each frontal slice
    Y_f_trans = zeros(n,r,k);
    X_f = zeros(m,r,k);
    T_omega_f_trans = zeros(n,m,k);
    omega_f_trans = zeros(n,m,k);
for i = 1: k
     Y_f_trans(:,:,i) = Y_f(:,:,i)';
     T_omega_f_trans(:,:,i) = T_omega_f(:,:,i)';
     omega_f_trans(:,:,i) = omega_f(:,:,i)';
end

iter=1;
while iter <=15
   [X_f_trans] = alter_min_LS_one_step(T_omega_f_trans, omega_f_trans * 1/k, Y_f_trans);
    for i =1:k
        X_f(:,:,i) = X_f_trans(:,:,i)';
    end
    
    %% Given X, do LS to get Y
   [Y_f] = alter_min_LS_one_step(T_omega_f, omega_f * 1/k, X_f);
    
    for i = 1: k
        Y_f_trans(:,:,i) = Y_f(:,:,i)';
    end
    iter = iter + 1;
end

%% The relative error:
X_est = ifft(X_f, [], 3); 
Y_est = ifft(Y_f, [], 3);
T_est = tprod(X_est, Y_est);
Tubal_Alt_Min_RSE =  norm(T_est(:) - T(:)) / norm(T(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       TNN                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalize              =        max(T(:))                    ;
Xn                     =        T/normalize                  ;
[n1,n2,n3]             =        size(Xn);                      
alpha                  =        1                             ;
rho                    =        0.08                          ;
myNorm                 =        'tSVD_1'                      ; % dont change for now
A                      =        diag(sparse(double(omega(:)))); % sampling operator
b                      =        A * Xn(:)                     ; % available data
bb                     =        reshape(b,[n1,n2,n3]);

%% ================ main process of completion =======================
[tnnT, tnnCurve]   =    tensor_cpl_admm( A , b , Xn, rho , alpha , ...
                                        [n1,n2,n3], 300, myNorm , 0 );
tnnT               =        tnnT * normalize                 ;
tnnT               =        reshape(tnnT,[n1,n2,n3])         ;
TNNRSE =  norm(tnnT(:) - T(:)) / norm(T(:));
fprintf('***********************TNNRSE = %d ***********\n',TNNRSE); 


%% figure
figure;
subplot(1,3,1);
SeisPlot(T(:,5, :),{'figure', 'old'});
xlabel('CMP x number');ylabel('Time(ms)')
subplot(1,3,2);
SeisPlot(T_est(:,5, :),{'figure', 'old'});
xlabel('CMP x number');ylabel('Time(ms)') 
subplot(1,3,3);
SeisPlot(tnnT(:,5,:),{'figure', 'old'});
xlabel('CMP x number');ylabel('Time(ms)')


