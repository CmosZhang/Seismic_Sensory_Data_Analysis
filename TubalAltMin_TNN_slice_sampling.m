clc; 
clear all; 
close all;


%% data loading
% load('T_synthetic_tuabl_rank_2.mat'); %加载我们自己合成的人工合成数据;
% T = T(1:60,1:60,1:100);
load('volume.mat');  %数据大小：326*431*531
volume = volume(1:300,1:180,1:80);


%% 数据预处理
original_tubalRank = LowTubalCDF(volume,1);
[U,S,V] = tSVDs(volume,18);
temp = tprod(U,S);
T_test = tprod(temp,V);
T_test_RSE = norm(T_test(:) - volume(:))/ norm(volume(:));

original_tubalRank1 = LowTubalCDF(T_test,1);
T1 = T_test;


%% 调节秩
r = 18;

%% Slice sampling
szT1 = size(T1);
samplingrate = 0.01;
MatOmega = randsample([0 1],szT1(3),true,[samplingrate, 1-samplingrate]);
omega = ones(szT1(1),szT1(2),szT1(3));
for k=1:szT1(3)
    if(MatOmega(k)==0)
       omega(:,:,k)=zeros(szT1(1),szT1(2));
    end
end


%% 5-th frontal slice missing  
%最终画图画第五个slice,第五个面缺失主要是为了画图方便
omega(:,:,5)=zeros(szT1(1),szT1(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Tubal_Alt_Min                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% observations
[m,n,k] = size(T1);
T_f = fft(T1, [], 3);

T_omega = omega .* T1;  
T_omega_f = fft(T_omega,[],3);
omega_f = fft(omega, [], 3);

%% Alternating Minimization
% X: m * r * k
% Y: r * n * k
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

% 考虑增加迭代次数。
maxiter = 25;
iter=1;
while iter <= maxiter
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
temp = 0;
X_est = ifft(X_f, [], 3); 
Y_est = ifft(Y_f, [], 3);
T_est = tprod(X_est, Y_est);
RSE =  norm(T_est(:) - T1(:)) / norm(T1(:));
fprintf('***********************Tubal_Alt_Min_RSE = %d ***********\n',RSE); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       TNN                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalize              =        max(T1(:))                    ;
Xn                     =        T1/normalize                  ;
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
TNNRSE =  norm(tnnT(:) - T1(:)) / norm(T1(:));
fprintf('***********************TNNRSE = %d ***********\n',TNNRSE); 
   

% figure 画第五个slice。
figure;
subplot(1,3,1);
SeisPlot(squeeze(T(:,:, 5))',{'figure', 'old'});
xlabel('CMP x number');ylabel('Time(ms)')
subplot(1,3,2);
SeisPlot(squeeze(tnnT(:,:, 5))',{'figure', 'old'});
xlabel('CMP x number');ylabel('Time(ms)') 
subplot(1,3,3);
SeisPlot(squeeze(T_est(:,:,5))',{'figure', 'old'});
xlabel('CMP x number');ylabel('Time(ms)')