clc; 
clear all; 
close all;


%% data loading
% load('T_synthetic_tuabl_rank_2.mat');   %加载我们自己合成的人工合成数据;
load('volume.mat');  %数据大小：326*431*531
T = volume(1:300,1:120,1:80);
% T = permute(T,[3,1,2]);

%% 数据预处理
original_tubalRank = LowTubalCDF(T,1);
[U,S,V] = tSVDs(T,15);
temp = tprod(U,S);
T_test = tprod(temp,V);
T_test_RSE = norm(T_test(:) - T(:))/ norm(T(:));

original_tubalRank1 = LowTubalCDF(T_test,1);
T1 = T_test;



srs =[0.02 : 0.02 : 0.2];
sampling_rse = zeros(2, length(srs));
O_sampling_rse = zeros(2, length(srs));

szT1 = size(T1);   

%% 调节tubal-rank
r = 15;


for i = 1 : 10
for loop = 1 : length(srs)
    samplingrate = srs(loop);
    MatOmega1 = randsample([0 1],szT1(3),true,[samplingrate, 1-samplingrate]);
    omega = ones(szT1(1),szT1(2),szT1(3));
    for k=1:szT1(3)
        if(MatOmega1(k)==0)
           omega(:,:,k)=zeros(szT1(1),szT1(2));
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Tubal-Alt-Min                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RSE = zeros(1,2); 
OriginalRSE = zeros(1,2);
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



   maxiter = 15;
   TubalAMCurve = zeros(maxiter, 1);
   iter=1;
   while iter <=maxiter
        [X_f_trans] = alter_min_LS_one_step(T_omega_f_trans, omega_f_trans * 1/k, Y_f_trans);
    
        for i =1:k
            X_f(:,:,i) = X_f_trans(:,:,i)';
        end
    
    %% Given X, do LS to get Y
        [Y_f] = alter_min_LS_one_step(T_omega_f, omega_f * 1/k, X_f);
    
        for i = 1: k
            Y_f_trans(:,:,i) = Y_f(:,:,i)'; 
        end
    %% The relative error:
       temp = 0;
       X_est = ifft(X_f, [], 3); 
       Y_est = ifft(Y_f, [], 3);
       T_est = tprod(X_est, Y_est);
       RSE(1) =  norm(T_est(:) - T1(:)) / norm(T1(:));
       OriginalRSE(1) =  norm(T_est(:) - T(:)) / norm(T(:));
       TubalAMCurve(iter) = RSE(1);
       iter = iter + 1;
   end
   fprintf('***********************RSE = %d ***********\n',RSE(1));  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       TNN                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normalize              =        max(T1(:))                     ;
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
                                            [n1,n2,n3] , 200, myNorm , 0 );
    tnnT               =        tnnT * normalize                 ;
    tnnT               =        reshape(tnnT,[n1,n2,n3])         ;
   TNNQ1 = -20*log10(norm(tnnT(:) - T1(:)) / norm(T1(:)));
   RSE(2) =  norm(tnnT(:) - T1(:)) / norm(T1(:));
   OriginalRSE(2) =  norm(tnnT(:) - T(:)) / norm(T(:));
   fprintf('***********************RSE = %d ***********\n',RSE(2)); 
   
   sampling_rse(:, loop) = sampling_rse(:, loop) + RSE(:);
   O_sampling_rse(:, loop) = O_sampling_rse(:, loop) + OriginalRSE(:);
end
end
sampling_rse = sampling_rse ./ 10;
O_sampling_rse = O_sampling_rse ./ 10;

% 绘制收敛速度
len = max([length(TubalAMCurve), length(tnnCurve)]);
figure;semilogy(TubalAMCurve(1:length(TubalAMCurve)), 'o-'); title(['Convergence Speed ']);
hold on; semilogy(tnnCurve(1:length(tnnCurve)), '*-');
legend('Tubal-Alt-Min', 'TNN'); 
xlabel('Iteration');ylabel('RSE in log-scale');
grid on;

% 重构误差
figure;plot([0.98: -0.02 : 0.8]*100, sampling_rse(1,:), '+-'); title(['Reconstruction Error']);
hold on; plot([0.98: -0.02 : 0.8]*100, sampling_rse(2, :), '*-');
legend( 'Tubal-Alt-Min','TNN'); 
xlabel('Slice Sampling Rate * 100 (%)');ylabel('RSE');
grid on;