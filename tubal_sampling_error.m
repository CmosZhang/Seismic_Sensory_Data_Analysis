clc; 
clear all; 
close all;

%% data loading
load('T_synthetic_tuabl_rank_2.mat');  %加载我们自己合成的人工合成数据;
% load('T_synthetic_tuabl_rank_3.mat'); %加载我们自己合成的人工合成数据;
% load('traces_100_100_1000.mat');     %加载真实地震数据


tubalRank = LowTubalCDF(T, 1);
r = tubalRank;

srs = [0.05 : 0.05 : 0.95];
sampling_rse = zeros(2, length(srs));

for i = 1 : 10
for loop = 1 : length(srs)
 numIter = 100;
 samplingRate = srs(loop);
  
  
 %% data sampling 
 szT = size(T);
 omega = repmat((rand(szT(1:2)) > samplingRate), [1, 1, szT(3)]);


  T1 = T;
  T1((omega)) = 0; 
  omega = abs(1 - omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Tubal-Alt-Min                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RSE = zeros(1,2); 
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
    X_est = ifft(X_f, [], 3); 
    Y_est = ifft(Y_f, [], 3);
    T_est = tprod(X_est, Y_est);
    RSE(1) =  norm(T_est(:) - T1(:)) / norm(T1(:));
    TubalAMCurve(iter) = RSE(1);
    iter = iter + 1;
    end
   fprintf('***********************Tubal_Alt_Min RSE = %d ***********\n',RSE(1));  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       TNN                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        normalize              =        max(T(:))                     ;
        Xn                     =        T/normalize                   ;
        [n1,n2,n3]             =        size(Xn);                      

        alpha                  =        1                             ;
        rho                    =        0.08                          ;
        myNorm                 =        'tSVD_1'                      ; % dont change for now

        A                      =        diag(sparse(double(omega(:)))); % sampling operator
        b                      =        A * Xn(:)                     ; % available data
        bb                     =        reshape(b,[n1,n2,n3]);

        %% ================ main process of completion =======================
        [tnnT, tnnCurve]   =    tensor_cpl_admm( A , b , Xn, rho , alpha , ...
                                                [n1,n2,n3] , numIter+200, myNorm , 0 );
        tnnT                   =        tnnT * normalize                 ;
        tnnT                   =        reshape(tnnT,[n1,n2,n3])         ;
        
       RSE(2) = norm(tnnT(:) - T(:)) / norm(T(:));
       fprintf('***********************TNN RSE = %d ***********\n', RSE(2));
       
       sampling_rse(:, loop) = sampling_rse(:, loop) + RSE(:);
end
end
sampling_rse = sampling_rse ./ 10;

% % 绘制收敛速度
% len = min([length(TubalAMCurve), length(tnnCurve)]);
% figure;semilogy(TubalAMCurve(1:len), 'o-'); title(['Convergence Speed ']);
% hold on; semilogy(tnnCurve(1:len), '*-');
% legend('Tubal-Alt-Min', 'TNN'); 
% xlabel('Iteration');ylabel('RSE in log-scale');
% grid on;

%% 重构误差
figure;semilogy([0.05:0.05:0.95]*100, sampling_rse(1,:), '+-'); title(['Reconstruction Error']);
hold on; semilogy([0.05:0.05:0.95]*100, sampling_rse(2, :), '*-');
legend( 'Tubal-Alt-Min', 'TNN'); 
xlabel('Sampling Rate %');ylabel('RSE in log-scale');