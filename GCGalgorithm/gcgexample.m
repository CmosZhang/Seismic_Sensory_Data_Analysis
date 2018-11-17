clc; 
clear all; 
close all;

%% data loading
load('T_synthetic_tuabl_rank_2.mat');
% load('T_synthetic_tuabl_rank_3.mat');
% load('T_synthetic_tuabl_rank_3_30_30_600.mat');
% load('T_synthetic_tuabl_rank_3_32_32_300.mat');
% load('traces_100_100_1000.mat'); 


szT = size(T);                %原始数据大小 
tubalRank = LowTubalCDF(T, 1);

%% add noise
% T1 = generatedata (T,10);  %加10dB的噪声
T1 = T;  %不加噪声


%% data sampling (tubal-sampling)
samplingRate = 0.5;
szT1 = size(T1);
Omega = repmat((rand(szT1(1:2)) > samplingRate), [1, 1, szT1(3)]);
T2 = T1;
T2(Omega) = 0; 
Omega = abs(1 - Omega);

%% DCT 变换
T2_d = tensordct2(T2);

%% GCG algorithm
lambda = 0.1; %目标函数中的参数
[~,~,c]=size(T2_d);
for i=1:c
    A = T2_d(:,:,i);
    [a,b]=size(A);
    
    % use boosting solver
    opts = init_opts();
    r = rank(T2_d(:,:,i));
    opts.init_rank = r;
    X0 = zeros(a, b);
    evalf = [];   % no functor for performance evaluation at each iteration
    [T3(:,:,i), opt_obj_boost_solver,iter, msg] = solve_trace_reg(@(X)my_obj(X,A), lambda, X0, evalf, opts);
end
  

%% 逆DCT变换
r_T = tensoridct2(T3);

%% figure
slice = 5;
SeisPlot(T(:,slice,:),{'title','Original Data'});
SeisPlot(T2(:,slice,:),{'title','Corrupted Data'});
SeisPlot(r_T(:,slice,:),{'title','Recovery Data'});

%% reconstruction error
GCGRSE =  norm(r_T(:) - T(:)) / norm(T(:));
S_N = -10*log(GCGRSE);
fprintf('***********************GCGRSE = %d ***********\n',GCGRSE);  