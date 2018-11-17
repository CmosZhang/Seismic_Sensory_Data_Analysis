clc; 
clear all; 
close all;



%% data loading
% load('T_synthetic_tuabl_rank_2.mat');% T = T(:, :, 1:256);   %加载我们自己合成的人工合成数据;
load('volume.mat');  %数据大小：326*431*531
volume = volume(1:300,1:120,1:80);


%% 数据预处理
original_tubalRank = LowTubalCDF(volume,1);
[U,S,V] = tSVDs(volume,15);
temp = tprod(U,S);
T_test = tprod(temp,V);
T_test_RSE = norm(T_test(:) - volume(:))/ norm(volume(:));

original_tubalRank1 = LowTubalCDF(T_test,1);
T1 = T_test;

srs = [0.01 : 0.01 : 0.1];
sampling_rse = zeros(1, length(srs));


%% 调节tubal-rank
r = 15;


szT1 = size(T1); 

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
temp = 0;
X_est = ifft(X_f, [], 3); 
Y_est = ifft(Y_f, [], 3);
T_est = tprod(X_est, Y_est);
RSE =  norm(T_est(:) - T(:)) / norm(T(:));
     
sampling_rse(:, loop) = sampling_rse(:, loop) + RSE(:);
end
end
sampling_rse = sampling_rse ./ 10;

% 重构误差
figure;semilogy([0.01 : 0.01 : 0.1]*100, sampling_rse(1,:), '+-'); title(['Reconstruction Error']);
legend( 'Tubal-Alt-Min'); 
xlabel('Slice Missing Rate %');ylabel('RSE');