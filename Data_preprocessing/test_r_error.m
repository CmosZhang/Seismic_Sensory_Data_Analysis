clc; 
clear all; 
close all;

%% data loading
load('volume.mat');  %���ݴ�С��326*431*531
T = volume(:,:,:);
% load('T_synthetic_tuabl_rank_2.mat');%T = T(:, :, 1:256);   %���������Լ��ϳɵ��˹��ϳ�����;

[tubalRank,CDF1] = LowTubalCDF(T, 1);

[U,S,V] = tSVDs(T,32); % ����tSVDs(�ض�tensor SVD�ֽ���ȣ��������)

temp = tprod(U,S);
T_test = tprod(temp,V);
RSE =  norm(T_test(:) - T(:)) / norm(T(:));

[tubalRank1,CDF2] = LowTubalCDF(T_test, 1);
CDF=[CDF1;CDF2];
PlotCDF1(CDF,[1,2]);

SeisPlot(squeeze(T(:,:, 5))', {'title', 'Original data'});
SeisPlot(squeeze(T_test(:,:,5))', {'title', 'Recovered data'});