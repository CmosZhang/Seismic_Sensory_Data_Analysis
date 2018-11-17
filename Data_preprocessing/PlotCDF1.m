function PlotCDF1(CDF, dimensional)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

str = {};
sz = length(CDF(1, :));
figure; 
markers = ['-o'; '-+'; '-*'; '-s'; '-d'; '->'; '-p'; '-h'];
for i = 1 : size(CDF, 1)
    plot((1:sz)/sz * 100, CDF(i, :), markers(i, :) ,'MarkerSize', 4);hold on; 
end
grid on; title('Empirical CDF'); xlabel('Top portion of singular values *100(%)'); ylabel('CDF');
legend('original data','processed data');
end
