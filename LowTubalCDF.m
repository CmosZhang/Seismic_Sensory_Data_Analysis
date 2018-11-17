function [r, CDF]=LowTubalCDF(T, sampleRate)
%LOWTUBALCDF 用于检查数据的低秩性;
% Input : 
%     T - 原始数据，3d tensor
%     sampleRate - 采样率向量, 采样率由小到大递增;
% Output :
%     r - T 的秩;
% Author ： Chirl
% Edit Date： 2017/4/25

if ~exist('sampleRate', 'var')
    sampleRate = [0.2 : 0.2 : 1];
else
    sampleRate = sort(sampleRate);
end

szT = size(T);
numSample = length(sampleRate);
sz = min(szT(1:2));
TubalRank = zeros(numSample, sz);
% numRepet = 40;
% for j = 1 : numRepet 
    for i = 1 : numSample
        Omega = rand(size(T)) < sampleRate(i);
        [~, S, ~] = tSVD(tPointProd(Omega, T));
        
        S1 = S(1:sz, 1:sz, :);
        for k = 1 : sz(1)
            TubalRank(i, k) = TubalRank(i, k) + norm(squeeze(S1(k, k, :)));
        end
    end
% end
% TubalRank = TubalRank / numRepet;

sumN = sum(TubalRank, 2);
CDF = zeros(size(TubalRank));
r = 0;
for i = 1 : sz
    CDF(:, i) = 100*sum(TubalRank(:, 1:i), 2) ./ sumN;
    if CDF(end, i) > 90 && r == 0
        r = i;
    end
end

str = {};
figure; 
markers = ['-o'; '-+'; '-*'; '-s'; '-d'; '->'; '-p'];
for i = 1 : numSample
    plot((1:sz)/sz * 100, CDF(i, :), markers(i, :) ,'MarkerSize', 4);hold on; 
    str = [str; [num2str(sampleRate(i)*100), '% sampling rate']];
end
grid on; title('Empirical CDF'); xlabel('Top portion of singular value %'); ylabel('CDF');
legend(str);



end

