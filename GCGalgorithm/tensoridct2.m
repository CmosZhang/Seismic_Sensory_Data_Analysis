function [T]=tensoridct2(dT)
[a,b,c]=size(dT);
T = zeros(a,b,c);
for i=1:a
    for j=1:b
        T(i,j,:)=idct(dT(i,j,:));
    end
end
end
