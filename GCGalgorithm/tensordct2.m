function [dT]=tensordct2(T)
[a,b,c]=size(T);
dT = zeros(a,b,c);
for i=1:a
    for j = 1:b
        dT(i,j,:)=dct(T(i,j,:));
    end
end
end
