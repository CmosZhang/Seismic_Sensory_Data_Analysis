function [Y_f] = alter_min_LS_one_step(T_omega_f, omega_f, X_f)
% the target dimension: r * s * k

[m,n,k] = size(T_omega_f);
[m,r,k] = size(X_f);

Y_f = zeros(r, n, k);
%% prepare variables
%X_f_new = zeros(m*k, r*k);
%for i=1:k
%    X_f_new((i-1)*m + 1 : i*m, (i-1)* r + 1 : i*r ) = X_f(:,:,i);
%end

X_f_new = [];
for i=1:k
    X_f_new = blkdiag(X_f_new, X_f(:,:,i));  % make it block diagonal
end


%% we recover the lateral slice one-by-one.
residual = 0;
tensor_V = zeros(k * m, 1);
temp_Y_f = zeros(r * k, 1);
for i = 1:n
    for j = 1:k
        tensor_V((j-1)*m + 1 : j*m) = squeeze(T_omega_f(:,i,j)); 
    end
    
    omega_f_3D = zeros(k,k,m);
    omega_f_new = zeros(k*m, k*m);
    for j = 1:m
        % omega_f_3D(:,:,j) = circ(squeeze(omega_f(j,i,:)));
        temp = gallery('circul',reshape(omega_f(j,i,:),1,k));
        omega_f_3D(:,:,j) = temp.';
    end
    
    for a=1:k
        for b=1:k
            for c=1:m
                row = (a-1)*m + c;
                col = (b-1)*m + c;
                omega_f_new(row,col) = omega_f_3D(a,b,c);
            end
        end
    end
    
    %[temp_Y_f, restnorm] = lsqlin(sparse(omega_f_new) * sparse(X_f_new), tensor_V);
    temp =  (omega_f_new) *  (X_f_new);   %%%%%%%%%%%%%% 
    temp_Y_f = pinv(temp'*temp)*temp'* tensor_V;
    %residual = residual + sqrt(restnorm)
    for j = 1:k
        Y_f(:,i,j) = temp_Y_f((j-1)*r + 1 : j*r);
    end 
end
