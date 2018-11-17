function [f, G] = my_obj(X,A)
  G = X-A;
  G(A==0)=0;
  f = 0.5*norm(G,'fro')^2;
end