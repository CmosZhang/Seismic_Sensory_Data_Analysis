function [solution, obj, iter, msg] = ...
          solve_trace_reg(objf, lambda, X0, evalf, opts)

        
% Partially corrective boosting
% Solve  Min_X  func(X) + lambda * ||X||_*
% Inputs:
%     objf: a function handle: [f, G] = objf(X), 
%     where f and G are objective value and gradient of loss L at X, resp.
%     Both X and G are matrices
%   evalf: a function handle that evaluates the performance of the solution
%           at each iteration (eg computes NMAE on the test set)
%   opts: parameters for the solver
% Outputs:
%   solution: the solution
%   obj: the objective value at the solution
%   iter: number of boosting iterations
%   msg: termination message

[m, n] = size(X0);
k = opts.init_rank;  % Initialize U and V with rank k

U = randn(m, k);
V = randn(k, n);

total_time = 0;     % Total time consumed by the algorithm
local_search_time = 0;  % Time spend on local search

global feval_logistic time_logistic

feval_logistic = 0;     
time_logistic = 0;

perf = 0;
for i = 1 : opts.max_iter
      
    t1 = cputime;
    k = size(U, 2);

  if opts.use_local     
    t2 = cputime;
    nel = k*(m+n);
    ub = Inf(nel, 1);
    
    % local search    
    [UV, obj, ~, ~, msg] = lbfgsb([U(:); V(:)], -ub, ub, @obj_UV, ...
                                              [], [], opts.lbfgsb_in);
    U = reshape(UV(1:k*m), [m, k]);
    V = reshape(UV(1+k*m:end), [k, n]);
    t2 = cputime - t2;    
    local_search_time = local_search_time + t2;    
  end
    
  norm_U = sqrt(sum(U.*U))';
  norm_V = sqrt(sum(V.*V, 2));
  norm_UV = sqrt(norm_U .* norm_V);
  idx = (norm_U > 1e-5) & (norm_V > 1e-5);
  nidx = sum(idx);
  if nidx < length(idx)
    U = U(:, idx);
    V = V(idx, :);
    norm_U = norm_U(idx);
    norm_V = norm_V(idx);
    norm_UV = norm_UV(idx);
  end
  
  if opts.use_local
     loss = obj - 0.5*lambda*(norm_U'*norm_U + norm_V'*norm_V);
     if isempty(loss)
     loss = obj;
     end
  end
  
  if opts.use_local || i > 1    
     fprintf('\niter=%d, loss = %g, obj=%g, perf=%g, k=%d, time=%g, ls_time=%g, msg=%s\n', ...
          i, loss, obj, perf, k, total_time, local_search_time, msg2str(msg, 'lbfgsb'));
  end
  
  if nidx > 0
     sk = 0.5*(norm_U'*norm_U + norm_V'*norm_V); % Xk = U*V;    
  else
     sk = 0; % Xk = zeros(m, n);
  end

  if i > 1 && abs(pre_obj-obj) / min(abs([pre_obj, obj])) < opts.rtol
     msg = 'Stop with small relative change';
     break;
  elseif total_time > opts.max_time
     msg = 'Stop by max_time';
    break;
  elseif i == opts.max_iter
     msg = 'Stop with max iteration';  
    break;
  else
     pre_obj = obj;
  end

  [~, G] = objf(U*V);    % G is a sparse matrix
  [u, ~, v] = svds(G, 1);
  v = -v;

  % line search
  [weights, obj, ~, ~, msg] = lbfgsb([1; 0.5], [0; 0], [inf; inf], ...
                                 @obj_ls, [], [], opts.lbfgsb_in);
  
  loss = obj - lambda*(sk*weights(1) + weights(2));
  weights = sqrt(weights);
  if nidx > 0
      U = [bsxfun(@times, U, weights(1)*(norm_UV./norm_U)'), weights(2) * u];
      V = [bsxfun(@times, V, weights(1)*norm_UV./norm_V); weights(2) * v']; 
  else
      U = weights(2) * u;
      V = weights(2) * v'; 
  end
  
  t1 = cputime - t1;  
  total_time = total_time + t1;
  
  if ~isempty(evalf)
     perf = evalf(U*V');
  end
  
end

  solution = U*V;
  iter = i;

if opts.verbose
   fprintf(1, '%s\n', msg);
end
   
  % local search objective
  function [f, g] = obj_UV(X)
    r = length(X) / (m+n);
    UU = reshape(X(1:r*m), [m, r]);
    VV = reshape(X(1+r*m:end), [r, n]);
        
    [f, G] = objf(UU*VV);
    f = f + 0.5*lambda*norm(X)^2;
    g = [vec(G*VV' + lambda*UU)
        vec(UU'*G +lambda*VV)];
  end
  
  % Line search objective
  function [f, g] = obj_ls(x)
    [f, G] = objf([sqrt(x(1))*U, sqrt(x(2))*u]*[sqrt(x(1))*V; sqrt(x(2))*v']);
    f = f + lambda*(sk*x(1) + x(2));
    g = [sum(sum(U .* (G*V')))+lambda*sk
         u'*(G*v)+lambda];
  end

  function out_v = vec(in_v)
           out_v = in_v(:);
  end
end
