function opts = init_opts(original)

  if nargin > 0
    opts = original;
  else
    opts = [];
  end
  
%   opts.init_rank = 30;
  opts.max_iter = 300;    % maximum number of iterations for boosting
  opts.rtol = 1e-5;       % terminate if relative difference falls below it
  opts.use_local = true;  % whether to use local optimization at each iteration
  
  % parameters related to testing 
  lbfgsb_in = [];
  lbfgsb_in.maxIter = 300;		 % max number of iterations
  lbfgsb_in.maxFnCall = 1000;	 % max number of calling the function
  lbfgsb_in.relCha = 1e-5;	 % tolerance of constraint satisfaction
  lbfgsb_in.tolPG = 1e-5;		 % final objective function accuracy parameter
  lbfgsb_in.m = 20;

  opts.lbfgsb_in = lbfgsb_in;
  
  opts.max_time = 3600*3;
  opts.verbose = 1;
  
end
