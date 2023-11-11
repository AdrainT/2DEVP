function [ muleigopt,lambdaleigopt,flag,i,info_eigopt,niter, tSub  ] = minimaxRay_leigopt( funAs,funBs,n,opts )
%   Use leigopt to solve the following minimax:
%      min_{x\neq 0} max{x^HAx/x^Hx, x^HBx/x^Hx}
%      transform to max\mu\in[0,1] \lambda_min(A-mu C), C = A-B;

i = 1;
d = 1;
info_eigopt = [];
flag = true;
pars.tol = opts.tol;
pars.gamma = 0;
pars.bounds.lb = 0;
pars.bounds.ub = 1;
if isfield(opts,'useCvx')
    pars.useCvx = opts.useCvx;
end
[lambdaleigopt,muleigopt,tSub,niter] = leigopt_max_MIMO(funAs,funBs,d,n,pars);

