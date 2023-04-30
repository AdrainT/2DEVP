function [ mueigopt,lambdaeigopt,flag,i,info_eigopt,niter ] = minimaxRay_eigopt( A,C,opts )
%   Use eigopt to solve the following minimax problem:
%      min_{x\neq 0} max{x^HAx/x^Hx, x^HBx/x^Hx}
%      transform to max\mu\in[0,1] \lambda_min(A-mu C), C = A-B;

i = 1;
info_eigopt = [];
flag = true;
parsin.A{1} = A;
parsin.A{2} = C;
parsin.tol = opts.tol;
parsin.gamma = 0;
parsin.itertol = 100000;
parsin.fname = 'funMIMOAC';
parsin.minmax = 1;
bounds.lb = 0;
bounds.ub = 1;
[lambdaeigopt,mueigopt,parsout] = eigopt('funMIMO',bounds,parsin);
niter = parsout.nfevals;
