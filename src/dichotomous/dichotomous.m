function [mu,lambda,X0,info] = dichotomous(A,C,opt)
% Find the maximum minium eigenvalue of A - mu C
a = opt.bd(1);
b = opt.bd(2);

epsilon = opt.tol/3;
func = @(mu) min(eig(A-mu*C));

info.len = zeros(1,opt.maxit); info.len(1) = b-a;
for ii=1:opt.maxit-1
    c  = (a+b)/2;
    cl = c-epsilon;
    cr = c+epsilon;
    lamcl = func(cl); lamcr = func(cr);
    if lamcl<lamcr
        a = cl;
    elseif lamcl>lamcr
        b = cr;
    else
        a = cl;
        b = cr;
    end
    info.len(ii+1) = b-a;
    if info.len(ii+1) <= opt.tol
        info.len(ii+2:end) = [];
        break;
    end
end
mu = (a+b)/2;
[x,lambda] = eig(A-mu*C,'vector');
[~,ind] = sort(lambda);
lambda = lambda(ind(1));
X0 = x(:,ind(1:2));