function [mu,lambda,X0,info] = dichotomous_largescale(funAs,funBs,opt)
% Find the maximum minium eigenvalue of A - mu C, where
% Ax = funAs(x), Cx = Ax - funBs(x)
a = opt.bd(1);
b = opt.bd(2);
n = opt.n;
epsilon = opt.tol/3;
opt_eig.issym = 1;
opt_eig.isreal = false;
func = @(mu) (eigs(@(u) (1-mu)*funAs(u)+mu*funBs(u),n,1,'SR',opt_eig));
info.len = zeros(1,opt.maxit); info.len(1) = b-a;
for ii=1:opt.maxit-1
    c  = (a+b)/2;
    cl = c-epsilon;
    cr = c+epsilon;
   lamcl = func(cl); lamcr = func(cr);
   lamcl = real(lamcl); lamcr = real(lamcr);
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
[X0,lambda] = eigs(@(u) (1-mu)*funAs(u)+mu*funBs(u),n,2,'SR',opt_eig);
if abs(lambda(1,1)-lambda(2,2))>abs(min(real(diag(lambda))))*1e-8
    [lambda,ind] = min(real(diag(lambda)));
    X0 = X0(:,ind)/norm(X0(:,ind));
else
    lambda = min(real(diag(lambda)));
    XCX = X0'*[funAs(X0(:,1))-funBs(X0(:,1)),funAs(X0(:,2))-funBs(X0(:,2))];
    [V0,e0]=eig(XCX,'vector');
    if e0(1)*e0(2)>0
        [~,ind] = min(abs(e0));
        X0 = X0*V0(:,ind);
        X0 = X0/norm(X0);
    else
        z0=1./sqrt(abs(e0));
        z0=z0/norm(z0);
        X0 = X0*(V0*z0);
    end
end
