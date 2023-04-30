function [x,mu, lambda, info]=GRQI(A,C,opt)
% use 2d Generalized Rayleigh Quotient Iteration method to solve equation
%            A-mu C-\lambda I=0,
%            x'*C*x=0,
%            x'*x=1.
% Inputs:
%            A,C: matrices in the 2d EVP
%            opt: additionally specifies options using a structure
%            opt.tol: the converged tolerance: backward error < opts.tol
%            opt.maxit: the maxit number of iteration steps
%            opt.x0,mu0,lambda0: the initial solution
% Outputs:
%            x,mu,lambda: solutions to the 2d EVP
%            info: iteration history using a structure
%            including info.x,info.xCx,info.backerror,info.resiA,
%            info.lambda,info.mu
%            info.converged: indicates whether the algorithm converges
% Comparing with the Newton method, GRQI can be used to solve
% complex (Hermitian) problem
% Select triples with distance
% (u-u0)^2+(lam-lam0)^2

n=size(A,1);
if nargin<3 || ~isfield(opt,'x0')
    x0=randn(n,1);
    x0=x0/norm(x0);
else
    x0=opt.x0;
end
if  nargin<3 || ~isfield(opt,'maxit')
    maxit=10;
else
    maxit=opt.maxit;
end

if  nargin<3 || ~isfield(opt,'tol')
    tol=n*eps;
else
    tol=opt.tol;
end
info.x=zeros(n,maxit+1);
x=x0;x=x/norm(x);
info.x(:,1)=x;
info.backerror=zeros(maxit,1);
Ax=A*x;
Cx=C*x;
if nargin<3 ||~isfield(opt,'mu0') || ~isfield(opt,'lambda0')
    temp=[Cx,x]\Ax;
    mu=temp(1);
    lambda=temp(2);
else
    mu=opt.mu0;
    lambda=opt.lambda0;
end
mu0 = mu; lambda0 = lambda;
info.converged=false;
normA=norm(A,1);
normC=norm(C,1);
I = speye(n);
for it=1:maxit
    resi=Ax-mu*Cx-lambda*x;
    info.resi(it)=norm(resi);
    info.xCx(it)=abs(x'*Cx);
    info.resiA(it)=abs(x'*Ax-lambda);
    info.mu(it)=mu;
    info.lambda(it)=lambda;
    
    %check convergence by backward error
    bkerror=[info.resi(it)/(normA+abs(mu)*normC),...
        info.xCx(it)/normC,...
        info.resiA(it)/normA];
    info.backerror(it)=max(bkerror);
    if info.backerror(it)<tol
        info.converged=true;
        info.numit=it;
        info.x(:,it+1:end)=[];
        info.backerror(it+1:end)=[];
        return;
    end

    H = [A-mu*C-lambda*I,-Cx,-x;-Cx',0,0;-x',0,0];
    [L_M,U_M,P_M]=lu(H);
    P_M=sparse(P_M);
    F=U_M\(L_M\(P_M*sparse([n+1;n+2],[1;2],[1;1],n+2,2) ));
    
    Q=orth(full(F(1:n,:)));
    %[Q,~] = qr(full(F(1:n,:)),0);
    C_half=C*Q;
    C_k=Q'*C_half;
    [V,e]=eig(C_k,'vector');
    if e(1)*e(2)>=0
        [~,ind] = min(abs(e));
        x   = Q*V(:,ind);
        Cx  = C_half*V(:,ind);
        Ax  = A*x;
        evls = [Cx,x]\Ax;
        evls = real(evls);
        mu = evls(1);
        lambda = evls(2);
    else
        z=1./sqrt(abs(e));
        z=z/norm(z);
        T0=Q*V;
        A_half=A*T0;
        A_k=T0'*A_half;
        C_half=C_half*V;
        C_k=V'*C_k*V;
        if A_k(2,1)~= 0
            alpha=A_k(2,1)/abs(A_k(2,1));
        else
            alpha = 1;
        end
        z1=[z(1);z(2)*alpha]; z2=[z(1);-z(2)*alpha];
        x1=T0*z1;x2=T0*z2;
        mu1=real(z1'*C_k*A_k*z1/norm(C_k*z1)^2); lam1=real(z1'*A_k*z1);
        mu2=real(z2'*C_k*A_k*z2/norm(C_k*z2)^2); lam2=real(z2'*A_k*z2);
        if abs(mu1-mu)+abs(lam1-lambda) <= abs(mu2-mu)+abs(lam2-lambda)
            x=x1;
            mu=mu1;
            lambda=lam1;
            Cx=C_half*z1;
            Ax=A_half*z1;
        else
            x=x2;
            mu=mu2;
            lambda=lam2;
            Cx=C_half*z2;
            Ax=A_half*z2;
        end
    end
    %update Ax and Cx
    info.x(:,it+1)=x;
end
end
