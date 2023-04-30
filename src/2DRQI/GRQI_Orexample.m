function [x,mu, lambda, info]=GRQI_Orexample(Ln,Bn,opt)
% use 2d Generalized Rayleigh Quotient Iteration method to solve equation
%            A-mu C-\lambda I=0,
%            x'*C*x=0,
%            x'*x=1.
% Inputs:
%            Ln,Bn: A=Ln\Bn
%            opt: additionally specifies options using a structure
%            opt.tol: the converged tolerance: backward error < opts.tol
%            opt.maxit: the maxit number of iteration steps
%            opt.x0,mu0,lambda0: the initial solution
%            opt.normA,opt.normC: the norm of A and C, it can be any norm
%            
% Outputs:
%            x,mu,lambda: solutions to the 2d EVP
%            info: iteration history using a structure
%            including info.x,info.xCx,info.backerror,info.resiA,
%            info.lambda,info.mu
%            info.converged: indicates whether the algorithm converges
% Comparing with the Newton method, GRQI can be used to solve
% complex (Hermitian) problem
if nargin<3
    error('In this version, opt must be specified for attributes n, normA and normC');
end
nn = size(Ln,2);
R = chol(-Ln);
funA = @(x) [-R\((R')\(Bn*x(nn+1:end,:)));-Bn'*(R\((R')\x(1:nn,:)))];
funC = @(x) [1i*x(nn+1:end,:);-1i*x(1:nn,:)];
n = 2*nn;
if ~isfield(opt,'x0')
    x0=randn(n,1);
    x0=x0/norm(x0);
else
    x0=opt.x0;
end
if ~isfield(opt,'maxit')
    maxit=15;
else
    maxit=opt.maxit;
end

if  ~isfield(opt,'tol')
    tol=n*eps;
else
    tol=opt.tol;
end
info.x=zeros(n,maxit+1);
x=x0;x=x/norm(x);
info.x(:,1)=x;
info.backerror=zeros(maxit,1);
Ax=funA(x);
Cx=funC(x);
pern = (1:(n+2)); pern(1:2:n) = 1:nn; pern(2:2:n) = (1:nn)+nn;

if ~isfield(opt,'mu0') || ~isfield(opt,'lambda0')
    temp=[Cx,x]\Ax;
    mu=temp(1);
    lambda=temp(2);
else
    mu=opt.mu0;
    lambda=opt.lambda0;
end
info.converged=false;
normA = opt.normA;
normC = opt.normC;
for it=1:maxit
    resi=Ax-mu*Cx-lambda*x;
    info.resi(it)=norm(resi);
    info.xCx(it)=abs(x'*Cx);
    info.resiA(it)=abs(x'*Ax-lambda);
    info.mu(it)=mu;
    info.lambda(it)=lambda;
    
    %check convergence by backward error
    bkerror=[sqrt(2)*info.resi(it)/normA,...
        abs(imag(x(1:nn)'*x(nn+1:end)))];
    info.backerror(it)=max(bkerror);
    if (info.backerror(it)<tol || (it>=3 && info.backerror(it) >= mean(info.backerror(it-2:it))))
        info.converged=true;
        info.numit=it;
        info.x(:,it+1:end)=[];
        info.backerror(it+1:end)=[];
        return;
    end
    
    %GRQI step
%     L13 = [Cx,x]; L13(1:nn,:) = Ln*L13(1:nn,:); L13 =-normA*L13;
%     leftmatrix = [-lambda*Ln^2, Bn-mu*1i*Ln,L13(1:nn,:);
%                   Bn'+mu*1i*Ln, -lambda*speye(nn),L13(nn+1:end,:);
%                   L13',sparse(2,2)];
   leftmatrix = [-lambda*Ln^2, Bn-mu*1i*Ln, -normA*Ln*Cx(1:nn),-normA*Ln*x(1:nn);
             Bn'+mu*1i*Ln, -lambda*speye(nn),-normA*Cx(nn+1:end),-normA*x(nn+1:end);
             -normA*Cx(1:nn)'*Ln,-normA*Cx(nn+1:end)',0,0;
             -normA*x(1:nn)'*Ln,-normA*x(nn+1:end)',0,0];
    %rightterms = zeros(n+2,2); rightterms(n+1,1) = 1; rightterms(n+2,2) = 1;
    rightterms = sparse([n+1;n+2],[1;2],[1;1],n+2,2);
    leftmatrix = leftmatrix(pern,pern);
%    rightterms = rightterms(pern,:);   

 

    Ainvc = leftmatrix(1:n,1:n)\leftmatrix(1:n,n+1:n+2);
    cAinvc = leftmatrix(n+1:n+2,1:n)*Ainvc;    
    F = Ainvc/cAinvc;
%    [LL,UU, pp, qq] = lu(leftmatrix,'vector');
%    F = UU \ (LL \ (rightterms(pp,:)));
%    F = F(qq,:);

%    [LL,UU, pp, qq] = lu(leftmatrix,[1,1]);
%    pp = sparse(pp); qq = sparse(qq);
%    F  = qq*((UU\(LL\(pp*rightterms))));


    F = [F(1:2:n,:);F(2:2:n,:)];
    F(1:nn,:) = Ln*F(1:nn,:);
    Q = orth(full(F(1:n,:)));

    C_half=funC(Q);
    C_k=Q'*C_half;
    [V,e]=eig(C_k,'vector');
    if e(1)*e(2)>=0
        disp('Different signs!');
        if abs(e(1))>=abs(e(2))
            z=V(:,2);
        else
            z=V(:,1);
        end
        x=Q*z;
        Cx=funC(x);
        Ax=funA(x);
        temp=[Cx,x]\Ax;
        mu=real(temp(1));
        lambda=real(temp(2));
    else
        z=1./sqrt(abs(e));
        z=z/norm(z);
        T0=Q*V;
        A_half=funA(T0);
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
        %if 1-abs(x'*x1) <= 1-abs(x'*x2) 
            x=x1;
            x(1:nn) = (x(1:nn)/norm(x(1:nn)))/sqrt(2);
            x(nn+1:end) = (x(nn+1:end)/norm(x(nn+1:end)))/sqrt(2);
            mu=mu1;
            lambda=lam1;
            Cx=funC(x);
            Ax=funA(x);
        else
            x=x2;
            x(1:nn) = (x(1:nn)/norm(x(1:nn)))/sqrt(2);
            x(nn+1:end) = (x(nn+1:end)/norm(x(nn+1:end)))/sqrt(2);
            mu=mu2;
            lambda=lam2;
            Cx=funC(x);
            Ax=funA(x);
        end
    end
    %update Ax and Cx
    info.x(:,it+1)=x;
end
end
