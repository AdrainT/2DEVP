function [ V,D,res ] = clay_Arnoldi2( Ln,Bn,k,opts )
%   Calculate the rightmost evl using Clay-Arnoldi method
if nargin <= 3
    maxit = 100;
    tol = 1e-6;
    m = 2*k+4;
else
    if isfield(opts,'maxit')
        maxit = opts.maxit;
    else
        maxit = 100;
    end
    if isfield(opts,'tol')
        tol = opts.tol;
    else
        tol = 1e-8;
    end
    if isfield(opts,'m')
        m = opts.m;
    else
        m = 2*k+4;
    end
end
n = size(Ln,1);
m = min(m,n);
v = randn(n,1);
v = v/norm(v);

Ln_R = chol(-Ln);
funLinv = @(x) -Ln_R\((Ln_R')\x);
funv = @(x) funLinv(Bn*x);

Q = zeros(n,m+1);
Q(:,1)=v;
H=zeros(m+1,m);
for j=1:m
    w=funv(Q(:,j));
    for i=1:j
        H(i,j)=Q(:,i)'*w;
        w=w-H(i,j)*Q(:,i);
    end
    for i=1:j              %reorthogonalization
        s=Q(:,i)'*w;
        H(i,j)=H(i,j)+s;
        w=w-s*Q(:,i);
    end
    H(j+1,j)=norm(w,2);
    if H(j+1,j)<eps
        H(j+1,j) = 0;
        v = randn(n,1);
        for i=1:j
            tmp = Q(:,i)'*v;
            v   = v - tmp*Q(:,i);
        end
        Q(:,j+1) = v/norm(v);
    else
        Q(:,j+1)=w/H(j+1,j);
    end
end

[V,D] = eig(H(1:m,1:m),'vector');
D = 1./D;
[~,ind] = sort(real(D),'descend');

v = randn(n,1); v = v/norm(v);
res = zeros(maxit,1);
for count=1:maxit
    Q(:,1) = v;
    % determine sigma and tau
    sigma = D(ind(1));
    tau = 2*D(ind(k+1)) - sigma;
%    if issparse(Bn)
%        [LL,UU,PP,QQ,DD] = lu(Bn-sigma*Ln);
%        funv = @(x) x + (sigma-tau)*QQ*(UU \ (LL \ (PP*(DD\(Ln*x)))));
%    else
        [LL,UU,PP] = lu(Bn-sigma*Ln);
        funv = @(x) x + (sigma-tau)*(UU \ (LL \ (PP*(Ln*x))));
%    end
    for j=1:m
        w=funv(Q(:,j));
        for i=1:j
            H(i,j)=Q(:,i)'*w;
            w=w-H(i,j)*Q(:,i);
        end
        for i=1:j              %reorthogonalization
            s=Q(:,i)'*w;
            H(i,j)=H(i,j)+s;
            w=w-s*Q(:,i);
        end
        H(j+1,j)=norm(w,2);
        if H(j+1,j)<eps
            H(j+1,j) = 0;
            v = randn(n,1);
            for i=1:j
                tmp = Q(:,i)'*v;
                v   = v - tmp*Q(:,i);
            end
            Q(:,j+1) = v/norm(v);
        else
            Q(:,j+1)=w/H(j+1,j);
        end
    end
    [V,D] = eig(H(1:m,1:m),'vector');
    D = (sigma-tau)./(D-1)+sigma;
    [~,ind] = sort(real(D),'descend');
    v = Q(:,1:m)*V(:,ind(1));
    v = v/norm(v);
    res(count) = norm(V(m,ind(1:k))*H(m+1,m));
    if res(count) < tol
        break;
    end
end
res = res(1:count);
V = Q(:,1:m)*V(:,ind(1:k));
D = D(ind(1:k));
end
