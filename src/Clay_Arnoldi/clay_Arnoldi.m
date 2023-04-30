function [ V,D,res ] = clay_Arnoldi( A,k,opts )
% An implementation of the Clay-Arnoldi method using complex shift for computing rightmost evl of the matrix A.
% Reference: Meerbergen K, and Roose D. Matrix transformations for computing rightmost eigenvalues of large sparse non-symmetric
% eigenvalue problems. IMA J Numer Anal. 1996
% k is the number of wanted rightmost evls.
if nargin <= 2
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
        tol = 1e-10;
    end
    if isfield(opts,'m')
        m = opts.m;
    else
        m = 2*k+4;
    end
end
n = size(A,2);
m = min(m,n);
v = randn(n,1);
v = v/norm(v);
if issparse(A)
    [L,U,PA,QA,DA] = lu(A);
    funv = @(x) QA*(U \ (L \ (PA*(DA\x))));
else
    [L,U,PA] = lu(A);
    PA = sparse(PA);
    funv = @(x) U\(L\(PA*(x)));
end

Q = zeros(n,m+1);
Q(:,1)=v;
H=zeros(m+1,m);
for j=1:m
    w=funv(Q(:,j));
    wnorm0 = norm(w);
    for i=1:j
        H(i,j)=Q(:,i)'*w;
        w=w-H(i,j)*Q(:,i);
    end
    wnorm = norm(w);
    if wnorm<0.717*wnorm0
        for i=1:j              %reorthogonalization
            s=Q(:,i)'*w;
            H(i,j)=H(i,j)+s;
            w=w-s*Q(:,i);
        end
        wnorm = norm(w);
    end
    H(j+1,j)=wnorm;
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
%   sigma = D(ind(1));
%   tau   = 2*D(ind(k+1))-sigma;
 %   sigma = mean(D(ind(1:k)));
 %   tau   = 2*D(ind(k+1))-sigma;
    sigma = 0.99*D(ind(1))+0.01*D(ind(k+1));
    for kk = 1:m-k
        tau = real(2*D(ind(k+kk)) - sigma)+1i*imag(sigma); 
        if abs(tau-sigma)>(abs(real(tau))+abs(real(sigma)))*1e-2
            break;
        end
    end
    if issparse(A)
        [L,U,PA,QA,DA] = lu(A-sigma*speye(n));
        funv = @(x) x + (sigma-tau)*(QA*(U \ (L \ (PA*(DA\x)))));
    else
        [L,U,PA] = lu(A-sigma*eye(n));
        PA = sparse(PA);
        funv = @(x) x + (sigma-tau)*(U\(L\(PA*x)));
    end
    for j=1:m
        w=funv(Q(:,j));
        wnorm0 = norm(w);
        for i=1:j
            H(i,j)=Q(:,i)'*w;
            w=w-H(i,j)*Q(:,i);
        end
        wnorm = norm(w);
        if wnorm<0.717*wnorm0
        for i=1:j              %reorthogonalization
            s=Q(:,i)'*w;
            H(i,j)=H(i,j)+s;
            w=w-s*Q(:,i);
        end
        wnorm = norm(w);
        end
        H(j+1,j)=wnorm;
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
