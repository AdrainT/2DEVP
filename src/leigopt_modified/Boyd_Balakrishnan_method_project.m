function [epsilon,omega]=Boyd_Balakrishnan_method_project(AV,V,opt)
%where
%       opt.maxit
%          .tol
% output epsilon is an estimation of distance of A
%reference: S. Boyd and V. Balakrishnan. A regularity result for the singular
% values of a transfer matrix and a quadratically convergent algorithm
% for computing its L-infty norm. Systems Control Lett., 15(1):, 1990.
%A is a stable matrix

Q  = orth([AV,V]);
AV = Q'*AV;
V  = Q'*V;
[m,n] = size(AV);
epsilon = min(svd(AV));
omega   = 0;


if  nargin<3||~isfield(opt,'tol')
    tol=1e-8;
else
    tol=opt.tol;
end

if  nargin<3||~isfield(opt,'maxit')
    maxit=50;
else
    maxit=opt.maxit;
end

iter=0;
err=1;
B = [V,zeros(m);zeros(n),V'];

VAAV = AV'*AV;
VAV  = V'*AV-AV'*V;
VV   = V'*V;

while err>tol && iter<maxit
    epsilon1=epsilon;
    H=[AV,-epsilon*eye(m);epsilon*eye(n),-AV'];
    lambda = eig(H,B);
    flag= abs(real(lambda))<1e-8;
    if isempty(find(flag, 1))
        break
    end
    lambda = sort(imag(lambda(flag)));      %find eigenvalues on imaginary axis, can use function eigenvalue_on_imaginary_axis
    numlam = length(lambda);
    % 找出区间
    if numlam==1
        return
    else
        sigma = zeros(numlam-1,1); index = 1;
        lambdamin = zeros(numlam-1,1);
        %        lambdamin = lambda; index = numlam + 1;
        for i=1:numlam-1
            lamtest = (lambda(i) + lambda(i+1))/2;
            %            [~,Apositive] = chol(VAAV+1i*lamtest*VAV+lamtest^2*VV-(epsilon^2)*(1+1e-8)*eye(n));
            sigmaTmp = min(svd(AV-1i*lamtest*V));
            
            lambdamin(index) = lamtest;
            sigma(index) = sigmaTmp;
            index = index + 1;
            
        end
    end
    if index==1
        disp(1);
    end
    [epsilon,ind]=min(sigma);
    omega = lambdamin(ind);
    err=abs(epsilon-epsilon1)/abs(epsilon1);
    iter=iter+1;
end