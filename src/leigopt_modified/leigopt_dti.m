function [beta,niter] = leigopt_dti(A,tol,mu0,maxit,~)
% An implementation of the subspace method by Kressner and Vandereycken for DTI of a stable 
% matrix A.
% Reference: Kressner D, and Vandereycken B. Subspace methods for computing the pseudospectral abscissa and the stability radius.
% SIAM J Matrix Anal Appl, 2014

if nargin<4
    pars.maxit = 100;
else
    pars.maxit = maxit;
end
% if nargin<3
%     claySolver = false;
% end
if nargin<2
    pars.tol	=	10^-12;
else
    pars.tol	=   tol;
end
% if strcmp(claySolver,'small')
%     [v0,e0] = eig(full(A),'vector');
%     [~,ind] = max(real(e0));
%     e0 = e0(ind);
% elseif ~claySolver
%     [v0,e0] = eigs(A, 1, 'lr');
% else
%     [v0,e0] = clay_Arnoldi(A, 1);
% end
% mu0 = imag(e0);


pars.bounds.lb	=	min(mu0-0.2*abs(mu0),mu0-5);
pars.bounds.ub	=	max(mu0+0.2*abs(mu0),mu0+5);
pars.z0	=	mu0;
pars.sq	=	[0	0	1	1	1]';
pars.gamma = -10;
n = size(A,1);
C{1} = A;
C{2} = speye(n);
C{3} = C{1}'*C{1};
C{4} = C{2}'*C{1};
C{5} = C{2}'*C{2};
[beta,mu,niter]	=	lsvdminopt_min_general2('distinstab_general',1,1,C,pars);