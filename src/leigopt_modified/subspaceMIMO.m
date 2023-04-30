function [z,f,tSub,teigs] = subspaceMIMO(funA,funB,n,tol)
pars.tol	=	tol;
pars.maxit	=	50;
pars.bounds.lb	=	0;
pars.bounds.ub	=	1;
pars.z0 = 0.5;
pars.isprint	=	1;	%	(set	this	1	to	print,	default	is	0	not	print)
pars.gamma	=	-1e-6;
[f,z,tSub,teigs]	=	leigopt_max_MIMO(funA,funB,1,n,n,pars);
end