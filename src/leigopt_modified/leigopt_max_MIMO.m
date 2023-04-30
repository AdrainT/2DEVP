function [f,z,tSub,iter] = leigopt_max_MIMO(funA,funB,d,n,pars)


%l = length(A);
%sp = issparse(A{1});

if isfield(pars,'bounds')
	bounds = pars.bounds;
else
	bounds.lb = -60*ones(d,1);
	bounds.ub = 60*ones(d,1);
end


if isfield(pars,'maxit')
	maxit = pars.maxit;
else
	maxit = round(sqrt(n));
end


if isfield(pars,'tol')
	tol = pars.tol;
else
	tol = 10^-8;
end


if isfield(pars,'gamma')
	gamma = pars.gamma;
else
	gamma = -2;
end

if isfield(pars,'z0')
	z = pars.z0;
else
	z = zeros(d,1);
end

if isfield(pars,'isprint')
	isprint = pars.isprint;
else
	isprint = 0;
end

opts.tol = 10^-2;
opts.issym = 1;
opts.maxit = 30000;
%opts.p = round(sqrt(n));
opts.isreal = false;
tSub  = 0;
teigs = 0;
V = [];
[~,ninit] = size(z);
for k = 1:ninit
    funM = @(x) (1-z)*funA(x) + z*funB(x);
%    tic;
    [Vt,D] = eigs(funM,n,1,'SR',opts); 		   
%    teigs = teigs + toc;
	V = [V Vt];
end

[P,~] = qr(V,0);


%opts.v0 = V(:,end);


iter = 1;
oldf = 0;
f = 1;
fU = inf;
fL = -inf;
z = z(:,1);


while (iter <= maxit) && (fU-fL>=tol) %((f - oldf) > tol)

	clear parsin;



	%%%%%%%%%%%
	% FORM THE REDUCED PROBLEM
	%%%%%%%%%%%
    parsin.A{1} = zeros(size(P,2));
    for ii=1:size(P,2)
        parsin.A{1}(:,ii) =P'*funA(P(:,ii));
    end    
    parsin.A{2} = zeros(size(P,2));
    for ii=1:size(P,2)
        parsin.A{2}(:,ii) =P'*funB(P(:,ii));
    end
    
    %%%%%%%%%%%
	%%%%%%%%%%%



	%%%%%%%%%%%
	% SOLVE THE REDUCED PROBLEM
	%%%%%%%%%%%
	parsin.tol = 0.1*tol;
	parsin.gamma = gamma;
	parsin.itertol = 100000;
	parsin.fname = 'funMIMOAB';
	parsin.minmax = 1;
	parsin.d = d;

	oldz = z;
	oldf = f;
    parsin.minmax = 1;
%    tic;
	[f,z,parsout] = eigopt('funMIMO',bounds,parsin);
%    tSub = tSub + toc;
    fU = min([fU,f]);		   
	if (iter == 1)
		oldf = f - 1;
	end
	%%%%%%%%%%%%
	%%%%%%%%%%%%
	if isprint
		fprintf('** iter: %d, f: %.14f z: %.14f ** \n', iter, f, z);
    end
	if (iter == 1) || abs(f - oldf) > tol
        funM = @(x) (1-z)*funA(x) + z*funB(x);
        				   
%		if abs(f - oldf) > 10^(-2)
%        tic;
		[V,D] = eigs(funM,n,1,'SR',opts);
%        teigs = teigs + toc;
%		else
%			[V,D] = eigs(funM,n,1,f+5*abs(f - oldf),opts);
%		end
		fL = max([fL,D]);
						   
%		if sp
			opts.v0 = V(:,end);
%		end

    	[P,~] = qr([V P],0);
	end



	if (iter > 2) && (abs(f - oldf) < 10^-2)
		opts.tol = max(abs(f - oldf)^2,0.1*tol);
	end


	iter = iter+1;

end



%if oldf > f
%	f = oldf;
%	z = oldz;
%end



return