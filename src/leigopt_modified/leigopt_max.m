function [f,z,iter] = leigopt_max(fname,d,j,A,pars)


l = length(A);
[n,m] = size(A{1});
sp = issparse(A{1});



if isfield(pars,'sp')
   sp = pars.sp; 
end

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
	gamma = -2*norm(full(A{1}));
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


if sp
	opts.tol = 10^-2;
	opts.issym = 1;
	opts.maxit = 30000;
	opts.p = round(sqrt(n));
end


V = [];
[~,ninit] = size(z);
for k = 1:ninit
	%[M,~] = feval(fname,z(:,k),A);
    M  = (A{1}*exp(1i*z) + A{1}'*exp(-1i*z))/2;
	if sp
		[Vt,D] = eigs(M,j,'LR',opts);
	else
		[Vt,D] = eig(M);
		[~,ind] = sort(diag(real(D)),'descend');
		Vt = Vt(:,ind(1:j));
	end

	V = [V Vt];
end

[P,~] = qr(V,0);

if sp
	opts.v0 = V(:,j);
end

iter = 1;
oldf = 0;
f = 1;
Id = eye(d);

z = z(:,1);


while (iter <= maxit) & ((f - oldf) > tol)

	clear parsin;



	%%%%%%%%%%%
	% FORM THE REDUCED PROBLEM
	%%%%%%%%%%%
	for k = 1:l 
		parsin.A{k} = full(P'*A{k}*P);
	end
	%%%%%%%%%%%
	%%%%%%%%%%%



	%%%%%%%%%%%
	% SOLVE THE REDUCED PROBLEM
	%%%%%%%%%%%
	parsin.tol = 0.1*tol;
	parsin.gamma = gamma;
	parsin.itertol = 100000;
	parsin.fname = fname;
	parsin.minmax = 1;
	parsin.j = j;
	parsin.d = d;

	oldz = z;
	oldf = f;
	[f,z,parsout] = eigopt('my_generic_fun',bounds,parsin);
						   
	if (iter == 1)
		oldf = f - 1;
	end
	%%%%%%%%%%%%
	%%%%%%%%%%%%
	



	if isprint
		fprintf('** iter: %d, f: %.14f z: %.14f ** \n', iter, f, z);
	end


	if (iter == 1) | (f - oldf > tol)
		[M,~] = feval(fname,z,A);


		if sp
			if isreal(M)
				opts.v0 = real(opts.v0);
			end
						   
			if abs(f - oldf) > 10^-2
				[V,D] = eigs(M,j,'LR',opts);
			else
				[V,D] = eigs(M,j,f+5*abs(f - oldf),opts);
			end
		else
			[Vt,D] = eig(M);
			[~,ind] = sort(diag(real(D)),'descend');
			V = Vt(:,ind(1:j));
		end

						   
		if sp
			opts.v0 = V(:,j);
		end



		%%%%%%%%%%%%%%
		% ALWAYS FRAMEWORK 2
		% (except when d=1, use framework 1 with last two iterations)
		%%%%%%%%%%%%%%
		h=norm(z-oldz);

		if (d > 1)
			for p=1:d
				ep=Id(:,p);

				for q=p:d
					eq=Id(:,q);
					if p==q
						epq=ep;
					else
						epq=(1/sqrt(2))*(ep+eq);
					end

					z1 = z + h*epq;


					[M,~] = feval(fname,z1,A);

						   
					if sp
						if abs(f - oldf) > 10^-2
						   [V1,D] = eigs(M,j,'LR',opts);
						else
						   [V1,D] = eigs(M,j,f+5*abs(f - oldf),opts);
						end
					else
						[Vt,D] = eig(M);
						[~,ind] = sort(diag(real(D)),'descend');
						V1 = Vt(:,ind(1:j));
					end
						   
					V = [V V1];

				end
			end
		end
		%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%


		%%%%%%%%%%%%%%%%%
		% FORGET THE PAST
		% AND FORM THE NEW SUBSPACE
		%%%%%%%%%%%%%%%%%
		if (d > 1)
			[P,~] = qr(V, 0);
		else
			[P,~] = qr([V P(:,1:j)],0);
		end
	end



	if sp & (iter > 2) & (abs(f - oldf) < 10^-2)
		opts.tol = max(abs(f - oldf)^2,0.1*tol);
	end


	iter = iter+1;

end



if oldf > f
	f = oldf;
	z = oldz;
end



return;