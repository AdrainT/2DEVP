function [f,z,iter,tsub] = lsvdminopt_min_general2(fname,d,j,A,pars)
l = length(A);
[n,m] = size(A{1});
p = min(n,m);
sp = issparse(A{1});

if isfield(pars,'sp')
    sp = pars.sp; 
end


subfun = isfield(pars,'subfun');

warning off;

if isfield(pars,'bounds')
	bounds = pars.bounds;
else
	bounds.lb = -60*ones(d,1);
	bounds.ub = 60*ones(d,1);
end

if isfield(pars,'eigsfun')
	eigsfun = 1;
	pars.eigsfun = str2func(pars.eigsfun);
else
	eigsfun = 0;
end


if isfield(pars,'sq')
	matrixsq = pars.sq;
else
	matrixsq = zeros(l,1);
end


if isfield(pars,'maxit')
	maxit = pars.maxit;
else
	maxit = round(sqrt(max(m,n)));
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


if sp | eigsfun
	opts.tol = 10^-2;
	opts.maxit = 30000;
	opts.isreal = false;
%	opts.p = round(sqrt(n));
end


V = [];
[~,ninit] = size(z);
for k = 1:ninit

    if eigsfun
        [ALU.L,ALU.U,ALU.P] = lu(A{1}-1i*z(:,k)*A{2});
        ALU.P = sparse(ALU.P);
        [Vt,D,Ut] = svds(@(u,flag) pars.eigsfun(u,flag,ALU,A),[n,n],j,'largest',opts);     
    else
        M = feval(fname,z(:,k),A);
		if sp
			[Ut,D,Vt] = svds(M,j,0,opts);
            if isempty(D)
                M = full(M);
                sp = false;
                [Ut,D,Vt] = svd(M);
                Vt = Vt(:,p-j+1:p);
                disp('Need full');
            end
                
		else
			[Ut,D,Vt] = svd(M);
			Vt = Vt(:,p-j+1:p);
		end
	end

	V = [V Vt];
end


[P,~] = qr(V,0); 

if sp | eigsfun
	opts.v0 = V(:,j);
end

iter = 1;
oldf = 2;
f = 1;
Id = eye(d);
z = z(:,1);
tsub = 0;

while ((iter <= maxit) &  (oldf - f > tol))

	clear parsin;


	%%%%%%%%%%%
	% FORM THE REDUCED PROBLEM
	%%%%%%%%%%%
	for k = 1:l 
		if subfun
			parsin.A = feval(pars.subfun,A,P);
		else
			if matrixsq(k)
				parsin.A{k} = full(P'*A{k}*P);
			else
				parsin.A{k} = full(A{k}*P);
			end
		end
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
	parsin.j = j;
	parsin.d = d;

	oldf = f;						
	oldz = z;
    %tic;
	[f,z] = Boyd_Balakrishnan_method_project(parsin.A{1},P);
    %tsub = tsub + toc;	
								   
	if (iter == 1)
		oldf = f + 1;
	end
	%%%%%%%%%%%%
	%%%%%%%%%%%% 

	if isprint
		fprintf('** iter: %d, f: %.14f ** \n', iter, f);
	end
						   
			   
	if (iter == 1) | (oldf - f > tol) %| 1				   

		if eigsfun					   		
            if abs(f - oldf) > 1e-2
                   [ALU.L,ALU.U,ALU.P] = lu(A{1}-1i*z(:,1)*A{2});
                   ALU.P = sparse(ALU.P);
                   [V,D,U] = svds(@(u,flag) pars.eigsfun(u,flag,ALU,A),[n,n],j,'largest',opts);
            else
                   [ALU.L,ALU.U,ALU.P] = lu(A{1}-1i*z(:,1)*A{2});
                   ALU.P = sparse(ALU.P);
                   [V,D,U] = svds(@(u,flag) pars.eigsfun(u,flag,ALU,A),[n,n],j+1,'largest',opts);

            end
		else
			
			[M,~,~] = feval(fname,z,A);
								   
			if sp
				if isreal(M)
					opts.v0 = real(opts.v0);
				end

				if abs(f - oldf) > 10^-2
					[U,D,V] = svds(M,j,0,opts);
				else
					% keyboard
					[U,D,V] = svds(M,j,f-5*abs(f - oldf),opts);
				end
			else
				[Ut,D,Vt] = svd(M);
				V = Vt(:,p-j+1:p);
			end
		end
				   
		if sp | eigsfun
			opts.v0 = V(:,j);
		end

		%%%%%%%%%%%%%%
		% ALWAYS FRAMEWORK 2
		% (except when d=1, use framework 1 with last two iterations)
		%%%%%%%%%%%%%%
		h=norm(z-oldz);

		%%%%%%%%%%%%%%%%%
		% FORGET THE PAST
		% AND FORM THE NEW SUBSPACE
		%%%%%%%%%%%%%%%%%
		if (d > 1)
			[P,~] = qr(V,0);
		else
			%[P,~] = qr([V P(:,1:j)],0);
            [P,~] = qr([V P],0);
		end
		
    end
						   
	if (sp | eigsfun) & (iter > 2) & (abs(f - oldf) < 10^-2)
		opts.tol = max(abs(f - oldf)^2,0.1*tol);
	end

	iter = iter+1;

end			   

if oldf < f
	f = oldf;
	z = oldz;
end

						   

return;