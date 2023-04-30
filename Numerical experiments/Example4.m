%% Add path
addpath ../utils
addpath ../src/2DRQI
addpath ../src/leigopt_modified
addpath ../src/Clay_Arnoldi

%% Clear and Format
clear; clc;
format compact; format longe;
%% Parameter setting
n         = 16000;  % Dimension of the test matrices
testtimes = 10;     % Number of test for timing per example/method
%% 2DRQI
[Ln,Bn]   = Orr_generator(n,1);

% Get initial point
tic;
for iii =1:testtimes
    Ln_R = chol(-Ln);
    funLinv = @(x) -Ln_R\((Ln_R')\x);
    funA = @(x) [  Ln*x(n+1:end,:)/100-1i*funLinv(  Un*(Ln*x(n+1:end,:)) + 2*x(n+1:end,:) )  ;
        Ln*x(1:n,:)/100+1i*( Ln*(Un*funLinv(x(1:n,:)))+2*funLinv(x(1:n,:)) )  ];
    funC = @(x) [1i*x(n+1:end,:);-1i*x(1:n,:)];    
    [ V,D ] = clay_Arnoldi2( Ln,Bn,1);
    [~,ind] = sort(real(D),'descend');
    D = D(ind(1));
    opt.mu0      = imag(D);
    opt.maxit    = 15;
    opt.n        = 2*n;
    opt.normC    = 1;
    
    opts.maxit = 100;opts.tol = 1e-6;opts.m = 6; 
    opts.isreal = false;
    [LL,UU,PP] = lu(Bn-opt.mu0*1i*Ln); PP = sparse(PP);
    funBLv = @(x) UU \ (LL \ (PP*(Ln*x)));
    PPv = PP';
    funBLTv = @(x) Ln*(PPv*( (LL')\((UU')\x) ));
    func = @(x,flag) svdfunc(funBLv,funBLTv,x,flag);
    [u,opt.lambda0,v] = svds(func,[n,n],1,'largest',opts);
    opt.lambda0 = 1/opt.lambda0;
    opt.x0    = [v;u]/sqrt(2);
end
t1 = toc;

% Run iterator
nGRQI = 0;
tic;
for iii = 1:testtimes
    % Estimate norm A
    normB        = sum(abs(Bn),1);
    [~,ind] = sort(normB,'descend');
    opt.normA    = 0;
    for ii = 1:5
        opt.normA = max(opt.normA,norm(funLinv(Bn(:,ind(ii))),1)  );
    end
    for ii = 1:5
        opt.normA = max(opt.normA,norm(funLinv(Bn(:,randi(n))),1)  );
    end
    opt.tol       = opt.normA*eps;    
    
    % Run
    [~,muGRQI, lamGRQI, infoGRQI]=GRQI_OrexampleN(Ln,Bn,opt);
    nGRQI = nGRQI + length(infoGRQI.backerror)-1;
end
t2 = toc;

fprintf('The resulf of RQI is %.10e\n',abs(lamGRQI));
fprintf('The niter of RQI is %d\n',nGRQI/testtimes);
fprintf('Timing of RQI for initial and iteration is %.3gs, %.3gs\n',t1/testtimes,t2/testtimes);
%% leigopt
C{1} = Bn;
C{2} = Ln;
%%
niterSubspace = 0;

pars.tol	=	10^-12;
pars.bounds.lb	=	-60;
pars.bounds.ub	=	80;
pars.z0	=	(pars.bounds.lb + pars.bounds.ub)/2;
pars.sq	=	[0	0	1	1	1]';
pars.subfun = @subfun_orrsommerfeld;
pars.eigsfun = 'eigsfun_orrsommerfeld2';
pars.gamma = 0; % Not be used.

tsubs = 0;
tic
for i=1:testtimes
    [f,z,niter,tsub]	=	lsvdminopt_min_general('distinstab_general2',1,1,C,pars);
    niterSubspace = niterSubspace + niter;
    tsubs = tsubs + tsub;
end
tTotal = toc;
fprintf('The resulf of leigopt is %.10e\n', f);
fprintf('The niter of leigopt is %d\n', niterSubspace/testtimes);
fprintf('Timing of leigopt for projected problem and total is %.3gs, %.3gs\n',tsubs/testtimes,tTotal/testtimes);

