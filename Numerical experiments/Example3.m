%% Add path
addpath ../utils
addpath ../src/2DRQI
addpath ../src/dichotomous

%% Clear and Format
clear; clc;
format compact; format longe;

%% Parameter setting
N         = 100;   % Dimension of the test matrices is n^2
testtimes = 5;     % Number of test for timing per example/method

%%
t_Ga   = 0;
t_GRQI1 = 0;
t_GRQI2 = 0;

niter_Ga = 0;
niter_GRQI = 0;
nstep = 0;
for jj=1:testtimes
    K = 2;
    n = N^2;
    
    % Generate the test matrix
    signoise = 10^(-10/10); signdk = 10^(-10/10);  SINR_require = 10^(3/10);
    h1  = rand(N,1)/sqrt(2)+1i*rand(N,1)/sqrt(2); h2 = rand(N,1)/sqrt(2)+1i*rand(N,1)/sqrt(2);
    g1  = rand(N,1)/sqrt(2)+1i*rand(N,1)/sqrt(2); g2 = rand(N,1)/sqrt(2)+1i*rand(N,1)/sqrt(2);
    G1  = g1*g1'; G2 = g2*g2';
    H11 = h1*g1'; H12 = h2*g1'; H21 = h1*g2'; H22 = h2*g2';
    a11 = H11';   a12 = H12';   a21 = H21';   a22 = H22';
    a11 = a11(:); a12 = a12(:); a21 = a21(:); a22 = a22(:);
    Htilde1 = h1*h1';             Htilde2 = h2*h2';
    
    Htmp = Htilde1.'+Htilde2.'+signoise*eye(N); Htmp = (Htmp+Htmp')/2;
    Ltilde = chol(Htmp);
    Lh2 = (Ltilde')\conj(h2); Lh1 = (Ltilde')\conj(h1); Ltildeinv = inv(Ltilde);
    P1_1 = ((Lh2*Lh2')+signoise*(Ltildeinv'*Ltildeinv)-Lh1*Lh1'/SINR_require)/(signdk^2);
    P1_2 = g1;
    
    P2_1 = (Lh1*Lh1'+signoise*(Ltildeinv'*Ltildeinv)-Lh2*Lh2'/SINR_require)/(signdk^2);
    P2_2 = g2;
        
    funA  = @(u) A*u;
    funAs = @(u) kron(P1_1*(P1_2'*reshape(u,N,N)).', P1_2);
    funBs = @(u) kron(P2_1*(P2_2'*reshape(u,N,N)).', P2_2);
    
    % Dichotomous method
    tic;
    mubd = [0,1];
    opt1.bd = mubd;
    opt1.tol = 1e-4;
    opt1.n = N^2;
    opt1.maxit = 5e2;
    [mustar,lamstar,~,info] = dichotomous_largescale(funAs,funBs,opt1);
    t_Ga = t_Ga + toc;
    niter_Ga = niter_Ga + length(info.len)-1;
    
    tic;
    opt2.maxit = 2e1;
    opt2.normA = max(abs(P1_2))*norm(P1_2,1)*norm(P1_1,1);
    opt2.normC = max(abs(P2_2))*norm(P2_2,1)*norm(P2_1,1)+opt2.normA;
    opt2.n = N^2;
    opt2.tol = opt2.normA*eps;
    opteigs.isreal = false;
    opteigs.issym = 1;
    [x0,lambda0] = eigs(funAs,n,1,'sr',opteigs);
    if real(lambda0) >= real(x0'*funBs(x0))
        muGRQI = 0;
        lambdaGRQI = real(lambda0);
        continue;
    end
    [x0,lambda0] = eigs(funBs,n,1,'sr',opteigs);
    if real(lambda0) >= real(x0'*funAs(x0))
        muGRQI = 1;
        lambdaGRQI = real(lambda0);
        continue;
    end
    t_GRQI1 = t_GRQI1 + toc;
    
    
    tic
    [ muGRQI,lambdaGRQI,flag,step,info_GRQI,niter ] = minimaxRay( funAs, funBs, n,opt2 );
    t_GRQI2 = t_GRQI2 + toc;
    nstep = nstep + step;
    niter_GRQI = niter_GRQI + niter;
    fprintf('%d: %.2e, %.2e\n' ,jj,abs(muGRQI-mustar),abs(muGRQI-mustar)/abs(mustar))
end

fprintf('Average timing for dichotomous: %.3fs\n',t_Ga/testtimes);
fprintf('Average iterative steps for dichotomous£º %.3f\n',niter_Ga/testtimes);
fprintf('Average timing for checking in GRQI is %.3fs, Average timing for GRQI without checking %.4fs\n',t_GRQI1/testtimes, t_GRQI2/testtimes);
fprintf('Average iterative steps for GRQI: %d\n',niter_GRQI/testtimes);

