function [ muGRQI,lambdaGRQI,flag,i,info_GRQI, niter] = minimaxRay( A, B, n,opts )
%   Use 2D RQI to solve the following minimax problem with large scale:
%      min_{x\neq 0} max{x^HAx/x^Hx, x^HBx/x^Hx}
%      transform to max\mu\in[0,1] \lambda_min(A-mu C), C = A-B;
if isfloat(A)
    funAs = @(x) A*x;
else
    funAs = A;
end
if isfloat(B)
    funBs = @(x) B*x;
else
    funBs = B;
end


interval = [0,1];
mu0 = (interval(1)+interval(2))/2;
%mu0 = -1000;
flag = 0;
niter = 0;
opteigs.isreal = false;
opteigs.issym = 1;
for i = 1:100
    % Prepare initial 2D evl
    opts.mu0     = mu0;
    funtmp = @(x) (1-mu0)*funAs(x) + mu0*funBs(x);
    [X0,lambda0] = eigs(funtmp,n,2,'sr',opteigs);
    opts.lambda0 = min(real(diag(lambda0)));
    % Prepare initial 2D eigenvector
    X0 = orth(X0);
    XCX = X0'*[funAs(X0(:,1))-funBs(X0(:,1)),funAs(X0(:,2))-funBs(X0(:,2))];
    [V0,e0]=eig(XCX,'vector');
    if e0(1)*e0(2)>0
        disp('Error!');
        [~,ind] = min(abs(e0));
        opts.x0 = X0*V0(:,ind);
        opts.x0 = opts.x0/norm(opts.x0);
    else
        z0=1./sqrt(abs(e0));
        z0=z0/norm(z0);
        T0 = X0*V0;
        A_k=T0'*[funAs(T0(:,1)),funAs(T0(:,2))];
        if A_k(2,1)~= 0
            alpha=A_k(2,1)/abs(A_k(2,1));
        else
            alpha = 1;
        end
        z1=[z0(1);z0(2)*alpha]; z2=[z0(1);-z0(2)*alpha];
        x1=T0*z1;x2=T0*z2;
        lam1=real(z1'*A_k*z1);
        lam2=real(z2'*A_k*z2);
        if lam1<lam2
            opts.x0 = x1;
        else
            opts.x0 = x2;  
        end
    end
    % Calculate
    [~,muGRQI, lambdaGRQI, info_GRQI] = GRQI_largescale_MIMO(funAs,funBs,opts);
     niter = niter + length(info_GRQI.mu)-1;
    % Check correctness
    if info_GRQI.backerror(end)<1e-8
        funtmp = @(x) (1-muGRQI)*funAs(x) + muGRQI*funBs(x);
        minlam = real(eigs(funtmp,n,1,'sr',opteigs));
        if abs(minlam-lambdaGRQI)<1e-8
            flag = 1;
            break;
        end
    end
    % Update interval and next mu0
    [~,ind] = min(real(diag(lambda0)));
    dlam = -real(X0(:,ind)'*(funAs(X0(:,ind))-funBs(X0(:,ind))));
    if dlam>0
        interval(1) = mu0;
    else
        interval(2) = mu0;
    end
    mu0 = (interval(1)+interval(2))/2;
end
