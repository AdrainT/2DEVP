function [ muGRQI,lambdaGRQI,flag,i,info_GRQI,niter ] = minimaxRay_small( A,C,opts )
%   Use 2D RQI to solve the following minimax problem with small scale:
%      min_{x\neq 0} max{x^HAx/x^Hx, x^HBx/x^Hx}
%      transform to max\mu\in[0,1] \lambda_min(A-mu C), C = A-B;



interval = [0,1];
mu0 = mean(interval);
flag = 0;
niter = 0;
opteigs.isreal = false;

for i = 1:10
    % Prepare initial 2D evl
    opts.mu0     = mu0;
    funtmp = A-mu0*C;
    [X0,D] = eig(funtmp,'vector');
    X00 = X0;
    [~,ind] = sort(real(D),'ascend');
    % Prepare initial 2D eigenvector
    X0 = orth(X0(:,ind(1:2)));
    XCX = X0'*C*X0;
    [V0,e0]=eig(XCX,'vector');
    if e0(1)*e0(2)>0
        disp('Error!');
    end
    z0=1./sqrt(abs(e0));
    z0=z0/norm(z0);
    T0 = X0*V0;
    A_k=T0'*A*T0;
    C_k=diag(e0);
    if A_k(2,1)~= 0
        alpha=A_k(2,1)/abs(A_k(2,1));
    else
        alpha = 1;
    end
    z1=[z0(1);z0(2)*alpha]; z2=[z0(1);-z0(2)*alpha];
    x1=T0*z1;x2=T0*z2;
    mu1=real(z1'*C_k*A_k*z1/norm(C_k*z1)^2);
    mu2=real(z2'*C_k*A_k*z2/norm(C_k*z2)^2);
    lam1=real(z1'*A_k*z1);
    lam2=real(z2'*A_k*z2);
    if lam1<lam2
        opts.x0 = x1;
    else
        opts.x0 = x2;
    end
    % Calculate
    [~,muGRQI, lambdaGRQI, info_GRQI] = GRQI(A,C,opts);
    niter = niter + length(info_GRQI.mu)-1;
    % Check correctness
    if info_GRQI.backerror(end)<1e-8
        funtmp = A-muGRQI*C;%(1-muGRQI)*A + muGRQI*B;
        minlam = min(real(eig(funtmp)));
        if abs(minlam-lambdaGRQI)<1e-8
            flag = 1;
            break;
        end
    end
    % Update interval and next mu0 
    dlam = -real(X00(:,ind(1))'*C*X00(:,ind(1)));
    if dlam>0
        interval(1) = mu0;
    else
        interval(2) = mu0;
    end
    mu0 = (interval(1)+interval(2))/2;
end
