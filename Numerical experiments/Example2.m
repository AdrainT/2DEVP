%% Add path
addpath ../utils
addpath ../src/2DRQI

%% Clear and Format
clear; clc;
format compact; format shorte;

%% Get matrix and plot Fig.3-Left
%----------- Problem data ----------------------
n           = 3;

A = [-0.7, 0.01, 0.2;
    0.01,  2,    0;
    0.2,   0,    0];

C = [0.3, 0.01, 0.2;
    0.01,  1,    0;
    0.2,   0,   -1];

nummu = 1000000;
muset = linspace(-1.5,1.5,nummu);
lamset = zeros(n,nummu);
for i=1:nummu
    lamset(:,i) = sort(real(eig(A-muset(i)*C)));
end
colorset0 = {'g-','b-','r-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-','k-'};
colorset2 = {'rx','gx','bx'};
for i=1:3
    colorset2 = [colorset2(:);colorset2(:)];
end

figure; hold on;
for i=1:n
    plot(muset,lamset(i,:),colorset0{i},'LineWidth',2.5);
end

% Calculate all 2D eigenpairs using brute force search
evlSet = [1;1];
% Find local extremes candidates
for i=1:n
    for j=2:nummu-1
        if ((lamset(i,j)<lamset(i,j+1)) && (lamset(i,j)<lamset(i,j-1))) || ((lamset(i,j)>lamset(i,j+1)) && (lamset(i,j)>lamset(i,j-1)))
            evlSet = [evlSet,[muset(j);lamset(i,j)]];
        end
    end
end
ind = [];
for i=1:size(evlSet,2)
    for j=i+1:size(evlSet,2)
        if norm(evlSet(:,i)-evlSet(:,j))<1e-5
            ind = [ind,j];
        end
    end
end
evlSet(:,ind) = [];
hold on
for i=1:size(evlSet,2)
    plot(evlSet(1,i),evlSet(2,i),colorset2{i},'MarkerSize',15,'LineWidth',2.5);
end
ylabel('\lambda');xlabel('\mu');
set(gca,'Fontsize',15)
grid on

evlSetP = ones(2,size(evlSet,2));
dmu = muset(2)-muset(1);
for i=2:size(evlSet,2)
    mucand  = evlSet(1,i);
    lamcand = evlSet(2,i);
    a = mucand - dmu;
    b = mucand + dmu;
    lams = sort(eig(A-mucand*C),'ascend');
    [~,ind0] = min(abs(lams-lamcand));
    e1 = sort(eig(A-(mucand-dmu)*C),'ascend');
    e2 = sort(eig(A-(mucand+dmu)*C),'ascend');
    e1 = e1(ind0);
    e2 = e2(ind0);
    if e1>lamcand && e2>lamcand
        flagMaximum = false;
    else
        flagMaximum = true;
    end
    
    for j=1:100
        mid = (a+b)/2;
        [V0,D0] = eig(A-mid*C,'vector');
        [~,ind] = sort(D0,'ascend');
        V0 = V0(:,ind(ind0));
        if -V0'*C*V0<0
            if flagMaximum
                b = mid;
            else
                a = mid;
            end
        else
            if flagMaximum
                a = mid;
            else
                b = mid;
            end
        end
    end
    mid = (a+b)/2;
    evlSetP(1,i) = mid;
    e = sort(eig(A-mid*C),'ascend');
    evlSetP(2,i) = e(ind0);
end
%
%% Plot Fig.3-Right
nummu = 100; numlam = 100;
muset2  = linspace(-1.5, 1.5, nummu);
lamset2 = linspace(-2, 2, numlam);
opt.maxit   = 50;
opt.tol     = norm(A)*eps;
I = eye(n);
ConvergeCase = zeros(nummu,numlam);
for i=1:nummu
    opt.mu0 = muset2(i);
    
    fprintf('Excuting i = %d ...\n', i);
    for j=1:numlam
        flagC = false;
        opt.lambda0 = lamset2(j);
        
        [X0,e] = eig(A-opt.mu0*C-opt.lambda0*I,'vector');
        [~,ind] = sort(abs(e),'ascend');
        X0 = orth(X0(:,ind(1:2)));
        XCX = X0'*C*X0;
        [V0,e0]=eig(XCX,'vector');
        if (e0(1)*e0(2)<0)
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
            if abs(lam1-opt.lambda0)+abs(mu1-opt.mu0) <= abs(lam2-opt.lambda0)+abs(mu2-opt.mu0)
                opt.x0 = x1;
            else
                opt.x0 = x2;
            end
        else
            [~,ind] = min(abs(e0));
            opt.x0 = X0*V0(:,ind);
        end
        [~,muGRQI, lambdaGRQI]=GRQI(A,C,opt);
        for k=1:size(evlSet,2)
            if (abs(muGRQI-evlSet(1,k))<=0.001) && (abs(lambdaGRQI-evlSet(2,k))<=0.001)
                flagC = true;
                break;
            end
        end
        if (flagC)
            ConvergeCase(i,j) = k;
        else
            ConvergeCase(i,j) = 0;
        end
    end
end
% Plot
figure; hold on;
colorset = [{'.k'};colorset2];
for i=1:nummu
    for j=1:numlam
        plot(muset2(i),lamset2(j),colorset{ConvergeCase(i,j)+1});
    end
end
hold on
for i=1:size(evlSet,2)
    plot(evlSet(1,i),evlSet(2,i),colorset2{i},'MarkerSize',10,'linewidth',2.5);
end
ylabel('\lambda');xlabel('\mu');
set(gca,'Fontsize',15)
grid on


%% Converge to (mu3, lam3)
mu3  = -0.145810069397437;
lam3 = -0.744080780565709;
opt.mu0 = 0.5;
opt.lambda0 = -0.5;
[V0,~] = eigs(A-opt.mu0*C-opt.lambda0*I,2,0);
Ck = V0'*C*V0;
[Q,e] = eig(Ck,'vector');
if (e(1)*e(2)<0)
    z0=1./sqrt(abs(e));
    z0=z0/norm(z0);
    T0 = V0*Q;
    A_k=T0'*A*T0;
    C_k=diag(e);
    if A_k(2,1)~= 0
        alpha=A_k(2,1)/abs(A_k(2,1));
    else
        alpha = 1;
    end
    z1_sub=[z0(1);z0(2)*alpha]; z2_sub=[z0(1);-z0(2)*alpha];
    x1_sub=T0*z1_sub;x2_sub=T0*z2_sub;
    mu1_sub=real(z1_sub'*C_k*A_k*z1_sub/norm(C_k*z1_sub)^2);
    mu2_sub=real(z2_sub'*C_k*A_k*z2_sub/norm(C_k*z2_sub)^2);
    lam1_sub=real(z1_sub'*A_k*z1_sub);
    lam2_sub=real(z2_sub'*A_k*z2_sub);
    if abs(lam1_sub-opt.lambda0)+abs(mu1_sub-opt.mu0) <= abs(lam2_sub-opt.lambda0)+abs(mu2_sub-opt.mu0)
        opt.x0 = x1_sub;
    else
        opt.x0 = x2_sub;
    end
else
    [~,ind] = min(abs(e));
    opt.x0 = V0*Q(:,ind);
end

%
[xGRQI,~, ~, info_GRQI]=GRQI_detailed(A,C,opt);

x_star  = xGRQI;
appinfpx = zeros(length(info_GRQI.mu),1);
appinfpmu = zeros(length(info_GRQI.mu),1);
appinfplam = zeros(length(info_GRQI.mu),1);
xend = x_star;
for ii=1:length(info_GRQI.mu)
    appinfpx(ii) = norm(info_GRQI.x(:,ii) - xend*(xend'*info_GRQI.x(:,ii)));
    appinfpmu(ii)=abs(info_GRQI.mu(ii)-mu3);
    appinfplam(ii)=abs(info_GRQI.lambda(ii)-lam3);
end
%
figure(3);
semilogy(0:(length(appinfpx)-1),appinfpx+eps,'-gp','Markersize',10,'LineWidth',2.5);hold on
semilogy(0:(length(appinfpx)-1),appinfpmu,'-rp','Markersize',10,'LineWidth',2.5);
semilogy(0:(length(appinfpx)-1),appinfplam,'-bp','Markersize',10,'LineWidth',2.5);
ylabel('Forward error');xlabel('Step');
set(gca,'Fontsize',15)
grid on
legend('x', '\mu','\lambda');
figure;
semilogy(0:length(info_GRQI.backerror)-1,info_GRQI.backerror,'-bp','Markersize',10,'LineWidth',2.5);
ylabel('Backward error');xlabel('Step');
set(gca,'Fontsize',15)
grid on
disp('a21:');
disp(info_GRQI.a21);
disp('C2x2:');
for i=1:length(info_GRQI.Cevl)
    disp(info_GRQI.Cevl{i});
end

fprintf('Convergence for mu3 = %d:\n',mu3);
disp(abs(info_GRQI.mu-mu3))
fprintf('Convergence for lambda3 = %d:\n',lam3);
disp(abs(info_GRQI.lambda-lam3))

%% Converge to (mu1, lam1)
mu1 = 1;
lam1 = 1;
opt.mu0 = 1.5;
opt.lambda0 = 2;

[V0,~] = eigs(A-opt.mu0*C-opt.lambda0*I,2,0);
Ck = V0'*C*V0;
[Q,e] = eig(Ck,'vector');
if (e(1)*e(2)<0)
    z0=1./sqrt(abs(e));
    z0=z0/norm(z0);
    T0 = V0*Q;
    A_k=T0'*A*T0;
    C_k=diag(e);
    if A_k(2,1)~= 0
        alpha=A_k(2,1)/abs(A_k(2,1));
    else
        alpha = 1;
    end
    z1_sub=[z0(1);z0(2)*alpha]; z2_sub=[z0(1);-z0(2)*alpha];
    x1_sub=T0*z1_sub;x2_sub=T0*z2_sub;
    mu1_sub=real(z1_sub'*C_k*A_k*z1_sub/norm(C_k*z1_sub)^2);
    mu2_sub=real(z2_sub'*C_k*A_k*z2_sub/norm(C_k*z2_sub)^2);
    lam1_sub=real(z1_sub'*A_k*z1_sub);
    lam2_sub=real(z2_sub'*A_k*z2_sub);
    if abs(lam1_sub-opt.lambda0)+abs(mu1_sub-opt.mu0) <= abs(lam2_sub-opt.lambda0)+abs(mu2_sub-opt.mu0)
        opt.x0 = x1_sub;
    else
        opt.x0 = x2_sub;
    end
else
    [~,ind] = min(abs(e));
    opt.x0 = V0*Q(:,ind);
end
[xGRQI,~, ~, info_GRQI]=GRQI_detailed(A,C,opt);

x1  = xGRQI;

appinfpx = zeros(length(info_GRQI.mu),1);
appinfpmu = zeros(length(info_GRQI.mu),1);
appinfplam = zeros(length(info_GRQI.mu),1);
xend = x1;
for ii=1:length(info_GRQI.mu)
    appinfpx(ii) = norm(info_GRQI.x(:,ii) - [0;sign(info_GRQI.x(2,ii))/sqrt(2);sign(info_GRQI.x(3,ii))/sqrt(2)]);
    appinfpmu(ii)=abs(info_GRQI.mu(ii)-mu1);
    appinfplam(ii)=abs(info_GRQI.lambda(ii)-lam1);
end
%
figure(3);
semilogy(0:(length(appinfpx)-1),appinfpx,'-gp','Markersize',10,'LineWidth',2.5);hold on
semilogy(0:(length(appinfpx)-1),appinfpmu,'-rp','Markersize',10,'LineWidth',2.5);
semilogy(0:(length(appinfpx)-1),appinfplam,'-bp','Markersize',10,'LineWidth',2.5);
ylabel('Forward error');xlabel('Step');
set(gca,'Fontsize',15)
grid on
legend('x', '\mu','\lambda');
figure;
semilogy(0:length(info_GRQI.backerror)-1,info_GRQI.backerror,'-bp','Markersize',10,'LineWidth',2.5);
ylabel('Backward error');xlabel('Step');
set(gca,'Fontsize',15)
grid on
disp('a21:');
disp(info_GRQI.a21);
disp('C2x2:');
for i=1:length(info_GRQI.Cevl)
    disp(info_GRQI.Cevl{i});
end

fprintf('Convergence for mu1 = %d:\n',mu1);
disp(abs(info_GRQI.mu-mu1))
fprintf('Convergence for lambda1 = %d:\n',lam1);
disp(abs(info_GRQI.lambda-lam1))