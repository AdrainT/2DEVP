%% Add path
addpath ../utils
addpath ../src/2DEVP
addpath ../src/leigopt_modified


%% Clear and Format
clear; clc;
format compact; format longe;

DataDir  = '../Data/';
%% Parameter setting for small size problems
testtimes = 10;     % Number of test for timing per example/method
MatricesDir = [DataDir,'testmatrices_small/'];

%% Example 4 -- Part I: Small matrices
tolCheck = 1e-4;
load([DataDir,'DTI_small.mat']);
warning('off');

fileList=dir(MatricesDir);
nList=length(fileList);
count = 0;
%One data per row for matrix size, timing that includes initials, average iteration steps, success counts
ResultOfRQI      = zeros(nList-2,4);
ResultOfLeigopt      = zeros(nList-2,4); %eigopt
MatrixNames     = cell(nList-2,1);

for i=1:nList
    if strcmp(fileList(i).name,'.')==1||strcmp(fileList(i).name,'..')==1
        continue;
    else
        count = count + 1;
        
        % Load data and preprocessing
        matName = fileList(i).name;
        fprintf('count = %d, %s\n',count, matName);
        MatrixNames{count,1} = matName;
        load([MatricesDir,matName]);
        
        % Find shiftUsed and exact DTI
        flagFind = false;
        for jj = 1:size(DTIInfo,1)
            if strcmp(DTIInfo{jj,1},matName)
                flagFind = true;
                break;
            end
        end
        if ~flagFind
            fprintf('Matrix %s not found in DTIInfo\n',matName);
            error('Exceptional matrices!');
        end
        shiftUsed = DTIInfo{jj,3};
        DTIexact = DTIInfo{jj,2};
        opt.maxit = 100;
        n = size(A,1);
        A = A - shiftUsed*speye(n);
        
        tRQIUsed1    = 0;   tRQIUsed2       = 0;     nRQISuccess  = 0;      nRQIIter = 0;
        tLeigoptUsed = 0;   nLeigoptSuccess = 0;     nLeigoptIter = 0;
        issp = (nnz(A)/numel(A)<0.1);
        if ~issp
            A = full(A);
        else
            A = sparse(A);
        end
        
        tic;
        for j=1:testtimes
            if ~issp || n<500
                D = eig(full(A),'vector'); [~,ind] = max(real(D)); D = D(ind(1));
            else
                D = eigs(A,1,'lr');
            end
            opt.mu0 = imag(D);
        end
        trightmost = toc;

        % RQI
        fprintf('Testing RQI\n');
        opt.tol   = n*eps;
        for j=1:testtimes
            fprintf('%d-th test...%d tests altogether.\n',j,testtimes);
            tic;
            if issp
                [ABig, CBig] = matrix2eigen(A,'s');
                [u,opt.lambda0,v]  = svds(A-opt.mu0*1i*speye(n),1,0);
            else
                [ABig, CBig] = matrix2eigen(A,'d');
                [uset,lamSet,vset]     = svd(A-opt.mu0*1i*speye(n));
                opt.lambda0 = lamSet(end,end);
                u = uset(:,end);
                v = vset(:,end);
            end
            opt.x0    = [u;v]/sqrt(2);
            t1 = toc;
            tic;
            [xRQI, muRQI, lambdaRQI, infoRQI]=GRQI_D(ABig, CBig, opt);
            t2 = toc;
            if abs(abs(lambdaRQI)-DTIexact)<tolCheck*DTIexact
                tRQIUsed1 = tRQIUsed1 + t1;
                tRQIUsed2 = tRQIUsed2 + t2;
                nRQIIter  = nRQIIter + length(infoRQI.mu)-1;
                nRQISuccess = nRQISuccess + 1;
            end
        end
        ResultOfRQI(count,1) = n;
        ResultOfRQI(count,2) = (tRQIUsed1 + tRQIUsed2+trightmost)/nRQISuccess;
        ResultOfRQI(count,3) = nRQIIter/nRQISuccess;
        ResultOfRQI(count,4) = nRQISuccess/testtimes;

        % leigopt
        fprintf('Testing Leigopt\n');
        if ~issp
            svdSolver = 'svd';
        else
            svdSolver = 'svds';
        end
        if ~issp || n<500
            flagE = 'small';
        else
            flagE = 'false';
        end
        
        opt.tol   = 1e-12;
        if strcmp(matName,'supg400.mat')
            A = full(A);
        end
        for j=1:testtimes
            fprintf('%d-th test...%d tests altogether.\n',j,testtimes);
            tic;
            [beta,niter] = leigopt_dti(A,opt.tol,opt.mu0,opt.maxit,svdSolver);
            t1 = toc;
            if abs(beta-DTIexact)<tolCheck*DTIexact
                tLeigoptUsed = tLeigoptUsed + t1;
                nLeigoptIter  = nLeigoptIter + niter;
                nLeigoptSuccess = nLeigoptSuccess + 1;
            end
        end
        ResultOfLeigopt(count,1) = n;
        ResultOfLeigopt(count,2) = (tLeigoptUsed+trightmost)/nLeigoptSuccess;
        ResultOfLeigopt(count,3) = nLeigoptIter/nLeigoptSuccess;
        ResultOfLeigopt(count,4) = nLeigoptSuccess/testtimes;

 
    end
    
end
save resultForExample4_partI ResultOfRQI ResultOfLeigopt DTIInfo MatrixNames
DisplayTable('Table_For_Example4_PartI.txt', MatrixNames, ResultOfRQI, ResultOfLeigopt, 'RQI method','leigopt', false, false, DTIInfo);


%% Parameter setting for large size problems
testtimes = 10;     % Number of test for timing per example/method
MatricesDir = [DataDir,'testmatrices_large/'];

%% Example 4 -- Part II: Large matrices
load([DataDir,'DTI_large.mat']);
warning('off');
fileList=dir(MatricesDir);
nList=length(fileList);
count = 0;
%One data per row for matrix size, timing that excludes initials, average iteration steps, success counts
ResultOfRQI     = zeros(nList-2,4);
ResultOfLeigopt = zeros(nList-2,4);
MatrixNames     = cell(nList-2,1);

for i=1:nList
    if strcmp(fileList(i).name,'.')==1||strcmp(fileList(i).name,'..')==1
        continue;
    else
        count = count + 1;        
        % Load data and preprocessing
        matName = fileList(i).name;
        fprintf('count = %d, %s\n',count, matName);
        MatrixNames{count,1} = matName;
        load([MatricesDir,matName]);
        
        % Find shiftUsed
        flagFind = false;
        for jj = 1:size(DTIInfo,1)
            if strcmp(DTIInfo{jj,1},matName)
                flagFind = true;
                break;
            end
        end
        if ~flagFind
            fprintf('Matrix %s not found in DTIInfo\n',matName);
            error('Exceptional matrices!');
        end
        shiftUsed = DTIInfo{jj,2};  
        n = size(A,1);
        A = A - shiftUsed*speye(n);
              
        tRQIUsed1    = 0;      tRQIUsed2    = 0;   nRQIIter = 0;
        tLeigoptUsed = 0;      nLeigoptIter = 0;
        
        fprintf('Calculating the rightmost evl');
        tic;
        for j=1:testtimes
            if strcmp('tols4000.mat',MatrixNames{count})
                [V,D] = clay_Arnoldi(A, 1);
            else
                D = eigs(A, 1, 'lr');
            end
            opt.mu0 = imag(D);
        end
        trightmost = toc;
        
        % RQI        
        fprintf('Testing RQI\n');
        opt.tol   = n*eps;
        opt.maxit = 30;
        for j=1:testtimes
            fprintf('%d-th test...%d tests altogether.\n',j,testtimes);
            tic
            [ABig, CBig] = matrix2eigen(A,'s');
            [u,opt.lambda0,v] = svds(A-opt.mu0*1i*speye(n),1,0);
            opt.tol   = n*eps;
            opt.x0    = [u;v]/sqrt(2);
            t1 = toc;
            
            tic;
            [xRQI, muRQI, lambdaRQI, infoRQI]=GRQI_D(ABig, CBig, opt);
            t2 = toc;
            tRQIUsed1 = tRQIUsed1 + t1;
            tRQIUsed2 = tRQIUsed2 + t2;
            nRQIIter  = nRQIIter + length(infoRQI.mu)-1;
        end
        ResultOfRQI(count,1) = n;
        ResultOfRQI(count,2) = (trightmost + tRQIUsed1 + tRQIUsed2)/testtimes;
        ResultOfRQI(count,3) = nRQIIter/testtimes;
        ResultOfRQI(count,4) = lambdaRQI;
        
        % Leigopt 
        fprintf('Testing Leigopt\n');
        opt.maxit = 30;
        opt.tol   = 1e-12;
        for j=1:testtimes
            fprintf('%d-th test...%d tests altogether.\n',j,testtimes);
            tic;
            [beta,niter] = leigopt_dti(A,opt.tol,opt.mu0,opt.maxit);
            t1 = toc;
            tLeigoptUsed = tLeigoptUsed + t1;
            nLeigoptIter  = nLeigoptIter + niter;
        end
        ResultOfLeigopt(count,1) = n;
        ResultOfLeigopt(count,2) = (trightmost + tLeigoptUsed)/testtimes;
        ResultOfLeigopt(count,3) = nLeigoptIter/testtimes;
        ResultOfLeigopt(count,4) = beta;
    end
end

save resultForExample4_partII ResultOfRQI ResultOfLeigopt DTIInfo MatrixNames
DisplayTable('Table_For_Example4_PartII.txt', MatrixNames, ResultOfRQI, ResultOfLeigopt, 'RQI method','leigopt', true, false, DTIInfo);