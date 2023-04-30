function [big_A,big_C]=matrix2eigen(A,type)
% Convert DTI problem of matrix A to 2DEVP formulation.
n=size(A,1);
if type == 's'
    big_A = [sparse(n,n),A;A',sparse(n,n)];
    big_C = [sparse(n,n),1i*speye(n);-1i*speye(n),sparse(n,n)];
else
    big_A=[zeros(n),A;A',zeros(n)];
    big_C=[zeros(n),1i*eye(n);-1i*eye(n),zeros(n)];
end