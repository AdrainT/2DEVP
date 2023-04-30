function [f,g] = funMIMO(z,pars)
[M,Md] = feval(pars.fname,z,pars.A);

[V,D] = eig(M);
[eigvals,ind] = sort(diag(real(D)),'descend');

f = real(eigvals(end));
eigvec = V(:,ind(end));

% keyboard
g = real(eigvec'*Md*eigvec);
return