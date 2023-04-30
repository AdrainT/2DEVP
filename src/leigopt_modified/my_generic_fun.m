function [f,g] = my_generic_fun(z,pars)

[M,Md] = feval(pars.fname,z,pars.A);
[n,~] = size(M);

[V,D] = eig(M);
[eigvals,ind] = sort(diag(real(D)),'descend');

f = real(eigvals(pars.j));
eigvec = V(:,ind(pars.j));

% keyboard
for k = 1:pars.d
	g(k) = real(eigvec'*Md(:,:,k)*eigvec);
end


return;