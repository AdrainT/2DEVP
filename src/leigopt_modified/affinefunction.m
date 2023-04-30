function [M,Md] = affinefunction(z,A)

M = A{1};

% keyboard

for k = 1:length(z)
	M = M + z(k)*A{k+1};
	Md(:,:,k) = full(A{k+1});
end


return;