function [M,Md,Msq] = distinstab_general(z,A)


% keyboard
M  = A{1} - i*z*A{2};
Md = -i*A{2}; 
Msq = A{3} + i*z*A{4} - i*z*A{4}' + z^2*A{5};
	 
return;