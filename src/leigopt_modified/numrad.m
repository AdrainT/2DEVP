function [M,Md] = numrad(z,A)

M  = (A{1}*exp(i*z) + A{1}'*exp(-i*z))/2;
Md = (i*A{1}*exp(i*z) - i*A{1}'*exp(-i*z))/2; 
	 
	 
return;