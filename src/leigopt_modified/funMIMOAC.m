function [M,Md] = funMIMOAC(z,A)

M  = A{1}-z*A{2};
Md = -A{2}; 

return