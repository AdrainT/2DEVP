function [M,Md] = funMIMOAB(z,A)

M  = (1-z)*A{1}+z*A{2};
Md = A{2} - A{1}; 

return