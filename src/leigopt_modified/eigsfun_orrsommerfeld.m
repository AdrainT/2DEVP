function y = eigsfun_orrsommerfeld(u,flag,A,z)
if strcmp(flag,'transp')
    y = A{2}'*(((A{1}-1i*z*A{2})')\u);
else
    y = (A{1}-1i*z*A{2})\(A{2}*u);
end
end
