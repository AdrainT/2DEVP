function y=svdfunc(funBLv,funBLTv,x,flag)
if strcmp(flag,'notransp')
    y = funBLv(x);
else
    y = funBLTv(x);
end