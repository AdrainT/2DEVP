function y = funABs(funAs,funBs,x,mu,FLAG,n)
if strcmp(FLAG,'dim')
    y= n;
elseif strcmp(FLAG,'real')
    y= 0;
else
    y= (1-mu)*funAs(x) + mu*funBs(x);
end
    
    