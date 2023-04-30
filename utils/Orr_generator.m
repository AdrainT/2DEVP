function [Ln,Bn] = Orr_generator(n,type)

if nargin<=1
    type = 1;
end

if type == 1    
    h = 2/(n+1); R = 1000;
    x = -1+h:h:1-h;
    Ln = sparse(1/h^2*( diag(-(2+h^2)*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) ));
    Un = sparse(diag(1-x.^2));
    Bn = 1/R*Ln^2-1i*(Un*Ln+2*speye(n));
elseif type == 2
    h = 2/(n+1); R = 1000;
    Ln = sparse(1/h^2*( diag(-(2+h^2)*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) ));
    x = -1 + (1:n)*h;
    Un = sparse(diag(1-x.^2));
    Bn = 1/R*Ln*Ln-1i*(Un*Ln-2*speye(n));
elseif type == 3
    h = 2/(n+1); R = 100;
    Ln = sparse(1/h^2*( diag(-(2+h^2)*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) ));
    x = 1 + (1:n)*h;
    Un = sparse(diag(1-x.^2));
    Bn = 1/R*Ln*Ln-1i*(Un*Ln+2*speye(n));
end

    
    

