function y = svdsfun_orrsommerfeld(u,flag,ALU,C)
if strcmp(flag,'transp')
    y = ALU.P'* ((ALU.L')\((ALU.U')\u));
    y = C{2}'*y;
else
    y = C{2}*u;
    y = ALU.U\(ALU.L\(ALU.P*y));
end
end
