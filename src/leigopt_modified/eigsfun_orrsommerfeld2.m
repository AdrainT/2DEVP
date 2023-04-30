function y = eigsfun_orrsommerfeld2(u,ALU,C)
y = ALU.P'* ((ALU.L')\((ALU.U')\u));
y = C{2}*(C{2}*y);
y = ALU.U\(ALU.L\(ALU.P*y));
end
