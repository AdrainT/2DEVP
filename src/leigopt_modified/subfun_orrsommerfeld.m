function Ar = subfun_orrsommerfeld(A,P)


Ar{1} = A{2}\(A{1}*P);
Ar{2} = P;

Ar{3} = Ar{1}'*Ar{1};
Ar{4} = Ar{2}'*Ar{1};
Ar{5} = Ar{2}'*Ar{2};
	 
return;