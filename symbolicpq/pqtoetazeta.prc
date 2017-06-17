#procedure pqtoetazeta()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    The canonical transform from the (p,q) to the (zeta,eta) variables
*   q -> (Hq/Hp)^(- 1/4))*(zeta + i_*eta)/Sqrt[2], 
*   p -> (Hq/Hp)^(1/4)*(zeta - i_*eta)*i_/Sqrt[2]
*==========================================================================
.sort
Symbols 
#do j=1,'NDIM'
	[C'j'],[D'j']
#enddo
;
repeat;
#do j=1,'NDIM'
		id q'j' = [C'j']*(zeta'j'+i_*eta'j');
		id p'j' = i_*[D'j']*(zeta'j'-i_*eta'j');
#enddo
endrepeat;

repeat;
#do j=1,'NDIM'
		id [C'j'] = Coef1'j';
		id [D'j'] = Coef2'j';
#enddo
endrepeat;
#call eidentities()
.sort
#endprocedure
