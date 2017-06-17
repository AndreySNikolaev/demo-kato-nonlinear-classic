#procedure etazetatopq()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    The canonical transform from the (zeta,eta) to (p,q) variables
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
		id zeta'j' = ([C'j']*q'j' - i_*[D'j']*p'j')/2;
		id  eta'j' = ( - i_*[C'j']*q'j' + [D'j']*p'j')/2;
#enddo
	endrepeat;
.sort
repeat;
#do j=1,'NDIM'
		id [C'j'] = 2*Coef2'j';
		id [D'j'] = 2*Coef1'j';
#enddo
endrepeat;
#call eidentities()
.sort
#endprocedure
