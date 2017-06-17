#procedure Lpq(G,F)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		calculate L_G F (or {F,G}) in the p,q variables 
*------------------------
.sort
#do j=1,'NDIM'
#if 'j'>1
	Skip Pb{'j'-1};
#endif
	L dFdq'j' = 'F';
	L dGdq'j' = 'G';
	id q'j'^k? = k*q'j'^(k-1);
.sort
	Skip dFdq'j',dGdq'j'
#if 'j'>1
	Pb{'j'-1}
#endif
	;
	L dFdp'j' = 'F';
	L dGdp'j' = 'G';
	id p'j'^k? = k*p'j'^(k-1);
.sort
	Drop dFdq'j',dGdq'j',dFdp'j',dGdp'j'
#if 'j'>1
	Pb{'j'-1}
#endif
	;
	L Pb'j' = dFdq'j'*dGdp'j'-dFdp'j'*dGdq'j'
#if 'j' > 1
	+ Pb{'j'-1}
#endif
	;
#call eidentities()
.sort
#enddo
.sort
#endprocedure;
