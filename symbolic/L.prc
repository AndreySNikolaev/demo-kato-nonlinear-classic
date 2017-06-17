#procedure L(G,F)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		calculate L_G F (or {F,G}) in (eta,zeta)  variables 
*------------------------
.sort
#do j=1,'NDIM'
#if 'j'>1
	Skip Pb{'j'-1};
#endif
	L dFdzeta'j' = 'F';
	L dGdzeta'j' = 'G';
	id zeta'j'^k? = k*zeta'j'^(k-1);
.sort
	Skip dFdzeta'j',dGdzeta'j'
#if 'j'>1
	Pb{'j'-1}
#endif
	;
	L dFdeta'j' = 'F';
	L dGdeta'j' = 'G';
	id eta'j'^k? = k*eta'j'^(k-1);
.sort
	Drop dFdzeta'j',dGdzeta'j',dFdeta'j',dGdeta'j'
#if 'j'>1
	Pb{'j'-1}
#endif
	;
	L Pb'j' = dFdzeta'j'*dGdeta'j'-dFdeta'j'*dGdzeta'j'
#if 'j' > 1
	+ Pb{'j'-1}
#endif
	;
#call eidentities()
.sort
#enddo
.sort
#endprocedure;
