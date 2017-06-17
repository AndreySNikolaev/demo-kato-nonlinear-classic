#procedure zRH0(am)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*      Calculate z R_H0  in the (zeta,eta) variables. 
*          P is the unperturbed averaging operator.
*              P F=Sum_((omega,m-n)=0) F_mn zeta^m eta^n  
*          S is the unperturbed integrating operator.
*              S F=Sum_((omega,m-n)!=0) 1/(i(omega,m-n))F_mn zeta^m eta^n  
*==========================================================================
.sort
cFunctions x;
Multiply x(0);
#do j=1,'NDIM'
        id zeta'j'^k?*eta'j'^m?*x(n?)=x(n+Omega'j'*k-Omega'j'*m)*zeta'j'^k*eta'j'^m;
#enddo
.sort
id x(0)=-1;
#if 'am'<{'MAXORDER'-1}
id z^m?*x(k?) =  sum_(n,m+1,{'MAXORDER'},z^n*(-i_/k)^(n-m));
#else
*==========================================================================
*      Computational optimization.
*==========================================================================
id z^'MAXORDER'*x(k?) =0;
id z^m?*x(k?) =   z^'MAXORDER'*(-i_/k)^('MAXORDER'-m);
#endif
.sort
#endprocedure
