#procedure PH0()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    Compute the P_H0 - the unperturbed averaging operator 
*==========================================================================
.sort
cFunctions xcount;
Multiply xcount(0);
#do j=1,'NDIM'
        id zeta'j'^k?*eta'j'^m?*xcount(n?)=xcount(n+Omega'j'*k-Omega'j'*m)*zeta'j'^k*eta'j'^m;
#enddo
id xcount(0)=1;
id xcount(k?)=0;
.sort
#endprocedure
