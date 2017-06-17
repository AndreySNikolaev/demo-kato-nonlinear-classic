#procedure RH()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*   Computation of the the perturbed resolvent
*
*-------------------- The unperturbed resolvent operator:
G RH0 = -P/z
#do n=0,{2*'MAXORDER'}
+ z^'n'*S^{'n'+1}
#enddo 
;
.store
*-------------------- The perturbation:
G LV=H-H[alpha^0];
#call toliouvillian()
.store
G F0= RH0;
.store
*-------------------- The Neumann series:
#do m=1,{'MAXORDER'}
L LF=LV*F{'m'-1};
.sort
G F'm'= -RH0*LF;
.store
#enddo
*-------------------- The perturbed resolvent:
write statistics;
G RH'MAXORDER'=
#do m=0,'MAXORDER'
+ F'm'
#enddo
;
Bracket z;
.store
#endprocedure
