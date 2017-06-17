#procedure HoriMersmanNormalizationpq(PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Normalization using Hori-Mersmann approach.
*       Algorithm by Mersman WA (1970), Celest. Mech. 3(1):81–89
*  PREFIX - for uniqueness. Op - initial Hamiltonian.
*  Names of all objects created inside the procedure will begin with this PREFIX.
*==========================================================================
CFunctions 'PREFIX'G;
nwrite statistics;
.global
*--------------------------------------------------
*         Generator
*--------------------------------------------------
Global  'PREFIX'G0= 
#do n = 0,'MAXORDER'-1
 +'PREFIX'G('n')*alpha^'n'
#enddo
;
G 'PREFIX'Op=Hpq;
B alpha;
.store

*--------------------------------------------------
*         Mersman algorithm:
*   exp(-alpha*L_G) H = \Sum_{n=0}^N alpha^n Ht_n
*    Ht_n = \Sum_{k=0}^n 1/k! F^(k)_{n-k}
*    F^(n)_k = - \Sum_m=0^k L_Gm F^(n-1)_{k-m}
*--------------------------------------------------
*         F^(0)_k - coefficient of alpha^n in Op:
*--------------------------------------------------
#do n = 0,'MAXORDER'
G 'PREFIX'F0o'n'='PREFIX'Op[alpha^'n'];
*Print;
.store
#enddo
*--------------------------------------------------
*         for each j compute all F^(n)_k such that n+k=j 
*--------------------------------------------------
#do j = 1,{'MAXORDER'}

#do n = 1,'j'
*--------------------------------------------------
G 'PREFIX'F'n'o{'j'-'n'}m0=0;
*--------------------------------------------------
.store
*    F^(n)_k = - \Sum_m=0^k L_Gm F^(n-1)_{k-m}
#do m= 0,{'j'-'n'}
#call Lpq('PREFIX'G{'j'-1}[alpha^'m'],'PREFIX'F{'n'-1}o{'j'-'n'-'m'})
.sort
*--------------------------------------------------
*     Workaround for FORM bug:
*--------------------------------------------------
L WA='PREFIX'F'n'o{'j'-'n'}m{'m'};
.sort
Drop WA;
G 'PREFIX'F'n'o{'j'-'n'}m{'m'+1}= WA  - Pb'NDIM';
*Print;
.store
#enddo
*--------------------------------------------------
*         F^(1)_ contain not known at this point L_G_n
*--------------------------------------------------
#if 'n' == 1 
G 'PREFIX'U'n'o{'j'-'n'}=
#else
G 'PREFIX'F'n'o{'j'-'n'}=
#endif
'PREFIX'F'n'o{'j'-'n'}m{'j'-'n'+1};
.store
#enddo

*--------------------
*   Compute the next order of Generator (G_n normalizes terms of alpha^{n+1})
*    This is FORM style pseudo-solution of homological equation 
*        L_H0 G_n = - (1-P)*Ht_(n+1)
*    here Ht_(n+1) is term of order of alpha^(n+1) in transformed Hamiltonian
*--------------------
Local GnO = 
#do k = 0,{'j'}
#if 'k' == 1 
-             'PREFIX'U'k'o{'j'-'k'}
#else
- 1/fac_('k')*'PREFIX'F'k'o{'j'-'k'}
#endif
#enddo
;
#call pqtoetazeta()
#call SH0()
#call etazetatopq()
.sort
*--------------------
*              compute L_{G_n} H0 
*--------------------
Global 'PREFIX'G{'j'}= 'PREFIX'G{'j'-1};
.sort
id 'PREFIX'G({'j'-1})=GnO;
B alpha;
*Print;
.store
#call Lpq('PREFIX'G{'j'}[alpha^{'j'-1}],H0)
.sort
*--------------------
*    substitute just computed G_n into partially normalized Hamiltonian
*--------------------
*     Workaround for FORM bug:
*--------------------------------------------------
L WA='PREFIX'U1o{'j'-1};
.sort
Global 'PREFIX'F1o{'j'-1}= WA - Pb'NDIM';
.store

G 'PREFIX'Hto{'j'} = 
#do k = 0,{'j'}
+ 1/fac_('k')*'PREFIX'F'k'o{'j'-'k'}
#enddo
;
.store

#enddo

G 'PREFIX'Ht{'MAXORDER'} = H0 + 
#do n = 1,'MAXORDER'
+ 'PREFIX'Hto'n'*alpha^'n'
#enddo
;
B alpha;
*Print;
.store
#endprocedure
