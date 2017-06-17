#procedure HoriMersmanInverse(G,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Compute inverse Hori transform exp(alpha*L_G) using Mersmann algorithm
*       Algorithm by Mersman WA (1970), Celest. Mech. 3(1):81–89
*  G=sum alpha^n Gn - Hori generator,  Op - operand.
*  PREFIX - for uniqueness. Op - initial Hamiltonian.
*  Names of all objects created inside the procedure will begin with this PREFIX.
*==========================================================================
*--------------------------------------------------
*         alpha Bracketing
*--------------------------------------------------
Global  'PREFIX'G='G';
G 'PREFIX'Op='Op';
B alpha;
*Print;
.store
*--------------------------------------------------
*         Mersman algorithm:
*   exp(alpha*L_G) F = \Sum_{n=0}^N alpha^n Ft_n
*    Ft_n = \Sum_{k=0}^n 1/k! F^(k)_{n-k}
*    F^(n)_k =  \Sum_m=0^k L_Gm F^(n-1)_{k-m}
*--------------------------------------------------
*         F^(0)_k - coefficient of alpha^k in Op:
*--------------------------------------------------
#do k = 0,'MAXORDER'
G 'PREFIX'F0o'k'='PREFIX'Op[alpha^'k'];
*Print;
.store
#enddo
*--------------------------------------------------
*         for each j compute all F^(n)_k such that n+k=j 
*--------------------------------------------------
#do j = 1,{'MAXORDER'}
#do n = 1,'j'
*--------------------------------------------------
*    F^(n)_k = \Sum_m=0^k L_Gm F^(n-1)_{k-m}
*--------------------------------------------------
#do m= 0,{'j'-'n'}
#call L('PREFIX'G[alpha^'m'],'PREFIX'F{'n'-1}o{'j'-'n'-'m'})
.sort
G 'PREFIX'F'n'o{'j'-'n'}m{'m'}= Pb'NDIM';
.store
#enddo
G 'PREFIX'F'n'o{'j'-'n'}=
#do m= 0,{'j'-'n'}
+'PREFIX'F'n'o{'j'-'n'}m{'m'}
#enddo
;
.store
#enddo

G 'PREFIX'Ft{'j'} = 
#do k = 0,{'j'}
+ 1/fac_('k')*'PREFIX'F'k'o{'j'-'k'}
#enddo
;
.store

#enddo
on statistics;
G 'PREFIX'HMinv{'MAXORDER'} = 'PREFIX'Op[alpha^0]  
#do n = 1,'MAXORDER'
+ 'PREFIX'Ft'n'*alpha^'n'
#enddo
;
B alpha;
*Print;
.store
#endprocedure
