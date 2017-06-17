#procedure HenrardInverse(W,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Compute the direct Depri transform U_W using the Henrard algorithm
*       Henrard, J., in "Recent advances in dynamical astronomy" ed. Szebehely, p. 248, 	
*  W - generator, Op - operand.
*  PREFIX - for uniqueness. Op - initial Hamiltonian.
*  Names of all objects created inside the procedure will begin with this PREFIX.
*==========================================================================
*         Henrard algorithm for inverse ordered exponent:
*    Instead of computation \tilde F= U_W F
*   algorithm computes \tilde F such that
*   U_W^-1 \tilde F = F
*   
*  \Sum_{n=0}^N alpha^n \tilde F_n
*       F_{n,k} == F{'n'}o{'k'}
*--------------------------------------------------
*         F_{0,k}  - coefficient of alpha^n in Op:
*--------------------------------------------------
L 'PREFIX'Op='Op';
L 'PREFIX'W='W';
B alpha;
.sort
#do n = 0,'MAXORDER'
G 'PREFIX'F'n'='PREFIX'Op[alpha^'n'];
G 'PREFIX'W'n'='PREFIX'W[alpha^'n'];
#enddo
.store
G 'PREFIX'F0o0='PREFIX'F0;
.store 
*--------------------------------------------------
*         Main loop
*--------------------------------------------------
#do n = 1,'MAXORDER'
#do k = 'n',1,-1
*--------------------------------------------------
*    F_{k,n-k} = -1/k \Sum_j=1^k L_W{j-1} F_{k-j,n-k}
*--------------------------------------------------
#do j=1,'k'
#call L('PREFIX'W{'j'-1},'PREFIX'F{'k'-'j'}o{'n'-'k'})
G 'PREFIX'F'k'o{'n'-'k'}m'j'= -Pb'NDIM';
.store
#enddo

G 'PREFIX'F'k'o{'n'-'k'}=
#do j=1,'k'
+ 'PREFIX'F'k'o{'n'-'k'}m'j'/'k'
#enddo
;
.store
#enddo

*--------------------------------------------------
*    \tilde F_n = F_{0,n}= F_n - \Sum_{k=1}^n F_{k,n-k}
*--------------------------------------------------
G 'PREFIX'F0o{'n'} = 'PREFIX'F'n'
#do k = 1,{'n'}
- 'PREFIX'F'k'o{'n'-'k'}
#enddo
;
.store

#enddo
*--------------------------------------------------
on statistics;
G 'PREFIX'Hinv'MAXORDER' = 'PREFIX'F0o0
#do n = 1,'MAXORDER'
+ alpha^{'n'}*'PREFIX'F0o{'n'}
#enddo
;
B alpha;
*Print;
.store
#endprocedure
