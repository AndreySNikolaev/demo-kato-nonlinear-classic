#procedure HenrardNormalization(PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Normalization using Henrard-Koseleff algorithm 
*       Henrard, J., in "Recent advances in dynamical astronomy" ed. Szebehely, p. 248, 	
*       Algorithm by Koseleff PV (1994), Cel. Mech. 58:17–36
*  PREFIX - for uniqueness. Names of all objects created inside the procedure will begin with this PREFIX.
*==========================================================================

*--------------------------------------------------
*         Henrard-Koseleff algorithm:
*   U^-1_W H = \Sum_{n=0}^N alpha^n \tilde H_n
*    \tilde H_n = \Sum_{k=0}^n F_{k,n-k}
*       F_{n,k} == F{'n'}o{'k'}
*--------------------------------------------------
*         F_{0,k}  - coefficient of alpha^n in Op:
*--------------------------------------------------
L 'PREFIX'Op=Hetazeta;
B alpha;
.sort
#do n = 0,'MAXORDER'
G 'PREFIX'F0o'n'='PREFIX'Op[alpha^'n'];
#enddo
.store
*--------------------------------------------------
*         Main loop
*--------------------------------------------------
#do n = 1,{'MAXORDER'}
#do k = 1,{'n'-1}
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
*     Part 1 of the homological equation for W{n-1}:
*     \Sum_{k=1}^{n-1}  L_W{k-1} F_{n-k,0} 
*--------------------------------------------------
#do k=1,{'n'-1}
#call L('PREFIX'W{'k'-1},'PREFIX'F{'n'-'k'}o0)
G 'PREFIX'1F'n'o0m'k'= Pb'NDIM';
.store
#enddo

.store
G 'PREFIX'1F'n'o0=
#do k=1,{'n'-1}
+ 'PREFIX'1F'n'o0m'k'/'n'
#enddo
;
.store

*--------------------------------------------------
*     Part 2 of the homological equation for G{n-1}:
*     \Sum_{k=0}^{n-1} F_{k,n-k}
*--------------------------------------------------
G 'PREFIX'2P{'n'} = 
#do k = 0,{'n'-1}
+ 'PREFIX'F'k'o{'n'-'k'}
#enddo
;
.store
*--------------------
*    The following is the FORM style pseudo-solution of the homological equation: 
*        L_H0 W_n = - (1-P)*n*Ht_n
*    here Ht_n is the term of order of alpha^n in the transformed Hamiltonian 
*--------------------
Global 'PREFIX'W{'n'-1} = 'n'*('PREFIX'1F'n'o0-'PREFIX'2P{'n'});
#call SH0()
.store
#call L('PREFIX'W{'n'-1},'PREFIX'F0o0)
G 'PREFIX'F'n'o0=-Pb'NDIM'/'n'-'PREFIX'1F'n'o0;
.store
G 'PREFIX'Hto'n'= 'PREFIX'2P{'n'} + 'PREFIX'F'n'o0;
.store
#enddo
*--------------------
G 'PREFIX'Ht{'MAXORDER'} = H0 + 
#do n = 1,'MAXORDER'
+ 'PREFIX'Hto'n'*alpha^'n'
#enddo
;
B alpha;
*Print;
.store
G 'PREFIX'W'MAXORDER' =
#do n = 0,{'MAXORDER'-1}
+'PREFIX'W'n'*alpha^'n'
#enddo
;
B alpha;
.store
#endprocedure
