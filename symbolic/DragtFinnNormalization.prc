#procedure DragtFinnNormalization(PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Normalization using Dragt-Finn chain of exponents:
*     \tilde H=exp(alpha^n L_G(n-1) )...exp(alpha L_G0) H
*  G=sum alpha^n Gn - pseudogenerator, Op - initial Hamiltonian.
*  PREFIX - for uniqueness. Names of all objects created inside this procedure 
*  will begin with PREFIX.
*       Algorithm by Koseleff PV (1994), Cel. Mech. 58:17–36
*==========================================================================
L 'PREFIX'Ht0=Hetazeta;
B alpha;
.sort
#do k=0,'MAXORDER'
*--------------------
*     F'n'k'k'm'm' = F^(n)_{k,m} =  1/m! L_Gn^m H_{k-m(n+1)}^{(n)}
*--------------------
G 'PREFIX'F0k'k'm0 = 'PREFIX'Ht0[alpha^'k'];
#enddo
.store
*--------------------
*   Main loop
*--------------------
#do n= 0,'MAXORDER'-1
*--------------------
*   Obtain the Generator for the next step. 
*   G_n block-diagonalizes the terms of alpha^(n+1)
*    The following is the FORM style pseudo-solution of the homological equation: 
*        L_H0 G_n =  (1-P)*H_n^(n)
*    here H_n^(n-1) is the term of order of alpha^n in the transformed Hamiltonian 
*--------------------
Global 'PREFIX'G'n' = 'PREFIX'F'n'k{'n'+1}m0;
*--------------------
#call SH0()
.store
#do k=0,'MAXORDER'
*--------------------
*   Sum_{m=1}^{[k/(n+1)] 1/m! L_G_n^m H_{k-m(n+1)}^(n)}
*--------------------
#do m=1,{'k'/('n'+1)}
#call L('PREFIX'G'n','PREFIX'F'n'k{'k'-('n'+1)}m{'m'-1})
G 'PREFIX'F'n'k'k'm'm'=Pb'NDIM'/'m';
.store
#enddo
G 'PREFIX'F{'n'+1}k'k'm0= 'PREFIX'F'n'k'k'm0
#do m=1,{'k'/('n'+1)}
+'PREFIX'F'n'k'k'm'm'
#enddo
;
.store
#enddo
#enddo
on statistics;
*--------------------
*	Pseudo-generator. Gn is the coefficient at alpha^n
*--------------------
G 'PREFIX'G{'MAXORDER'}= 
#do n = 0,'MAXORDER'-1
 + alpha^'n'*'PREFIX'G'n'
#enddo
;
G 'PREFIX'Ht{'MAXORDER'}= 
#do n = 0,'MAXORDER'
 + alpha^'n'*'PREFIX'F{'MAXORDER'}k'n'm0
#enddo
;
B alpha;
.store
#endprocedure
