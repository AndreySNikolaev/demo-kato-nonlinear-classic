#procedure DragtFinnInverse(G,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Calculation of Dragt-Finn inverse transform as a product of exponents
*    exp(-alpha L_G0) ...  exp(-alpha^n L_G(n-1) ) Op
*  G=sum alpha^n Gn - pseudogenerator , Op - operand.
*  PREFIX - for uniqueness. Names of all objects created inside this procedure 
*  will begin with PREFIX.
*       Algorithm by Koseleff PV (1994), Cel. Mech. 58:17–36
*==========================================================================
*--------------------
*	Pseudo-generator. Gn is the coefficient at alpha^n
*--------------------
L 'PREFIX'G = 'G';
L 'PREFIX'It'MAXORDER' = 'Op';
B alpha;
.sort
#do k=0,'MAXORDER'
*--------------------
*     F'n'k'k'm'm' = F^(n)_{k,m} =  (-1)^m/m! L_Gn^m I_{k-m(n+1)}^{(n)}
*--------------------
G 'PREFIX'F'MAXORDER'k'k'm0 = 'PREFIX'It'MAXORDER'[alpha^'k'];
G 'PREFIX'G'k' = -'PREFIX'G[alpha^'k'];
#enddo
.store
*--------------------
*   Main loop
*--------------------
#do n= 'MAXORDER'-1,0,-1
#do k=0,'MAXORDER'
*--------------------
*   Sum_{m=1}^{[k/(n+1)] 1/m! L_G_n^m I_{k-m(n+1)}^(n)}
*--------------------
#do m=1,{'k'/('n'+1)}
#call L('PREFIX'G'n','PREFIX'F{'n'+1}k{'k'-('n'+1)}m{'m'-1})
G 'PREFIX'F{'n'+1}k'k'm'm'= Pb'NDIM'/'m';
.store
#enddo
G 'PREFIX'F'n'k'k'm0 = 
+ 'PREFIX'F{'n'+1}k'k'm0
#do m=1,{'k'/('n'+1)}
+'PREFIX'F{'n'+1}k'k'm'm'
#enddo
;
.store
#enddo
#enddo
on statistics;
G 'PREFIX'DFinv'MAXORDER'= 
#do n = 0,'MAXORDER'
 + alpha^'n'*'PREFIX'F0k'n'm0
#enddo
;
B alpha;
.store
#endprocedure
