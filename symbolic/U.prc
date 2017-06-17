#procedure U(W,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*			Calculation of the direct Deprit transform U_'W' 
*   U_n =  1/n Sum_{m=0}^{n-1} U_m L_{W_{n-m-1}}
*   W - generator
*   PREFIX - for uniqueness. 
*   Names of all objects created inside this procedure 
*   will begin with PREFIX.
*==========================================================================
L 'PREFIX'W='W';
L 'PREFIX'Op='Op';
B alpha;
.sort
#do n = 0,'MAXORDER'
G 'PREFIX'W'n'='PREFIX'W[alpha^'n'];
G 'PREFIX'Op'n'='PREFIX'Op[alpha^'n'];
#enddo
.store
*--------------------
*	The seed expressions for the ordered exponential: Sum_{n=0}^N \hat V(n) F^(N)_n
*--------------------
#do n = 0,'MAXORDER'
G 'PREFIX'F'MAXORDER'o'n' = 
#do k= 'n','MAXORDER'
+ alpha^'k'*'PREFIX'Op{'k'-'n'}
#enddo
;
#enddo
.store
*--------------------
*    Main loop 
*--------------------
#do n = 'MAXORDER'-1,0,-1
*--------------------
*   F^(n)_k = F^(n+1)_k + 1/(n+1) L_W(n-k) F^(n+1)_n
*--------------------
#do k = 0,'n'
*                   L_Wm Op  
#call L('PREFIX'W{'n'-'k'},  'PREFIX'F{'n'+1}o{'n'+1} )
.sort
Global 'PREFIX'F'n'o{'k'} = 'PREFIX'F{'n'+1}o{'k'} + Pb'NDIM' /{'n'+1}; 
.store
#enddo
#enddo
On statistics;
Global 'PREFIX'U'MAXORDER' = 'PREFIX'F0o0;
B alpha;
*Print;
.store
#endprocedure
