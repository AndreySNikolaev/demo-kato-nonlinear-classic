#procedure U(W,PREFIX)
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
*--------------------
*     1. Calculate L_W. 
*     Warning: Bracketing alpha is mandatory here.
*--------------------
Global 'PREFIX'LW = 'W';
#call toliouvillian()
B alpha;
*Print;
.store
#do n = 0,'MAXORDER'
G 'PREFIX'F'MAXORDER'o'n' = alpha^'n';
#enddo
.store
*--------------------
*    2. Main loop 
*--------------------
#do n = 'MAXORDER'-1,0,-1
*--------------------
*   F^(n)_k = F^(n+1)_k + 1/(n+1) L_W(n-k) F^(n+1)_n
*--------------------
#do k = 0,'n'
Global 'PREFIX'F'n'o{'k'} = 'PREFIX'F{'n'+1}o{'k'} + 1/{'n'+1}*'PREFIX'LW[alpha^{'n'-'k'}]*'PREFIX'F{'n'+1}o{'n'+1}; 
.store
#enddo
#enddo
On statistics;
Global 'PREFIX'U = 'PREFIX'F0o0;
B alpha;
*Print;
.store
#endprocedure
