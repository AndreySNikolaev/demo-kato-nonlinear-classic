#procedure Deprit3transform(W,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Calculation of the direct Deprit transform U_W Op
*   	Triangular Deprit algorithm
*       Deprit A (1969), Celestial Mech 1:12–30
*  W - generator, Op - operand.
*  PREFIX - for uniqueness. Names of all objects inside the procedure will begin with this PREFIX.
*==========================================================================
Global 'PREFIX'W = 'W';
Global 'PREFIX'Op = 'Op';
*--------------------
*     Warning: Bracketing alpha is mandatory here.
*--------------------
B alpha;
.store
#do n = 0,'MAXORDER'
G 'PREFIX'F0o'n'='PREFIX'Op[alpha^'n'];
*Print;
.store
#enddo

*--------------------
*	F^[n+1]_{k-1} = k F^[n]_k + \sum_{m=1}^k L_{m-1} F^[n]_{k-m} 
*--------------------
#do n = 0,{'MAXORDER'-1}
#do k = 1,{'MAXORDER'-'n'}
L 'PREFIX'F'n'o'k'm0=0;
.sort
Hide 'PREFIX'F'n'o'k'm0;
#do m = 1,'k'
#call L('PREFIX'W[alpha^{'m'-1}],'PREFIX'F'n'o{'k'-'m'})
Drop 'PREFIX'F'n'o'k'm{'m'-1},Pb'NDIM';
Local 'PREFIX'F'n'o'k'm'm' = 'PREFIX'F'n'o'k'm{'m'-1} + Pb'NDIM';
.sort
Hide 'PREFIX'F'n'o'k'm'm';
#enddo
G 'PREFIX'F{'n'+1}o{'k'-1}= 'k'*'PREFIX'F'n'o'k'+'PREFIX'F'n'o'k'm'k';
.store
#enddo
#enddo
G 'PREFIX'U'MAXORDER' =
#do n = 0,'MAXORDER'
+ alpha^{'n'}/fac_('n')*'PREFIX'F'n'o0
#enddo
;
B alpha;
.store
#endprocedure
