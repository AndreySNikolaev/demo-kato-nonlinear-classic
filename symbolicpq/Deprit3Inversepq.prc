#procedure Deprit3Inversepq(W,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Calculation of the inverse Deprit transform U^-1_W Op
*   	Inverse Triangular Deprit algorithm in p,q variables
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
G 'PREFIX'F'n'o0=fac_('n')*'PREFIX'Op[alpha^'n'];
*Print;
.store
#enddo

*--------------------
*	F^{[n]}_{k} = \frac{1}{k} (F^{[n+1]}_{k-1} - \sum_{m=1}^{k} \Lo_{G_{m-1}} F^{[n]}_{k-m})
*--------------------
#do k = 1,{'MAXORDER'}
#do n = 0,{'MAXORDER'-'k'}
L 'PREFIX'F'n'o'k'm0=0;
.sort
Hide 'PREFIX'F'n'o'k'm0;
#do m = 1,'k'
#call Lpq('PREFIX'W[alpha^{'m'-1}],'PREFIX'F'n'o{'k'-'m'})
Drop 'PREFIX'F'n'o'k'm{'m'-1},Pb'NDIM';
Local 'PREFIX'F'n'o'k'm'm' = 'PREFIX'F'n'o'k'm{'m'-1} + Pb'NDIM';
.sort
Hide 'PREFIX'F'n'o'k'm'm';
#enddo
G 'PREFIX'F'n'o'k'= 1/'k'*('PREFIX'F{'n'+1}o{'k'-1}-'PREFIX'F'n'o'k'm'k');
.store
#enddo
#enddo
on statistics;
G 'PREFIX'Uinv'MAXORDER'= 
#do k = 0,'MAXORDER'
 + alpha^'k'*'PREFIX'F0o'k'
#enddo
;
B alpha;
.store
#endprocedure
