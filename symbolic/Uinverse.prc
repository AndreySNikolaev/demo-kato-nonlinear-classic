#procedure Uinverse(W,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*	Calculation of the inverse Deprit transform U_'W'^-1
*  U^-1_n = - 1/n Sum_{m=0}^{n-1} L_{W_{n-m-1}} U^-1_m 
*  using Koseleff-Henrard algorithm  	
*	Koseleff, P.-V., Celestial Mechanics & Dynamical Astronomy, Volume 58, Issue 1, pp.17
*   W - generator
*   PREFIX - for uniqueness. 
*   Names of all objects created inside this procedure 
*   will begin with PREFIX.
*==========================================================================
G 'PREFIX'W ='W';
G 'PREFIX'F0 ='Op';
B alpha;
.sort
#do n = 0,'MAXORDER'
G 'PREFIX'W'n'='PREFIX'W[alpha^'n'];
#enddo
.store
#do n = 1,'MAXORDER'
*--------------------
*			Calculate F_n  = -Sum_{m=0}^{n-1} L_{W_{n-m-1}} F_m:
*--------------------
#do m = 0,'n'-1
*
#call L('PREFIX'W{'n'-'m'-1},  'PREFIX'F'm' )
.sort
Global 'PREFIX'{'n'}F'm' = Pb'NDIM'  
#if 'm'>0
+ 'PREFIX'{'n'}F{'m'-1}
#endif
;
.store
#enddo
G 'PREFIX'F{'n'} = - 'PREFIX'{'n'}F{'n'-1}/'n';
.store
#enddo
On statistics;
G 'PREFIX'Uinv'MAXORDER' =
#do n = 0,'MAXORDER'
+ alpha^'n'*'PREFIX'F{'n'}
#enddo
;
B alpha;
*Print;
.store
#endprocedure
