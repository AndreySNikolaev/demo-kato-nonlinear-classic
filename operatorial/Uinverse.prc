#procedure Uinverse(W,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*			Operatorial 
*	Calculation of the inverse Deprit transform U_'W'^-1
*  U^-1_n = - 1/n Sum_{m=0}^{n-1} L_{W_{n-m-1}} U^-1_m 
*  G - generator
*  PREFIX - for uniqueness. Names of all objects created inside this procedure 
*  will begin with PREFIX.
*==========================================================================
*	Calculate L_G. 
*     Warning: Bracketing alpha is mandatory here.
*--------------------
Global 'PREFIX'LW = 'W';
#call toliouvillian()
B alpha;
*Print;
.store
Global 'PREFIX'0Ui = 1;
.store
#do n = 1,'MAXORDER'
*--------------------
*			Calculate Ui_n = Sum_{m=0}^{n-1} L_{W_{n-m-1}} Ui_m:
*--------------------
#do m = 0,'n'-1
*
Global 'PREFIX'{'n'}Ui'm' = 'PREFIX'LW[alpha^{'n'-'m'-1}]*'PREFIX'{'m'}Ui
#if 'm'>0
+ 'PREFIX'{'n'}Ui{'m'-1}
#endif
;
.store 
#enddo
G 'PREFIX'{'n'}Ui = - 'PREFIX'{'n'}Ui{'n'-1}/'n';
.store
#enddo
write statistics;
G 'PREFIX'Ui =
#do n = 0,'MAXORDER'
+ alpha^'n'*'PREFIX'{'n'}Ui
#enddo
;
B alpha;
*Print;
.store
#endprocedure
