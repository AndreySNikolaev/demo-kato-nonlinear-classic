#procedure Deprit3Normalizationpq(PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Normalization using the triangular Deprit algorithm in p,q variables
*       Deprit A (1969), Celestial Mech 1:12–30
*  PREFIX - for uniqueness. Names of all objects created inside the procedure will begin with this PREFIX.
*==========================================================================
L 'PREFIX'H = Hpq;
*--------------------
*     Warning: Bracketing alpha is mandatory here.
*--------------------
B alpha;
.sort
#do n = 0,'MAXORDER'
G 'PREFIX'H0o'n'='PREFIX'H[alpha^'n'];
#enddo
.store
*--------------------
*	Fill the Deprit triangle 
*--------------------

#redefine k "0"
*--------------------
*	main loop. 'nr' is the row number.
*--------------------
#do nr = 1,'MAXORDER'
*--------------------
*	First pass. 
*        These are Htilde{'n'}o{'k'} without the terms containing W{n-1}.  
*--------------------
L 'PREFIX'Ht0o'nr'='PREFIX'H0o'nr';
L 'PREFIX'W{'nr'-1}=0;
.sort
Hide 'PREFIX'Ht0o'nr','PREFIX'W{'nr'-1};
.sort
#do n = 1,'nr'
#redefine k "{'nr'-'n'}"
*--------------------
*	Since W{nr-1} have not computed yet, it is effectively treated as 0  
*         in the following Poisson brackets:
*--------------------
#do m = 0,'k'
#call Lpq('PREFIX'W'm','PREFIX'H{'n'-1}o{'k'-'m'} )
L 'PREFIX'Ht{'n'}o{'k'}A'm' = Pb'NDIM';
.sort
Drop Pb'NDIM';
Hide 'PREFIX'Ht{'n'}o{'k'}A'm';
#enddo
.sort
Drop
#do m = 0,'k'
 'PREFIX'Ht{'n'}o{'k'}A'm'
#enddo
;
*--------------------
*	F^[n]_{k} = (k+1) F^[n-1]_{k+1} + \sum_{m=0}^k L_{m} F^[n-1]_{k-m} 
*--------------------
L 'PREFIX'Ht{'n'}o{'k'}= {'k'+1}*'PREFIX'Ht{'n'-1}o{'k'+1}
#do m = 0,'k'
+ 'PREFIX'Ht{'n'}o{'k'}A'm'
#enddo
;
.sort
Hide 'PREFIX'Ht{'n'}o{'k'};
.sort
#enddo
*--------------------
*   Obtain the next order of W
*    This is FORM style pseudo-solution of homological equation 
*        L_H0 W_{n-1} =  1/(n-1)!*(1-P)*Ht^[n]_0
*--------------------
Local 'PREFIX'Wn='PREFIX'Ht{'nr'}o0/fac_({'nr'-1});
#call pqtoetazeta()
#call SH0()
#call etazetatopq()
.sort
Hide 'PREFIX'Wn;
*--------------------
*              Compute (j-1)! L_{G_{j-1}} H0
*--------------------
#call Lpq('PREFIX'Wn,'PREFIX'H0o0 )
Drop Pb'NDIM';
Drop 'PREFIX'W{'nr'-1};
L 'PREFIX'LW{'nr'-1}H0 = fac_({'nr'-1})*Pb'NDIM';
.sort
*--------------------
*	Second pass. Add previously omitted terms with W{n-1}
*--------------------
Global 'PREFIX'W{'nr'-1}='PREFIX'Wn;
#do n = 1,'nr'
#redefine k "{'nr'-'n'}"
G 'PREFIX'H{'n'}o{'k'}= 'PREFIX'Ht{'n'}o{'k'} + 'PREFIX'LW{'nr'-1}H0/fac_('k');
#enddo
.store
#enddo
*--------------------
*	Ht= \sum alpha^n/n! F{n}o0:
*--------------------
write statistics;
Global  'PREFIX'W =
#do n = 0,'MAXORDER'-1
+ alpha^{'n'}*'PREFIX'W{'n'}
#enddo
;
Global  'PREFIX'Ht{'MAXORDER'} =
#do n = 0,'MAXORDER'
+ alpha^{'n'}/fac_('n')*'PREFIX'H{'n'}o0
#enddo
;
B alpha;
.store
#endprocedure
