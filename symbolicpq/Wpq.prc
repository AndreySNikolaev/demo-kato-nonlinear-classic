#procedure Wpq(PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*--------------------------------------------------------------------------
*				Initialization: F_0= z R_H0(z) dH/dalpha
*    upto MAXORDER-1  in P,Q variables
*--------------------------------------------------------------------------
Symbol 'PREFIX'alpha(:{'MAXORDER'-1});
.global
Global 'PREFIX'V=Hpq-Hpq[alpha^0];
id alpha='PREFIX'alpha;
Bracket 'PREFIX'alpha;
.store
Local 'PREFIX'dHdalpha=Hpq;
id alpha^n? = n*'PREFIX'alpha^(n-1);
Bracket 'PREFIX'alpha;
.sort
Drop 'PREFIX'dHdalpha;
Global 'PREFIX'F1= 'PREFIX'dHdalpha;
#call pqtoetazeta()
#call zRH0(0)
#call etazetatopq()
Bracket z;
.store
*--------------------------------------------------------------------------
*				Main loop.  F_m =  z R_H0(z) L_V F_{m-1}
*--------------------------------------------------------------------------
#do m=2,{'MAXORDER'}
Local 'PREFIX'Ltmp=0;
Local 'PREFIX'Fa='PREFIX'F{'m'-1};
Local 'PREFIX'Va='PREFIX'V;
Bracket 'PREFIX'alpha;
.sort
Hide 'PREFIX'Fa,'PREFIX'Va;
.sort
#do n=0,{'MAXORDER'-1}
#do k=0,'n'
#call Lpq('PREFIX'Va['PREFIX'alpha^'k'],'PREFIX'Fa['PREFIX'alpha^{'n'-'k'}])
Drop Pb'NDIM';
Drop 'PREFIX'Ltmp;
Local 'PREFIX'Ltmp1= 'PREFIX'Ltmp-'PREFIX'alpha^'n'*Pb'NDIM';
.sort
Drop 'PREFIX'Ltmp1;
Local 'PREFIX'Ltmp='PREFIX'Ltmp1;
.sort
Hide 'PREFIX'Ltmp;
.sort
#enddo
#enddo
Global 'PREFIX'F'm'='PREFIX'Ltmp;
#call pqtoetazeta()
#call zRH0({'m'-1})
#call etazetatopq()
Bracket z;
.store
#enddo
write statistics;
*    Compute W
Global 'PREFIX'SH=
#do n=1,'MAXORDER'
+ 'PREFIX'F'n'[z^{'n'}]
#enddo
;
id 'PREFIX'alpha=alpha;
Bracket alpha;
*Print;
.store
#endprocedure
