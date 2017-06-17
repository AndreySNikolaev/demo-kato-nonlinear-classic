#procedure polysubst(var,subst,MAXJ)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		Efficient substitution of series var=subst into all active expressions 
*------------------------
.sort
PushHide;
.sort
L 'PREFIX'var1='subst';
B alpha;
.sort
#do j=1,'MAXJ'
*------------------------
*        Compute powers of 2 of 'subst' expression
*------------------------
#ifndef 'FASTPOWER'
L 'PREFIX'var{2^'j'}=('PREFIX'var{2^('j'-1)})^2;
B alpha;
.sort
#else
L 'PREFIX'Ltmp=0;
.sort
#do n=0,'MAXORDER'
#do k=0,'n'
Drop 'PREFIX'Ltmp;
L 'PREFIX'Ltmp1= 'PREFIX'Ltmp
+alpha^'n'*'PREFIX'var{2^('j'-1)}[alpha^'k']*'PREFIX'var{2^('j'-1)}[alpha^{'n'-'k'}];
.sort
Drop 'PREFIX'Ltmp1;
Local 'PREFIX'Ltmp='PREFIX'Ltmp1;
.sort
#enddo
#enddo
Drop 'PREFIX'Ltmp;
Local 'PREFIX'var{2^'j'}='PREFIX'Ltmp;
B alpha;
.sort
#endif
*------------------------
#enddo
PopHide;
.sort
#do j='MAXJ',0,-1
#do jj='j',0,-1
Skip 'PREFIX'var{2^'jj'};
#enddo
Drop 'PREFIX'var{2^'j'};
#ifndef 'FASTPOWER'
id ('var')^{2^'j'}='PREFIX'var{2^'j'};
#else
id alpha^n?*('var')^{2^'j'}=sum_(m,0,{'MAXORDER'}-n,alpha^(m+n)*('PREFIX'var{2^'j'}[alpha^m]));
#endif
.sort
#enddo
Drop 'PREFIX'var1;
#call eidentities()
*Print;
.sort
#endprocedure;
