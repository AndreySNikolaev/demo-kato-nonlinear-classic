#procedure GustavsonInverse(W,Op,PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Calculation of the Gustavson inverse transfrom as a product of exponents
*  W=sum alpha^n Wn - pseudo generating function , Op - operand.
*  PREFIX - for uniqueness. Names of all objects created inside this procedure 
*  will begin with PREFIX.
*    Algorithm by Koseleff, P.-V., Cel. Mech. 58, 1994
*==========================================================================
S
#do m=1,'NDIM'
  etat'm',zetat'm'
#enddo
;
.global
*--------------------
*	Pseudo-generator. Gn is the coefficient at alpha^n
*--------------------
L 'PREFIX'W = 'W';
B alpha;
.sort
#do n=1,'MAXORDER'
G 'PREFIX'W'n' = 'PREFIX'W[alpha^'n'];
#enddo
#do m=1,'NDIM'
id  eta'm'=etat'm';
#enddo
.store
*--------------------
*   Main loop
*--------------------
#do n= {'MAXORDER'},1,-1
*---------------------------
*  The derivatives of generating function
*---------------------------
#do m=1,'NDIM'
L 'PREFIX'W'n'eta'm'='PREFIX'W'n';
id etat'm'^k? = k*etat'm'^(k-1);
.sort
Hide 'PREFIX'W'n'eta'm';
.sort
L 'PREFIX'W'n'zeta'm'='PREFIX'W'n';
id zeta'm'^k? = k*zeta'm'^(k-1);
.sort
Hide 'PREFIX'W'n'zeta'm';
*---------------------------
*       Reversion of series for \eta: first iteration
*---------------------------
L 'PREFIX'eta0t'm'=eta'm';
.sort
Hide 'PREFIX'eta0t'm';
#enddo
*---------------------------
*   Iterations
*---------------------------
#do k=1,{'MAXORDER'/'n'}
#do m=1,'NDIM'
Drop  'PREFIX'eta{'k'-1}t'm';
#enddo
#do m=1,'NDIM'
L 'PREFIX'eta'k't'm'=eta'm'-alpha^'n'*'PREFIX'W'n'zeta'm';
#enddo
#do m=1,'NDIM'
*     efficient substitution
while ( match(etat'm') );
id,once etat'm'='PREFIX'eta{'k'-1}t'm';
endwhile;
#enddo
.sort
#enddo
#do m=1,'NDIM'
Hide 'PREFIX'eta{'MAXORDER'/'n'}t'm';
#enddo
.sort
*---------------------------
#do m=1,'NDIM'
L 'PREFIX'zeta{'MAXORDER'/'n'}t'm'=zeta'm'+alpha^'n'*'PREFIX'W'n'eta'm';
#enddo
#do mm=1,'NDIM'
#ifndef 'FASTPOWER'
while ( match(etat'mm') );
id,once etat'mm'='PREFIX'eta{'MAXORDER'/'n'}t'mm';
endwhile;
#else
#call polysubst(etat'mm','PREFIX'eta{'MAXORDER'/'n'}t'mm',{'MAXORDER'^%})
#endif
.sort 
#enddo
.sort
#do m=1,'NDIM'
Hide 'PREFIX'zeta{'MAXORDER'/'n'}t'm';
#enddo
.sort
*---------------------------
*  Inverted series \tilde \eta^{(n)}_i
*---------------------------
#do m=1,'NDIM'
Drop 'PREFIX'eta{'MAXORDER'/'n'}t'm','PREFIX'zeta{'MAXORDER'/'n'}t'm';
#enddo
#do m=1,'NDIM'
L 'PREFIX'etat'n'o'm'='PREFIX'eta{'MAXORDER'/'n'}t'm';
L 'PREFIX'zetat'n'o'm'='PREFIX'zeta{'MAXORDER'/'n'}t'm';
#enddo
.sort
Hide;
.sort
#enddo
*---------------------------
*  Gustavson integral
*---------------------------
*on stat;
.global
L 'PREFIX'It'MAXORDER' = 'Op';
.sort
Hide 'PREFIX'It'MAXORDER';
#do n={'MAXORDER'-1},0,-1
L 'PREFIX'It'n' = 'PREFIX'It{'n'+1};
#do m=1,'NDIM'
id  eta'm'=etat'm';
id zeta'm'=zetat'm';
#enddo
*Print;
.sort
#do m=1,'NDIM'
#call polysubst(etat'm','PREFIX'etat{'n'+1}o'm',{'MAXORDER'^%})
#call polysubst(zetat'm','PREFIX'zetat{'n'+1}o'm',{'MAXORDER'^%})
.sort
#enddo
Hide 'PREFIX'It'n';
#enddo
.sort
on statistics;
G 'PREFIX'Ginv'MAXORDER'= 'PREFIX'It0;
B alpha;
.store
#endprocedure
