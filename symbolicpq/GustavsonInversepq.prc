#procedure GustavsonInversepq(W,Op,PREFIX)
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
  pt'm',qt'm'
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
id  p'm'=pt'm';
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
L 'PREFIX'W'n'p'm'='PREFIX'W'n';
id pt'm'^k? = k*pt'm'^(k-1);
.sort
Hide 'PREFIX'W'n'p'm';
.sort
L 'PREFIX'W'n'q'm'='PREFIX'W'n';
id q'm'^k? = k*q'm'^(k-1);
.sort
Hide 'PREFIX'W'n'q'm';
*---------------------------
*       Reversion of series for q: first iteration
*---------------------------
L 'PREFIX'p0t'm'=p'm';
.sort
Hide 'PREFIX'p0t'm';
#enddo
*---------------------------
*   Iterations
*---------------------------
#do k=1,{'MAXORDER'/'n'}
#do m=1,'NDIM'
Drop  'PREFIX'p{'k'-1}t'm';
#enddo
#do m=1,'NDIM'
L 'PREFIX'p'k't'm'=p'm'-alpha^'n'*'PREFIX'W'n'q'm';
#enddo
#do m=1,'NDIM'
*     efficient substitution
while ( match(pt'm') );
id,once pt'm'='PREFIX'p{'k'-1}t'm';
endwhile;
#enddo
.sort
#enddo
#do m=1,'NDIM'
Hide 'PREFIX'p{'MAXORDER'/'n'}t'm';
#enddo
.sort
*---------------------------
#do m=1,'NDIM'
L 'PREFIX'q{'MAXORDER'/'n'}t'm'=q'm'+alpha^'n'*'PREFIX'W'n'p'm';
#enddo
#do mm=1,'NDIM'
#ifndef 'FASTPOWER'
while ( match(pt'mm') );
id,once pt'mm'='PREFIX'p{'MAXORDER'/'n'}t'mm';
endwhile;
#else
#call polysubst(pt'mm','PREFIX'p{'MAXORDER'/'n'}t'mm',{'MAXORDER'^%})
#endif
.sort 
#enddo
.sort
#do m=1,'NDIM'
Hide 'PREFIX'q{'MAXORDER'/'n'}t'm';
#enddo
.sort
*---------------------------
*  Inverted series \tilde p^{(n)}_i
*---------------------------
#do m=1,'NDIM'
Drop 'PREFIX'p{'MAXORDER'/'n'}t'm','PREFIX'q{'MAXORDER'/'n'}t'm';
#enddo
#do m=1,'NDIM'
L 'PREFIX'pt'n'o'm'='PREFIX'p{'MAXORDER'/'n'}t'm';
L 'PREFIX'qt'n'o'm'='PREFIX'q{'MAXORDER'/'n'}t'm';
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
*Print;
.sort
Hide 'PREFIX'It'MAXORDER';
#do n={'MAXORDER'-1},0,-1
L 'PREFIX'It'n' = 'PREFIX'It{'n'+1};
#do m=1,'NDIM'
id  p'm'=pt'm';
id q'm'=qt'm';
#enddo
*Print;
.sort
#do m=1,'NDIM'
#call polysubst(pt'm','PREFIX'pt{'n'+1}o'm',{'MAXORDER'^%})
#call polysubst(qt'm','PREFIX'qt{'n'+1}o'm',{'MAXORDER'^%})
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
