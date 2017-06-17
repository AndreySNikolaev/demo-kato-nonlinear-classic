#procedure GustavsonNormalization2d(PREFIX)
*==========================================================================
*   This file is part of demos for 
*    "Kato perturbation expansion in classical mechanics ..."
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*    arXiv:1307.3368 [math.DS], 2013, submitted to J. Phys. A.
*==========================================================================
*  Normalization using Birkhoff-Gustavson algorithm
*       Gustavson FG (1966), Astron J 71:670
*  PREFIX - for uniqueness. Op - initial Hamiltonian.
*  Names of all objects created inside the procedure will begin with this PREFIX.
*==========================================================================
*on stat;
*.global
Global 'PREFIX'H = Hetazeta;
*--------------------
*     Warning: Bracketing alpha is mandatory here.
*--------------------
B alpha;
.store
#do k = 0,'MAXORDER'
G 'PREFIX'Gn0k'k'j0m0='PREFIX'H[alpha^'k'];
#enddo
.store
*--------------------
*	main loop 
*--------------------
#do n = 0,{'MAXORDER'-1}
#do k = 0,'MAXORDER'
G 'PREFIX'Fn'n'k'k'j0m0='PREFIX'Gn'n'k'k'j0m0;
#enddo
.store
*----------------------------------
*  Solution of the homological equation:
*----------------------------------
G 'PREFIX'W{'n'+1}=-'PREFIX'Fn'n'k{'n'+1}j0m0;
.sort
#call SH0()
.store
#do k = 0,'n'
L 'PREFIX'Gn{'n'+1}k'k'j0m0='PREFIX'Fn'n'k'k'j0m0;
.sort
Hide 'PREFIX'Gn{'n'+1}k'k'j0m0;
#enddo
*---------------------------
*  The derivatives of generating function
*---------------------------
#do m=1,'NDIM'
L 'PREFIX'W{'n'+1}eta'm'='PREFIX'W{'n'+1};
id eta'm'^k? = k*eta'm'^(k-1);
.sort
Hide 'PREFIX'W{'n'+1}eta'm';
L 'PREFIX'W{'n'+1}zeta'm'='PREFIX'W{'n'+1};
id zeta'm'^k? = k*zeta'm'^(k-1);
.sort
Hide 'PREFIX'W{'n'+1}zeta'm';
#enddo
L 'PREFIX'W{'n'+1}zetaj0m0=1;
L 'PREFIX'W{'n'+1}etaj0m0=1;
.sort
Hide 'PREFIX'W{'n'+1}zetaj0m0,'PREFIX'W{'n'+1}etaj0m0;
*---------------------------
*  Gustavson's formula:
*---------------------------
#do k = {'n'+1},'MAXORDER'
#do j = 1,{'k'/('n'+1)}
*---------------------------
*  m=0
*  F(n)(k)j0=1/j d F(n)(k-(n+1))j0/d eta2, 
*  G(n+1)(k)j0=1/j d G(n+1)(k-(n+1))j0/d zeta2: 
*---------------------------
L 'PREFIX'Fn'n'k'k'j'j'm0='PREFIX'Fn'n'k{'k'-('n'+1)}j{'j'-1}m0/'j';
id eta2^k? = k*eta2^(k-1);
.sort
Hide 'PREFIX'Fn'n'k'k'j'j'm0;
L 'PREFIX'Gn{'n'+1}k'k'j'j'm0='PREFIX'Gn{'n'+1}k{'k'-('n'+1)}j{'j'-1}m0/'j';
id zeta2^k? = k*zeta2^(k-1);
.sort
Hide 'PREFIX'Gn{'n'+1}k'k'j'j'm0;
L 'PREFIX'W{'n'+1}etaj'j'm0='PREFIX'W{'n'+1}etaj{'j'-1}m0*'PREFIX'W{'n'+1}eta2;
L 'PREFIX'W{'n'+1}zetaj'j'm0='PREFIX'W{'n'+1}zetaj{'j'-1}m0*'PREFIX'W{'n'+1}zeta2;
.sort
Hide 'PREFIX'W{'n'+1}etaj'j'm0,'PREFIX'W{'n'+1}zetaj'j'm0;
L 'PREFIX'F0n{'n'+1}k'k'j'j'm0= 
'PREFIX'Fn'n'k'k'j'j'm0*'PREFIX'W{'n'+1}zetaj'j'm0 
-'PREFIX'Gn{'n'+1}k'k'j'j'm0*'PREFIX'W{'n'+1}etaj'j'm0;
#call eidentities()
.sort
Hide 'PREFIX'F0n{'n'+1}k'k'j'j'm0;
#do m=1,'j'
*---------------------------
*  F(n)(k)jm=1/m d F(n)(k-(n+1))(j-1)(m-1)/d eta1, 
*  G(n+1)(k)jm=1/m d G(n+1)(k-(n+1))(j-1)(m-1)/d zeta1 
*---------------------------
L 'PREFIX'Fn'n'k'k'j'j'm'm'='PREFIX'Fn'n'k{'k'-('n'+1)}j{'j'-1}m{'m'-1}/'m';
id eta1^k? = k*eta1^(k-1);
.sort
Hide 'PREFIX'Fn'n'k'k'j'j'm'm';
L 'PREFIX'Gn{'n'+1}k'k'j'j'm'm'='PREFIX'Gn{'n'+1}k{'k'-('n'+1)}j{'j'-1}m{'m'-1}/'m';
id zeta1^k? = k*zeta1^(k-1);
.sort
Hide 'PREFIX'Gn{'n'+1}k'k'j'j'm'm';
L 'PREFIX'W{'n'+1}etaj'j'm'm'='PREFIX'W{'n'+1}etaj{'j'-1}m{'m'-1}*'PREFIX'W{'n'+1}eta1;
L 'PREFIX'W{'n'+1}zetaj'j'm'm'='PREFIX'W{'n'+1}zetaj{'j'-1}m{'m'-1}*'PREFIX'W{'n'+1}zeta1;
.sort
Hide 'PREFIX'W{'n'+1}etaj'j'm'm','PREFIX'W{'n'+1}zetaj'j'm'm';
L 'PREFIX'F0n{'n'+1}k'k'j'j'm'm'= 
     'PREFIX'Fn'n'k'k'j'j'm'm'*'PREFIX'W{'n'+1}zetaj'j'm'm'
-'PREFIX'Gn{'n'+1}k'k'j'j'm'm'*'PREFIX'W{'n'+1}etaj'j'm'm'
;
#call eidentities()
.sort
Hide 'PREFIX'F0n{'n'+1}k'k'j'j'm'm';
#enddo
#enddo
L 'PREFIX'Gn{'n'+1}k'k'j0m0= 
'PREFIX'Fn'n'k'k'j0m0
#do j = 1,{'k'/('n'+1)}
#do m=0,'j'
 + 'PREFIX'F0n{'n'+1}k'k'j'j'm'm'
#enddo
#enddo
;
.sort
Hide 'PREFIX'Gn{'n'+1}k'k'j0m0;
#enddo
#enddo
*--------------------
*	Ht and pseudo generating function W
*--------------------
.sort
write statistics;
Global  'PREFIX'G{'MAXORDER'}= 
#do n = 1,'MAXORDER'
+ alpha^'n'*'PREFIX'W'n'
#enddo
;
Global  'PREFIX'Ht{'MAXORDER'} =
#do n = 0,'MAXORDER'
+ alpha^'n'*'PREFIX'Gn'MAXORDER'k'n'j0m0
#enddo
;
B alpha;
.store
#endprocedure
