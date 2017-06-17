*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*	Comparison of the normal forms, generators, and the Gustavson integrals
*       for the Deprit triangular and the explicit algorithms.
*       Toda2d System.
*       Deprit A (1969), Celestial Mech 1:12–30
*==========================================================================
*          MAXORDER - maximum order in alpha,
*          NDIM - number of dimensions
*-------------------- 
#define MAXORDER	"10"
#define NDIM    	"2"
CFunctions sqrt,qurt,asqrt,aqurt,exp; 
Symbols l,m,n,k,alpha(:{'MAXORDER'}),z
*-------------------- 
*			The following are the canonical variables:
*-------------------- 
#do j=1,'NDIM'
	zeta'j',eta'j',P'j',Q'j',p'j',q'j'
#enddo
;
Off statistics;
.global
*-------------------- 
*	1. The Hamiltonian in P,Q and eta,zeta variables is:
*-------------------- 
#define NAME "Toda2d"
Global HPQ=1/2*(P1^2 + P2^2)+ 1/24*(exp(2*Q2+2*sqrt(3)*Q1)+exp(2*Q2-2*sqrt(3)*Q1)+exp(-4*Q2))-1/8;
B alpha;
Print;
.sort
id exp(m?)=sum_(n,0,{'MAXORDER'+2},m^n*invfac_(n));
#call eidentities()
B alpha;
Print;
.store
*-------------------- 
*	2. A scaling transform introduces a small parameter:
*-------------------- 
G H=HPQ/alpha^2;
#do j=1,'NDIM'
id P'j'=alpha*p'j';
id Q'j'=alpha*q'j';
#enddo
.sort
#call eidentities()
B alpha;
Print;
.store
#call frequencies()
Global Hetazeta = H;
#call pqtoetazeta()
B alpha;
Print;
.store
*-------------------- 
*	The unperturbed part is as follows:
*-------------------- 
Global H0 = Hetazeta[alpha^0];
Print;
.store
*------------------------------------------------------------------------------ 
*     3. The Deprit triangular normalisation 
*        and the corresponding Gustavson-Hori integral: 
*------------------------------------------------------------------------------ 
#call Deprit3Normalization(dt)
Global W3=dtW;
Global d3Ht=dtHt'MAXORDER';
.store
#call Deprit3Inverse(dtW,H0,dti)
Global d3IGH0=dtiUinv'MAXORDER'; 
.store
*------------------------------------------------------------------------------ 
*     4. The explicit expression for the generator: 
*------------------------------------------------------------------------------ 
#call W(e)
Global W = eSH;
.store
*--------------------
*    The transformed Hamiltonian U_W H:
*--------------------
#call U(W,Hetazeta,e)
.store
Global eHt=eU'MAXORDER';
.store
*--------------------
*  The Gustavson-Hori integral U_G^-1 H0 :
*--------------------
#call Uinverse(W,H0,ei)
G eIGH0= eiUinv'MAXORDER';
.store
*------------------------------------------------------------------------------ 
*   5. The comparison between nonsecular and natural uniqueness conditions. 
*     1) The Gustavson integrals are identical:
*------------------------------------------------------------------------------ 
write statistics;
L deltaGustavsonIntegrals=eIGH0-d3IGH0;
B alpha;
Print;
.store
*-------------------- 
*     2) The generators differ starting from the 5th order:
*-------------------- 
write statistics;
G deltaGenerators=W-W3;
B alpha;
Print;
.store
*-------------------- 
*     3) The normalised Hamiltonians differ starting from the 8th order:
*-------------------- 
L deltaHt=eHt-d3Ht;
B alpha;
Print;
.store
*-------------------- 
*    4) The normalized Hamiltonians are connected by the Deprit transform
*    with the generator:
*
*    W_21=U_W deltaGenerators
*-------------------- 
#call U(W,deltaGenerators,dw)
Global W21=dwU'MAXORDER';
id alpha^'MAXORDER'=0;
B alpha;
Print;
.store
*--------------------
*    This generator W_21 is secular (1-P_H0)W_21=0 :
*-------------------- 
Local PH0W21=W21;
#call PH0()
.sort
Drop PH0W21;
Local delta=W21-PH0W21;
B alpha;
Print;
.store
*-------------------- 
*    The normalised Hamiltonians are connected by 
*    the Deprit exponent with the generator W_21:
*-------------------- 
#call Uinverse(W21,eHt,h21)
write statistics;
L delta=d3Ht-h21Uinv'MAXORDER';
B alpha;
Print;
.end
