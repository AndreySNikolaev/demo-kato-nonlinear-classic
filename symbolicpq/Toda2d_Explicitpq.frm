#:TermsInSmall   200000
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*	The normal form and the Gustavson integral
*       for the Toda2d Hamiltonian (Ford 1973),
*       using the explicit formula for Deprit generator.
*       Normalisation in p,q variables.
*===================================================================
*          MAXORDER - maximum order in alpha,
*          NDIM - number of dimensions
*-------------------- 
#define MAXORDER	"10"
#define NDIM    	"2"
CFunctions sqrt,qurt,asqrt,aqurt,exp; 
Symbols l,m,n,k,alpha(:{'MAXORDER'}), z,
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
*	1. The Hamiltonian in P,Q  variables is:
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
Global Hpq = H;
B alpha;
.store
*-------------------- 
*	The unperturbed part is as follows:
*-------------------- 
Global H0 = Hpq[alpha^0];
Print;
.store
*-------------------- 
*     3. The explicit expression for generator is as follows: 
*-------------------- 
#call Wpq(e)
Global W = eSH;
B alpha;
Print;
.store
*--------------------
*   4. The transformed Hamiltonian U_W H
*        in p,q variables is as follows:
*--------------------
#call Upq(W,Hpq,e)
.store
G Htpq=eU'MAXORDER';
B alpha;
Print;
.store
*--------------------
*   5. The transformed Hamiltonian
*      in eta,zeta variables is as follows:
*--------------------
G Htetazeta=Htpq;
#call pqtoetazeta()
B alpha;
Print;
.store
*--------------------
*   6. The normalised Hamiltonian is commutative with H0:
*--------------------
#call Lpq(H0,Htpq)
Drop Pb'NDIM';
L delta=Pb'NDIM';
B alpha;
Print;
.store
*--------------------
*  7. The Gustavson-Hori first integral I_G=U_G^-1 H0
*--------------------
G IGH=H0;
B alpha;
*Print;
.store
#call Uinversepq(W,IGH,h)
G IGH0= hUinv'MAXORDER';
B alpha;
Print;
.store
*--------------------
*   8. It is commutative with the perturbed Hamiltonian 
*--------------------
#call Lpq(Hpq,IGH0)
B alpha;
Print;
.store
*--------------------
*   9. The formal integral in  eta,zeta variables is:
*--------------------
G IGH0etazeta=IGH0;
#call pqtoetazeta()
B alpha;
Print;
.store
*--------------------
*   10. Its nontrivial part 
*--------------------
G Ipq=(H-IGH0)/alpha^2;
B alpha;
Print;
.end
