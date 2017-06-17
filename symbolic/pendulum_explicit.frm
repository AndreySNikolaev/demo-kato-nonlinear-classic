*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*	The normal form of the Pendulum Hamiltonian up to 10th order 
*	 using the explicit formula  for Deprit generator.
*===================================================================
*          MAXORDER - maximum order in alpha,
*          NDIM - number of dimensions
*-------------------- 
#define MAXORDER	"10"
#define NDIM    	"1"
CFunctions sqrt,qurt,asqrt,aqurt; 
Symbols l,m,n,k,alpha(:{'MAXORDER'}),z, J, cosphi, sinphi,
*-------------------- 
*			The canonical variables are:
*-------------------- 
#do j=1,'NDIM'
	zeta'j',eta'j',P'j',Q'j',p'j',q'j'
#enddo
;
nwrite statistics;
.global
*-------------------- 
*	1. The Hamiltonian in P,Q  variables is:
*            H = 1/2*P1^2+(1-cos(Q1))
*-------------------- 
#define NAME "Pendulum" 
Global HPQ = 1/2*P1^2 +1 
#do k=0,{'MAXORDER'+1}
-(-1)^{'k'}*(Q1)^({2*'k'})/fac_({2*'k'})
#enddo
;
Bracket Q1;
Print;
.store
*-------------------- 
*	2. A scaling transform introduces a small parameter:
*-------------------- 
G H=HPQ/alpha;
#do j=1,'NDIM'
id P'j'=sqrt(alpha)*p'j';
id Q'j'=sqrt(alpha)*q'j';
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
*	2. Unperturbed part and perturbation 
*-------------------- 
Global H0 = Hetazeta[alpha^0];
Print;
.store
*-------------------- 
*     3. The explicit expression for generator is as follows: 
*-------------------- 
#call W(e)
Global W =  eSH;
B alpha;
*Print;
.store
*--------------------
*   4. The transformed Hamiltonian U_W H
*        in eta,zeta variables is as follows:
*--------------------
#call U(W,Hetazeta,e)
.store
G Htetazeta=eU'MAXORDER';
B alpha;
Print;
.store
*--------------------
*   5. The transformed Hamiltonian
*      in p,q variables are as follows:
*--------------------
G Htpq=Htetazeta;
G Wt=W;
#call etazetatopq()
B alpha;
Print;
.store
*--------------------
*   6. The normalised Hamiltonian is commutative with H0:

*--------------------
#call L(H0,Htetazeta)
Drop Pb'NDIM';
L delta=Pb'NDIM';
B alpha;
Print;
.store
#if 'NDIM' == 1
*--------------------
*   7. The Generator and the transformed Hamiltonian
*      in Action-Angle variables are as follows:
*--------------------
G WAA=Wt;
id p1=sqrt(2)*sqrt(J)*cosphi;
id q1=sqrt(2)*sqrt(J)*sinphi;
#call eidentities()
B alpha;
Print;
.store
G Htaa=Htpq;
id p1^2=2*J-q1^2;
B alpha;
Print;
.store
*--------------------
*   8. Birkhoff series for the transformed Hamiltonian 
*--------------------
Symbol [P^2+Q^2];
G HtPQ=Htaa*alpha;
id J=[P^2+Q^2]/2/alpha;
B [P^2+Q^2];
Print;
.store
#endif
Symbol [p^2+p^2];
G Hprint=Htpq;
id q1^2=[p^2+p^2]-p1^2;
B alpha;
Print;
.end