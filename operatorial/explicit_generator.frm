*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
#define MAXORDER	"3"
#define MAXREPS "10"
Functions L,S,P,Q,I,Op,aL
#do k = 0,'MAXORDER'+1
 ,H'k'
#enddo
;
Set op:L,S,P,Q,I,Op;
.global
nwrite statistics;
Symbols l,m,n,k,alpha(:{'MAXORDER'}),z(:{'MAXORDER'});
.global
Global H  = H0 + 
#do k=1,'MAXORDER'
 + alpha^'k'*H'k'
#enddo
;
B alpha;
Print;
.store
*-------------------- 
*     1. The perturbed integrating operator is:
*-------------------- 
#call RH()
G Sh=RH'MAXORDER'[z^0];
.store
*-------------------- 
*     2. dH/dalpha is:
*-------------------- 
G dHdalpha=H;
id alpha^n? = n*alpha^(n-1);
.store
*-------------------- 
*     3. The explicit expression for generator 
*         W=S_H dH/dalpha
*-------------------- 
G W=Sh*dHdalpha;
B alpha;
Print;
.store
*--------------------
*     4. U_G up to alpha^MAXORDER:
*--------------------
#call U(W,demo1)
Global U =demo1U; 
B alpha;
Print;
.store
*--------------------
*     5. The transformed Hamiltonian:
*--------------------
Global HH= U*H;
B alpha;
Print;
.sort
#call identities4()
.sort
#call identities4()
.sort
.sort
#call identities3()
write statistics;
id P*L(H1)*S*S*H1*I=0;
id P*L(H2)*S*S*H2*I=0;
id P*L(H2)*S*H1*I=P*L(H1)*S*H2*I;
id P*L(H3)*S*H1*I=P*L(H1)*S*H3*I;
id I=1;
B alpha;
Print;
.store
*--------------------
*   4. U_G^-1 up to alpha^MAXORDER:
*--------------------
#call Uinverse(W,demo2)
Global Uinverse =demo2Ui; 
B alpha;
Print;
.store
*--------------------
*  5. The Gustavson-Hori first integral I_G=U_G^-1 H0
*     and its nontrivial part  
*--------------------
Global IH0= Uinverse*H0;
B alpha;
*Print;
.sort
Global IG= (H-IH0)/alpha;
B alpha;
*Print;
.sort
*--------------------
*  6. Canonical simplifications of Hori integral
*--------------------
#call identities4()
B alpha;
*Print;
.sort
#call identities4()
.sort
id P*L(H2)*S*H1*I=P*L(H1)*S*H2*I;
id I=1;
B alpha;
Print;
.store
*--------------------
*  7. It is commutative with the perturbed Hamiltonian
*--------------------
G LH=H;
#call toliouvillian()
B alpha;
*Print;
.store
on statistics;
L delta=LH*IH0;
B alpha;
Print;
.sort
#call identities7()
B alpha;
Print;
.end
