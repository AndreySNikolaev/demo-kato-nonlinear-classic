*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Properties of the Deprit U_G and U_G^-1
*==========================================================================
#define MAXORDER	"6"

*-----------------------
* 1. Noncommutative functions represent operators
*    L(G0) means Liouvillian operator L_{G0} = [.,G0] 
*-----------------------
Functions L,I,
#do n = 0,'MAXORDER'-1
 ,W'n'
#enddo
;
.global
off statistics;
Symbols k,alpha(:{'MAXORDER'});
.global
*--------------------
*   General Generator: 
*--------------------
Global W = 
#do n = 0,'MAXORDER'-1
+ alpha^'n'*W'n'
#enddo
;
B alpha;
Print;
.store
*--------------------
*  2. Compute U_W up to alpha^MAXORDER. 
*       Arguments are: generator and an unique PREFIX
*--------------------
#call U(W,demo1)
.store
Global U =demo1U; 
B alpha;
Print;
.store
*--------------------
*  3. Compute U_W^-1 up to alpha^MAXORDER.
*--------------------
#call Uinverse(W,demo2)
.store
Global Uinv =demo2Ui; 
B alpha;
Print;
.store
on statistics;
.global
*--------------------
*  4. Demos for properties of Deprit exponents
*    1) They are mutually inverse:
*--------------------
L EEinv= U*Uinv;
L EinvE= Uinv*U;
B alpha;
Print;
.store
*--------------------
*    2) d(U)/d(alpha) = U L_W :
*--------------------
G LW=W;
*--------------------
*	  Calculate L_W. 
*--------------------
#call toliouvillian()
B alpha;
Print;
.store
*--------------------
*	  Compute d(U)/d(alpha) up to alpha^(MAXORDER-1)
*--------------------
L DUDalpha= U;
id alpha^k? = k*alpha^(k-1);
B alpha;
Print;
.sort
*--------------------
*	Because we computed the derivative only up to alpha^(MAXORDER-1), 
*       the identitiy holds up to the same order
*--------------------
L Delta=DUDalpha - U*LW;
id alpha^'MAXORDER' =0;
B alpha;
Print;
.store
*--------------------
*    3) d(U^-1)/d(alpha) = -L_W U:
*--------------------
L DUinvDalpha= Uinv;
id alpha^k? = k*alpha^(k-1);
B alpha;
Print;
.sort
L Delta=DUinvDalpha + LW*Uinv;
id alpha^'MAXORDER' =0;
B alpha;
Print;
.end
