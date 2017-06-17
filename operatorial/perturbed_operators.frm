*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		Kato series for perturbed operators S_H, P_H and D_H
*==========================================================================
#define MAXORDER	"4"
Functions L,S,P,I,Hi
#do n = 0,'MAXORDER'
 ,H'n'
#enddo
;
Symbols alpha(:{'MAXORDER'}),z(:{'MAXORDER'});
off statistics;
.global
Global H  = H0 + 
#do k=1,'MAXORDER'
 + alpha^'k'*H'k'
#enddo
;
B alpha;
Print;
.store
G LH=H;
#call toliouvillian()
B alpha;
Print;
.store
*-------------------- The unperturbed resolvent operator:
G RH0 = -P/z
#do n=0,{'MAXORDER'}
+ z^'n'*S^{'n'+1}
#enddo 
;
.store
*-------------------- The Neumann series:
G LV=LH-L(H0);
.store
Global F0= RH0;
.store
#do m=1,{'MAXORDER'}
Global F'm'= -RH0*LV*F{'m'-1};
.store
#enddo
*-------------------- The perturbed resolvent:
write statistics;
*    Compute W
Global RH=
#do m=0,'MAXORDER'
+ F'm'
#enddo
;
Bracket z;
.store
*-------------------- Perturbed Integrating, Averaging and  Quasi-nilpotent operators:
Global Sh = RH[z^0];
Global Ph = -RH[z^-1];
Global Dh = -RH[z^-2];
B alpha;
Print;
.store
*=======================================================
*                     Properties of the perturbed Kato operators
on statistics;
.global
*-------------------- 
*    1)  P_H S_H = 0
*-------------------- 
Local  PhSh=Ph*Sh;
Local  ShPh=Sh*Ph;
id P*S=0;
id S*P=0;
id P*P=P;
B alpha;
Print;
.store
*-------------------- 
*    2)  P_H P_H = P_H
*-------------------- 
Local  delta= Ph*Ph - Ph;
id P*S=0;
id S*P=0;
id P*P=P;
B alpha;
Print;
.store
*-------------------- 
*    3)  S_H D_H = 0
*-------------------- 
Local  ShDh=Sh*Dh;
Local  DhSh=Dh*Sh;
id P*S=0;
id S*P=0;
id P*P=P;
B alpha;
Print;
.store
*-------------------- 
*    4)  P_H D_H = D_H
*-------------------- 
Local  delta=Ph*Dh-Dh;
id P*S=0;
id S*P=0;
id P*P=P;
B alpha;
Print;
.store
*-------------------- 
*    5)  P_H H = H
*-------------------- 
Local  PhH=Ph*H;
*
*     You can uncomment three following lines to take a look on unsimplified expression:
*
*B alpha;
*Print;
*.end
.sort
*-------------------- 
*    		The following procedure implements simple canonical identities:
*-------------------- 
#call identities1()
id I=1;
B alpha;
Print;
.store
*-------------------- 
*    6)  S_H H = 0
*-------------------- 
Local  ShH=Sh*H;
#call identities1()
id I=1;
B alpha;
Print;
.store
*-------------------- 
*    7)  D_H H = 0
*-------------------- 
Local  DhH=Dh*H;
#call identities1()
id I=1;
B alpha;
Print;
.store

*-------------------- 
*    8)  L_H P_H = D_H
*-------------------- 
Local  delta= LH*Ph - Dh;
.sort
#call identities1()
id I=1;
B alpha;
Print;
.store
*-------------------- 
*    9)  L_H S_H = 1-P_H
*-------------------- 
Local  delta= LH*Sh -(1-Ph);
.sort
#call identities1()
id I=1;
B alpha;
Print;
.store
*-------------------- 
*    10)  L_H D_H = Dh^2
*-------------------- 
Local  delta= LH*Dh - Dh*Dh;
.sort
#call identities1()
id I=1;
B alpha;
Print;
.end
