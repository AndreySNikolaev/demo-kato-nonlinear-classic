*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Comparison of the Gustavson [1], Hori [2], Deprit [3], Henrard [4],
*    Dragt-Finn [5] and the explicit [6]  algorihms for the Toda 2D system
*  1. Gustavson FG (1966), Astron J 71:670
*  2. Hori G (1966), Publ Astron Soc Japan 18:287
*  3. Deprit A (1969), Celestial Mech 1:12–30
*  4. Henrard, J., in "Recent advances in dynamical astronomy" ed. Szebehely, p. 248, 	
*  5. Dragt AJ, Finn JM (1976), J. Math. Phys. 17:2215–2227
*  6. Nikolaev A.S., arXiv:1612.05207
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
*	The Hamiltonian in P,Q  variables is:
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
*	A scaling transform introduces a small parameter:
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
*     1. The explicit expression for the generator: 
*        Nikolaev A.S., arXiv:1612.05207
*------------------------------------------------------------------------------ 
#call W(e)
Global W = eSH;
.store
*--------------------
*         The transformed Hamiltonian U_W H:
*--------------------
#call U(W,Hetazeta,e)
.store
Global HtExplicit=eU'MAXORDER';
.store
*--------------------
*        The Gustavson-Hori integral U_G^-1 H0 :
*--------------------
#call Uinverse(W,H0,ei)
G IGH0Explicit= eiUinv'MAXORDER';
.store
*------------------------------------------------------------------------------ 
*     2. The Henrard normalisation.
*     Henrard, J., in "Recent advances in dynamical astronomy" ed. Szebehely, p. 248	
*     Algorithm by Koseleff PV (1994), Cel. Mech. 58:17-36
*------------------------------------------------------------------------------ 
#call HenrardNormalization(hk)
Global HtHenrard=hkHt'MAXORDER';
.store
#call HenrardInverse(hkW'MAXORDER',H0,hki)
Global IGH0Henrard= hkiHinv'MAXORDER';
.store
*------------------------------------------------------------------------------ 
*     The comparison of the Explicit and Henrard normalisations.
*------------------------------------------------------------------------------ 
*      1) The Gustavson integrals are identical:
*-----------------
write statistics;
Local dGustavsonIntegrals=IGH0Henrard-IGH0Explicit;
B alpha;
Print;
.store
*      2) The Hamiltonians normalised by the Henrard and Explicit methods are identical:
*-----------------
write statistics;
Local dHtExplicitvsHenrard=HtHenrard-HtExplicit;
B alpha;
Print;
.store
*------------------------------------------------------------------------------ 
*     3. The Deprit triangular normalisation. 
*        Deprit A (1969), Celestial Mech 1:12–30 
*------------------------------------------------------------------------------ 
#call Deprit3Normalization(dt)
Global HtDeprit3=dtHt'MAXORDER';
.store
#call Deprit3Inverse(dtW,H0,dti)
Global IGH0Deprit3=dtiUinv'MAXORDER'; 
.store
*      1) The Gustavson integrals are identical:
*-----------------
write statistics;
Local dGustavsonIntegrals=IGH0Henrard-IGH0Deprit3;
B alpha;
Print;
.store
*-------------------- 
*     2) The normalised Hamiltonians differ starting from the 8th order:
*-------------------- 
Local dHtDeprit3vsHenrard=HtHenrard-HtDeprit3;
B alpha;
Print;
.store
*------------------------------------------------------------------------------ 
*     4. The Hori normalisation. 
*       Hori G (1966), Publ Astron Soc Japan 18:287
*       Algorithm by Mersman W.A., Celest. Mech. 3 (1970) 81–89.
*------------------------------------------------------------------------------ 
#call HoriMersmanNormalization(hm)
Global HtHori=hmHt'MAXORDER';
.store
#call HoriMersmanInverse(hmG'MAXORDER',H0,h)
Global IGH0Hori=hHMinv'MAXORDER'; 
.store
*      1) The Gustavson integrals are identical:
*-----------------
write statistics;
Local dGustavsonIntegrals=IGH0Hori-IGH0Deprit3;
B alpha;
Print;
.store
*-------------------- 
*     3) The normalised Hamiltonians differ starting from the 8th order:
*-------------------- 
Local dHtHorivsDeprit3=HtHori-HtDeprit3;
Local dHtHorivsHenrard=HtHori-HtHenrard;
B alpha;
Print;
.store
*------------------------------------------------------------------------------ 
*     5. The Dragt-Finn normalisation. 
*       Dragt AJ, Finn JM (1976), J. Math. Phys. 17:2215–2227
*       Algorithm by Koseleff PV (1994), Cel. Mech. 58:17–36.
*------------------------------------------------------------------------------ 
#call DragtFinnNormalization(df)
Global HtDragtFinn=dfHt'MAXORDER';
.store
#call DragtFinnInverse(dfG'MAXORDER',H0,dfi)
Global IGH0DragtFinn=dfiDFinv'MAXORDER';
.store
*      1) The Gustavson integrals are identical:
*-----------------
write statistics;
Local dGustavsonIntegrals=IGH0DragtFinn-IGH0Deprit3;
B alpha;
Print;
.store
*-------------------- 
*     3) The normalised Hamiltonians differ starting from the 8th order:
*-------------------- 
Local dHtDragtFinnvsDeprit3=HtDragtFinn-HtDeprit3;
Local dHtDragtFinnvsHenrard=HtDragtFinn-HtHenrard;
Local dHtDragtFinnvsHori   =HtDragtFinn-HtHori;
B alpha;
Print;
.store
*------------------------------------------------------------------------------ 
*     6. The Gustavson normalisation. 
*       Gustavson FG (1966), Astron J 71:670
*------------------------------------------------------------------------------ 
#call GustavsonNormalization2d(g)
Global HtGustavson=gHt'MAXORDER';
.store
#call GustavsonInverse(gG'MAXORDER',H0,gh)
Global IGH0Gustavson=ghGinv'MAXORDER';
.store
*      1) The Gustavson integrals are identical:
*-----------------
write statistics;
Local dGustavsonIntegrals=IGH0Gustavson-IGH0Deprit3;
B alpha;
Print;
.store
*-------------------- 
*     3) The normalised Hamiltonians differ starting from the 8th order:
*-------------------- 
Local dHtGustavsonvsDeprit3  =HtGustavson-HtDeprit3;
Local dHtGustavsonvsHenrard  =HtGustavson-HtHenrard;
Local dHtGustavsonvsHori     =HtGustavson-HtHori;
Local dHtGustavsonvsDragtFinn=HtGustavson-HtDragtFinn;
B alpha;
Print;
.store
.end
