#procedure identities4(a)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    Reduction using canonical identities
*==========================================================================
.sort
Functions aL;
multiply,right,I;
repeat;
repeat;
	id I*I = I;
*---------------------- 
*        aL,I - ancillary operators used for implementation of identity 
*           L_(...) H0 = - L_H0 (...) 
*---------------------- 
	id L(?a)*H0*I = -L(H0)*aL(?a)*I;
	id P*L(H0) = 0;
	id L(H0)*P = 0;
	id S*L(H0) = (1-P);
	id L(H0)*S = (1-P);
endrepeat;
#do n = 0,'MAXORDER'+1
	id L(?a)*H'n'*I = -L(H'n')*aL(?a)*I;
#enddo
repeat;
      id aL(Op?,?a) = Op*aL(?a);
	id aL(Op?(?e),?a) = Op(?e)*aL(?a); 
endrepeat;
endrepeat;
id aL = 1;
.sort
#do jj=1,'MAXREPS'
repeat;
id S*L(P,?a) = L(P,?a)*S;
endrepeat;
id L(P,?a)*S*S = (S*L(?a)*S-P*L(?a)*S*S-S*S*L(?a)*P-L(S,?a)*S+S*L(S,?a));
#call identities5()

repeat;
id S*L(P,?a) = L(P,?a)*S;
endrepeat;
.sort
id  L(S,?a)*S  = (S*L(?a)*S-P*L(?a)*S*S-S*S*L(?a)*P +S*L(S,?a)) - L(P,?a)*S*S;
#call identities5()

repeat;
id L(P,?a)*S = S*L(P,?a);
endrepeat;
.sort
id S*S*L(P,?a) =(S*L(?a)*S-P*L(?a)*S*S-S*S*L(?a)*P-L(S,?a)*S+S*L(S,?a));
#call identities5()

repeat;
id L(P,?a)*S = S*L(P,?a);
endrepeat;
.sort
#enddo
#endprocedure
