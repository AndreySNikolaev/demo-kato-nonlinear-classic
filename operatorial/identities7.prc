#procedure identities7()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		Set of canonical identities 
*==========================================================================
.sort
Functions aL;
multiply,right,I;
repeat;
repeat;
	id I*I = I;
#do n = 0,'MAXORDER'+1
      id L(H'n')*H'n'*I=0;
      id L(H'n')*H0*I = -L(H0)*H'n'*I;
      id P*L(H'n')*P*H'n'*I=0;
#enddo
      id P*L(H2)*P*H1*I = - P*L(H1)*P*H2*I;
      id P*L(H3)*P*H1*I = - P*L(H1)*P*H3*I;
      id P*L(H1)*S*L(H1)*P*H1*I = - 1/2* P*L(H1)*P*L(H1)*S*H1*I;
*---------------------- 
*        aL,I - ancillary operators used for implementation of identity 
*           L_(...) H0 = - L_H0 (...) 
*---------------------- 
*	id L(?a)*H0*I = -L(H0)*aL(?a)*I;
	id P*L(H0) = 0;
	id L(H0)*P = 0;
	id S*L(H0) = (1-P);
	id L(H0)*S = (1-P);
endrepeat;
	id L(?a)*H1*I = -L(H1)*aL(?a)*I;
repeat;
      id aL(Op?,?a) = Op*aL(?a);
	id aL(Op?(?e),?a) = Op(?e)*aL(?a); 
endrepeat;
endrepeat;
id aL = 1;
.sort
#endprocedure
