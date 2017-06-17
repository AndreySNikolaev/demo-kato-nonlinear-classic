#procedure identities1(a)
.sort
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    Reduction using canonical identities
*==========================================================================
*      Ancillary functions
*----------------------------
Functions aL, aOP,I;
multiply,right,I; 

repeat;
repeat;
	id P*P = P;
	id S*P = 0;
	id P*S = 0;
	id I*I = I;
	id L(?a)*H0*I = -L(H0)*aL(?a)*I;
	id P*L(H0) = 0;
	id P*H0*I = H0*I;
	id S*H0*I = 0;
	id L(H0)*P = 0;
	id S*L(H0) = (1-P);
	id L(H0)*S = (1-P);
#do n = 1,'MAXORDER'
	id L(P,H'n')*P*H'n'*I = 0;
	id L(?a)*P*H'n'*I = -L(P,H'n')*aL(?a)*I;
#enddo
repeat;
        id aL(aOP?,?a) = aOP*aL(?a);
	id aL(aOP?(?e),?a) = aOP(?e)*aL(?a); 
endrepeat;
endrepeat;
   id aL*I = I;
endrepeat;
.sort
#do m = 1,'MAXORDER'
repeat;
repeat;
#do n = 1,'MAXORDER'
	id L(H'n')*H'n'*I =0;
#enddo
	id L(?a)*H'm'*I = -L(H'm')*aL(?a)*I;
repeat;
        id aL(aOP?,?a) = aOP*aL(?a);
	id aL(aOP?(?e),?a) = aOP(?e)*aL(?a); 
endrepeat;
endrepeat;
   id aL*I = I;
endrepeat;
.sort
#enddo
#endprocedure
