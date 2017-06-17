#procedure identities5(a)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    Reduction using canonical identities
*==========================================================================
repeat;
repeat;
	id P*P = P;
	id S*P = 0;
	id P*S = 0;
	id L(P,S,?a) = 0;
	id L(S,P,?a) = 0;
	id L(P,P,?a) = L(P,?a);
	id L(?a,P,P,?s) = L(?a,P,?s);
	id L(?a,P,S,?s) = 0;
	id L(?a,S,P,?s) = 0;
	id L(?a,P,P,?s) = L(?a,P,?s);
	id P*L(H1)*P*H1*I = 0;
	id P*L(H2)*P*H2*I = 0;
endrepeat;
repeat;
	id P*L(S,?a) = -P*L(?a)*S;
	id L(S,?a)*P = S*L(?a)*P;
endrepeat;
repeat;
	id L(P,?a)*P = P*L(?a)*P;
	id P*L(P,?a) = P*L(?a)*P;
	id L(L(?H1),?a) = L(?H1)*L(?a)-L(?a)*L(?H1);
endrepeat;
endrepeat;
.sort
#endprocedure
