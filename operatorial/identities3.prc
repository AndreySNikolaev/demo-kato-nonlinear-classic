#procedure identities3(a)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*    Reduction using canonical identities
*==========================================================================
repeat;
	id S*L(P,?a) = L(P,?a)*S;
endrepeat;
*id  S*L(S,?a) = -S*L(?a)*S+P*L(?a)*S*S+S*S*L(?a)*P+L(S,?a)*S+L(P,?a)*S*S;
id S*L(S,?a) = -(S*L(?a)*S-P*L(?a)*S*S-S*S*L(?a)*P-L(S,?a)*S)+L(P,?a)*S*S ;
#call identities5()
repeat;
	id L(P,?a)*S = S*L(P,?a);
endrepeat;
.sort
repeat;
	id S*L(P,?a) = L(P,?a)*S;
endrepeat;
id  S*L(S,?a) = -S*L(?a)*S+P*L(?a)*S*S+S*S*L(?a)*P+L(S,?a)*S+L(P,?a)*S*S;
#call identities5()
repeat;
	id L(P,?a)*S = S*L(P,?a);
endrepeat;
.sort
#endprocedure
