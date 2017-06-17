#procedure toliouvillian()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*   Convert active expressions into Liouville operators 
*==========================================================================
.sort
*-------------------------
* 		Ancillary noncommutative operators 
*-------------------------
Functions aL,aOp;
multiply,left,aL;
multiply,right,I;

repeat;
	id aL(?a)*aOp? = aL(?a,aOp); 
	id aL(?a)*aOp?(?e) = aL(?a,aOp(?e)); 
endrepeat;
id aL(?a,I) = L(?a);
*Print;
.sort
#endprocedure
