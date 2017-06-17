#procedure eidentities()
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		Identities for elementary functions
*==========================================================================
repeat;
	id asqrt(1) = 1;
	id qurt(1) = 1;
	id aqurt(1) = 1;
	id sqrt(1) = 1;
	id asqrt(1) = 1;
	id sqrt(n?)*asqrt(n?) = 1;
	id qurt(n?)*aqurt(n?) = 1;
	id sqrt(n?)*sqrt(n?) = n;
	id asqrt(n?)*asqrt(n?) = 1/n;
	id qurt(n?)*qurt(n?) = sqrt(n);
	id aqurt(n?)*aqurt(n?) = asqrt(n);

	id sqrt(4) = 2;
	id asqrt(4) = 1/2;
	id qurt(4) = sqrt(2);
	id aqurt(4) = asqrt(2);

	id sqrt(16) = 4;
	id asqrt(16) = 1/4;
	id qurt(16) = 2;
	id aqurt(16) = 1/2;
*	id sqrt(25) = 5;
*	id asqrt(25)=1/5;
*	id sqrt(49)=7;
*	id asqrt(49)=1/7;
*	id qurt(100) = sqrt(10);
*	id aqurt(100) = asqrt(10);
*	id qurt(196) = sqrt(14);
*	id aqurt(196) = asqrt(14);
endrepeat;
id asqrt(n?)=sqrt(n)/n;
#endprocedure
