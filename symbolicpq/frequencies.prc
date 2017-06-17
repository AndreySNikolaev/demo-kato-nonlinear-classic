#procedure frequencies(a)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*    Andrey S. Nikolaev (Andrey.Nikolaev@rdtex.ru), 
*==========================================================================
*		Calculation of the frequencies omega_i.
*           H0 assumed to be already diagonalized H0=Sum (H0p_i*p_i^2+ H0q_i*q_i^2)
*==========================================================================		
#do i=1,'NDIM'
Local HO = H;
id alpha = 0;
Bracket p'i',q'i';
.sort
Local Tm1 = HO[p'i'^2]*HO[q'i'^2];
Local Tm2 = HO[q'i'^2]/HO[p'i'^2];
#call eidentities()
.sort
G Omega'i' = sqrt(4*Tm1);
G Coef1'i' = aqurt(4*Tm2);
G Coef2'i' = qurt(4*Tm2)/2;
#call eidentities()
Print Omega'i';
*,Coef1'i',Coef2'i';
.store
#enddo
#endprocedure
