# demo-kato-nonlinear-classic
The supplementary demonstration files for:   "Generalisation of the explicit expression for the Deprit generator     to Hamiltonians nonlinearly dependent on small parameter",       arXiv:1612.05207    Andrey Nikolaev, ICPT, RDTeX, Provino, Russia, 2016 

----------------------------------------
There are two sets of demonstrations:

 "Operatorial" subdirectory contains the scripts illustrating: 
  The algebraic properties of P and S operators; Deprit transforms;
  Kato expansions; the explicit formula for Deprit generator.

 "Symbolic" subdirectory contains the examples of perturbative expansions 
  for particular systems (the Pendulum and the Toda 2D system).

 "Symbolicpq" subdirectory contains the examples of normalisation in p,q variables.

I. Prerequisites:

 The symbolic manipulation system FORM [2] can be downloaded from
 http://www.nikhef.nl/~form/maindir/binaries/binaries.html

 The Windows version of Form is available at:
 http://www.nikhef.nl/~form/oldsite/maindir/binaries/windows/windows.html

 Place the FORM executable "form.exe" into the "operatorial" and both "symbolic" subdirectories.


II. Notations:

 Operators  are represented in the FORM system by noncommutative functions.

 I is the identity,  L(G) means Liouvillian L_G, 
 P and S is unperturbed averaging and integrating operators P_H0 and S_H0.
 L(S,F) means L_{SF} and so on.

 * Lines beginning with the asterisk  are comments.

 In order to make the output compact, demonstrations usually were performed  up to o(alpha^10). 
 If needed, the MAXORDER variable can be set to the desired value.

II. Demonstrations:

 1. "Operatorial":

  1) Properties of the Deprit transforms U_G and U_G^-1:

   Usage:
     form -l Deprit_transforms.frm
   Output:
     Deprit_transforms.log

  2) Kato series for perturbed operators:

   Usage:
     form -l perturbed_operators.frm
   Output:
     perturbed_operators.log

  4) The explicit expression for the Deprit generator:  

   Usage:
     form -l explicit_generator.frm
   Output:
     explicit_generator.log	

  5) The Hamiltonian normalized to o(alpha^6): 

   Usage:
     form -l Normalized_Hamiltonian6.frm
   Output:
     Normalized_Hamiltonian6.log  -    


II. "Symbolic". The normalisation in eta,zeta variables:

  1) The normal form for the Pendulum Hamiltonian 
     using the explicit formula for Deprit generator [1].

   Usage:
     form -l pendulum_explicit.frm
   Output:
     pendulum_explicit.log	

  2) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the explicit formula for Deprit generator [1].

   Usage:
     form -l Toda2d_explicit.frm
   Output:
     Toda2d_explicit.log

  3) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Deprit triangular algorithm [3].

   Usage:
     form -l Toda2d_Deprit3.frm
   Output:
     Toda2d_Deprit3.log

  4) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Dragt-Finn method [4].

   Usage:
     form -l Toda2d_DragtFinn.frm
   Output:
     Toda2d_DragtFinn.log

  5) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Henrard method [5]

   Usage:
     form -l Toda2d_Henrard.frm
   Output:
     Toda2d_Henrard.log

  6) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Hori-Mersman algorithm [6].

   Usage:
     form -l Toda2d_HoriMersman.frm
   Output:
     Toda2d_HoriMersman.log

  7) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Gustavson algorithm [7].

   Usage:
     form -l Toda2d_Gustavson.frm
   Output:
     Toda2d_Gustavson.log

 8) Comparison of the normal forms, generators, and the Gustavson integrals
     for the Deprit triangular and the explicit algorithms:

   Usage:
     form -l Deprit3_vs_explicit.frm
   Output:
     Deprit3_vs_explicit.log


  9) Comparison of the explicit [1], Deprit [3], Dragt-Finn [4], Henrard [5],
     Hori [6] and Gustavson [7] algorihms
     for the Toda 2D system

   Usage:
     form -l Toda2d_comparison.frm
   Output:
     Toda2d_comparison.log

  10) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Dragt-Finn-Lunter algorithm

   Usage:
     form -l Toda2d_DragtFinnlunter.frm
   Output:
     Toda2d_DragtFinnlunter.log

III. "Symbolicpq". The normalisation in p,q variables: 

  1) The normal form for the Pendulum Hamiltonian 
     using the explicit formula for Deprit generator
     in the p,q variables.

   Usage:
     form -l pendulum_explicitpq.frm
   Output:
     pendulum_explicitpq.log	

  2) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the explicit formula for Deprit generator 
     in the p,q variables.

   Usage:
     form -l Toda2d_explicitpq.frm
   Output:
     Toda2d_explicitpq.log

  3) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Deprit triangular algorithm
     in the p,q variables.

   Usage:
     form -l Toda2d_Deprit3pq.frm
   Output:
     Toda2d_Deprit3pq.log

  4) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Dragt-Finn algorithm
     in the p,q variables.

   Usage:
     form -l Toda2d_DragtFinnpq.frm
   Output:
     Toda2d_DragtFinnpq.log

  5) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Henrard algorithm
     in the p,q variables.

   Usage:
     form -l Toda2d_Henrardpq.frm
   Output:
     Toda2d_Henrardpq.log

  6) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Hori-Mersman algorithm
     in the p,q variables.

   Usage:
     form -l Toda2d_HoriMersmanpq.frm
   Output:
     Toda2d_HoriMersmanpq.log

  7) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Gustavson algorithm
     in the p,q variables.

   Usage:
     form -l Toda2d_Gustavsonpq.frm
   Output:
     Toda2d_Gustavsonpq.log

  8) Comparison of the normal forms, generators, and the Gustavson integrals
     for the Deprit triangular and the explicit algorithms
     in the p,q variables:

   Usage:
     form -l Deprit3_vs_explicit.frm
   Output:
     Deprit3_vs_explicit.log


  9) Comparison of the explicit [1], Deprit [3], Dragt-Finn [4], Henrard [5],
     Hori [6] and Gustavson [7] algorihms
     for the Toda 2D system
     in the p,q variables:

   Usage:
     form -l Toda2d_comparison.frm
   Output:
     Toda2d_comparison.log

  10) The normal form and the Gustavson integral
     for the Toda 2D Hamiltonian 
     using the Dragt-Finn-Lunter algorithm
   Usage:
     form -l Toda2d_DragtFinnlunterpq.frm
   Output:
     Toda2d_DragtFinnlunterpq.log

----------------------------------------
References:
1. Nikolaev A., arXiv:1612.05207 [math.DS]
2. New features of FORM - Vermaseren, J.A.M. math-ph/0010025 
3. Deprit A (1969), Celestial Mech 1:12–30
4. Dragt AJ, Finn JM (1976), J. Math. Phys. 17:2215–2227
   Algorithm by Koseleff PV (1994), Cel. Mech. 58:17–36
5. Henrard, J., in "Recent advances in dynamical astronomy" ed. Szebehely, p. 248, 	
   Algorithm by Koseleff PV (1994), Cel. Mech. 58:17–36
6. Hori G (1966), Publ Astron Soc Japan 18:287
   Algorithm by Mersman W.A., Celest. Mech. 3 (1970) 81–89.
7. Gustavson FG (1966), Astron J 71:670
8. Broer, H.; Hoveijn, I.; Lunter, G.; Vegter, G.,"Birkhoff normalization"

