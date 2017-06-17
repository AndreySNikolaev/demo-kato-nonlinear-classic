#procedure DragtFinnLunterNormalization(PREFIX)
*==========================================================================
*   This file is part of the demos for 
*    "Generalization of the explicit expression for the Deprit generator 
*       to  Hamiltonians nonlinearly dependent on small parameter"
*     Andrey Nikolaev, ICPT, RDTeX, Russia, 2016, Andrey.Nikolaev@rdtex.ru
*==========================================================================
*  Normalization using Dragt-Finn chain of exponents.
*    Algorithm by Lunter G.
*   Birkhoff normalization. / Broer, H.; Hoveijn, I.; Lunter, G.; Vegter, G.
*   EPRINTS-BOOK-TITLE. University of Groningen, Johann Bernoulli Institute 
*   for Mathematics and Computer Science, 2003.
*   http://www.rug.nl/research/portal/publications/birkhoff-normalization(6068d4bc-328a-4f2c-a6b1-79ed414b95fe).html
*==========================================================================
*   We use here the original Lunter notation
*     \tilde H=exp(alpha^n L_Fn )...exp(alpha L_F1) H
*  F=sum alpha^n Fn - pseudogenerator, Op - initial Hamiltonian.
*  PREFIX - for uniqueness. Names of all objects created inside this procedure 
*  will begin with PREFIX.
*==========================================================================
*      An ancillary variable that will be used for a summation:
Symbol 'PREFIX'A;
*--------------------
L 'PREFIX'K=Hetazeta;
B alpha;
.sort

*--------------------
*     	K_i <- H_i
*--------------------
Drop 'PREFIX'K;
#do i=0,'MAXORDER'
L 'PREFIX'K'i' = 'PREFIX'K[alpha^'i']+'PREFIX'A;
#enddo
.sort
#do i=0,'MAXORDER'
Hide 'PREFIX'K'i';
#enddo
.sort
*--------------------
*   Main loop:
*--------------------
#do i= 1,'MAXORDER'
L 'PREFIX'F'i' = 'PREFIX'K'i';
id 'PREFIX'A=0;
#call SH0()
.sort
Hide 'PREFIX'F'i';
.sort
#do j='MAXORDER'-'i',0,-1
*--------------------
*	F <- adKj (Fi)
*	Ki+j <- Ki+j + F
*	q <- 2; While qi + j <= k do the following:
*--------------------
L 'PREFIX'F='PREFIX'K'j';
id 'PREFIX'A=0;
.sort
Hide 'PREFIX'F;
#do q=1,{('MAXORDER'-'j')/'i'}
*--------------------
*	F <- adFi (F)/q
*--------------------
#call L('PREFIX'F'i','PREFIX'F)
Drop 'PREFIX'F;
.sort
Drop Pb'NDIM';
L 'PREFIX'F=Pb'NDIM'/'q';
.sort
Hide 'PREFIX'F;
*--------------------
*	Kqi+j <- Kqi+j + F
*--------------------
Unhide 'PREFIX'K{'q'*'i'+'j'};
id,once 'PREFIX'A='PREFIX'A+'PREFIX'F;
.sort
Hide 'PREFIX'K{'q'*'i'+'j'};
.sort
#enddo
Drop 'PREFIX'F;
.sort
#enddo
#enddo
.sort
on statistics;
*--------------------
*	Pseudo-generator. Gn is the coefficient at alpha^n
*--------------------
G 'PREFIX'G{'MAXORDER'}= 
#do n = 1,'MAXORDER'
 + alpha^'n'*'PREFIX'F'n'
#enddo
;
G 'PREFIX'Ht{'MAXORDER'}= 
#do n = 0,'MAXORDER'
 + alpha^'n'*'PREFIX'K'n'
#enddo
;
id 'PREFIX'A=0;
B alpha;
.store
#endprocedure
