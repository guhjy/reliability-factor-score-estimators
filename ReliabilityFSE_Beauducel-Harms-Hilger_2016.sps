* ' This function computes and returns reliability estimates for three commonly used
  ' Factor Score Estimators in Factor Analyses,
  ' 
  ' Explanations of the algebraic formulas are presented in the manuscript 
  '
  ' André Beauducel (\email{beauducel@uni-bonn.de})
  ' Christopher Harms (\email{christopher.harms@uni-bonn.de})
  ' Norbert Hilger (\email{nhilger@uni-bonn.de})
/*.
MATRIX.
    * Users may enter their respective numbers into the loading matrix:.
    compute L={
        0.50,-0.10, 0.10;
        0.50, 0.10, 0.10;
        0.50, 0.10,-0.10;
       -0.10, 0.50, 0.15;
        0.15, 0.50, 0.10;
       -0.15, 0.50, 0.10;
        0.10, 0.10, 0.60;
        0.10,-0.10, 0.60;
        0.10, 0.10, 0.60
    }.
    print L/format=F5.2.

    * Enter respective numbers into factor inter-correlations.
    compute Phi={
        1.00, 0.30, 0.20;
        0.30, 1.00, 0.10;
        0.20, 0.10, 1.00
    }.
    print Phi/format=F5.2.
     
    * Reproduce the observed covariances from the parameters of the factor model.
    compute Sig=L*Phi*T(L).
    compute Sig=Sig-Mdiag(diag(Sig))+ident(nrow(L),nrow(L)).

    * Spezifität/Uniqueness/Error der Items berechnen.
    compute Psi=Mdiag(diag(Sig-L*Phi*T(L)))&**0.5.

    * Equation 6.
    compute Rtt_r = INV( Mdiag(diag( Phi*T(L)*INV(Sig)*L*Phi )) )&**0.5 *
    Mdiag(diag(Phi*T(L)*INV(Sig)*L*Phi*T(L)*INV(Sig)*L*Phi)) *
    INV(Mdiag(diag(Phi*T(L)*INV(Sig)*L*Phi)))&**0.5 .

    * Equation 11.
    compute Rtt_b=INV( Mdiag(diag(INV(T(L)*INV(Psi)&**2*L) + Phi)) ).

    * Equation 12.
    CALL svd(phi, QQ, eig, QQQ).
    compute N=QQ*abs(eig)&**0.5.
    compute help=T(N)*T(L)*INV(Psi)&**2*Sig*INV(Psi)&**2*L*N.
    CALL svd(help, QQ, eig, QQQ).
    compute help12=QQ*((eig)&**0.5)*T(QQ).

    compute Rtt_m=Mdiag(diag(
    INV(help12)*T(N)*T(L)*INV(Psi)&**2*L*Phi*T(L)*INV(Psi)&**2*L*N*INV(help12)
    )).
    print/Title "Reliabilities for Regression factor score estimators:".
    print {T(diag(rtt_r))}/Format=F6.3.
    print/Title "Reliabilities for Bartlett factor score estimators:".
    print {T(diag(rtt_b))}/Format=F6.3.
    print/Title "Reliabilities for McDonald factor score estimators:".
    print {T(diag(rtt_m))}/Format=F6.3.
END MATRIX.
