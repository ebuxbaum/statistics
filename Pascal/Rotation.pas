UNIT Rotation;

{ Rotation of loading matrices by gradient projection according to
  C.A. Bernaards & R.I. Jennrich: Gradient Projection Algorithms and Software
  for Arbitrary Rotation Criteria in Factor Analysis,
  Educ. Psychol. Meas. 65:5 (2005) 676-696, doi:10.1177/0013164404272507.
  Some code in R and SAS is deposited on http://www.stat.ucla.edu/research/gpa/ }

INTERFACE

USES MathFunc, Vector, Matrix, SingularValue;

CONST
  Convergence = -5; // log OF convergence criterium considered success

TYPE
  RotType = (Varimax, Quartimax, Equamax, Parsimax, Parsimony, Quartimin,
    Biquartimin, Covarimin, Entropy, Tandem1, Tandem2);

PROCEDURE GradProjAlgOrth(VAR Loading, Rotation: MatrixTyp; RotAlg: RotType);
{ Gradient Projection Algorithm for orthogonal rotation. Rotation algorithm used
  is selected by RotAlg.  }

PROCEDURE CalculateSAL(CONST Loading: MatrixTyp; VAR SAL: VectorTyp);
{ Sorted absolute loadings (SAL) of a loading matrix. }


IMPLEMENTATION


PROCEDURE FVarimax(CONST ActLoading: MatrixTyp; VAR Gradient: MatrixTyp;
  VAR Quality: double);

VAR
  cm, L2, QL, Zwischen: MatrixTyp;
  n: WORD;

BEGIN
  HadamardSchurProduct(ActLoading, ActLoading, L2);            // L^2
  n := MatrixRows(L2);
  CreateMatrix(cm, n, n, 1.0 / n);
  MatrixInnerProduct(cm, L2, Zwischen);
  DestroyMatrix(cm);
  NegativeMatrix(Zwischen);
  MatrixAdd(L2, Zwischen, QL);
  DestroyMatrix(L2);
  DestroyMatrix(Zwischen);
  Quality := -FrobeniusNorm(QL) / 4;
  CopyMatrix(ActLoading, Zwischen);
  NegativeMatrix(Zwischen);
  HadamardSchurProduct(Zwischen, QL, Gradient);
  DestroyMatrix(Zwischen);
  DestroyMatrix(QL);
END;


PROCEDURE FOblimin(CONST ActLoading: MatrixTyp; gamma: double;
  VAR Gradient: MatrixTyp; VAR Quality: double);
{ gamma determines type:
  0.0 quartimin (most oblique)
  0.5 Bi-quartimin
  1.0 covarimin }

VAR
  Rows, Columns, i, j: WORD;
  Ident, N, M, Zwischen, Zwischen2, L2: MatrixTyp;

BEGIN
  Rows := MatrixRows(ActLoading);
  Columns := MatrixColumns(ActLoading);
  CreateIdentityMatrix(Ident, Columns);
  NegativeMatrix(Ident);
  CreateMatrix(Zwischen, Columns, Columns, 1.0);
  MatrixAdd(Zwischen, Ident, N);
  DestroyMatrix(Zwischen);
  DestroyMatrix(Ident);
  CreateIdentityMatrix(Ident, Rows);
  FOR i := 1 TO Rows DO
    FOR j := 1 TO Rows DO
      SetMatrixElement(Ident, i, j, GetMatrixElement(Ident, i, j) - gamma);
  CreateMatrix(M, Rows, Rows, 1 / Rows);
  HadamardSchurProduct(ActLoading, ActLoading, L2);            // L^2
  MatrixInnerProduct(L2, N, Zwischen);
  HadamardSchurProduct(Ident, M, Zwischen2);
  DestroyMatrix(Ident);
  DestroyMatrix(M);
  DestroyMatrix(N);
  MatrixInnerProduct(Zwischen2, Zwischen, N);
  // [(I(p)-gamma \cdot M) * L2 * N]
  DestroyMatrix(Zwischen);
  DestroyMatrix(Zwischen2);
  Quality := FrobeniusSkalarProduct(L2, N) / 4;
  // <L2, [(I(p)-gamma \cdot M) * L2 * N]>
  HadamardSchurProduct(ActLoading, N, Gradient);
  // L \cdot [(I(p)-gamma \cdot M) * L2 * N]
  DestroyMatrix(N);
  DestroyMatrix(L2);
END;


PROCEDURE FEntropy(CONST ActLoading: MatrixTyp; VAR Gradient: MatrixTyp;
  VAR Quality: double);
{ Bernaards C.A., Jennrich R.I.: Gradient Projection Algorithms and Software for
  Arbitrary Rotation Criteria in Factor Analysis, Educ. Psychol. Meas. 65 (2005)
  676–696. Rotation very good for orthogonal, but poor for oblique rotations. }

VAR
  L2, lnL2, Zwischen: MatrixTyp;
  i, j: WORD;

BEGIN
  HadamardSchurProduct(ActLoading, ActLoading, L2);            // L^2
  CopyMatrix(L2, lnL2);
  FOR i := 1 TO MatrixRows(lnL2) DO
    FOR j := 1 TO MatrixColumns(lnL2) DO
      SetMatrixElement(lnL2, i, j, Ln(GetMatrixElement(lnL2, i, j)));
  Quality := -FrobeniusSkalarProduct(L2, lnL2) / 2;          // -<L2, log(L2)>/2
  DestroyMatrix(L2);
  HadamardSchurProduct(ActLoading, lnL2, Zwischen);
  DestroyMatrix(lnL2);
  MatrixAdd(Zwischen, ActLoading, Gradient);            // L \cdot Ln(L2) - L
  NegativeMatrix(Gradient);
  DestroyMatrix(Zwischen);
END;


PROCEDURE FTandem1(CONST ActLoading: MatrixTyp; VAR Gradient: MatrixTyp;
  VAR Quality: double);
{ A.L. Comrey: Tandem criteria for analytic rotation in factor analysis,
  Psychometrika 32(1967) 143-54. Rotation with criterium I is employed to
  determine the number of factors to be retained, the resulting loadings
  are then rotated by criterium II. Works with orthogonal, but not oblique,
  rotation. }

VAR
  TL, TL2, LL, LL2, L2, gq1, gq2, Zwischen, Zwischen1: MatrixTyp;

BEGIN
  MatrixTranspose(ActLoading, TL);
  MatrixInnerProduct(ActLoading, TL, LL);
  DestroyMatrix(TL);
  HadamardSchurProduct(LL, LL, LL2);
  HadamardSchurProduct(ActLoading, ActLoading, L2);
  MatrixTranspose(L2, TL2);
  MatrixInnerProduct(LL2, L2, Zwischen);
  DestroyMatrix(LL2);
  MatrixInnerProduct(TL2, Zwischen, Zwischen1);
  Quality := -MatrixTrace(Zwischen1);
  DestroyMatrix(Zwischen1);
  HadamardSchurProduct(ActLoading, Zwischen, gq1);
  SkalarMultiplikation(gq1, 4);
  DestroyMatrix(Zwischen);
  MatrixInnerProduct(L2, TL2, Zwischen);
  DestroyMatrix(TL2);
  DestroyMatrix(L2);
  HadamardSchurProduct(LL, Zwischen, Zwischen1);
  DestroyMatrix(LL);
  DestroyMatrix(Zwischen);
  MatrixInnerProduct(Zwischen1, ActLoading, gq2);
  SkalarMultiplikation(gq2, 4);
  DestroyMatrix(Zwischen1);
  NegativeMatrix(gq1);
  NegativeMatrix(gq2);
  MatrixAdd(gq1, gq2, Gradient);
  DestroyMatrix(gq1);
  DestroyMatrix(gq2);
END;


PROCEDURE FTandem2(CONST ActLoading: MatrixTyp; VAR Gradient: MatrixTyp;
  VAR Quality: double);

VAR
  TL, TL2, LL, LL2, L2, U, TU, UU, gq1, gq2,
  Zwischen, Zwischen1, Copy                       : MatrixTyp;
  p: WORD;

BEGIN
  p := MatrixRows(ActLoading);
  MatrixTranspose(ActLoading, TL);
  MatrixInnerProduct(ActLoading, TL, LL);
  DestroyMatrix(TL);
  HadamardSchurProduct(LL, LL, LL2);
  NegativeMatrix(LL2);
  DestroyMatrix(LL);
  CopyMatrix(ActLoading, Copy);
  HadamardSchurProduct(Copy, Copy, L2);
  CreateMatrix(U, p, 1, 1.0);
  CreateMatrix(TU, 1, p, 1.0);
  MatrixInnerProduct(U, TU, UU);
  DestroyMatrix(U);
  DestroyMatrix(TU);
  MatrixAdd(UU, LL2, Zwischen);
  DestroyMatrix(UU);
  MatrixInnerProduct(Zwischen, L2, Zwischen1);
  Quality := FrobeniusSkalarProduct(L2, Zwischen1);
  SkalarMultiplikation(Copy, 4);
  HadamardSchurProduct(Copy, Zwischen1, gq1);   {******}
  DestroyMatrix(Copy);
  DestroyMatrix(Zwischen);
  DestroyMatrix(Zwischen1);
  MatrixTranspose(L2, TL2);
  MatrixInnerProduct(L2, TL2, Zwischen);
  DestroyMatrix(L2);
  DestroyMatrix(TL2);
  HadamardSchurProduct(LL2, Zwischen, Zwischen1);
  DestroyMatrix(Zwischen);
  DestroyMatrix(LL2);
  MatrixInnerProduct(Zwischen1, Copy, gq2);
  DestroyMatrix(Zwischen1);
  SkalarMultiplikation(gq2, -4);
  MatrixAdd(gq1, gq2, Gradient);
  DestroyMatrix(gq1);
  DestroyMatrix(gq2);
END;


PROCEDURE FQuartimax(CONST ActLoading: MatrixTyp; VAR Gradient: MatrixTyp;
  VAR Quality: double);

VAR
  SqrLoading: MatrixTyp;

BEGIN
  HadamardSchurProduct(ActLoading, ActLoading, SqrLoading);            // L^2
  Quality := -0.25 * Sqr(FrobeniusNorm(SqrLoading));                   // -1/4 <L^2>^2
  HadamardSchurProduct(SqrLoading, ActLoading, Gradient);              // L^3
  NegativeMatrix(Gradient);                                            // -L^3
  DestroyMatrix(SqrLoading);
END;


PROCEDURE FCrawfordFerguson(CONST ActLoading: MatrixTyp; kappa: double;
  VAR Gradient: MatrixTyp; VAR Quality: double);
{ Crawford C.B., Ferguson G.A.: A General Rotation Criterion and its Use in
  Orthogonal Rotation, Psychometrika 35 (1970) 321–332.
  kappa determines rotation type (kappa =
  0 Quartimax,
  1/p  Varimax
  q/(2*p) Equamax.
  (q-1)/(p+q-2) Parsimax
  1  Factor parsimony. }

VAR
  Rows, Columns: WORD;
  Ident, N, M, Zwischen, L2, g1, g2: MatrixTyp;
  f1, f2: double;

BEGIN
  Rows := MatrixRows(ActLoading);
  Columns := MatrixColumns(ActLoading);
  CreateIdentityMatrix(Ident, Columns);
  NegativeMatrix(Ident);
  CreateMatrix(Zwischen, Columns, Columns, 1.0);
  MatrixAdd(Zwischen, Ident, N);
  DestroyMatrix(Ident);
  DestroyMatrix(Zwischen);
  CreateIdentityMatrix(Ident, Rows);
  NegativeMatrix(Ident);
  CreateMatrix(Zwischen, Rows, Rows, 1.0);
  MatrixAdd(Zwischen, Ident, M);
  DestroyMatrix(Ident);
  DestroyMatrix(Zwischen);
  HadamardSchurProduct(ActLoading, ActLoading, L2);            // L^2
  MatrixInnerProduct(L2, N, Zwischen);
  f1 := (1 - kappa) * FrobeniusSkalarProduct(L2, Zwischen) / 4;
  // (1-kappa) <L2, (L2*N)> /4
  HadamardSchurProduct(ActLoading, Zwischen, g1);             // L \cdot (L2 * N)
  DestroyMatrix(Zwischen);
  DestroyMatrix(N);
  MatrixInnerProduct(M, L2, Zwischen);
  f2 := kappa * FrobeniusSkalarProduct(L2, Zwischen) / 4;      // kappa * <L2, (M*L2)> /4
  HadamardSchurProduct(ActLoading, Zwischen, g2);             // L \cdot (M * L2)
  DestroyMatrix(Zwischen);
  DestroyMatrix(M);
  Quality := f1 + f2;
  SkalarMultiplikation(g1, 1 - kappa);
  SkalarMultiplikation(g2, kappa);
  MatrixAdd(g1, g2, Gradient);
  DestroyMatrix(L2);
  DestroyMatrix(g1);
  DestroyMatrix(g2);
END;


PROCEDURE FBentler(CONST ActLoading: MatrixTyp; VAR Gradient: MatrixTyp;
  // funktioniert nicht
  VAR Quality: double);
{ Bentler's invariant pattern simplicity, useful for oblique and orthogonal
  rotation. P.M. Bentler: Factor simplicity index and transformations,
  Psychometrika 42 (1977) 277-295 }

VAR
  L2, TL2, M, D, Zwischen, Zwischen1: MatrixTyp;
  Det1, Det2: double;

BEGIN
  HadamardSchurProduct(ActLoading, ActLoading, L2);
  MatrixTranspose(L2, TL2);
  MatrixInnerProduct(TL2, L2, M);
  DestroyMatrix(TL2);
  CopyMatrix(M, D);
  Det1 := Determinante(D);
  Diag(D);
  Det2 := Determinante(D);
  Quality := -Ln(Det1 / Det2);
  InverseMatrix(M);
  InverseMatrix(D);
  MatrixInnerProduct(M, D, Zwischen);
  DestroyMatrix(M);
  DestroyMatrix(D);
  MatrixInnerProduct(L2, Zwischen, Zwischen1);
  DestroyMatrix(L2);
  DestroyMatrix(Zwischen);
  CopyMatrix(Zwischen1, Zwischen);
  SkalarMultiplikation(Zwischen, -4);
  HadamardSchurProduct(Zwischen, Zwischen1, Gradient);
  DestroyMatrix(Zwischen);
  DestroyMatrix(Zwischen1);
END;


PROCEDURE GradProjAlgOrth(VAR Loading, Rotation: MatrixTyp; RotAlg: RotType);

VAR
  alpha, s, s1, Q, Qold: double;
  Columns, Rows, j, Iter: WORD;
  NewLoading, LoadingTrans, RotationTrans, GradQ, Grad, GradP,
  Manifold, ManifoldTrans, Zwischen, Skm2, X, V, VT: MatrixTyp;
  Delta: VectorTyp;

BEGIN
  alpha := 1.0;
  Columns := MatrixColumns(Loading);
  Rows := MatrixRows(Loading);
  Writeln;
  Writeln('Rotation:');
  Writeln('Iter Q      s           Log10(s)  alpha ');
  MatrixInnerProduct(Loading, Rotation, NewLoading);
  CreateMatrix(GradQ, Columns, Columns, 0.0);
  CASE RotAlg OF
    Varimax: FCrawfordFerguson(NewLoading, 1 / Rows, GradQ, Q);
    Quartimax: FCrawfordFerguson(NewLoading, 0.0, GradQ, Q);
    Equamax: FCrawfordFerguson(NewLoading, Columns / (2 * Rows), GradQ, Q);
    Parsimax: FCrawfordFerguson(NewLoading, Pred(Columns) /
        (Rows + Columns - 2), GradQ, Q);
    Parsimony: FCrawfordFerguson(NewLoading, 1.0, GradQ, Q);
    Quartimin: FOblimin(NewLoading, 0.0, GradQ, Q);
    Biquartimin: FOblimin(NewLoading, 0.5, GradQ, Q);
    Covarimin: FOblimin(NewLoading, 1.0, GradQ, Q);
    Entropy: FEntropy(NewLoading, GradQ, Q);
    Tandem1: FTandem1(NewLoading, GradQ, Q);
    Tandem2: FTandem2(NewLoading, GradQ, Q);
    //    Bentler    : FBentler(NewLoading, GradQ, Q);
  END;
  Qold := Q;
  MatrixTranspose(Loading, LoadingTrans);
  MatrixInnerProduct(LoadingTrans, GradQ, Grad);    // calculate gradient
  Iter := 0;
  REPEAT
    MatrixTranspose(Rotation, RotationTrans);
    MatrixInnerProduct(RotationTrans, Grad, Manifold);        // M = T' * G
    MatrixTranspose(Manifold, ManifoldTrans);
    MatrixAdd(Manifold, ManifoldTrans, Skm2);
    SkalarMultiplikation(Skm2, 0.5);                      // S = (M+M')/2
    MatrixInnerProduct(Rotation, Skm2, Zwischen);
    DestroyMatrix(Skm2);
    NegativeMatrix(Zwischen);
    MatrixAdd(Grad, Zwischen, GradP);                     // Gp = G - T*S
    DestroyMatrix(Zwischen);
    s := FrobeniusNorm(GradP);
    s1 := log(s, 10);
    Writeln(Iter: 3, '  ', Q: 2: 3, ' ', s: 2: 8, ' ', s1: 2: 3, '     ', alpha: 2: 2);
    IF (s1 > Convergence)
      THEN
        BEGIN
          alpha := 2 * alpha;
          j := 0;
          REPEAT                                             // ensure that Q IS decreased
            DestroyMatrix(NewLoading);
            SkalarMultiplikation(GradP, -alpha);
            MatrixAdd(Rotation, GradP, X);          // X = T - a1*Gp
            SVDCmp(X, Delta, V);                             // X now contains U
            MatrixTranspose(V, Vt);
            DestroyMatrix(V);
            MatrixInnerProduct(X, Vt, Zwischen);               // Tt = U T'
            MatrixInnerProduct(Loading, Zwischen, NewLoading); // L = A * Tt
            CASE RotAlg OF
              Quartimax: FCrawfordFerguson(NewLoading, 0.0, GradQ, Q);
              Varimax: FCrawfordFerguson(NewLoading, 1 / Rows, GradQ, Q);
              Equamax: FCrawfordFerguson(NewLoading, Columns / (2 * Rows), GradQ, Q);
              Parsimax: FCrawfordFerguson(NewLoading, pred(Columns) /
                  (Rows + Columns - 2), GradQ, Q);
              Parsimony: FCrawfordFerguson(NewLoading, 1.0, GradQ, Q);
              Quartimin: FOblimin(NewLoading, 0.0, GradQ, Q);
              Biquartimin: FOblimin(NewLoading, 0.5, GradQ, Q);
              Covarimin: FOblimin(NewLoading, 1.0, GradQ, Q);
              Entropy: FEntropy(NewLoading, GradQ, Q);
              Tandem1: FTandem1(NewLoading, GradQ, Q);
              Tandem2: FTandem2(NewLoading, GradQ, Q);
              //              Bentler    : FBentler(NewLoading, GradQ, Q);
            END;
            IF (Q > (Qold - 0.5 * s * s * alpha))
              THEN alpha := alpha / 2;
            Inc(j);
            DestroyMatrix(X);
            DestroyMatrix(Vt);
            DestroyVector(Delta);
          UNTIL ((j = 10) OR (Q < (Qold - 0.5 * s * s * alpha)));
          CopyMatrix(Zwischen, Rotation);
          DestroyMatrix(Zwischen);
        END;
    Qold := Q;
    MatrixInnerProduct(LoadingTrans, GradQ, Grad);
    Inc(Iter);
    DestroyMatrix(RotationTrans);
    DestroyMatrix(Manifold);
    DestroyMatrix(ManifoldTrans);
  UNTIL ((s1 < Convergence) OR (Iter > MaxIter));
  DestroyMatrix(NewLoading);
  MatrixInnerProduct(Loading, Rotation, NewLoading);
  DestroyMatrix(Loading);
  CopyMatrix(NewLoading, Loading);
  DestroyMatrix(NewLoading);
  DestroyMatrix(LoadingTrans);
  DestroyMatrix(GradQ);
  DestroyMatrix(GradP);
  DestroyMatrix(Grad);
END;


PROCEDURE CalculateSAL(CONST Loading: MatrixTyp; VAR SAL: VectorTyp);

VAR
  i, j, Rows, Columns, Length: word;

BEGIN
  Rows := MatrixRows(Loading);
  Columns := MatrixColumns(Loading);
  Length := Rows * Columns;
  CreateVector(SAL, Length, 0.0);
  FOR i := 1 TO Rows DO
    FOR j := 1 TO Columns DO
      SetVectorElement(SAL, i * j, abs(GetMatrixElement(Loading, i, j)));
  Shellsort(SAL);
END;

END.



