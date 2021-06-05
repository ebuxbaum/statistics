UNIT PerformShrinkage;

INTERFACE

USES Math, MathFunc, Vector, Matrix, EigenValues;

PROCEDURE ErzeugeLambda(CONST Eigenvalues: VectorTyp; VAR Lambda: MatrixTyp);
{ Turn vector of eigenvalues into a diagonal matrix }

PROCEDURE Average(VAR R: MatrixTyp; omega: double);
{ Shrink the correlation matrix R by calculating the weighted average between R
  and the row means of R for each i and j. omega (0..1) is the relative weight
  of of the average, the relative weight of R is (1-omega) }

PROCEDURE Shrinkage(VAR R: MatrixTyp; omega: double);
{ Shrink the correlation matrix R by calculating the weighted average between R
  and the identity matrix. omega (0..1) is the relative weight of I, the relative
  weight of R is (1-omega) }

PROCEDURE NormaliseCorrelations(VAR R: MatrixTyp);
{ make an indefinite correlation matrix positive definite for eigenanalysis,
  by setting all eigenvalues < 0 to zero and re-calculating R = E Lambda E^{-1}. }

IMPLEMENTATION

PROCEDURE ErzeugeLambda(CONST Eigenvalues: VectorTyp; VAR Lambda: MatrixTyp);

VAR
  p, i, j: WORD;

BEGIN
  p := VectorLength(Eigenvalues);
  CreateIdentityMatrix(Lambda, p);
  FOR j := 1 TO p DO
    SetMatrixElement(Lambda, j, j, GetVectorElement(Eigenvalues, j));
END;


PROCEDURE Average(VAR R: MatrixTyp; omega: double);

VAR
  p, i, j: WORD;
  Averages: VectorTyp;
  Sum, w: double;

BEGIN
  IF (omega < 0) OR (omega > 1)
    THEN
      BEGIN
        Writeln('Shrinkage: omega not in [0..1]');
        HALT;
      END;
  p := MatrixRows(R);
  CreateVector(Averages, p, 0);                // calculate row-averages OF R
  FOR i := 1 TO p DO
    BEGIN
      Sum := 0;
      FOR j := 1 TO p DO
        IF (i = j)
          THEN  // ignore correlation WITH self
          ELSE  Sum := Sum + GetMatrixElement(R, i, j);
      SetVectorElement(Averages, i, Sum / Pred(p));
    END;
  w := 1 - omega;
  omega := omega / 2;                               // weight FOR each r_i AND r_j
  FOR i := 1 TO p DO                              // calculate NEW elements OF R
    FOR j := Succ(i) TO p DO
      BEGIN
        Sum := (W * GetMatrixElement(R, i, j) + omega *
          GetVectorElement(Averages, i) + omega * GetVectorElement(Averages, j));
        SetMatrixElement(R, i, j, Sum);
        SetMatrixElement(R, j, i, Sum);
      END;
  DestroyVector(Averages);
END;


PROCEDURE Shrinkage(VAR R: MatrixTyp; omega: double);

VAR
  p, i, j: WORD;
  Sum, w: double;

BEGIN
  IF (omega < 0) OR (omega > 1)
    THEN
      BEGIN
        Writeln('Shrinkage: omega not in [0..1]');
        HALT;
      END;
  w := 1 - omega;
  p := MatrixRows(R);
  FOR i := 1 TO p DO                              // calculate NEW elements OF R
    FOR j := Succ(i) TO p DO
      BEGIN
        Sum := (W * GetMatrixElement(R, i, j)); // non-diagonal elements OF I are 0
        SetMatrixElement(R, i, j, Sum);        // AND can be ignored
        SetMatrixElement(R, j, i, Sum);
      END;
END;


PROCEDURE NormaliseCorrelations(VAR R: MatrixTyp);

VAR
  Eigenvalues: VectorTyp;
  Eigenvectors, Hilfs, Lambda, RecipEV: MatrixTyp;
  i, iter, n, NoNegs, j: WORD;
  Rneu: double;

BEGIN
  n := MatrixRows(R);
  Writeln('Initial calculation of eigenvalues');
  iter := 1;
  REPEAT
    INC(iter);
    i := Jacobi(R, Eigenvalues, Eigenvectors, j);  // eigenanalysis
    Write(iter: 2, ': ', j: 3, ' iterations, Result = ', i: 1);
    CASE i OF
      0    : Write(' ok             ');
      5    : Write(' no convergence ');
      ELSE   Write(' error          ');
    END;
    NoNegs := 0;
    FOR i := 1 TO n DO
    BEGIN
      IF (GetVectorElement(Eigenvalues, i) < 0)
        THEN
          BEGIN
            SetVectorElement(Eigenvalues, i, 0);
            INC(NoNegs);
          END;
    END;
    Writeln(NoNegs: 3);
    DestroyMatrix(R);
    ErzeugeLambda(Eigenvalues, Lambda);
    CopyMatrix(EigenVectors, RecipEV);
    InverseMatrix(RecipEV);
    MatrixInnerProduct(EigenVectors, Lambda, Hilfs);    // calculate R = E Lambda E^{-1}
    MatrixInnerProduct(Hilfs, RecipEV, R);
    DestroyMatrix(EigenVectors);
    DestroyMatrix(RecipEV);
    DestroyMatrix(Hilfs);
    DestroyMatrix(Lambda);
    DestroyVector(EigenValues);
    FOR i := 1 TO n DO
      BEGIN
        FOR j := 1 TO n DO
          BEGIN
            Rneu := GetMatrixElement(R, i, j);
            IF (Abs(Rneu) > 1.0)
              THEN SetMatrixElement(R, i, j, sign(Rneu));  // max + OR - 1
          END;
        SetMatrixElement(R, i, i, 1.0);                 // diagonal elements 1
      END;
  UNTIL (iter >= 100) OR (NoNegs = 0);
END;

END.

