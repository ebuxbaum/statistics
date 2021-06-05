UNIT EigenValues;
{ Literatur:
  G. Engelen-Müllges & F. Reuter: Formelsammlung zur Numerischen Mathematik mit
    Turbo-Pascal Programmen. 2. Aufl. Mannheim (Wissenschaftsverlag) 1987
  W.H. Press et. al.: Mumerical recipes in Pascal: The art of scientific computation
    Cambridge (Cambridge University Press) 1989
  G.W. Steward: Introduction to matrix computation. New York (Academic Press) 1973
  R. Sedgewick: Algorithms. Reading (Mass.) (Addison-Wesley) 1983
  L.V. Atkinson & P.J. Harley: An introduction to numerical methods with Pascal,
    London (Addison-Wesley) 1983
  Borland: Turbo-Pascal Mathe-Toolboox version 4.0, Borlant Int. Inc., 1986
}

INTERFACE

USES MathFunc, Vector, Matrix, SystemSolve;

CONST EigenError: StrArrayType = ('ok', 'dimension error', 'tolerance < 0',
    'MaxIter < 1', 'no convergence', 'matrix singular', 'matrix not square',
    'matrix not symmetrical', 'matrix not positive definite',
    'last two roots not real',
    'not enough memory', 'Null-matrix');

TYPE  IntVec = ARRAY[1..MaxIter] OF WORD;  // iterations FOR eigenvalue calculations


  { **************************** EigenValues ********************************* }

FUNCTION DominantEigenValue(CONST Mat: MatrixTyp; VAR EigenVector: VectorTyp;
    VAR EigenValue: double): BYTE;

FUNCTION NextEigenValue(CONST Mat: MatrixTyp; VAR EigenVector: VectorTyp;
    VAR EigenValue: double): BYTE;

FUNCTION Jacobi(VAR Mat : MatrixTyp; VAR Eigenvalues : VectorTyp;
    VAR Eigenvectors : MatrixTyp; VAR Iter : WORD) : BYTE;

PROCEDURE SortEigenValues(VAR EigenValues: VectorTyp; VAR EigenVectors: MatrixTyp);
{ Sorts eigenvalues largest first and puts eigenvectors in the same order }

IMPLEMENTATION

VAR Iter: WORD;
    CH  : CHAR;


FUNCTION Convergenz(OldApprox, NewApprox : VectorTyp): BOOLEAN;
{ wird true, wenn die Differenz aller Elemente in NewApprox und OldApprox
  kleiner als MaxError wird }

VAR Index : WORD;
    Found : BOOLEAN;

BEGIN
  Index := 0;
  Found := TRUE;
  WHILE Found AND (Index < VectorLength(OldApprox)) DO
    BEGIN
      INC(Index);
      IF Abs(GetVectorElement(OldApprox,Index) - GetVectorElement(NewApprox, Index)) > MaxError
        THEN found := FALSE;
    END;
  Result := found;
END;


FUNCTION TestDataAndInit(Mat: MatrixTyp; VAR EigenVector, NewApprox, OldApprox: VectorTyp;
  VAR EigenValue: double): BYTE;

VAR Error: BYTE;

BEGIN
  Error := 0;
  IF MatrixRows(Mat) < 1 THEN Error := 1;
  IF MaxIter < 1 THEN Error := 3;
  LoadConstant(EigenVector, 1);
  IF Error = 0
    THEN
      BEGIN
        Iter := 0;
        CopyVector(EigenVector, OldApprox);
        CopyVector(EigenVector, NewApprox);
      END;
  IF MatrixRows(Mat) = 1
    THEN
      BEGIN
        EigenValue := GetMatrixElement(Mat, 1, 1);
        SetVectorElement(EigenVector, 1, 1);
      END;
  Result := Error;
END; { TestDataAndInit }


FUNCTION DominantEigenValue(CONST Mat : MatrixTyp; VAR EigenVector : VectorTyp;
  VAR EigenValue : double) : BYTE;

TYPE Hilfsmatrix = ARRAY [0..2, 1..MaxVector] OF double;  { Spart Speicherplatz }

VAR Remainder, Index, n     : WORD;
    Error                   : BYTE;
    ApproxEigenValue, Denom : double;
    AitkenVector            : HilfsMatrix;
    Gefunden                : BOOLEAN;
    NewApprox, OldApprox    : VectorTyp;

BEGIN
  n := MatrixRows(Mat);
  CreateVector(NewApprox, n, 0.0);
  CreateVector(OldApprox, n, 0.0);
  IF VectorError
    THEN
      BEGIN
        VectorError := FALSE;
        Result := 10;
        EXIT;
      END;
  Error := TestDataAndInit(Mat, EigenVector, NewApprox, OldApprox, ApproxEigenValue);
  IF Error <> 0
    THEN
      BEGIN
        Result := Error;
        EXIT;
      END;
  IF n = 1
    THEN
      BEGIN
        Result := 0;
        EigenValue := ApproxEigenValue;
        EXIT;
      END;
  Gefunden := FALSE;
  FillChar(AitkenVector, SizeOf(AitkenVector), 0);
  ApproxEigenValue := FindLargest(OldApprox);
  DivConstant(OldApprox, ApproxEigenValue);
  WHILE (Iter < MaxIter) AND NOT Gefunden DO
    BEGIN
      INC(Iter);
      Remainder := Iter MOD 3;
      IF Remainder = 0
        THEN                    { Aitken's Beschleunigung }
          BEGIN
            FOR Index := 1 TO n DO
              SetVectorElement(OldApprox, Index, AitkenVector[0, Index]);
            FOR Index := 1 TO n DO
              BEGIN
                Denom := AitkenVector[2, Index] - 2 *
                AitkenVector[1, Index] + AitkenVector[0, Index];
                IF Abs(Denom) > MaxError
                  THEN
                    SetVectorElement(OldApprox, Index, AitkenVector[0, Index] -
                      Sqr(AitkenVector[1, Index] - AitkenVector[0, Index]) / Denom);
              END; { for }
          END; {if Remainder }
      MultMatrixVector(Mat, OldApprox, NewApprox);
      ApproxEigenValue := FindLargest(NewApprox);
      IF Abs(ApproxEigenValue) < MaxError
        THEN
          BEGIN
            ApproxEigenValue := 0;
            Gefunden := TRUE;
          END
        ELSE
          BEGIN
            DivConstant(NewApprox, ApproxEigenValue);
            Gefunden := Convergenz(OldApprox, NewApprox);
            CopyVector(NewApprox, OldApprox);
          END;
      FOR Index := 1 TO n DO
        AitkenVector[Remainder, Index] := GetVectorElement(NewApprox, Index);
    END; { while }
  EigenValue := ApproxEigenValue;
  CopyVector(OldApprox, EigenVector);
  IF Iter >= MaxIter THEN Error := 4;
  Result := Error;
  DestroyVector(NewApprox);
  DestroyVector(OldApprox);
END; { DominantEigenValue }


FUNCTION NextEigenValue(CONST Mat : MatrixTyp; VAR EigenVector : VectorTyp;
  VAR EigenValue : double) : BYTE;

VAR i, n                 : WORD;
    Error                : BYTE;
    ok, Gefunden         : BOOLEAN;
    mu                   : double;
    NewApprox, OldApprox : VectorTyp;
    A, Decomp, Permute   : MatrixTyp;

BEGIN
  n := MatrixRows(A);
  CopyMatrix(Mat, A);
  CreateMatrix(Decomp, n, n, 0.0);
  CreateMatrix(Permute, n, n, 0.0);
  CreateVector(NewApprox, n, 0.0);
  CreateVector(OldApprox, n, 0.0);
  IF VectorError OR MatrixError
    THEN
      BEGIN
        VectorError := FALSE;
        MatrixError := FALSE;
        Result := 8;
        EXIT;
      END;
  Error := TestDataAndInit(A, EigenVector, NewApprox, OldApprox, EigenValue);
  IF Error <> 0
    THEN
      BEGIN
        Result := Error;
        EXIT;
      END;
  IF n = 1
    THEN
      BEGIN
        Result := 0;
        EXIT;
      END;
  FOR i := 1 TO n DO
    BEGIN
      SetMatrixElement(A, i, i, GetMatrixElement(A, i, i) - EigenValue);
      SetVectorElement(OldApprox, i, 1);
    END;
  IF (LU_Decompose(A, Decomp, Permute) <> 0) OR
    (LU_Solve(Decomp, Permute, OldApprox, NewApprox) <> 0)
  { AusValueung kurzgeschlossen! }
    THEN
      BEGIN
        Result := 5;
        EXIT;
      END;
  REPEAT
    INC(Iter);
    mu := 0;
    FOR i := 1 TO n DO
      IF (Abs(GetVectorElement(NewApprox, i)) > Abs(mu))
        THEN mu := GetVectorElement(NewApprox, i);
    Gefunden := TRUE;
    DivConstant(NewApprox, mu);             { y normalisieren }
    Gefunden := Convergenz(NewApprox, OldApprox);
    CopyVector(NewApprox, OldApprox);
    IF NOT Gefunden
      THEN
        IF Iter >= MaxIter
          THEN Error := 4
          ELSE ok := LU_Solve(Decomp, Permute, OldApprox, NewApprox) = 0;
  UNTIL Gefunden OR (Error <> 0);
  IF Gefunden
    THEN Result := 0
    ELSE Result := Error;
  EigenValue := EigenValue + 1 / mu;
  CopyVector(NewApprox, EigenVector);
  DestroyVector(NewApprox);
  DestroyVector(OldApprox);
END; { NextEigenValue }


FUNCTION Jacobi(VAR Mat : MatrixTyp; VAR Eigenvalues : VectorTyp;
           VAR Eigenvectors : MatrixTyp; VAR Iter : WORD) : BYTE;

VAR Row, Column, Diag, Dimen          : WORD;
    SinTheta, CosTheta, SumSquareDiag : double;
    Done                              : BOOLEAN;

  FUNCTION TestData(Dimen: WORD; VAR Mat: MatrixTyp; MaxIter: WORD): BYTE;
    {- This procedure tests the input data for errors. -}

  BEGIN
    IF (Dimen < 1)
      THEN
        BEGIN
          Result := 1;
          EXIT;
        END;
    IF NOT (MatrixSquare(Mat))
      THEN
        BEGIN
          Result := 6;
          EXIT;
        END;
    IF MaxIter < 1
      THEN
        BEGIN
          Result := 3;
          EXIT;
        END;
    IF NOT (MatrixSymmetric(Mat))
      THEN
        BEGIN
          Result := 7;
          EXIT;
        END;
    Result := 0;
  END; { procedure TestData }


  PROCEDURE CalculateRotation(RowRow, RowCol, ColCol : double; VAR SinTheta, CosTheta: double);

  { This procedure calculates the sine and cosine of the angle Theta through which
    to rotate the matrix Mat. Given the tangent of 2-Theta, the tangent of Theta
    can be calculated with the quadratic formula.  The cosine and sine are easily
    calculable from the tangent. The rotation must be such that the Row, Column
    element is MaxError. RowRow is the Row,Row element; RowCol is the Row,Column element;
    ColCol is the Column,Column element  of Mat. }

  VAR TangentTwoTheta, TangentTheta, Dummy: double;

  BEGIN
    IF Abs(RowRow - ColCol) > MaxError
      THEN
        BEGIN
          TangentTwoTheta := (RowRow - ColCol) / (2 * RowCol);
          Dummy := Sqrt(Sqr(TangentTwoTheta) + 1);
          IF TangentTwoTheta < 0
            THEN TangentTheta := -TangentTwoTheta - Dummy
            ELSE TangentTheta := -TangentTwoTheta + Dummy;
          CosTheta := 1 / Sqrt(1 + Sqr(TangentTheta));
          SinTheta := CosTheta * TangentTheta;
        END
      ELSE
        BEGIN
          CosTheta := Sqrt(1 / 2);
          IF RowCol < 0
            THEN SinTheta := -Sqrt(1 / 2)
            ELSE SinTheta := Sqrt(1 / 2);
        END;
  END; { procedure CalculateRotation }

  PROCEDURE RotateMatrix(SinTheta, CosTheta: double; Row: WORD; Col: WORD;
    VAR Mat: MatrixTyp);
  { This procedure rotates the matrix Mat through an angle Theta.  The rotation
    matrix is the identity matrix execept for the Row,Row; Row,Col; Col,Col; and
    Col,Row elements. The rotation will make the Row,Col element of Mat to be
    MaxError. }

  VAR CosSqr, SinSqr, SinCos, MatRowRow, MatColCol,
      MatRowCol, MatRowIndex, MatColIndex          : double;
      Index                                        : WORD;

  BEGIN
    CosSqr := Sqr(CosTheta);
    SinSqr := Sqr(SinTheta);
    SinCos := SinTheta * CosTheta;
    MatRowRow := GetMatrixElement(Mat, Row, Row) * CosSqr + 2 *
      GetMatrixElement(Mat, Row, Col) * SinCos + GetMatrixElement(Mat, Col, Col) * SinSqr;
    MatColCol := GetMatrixElement(Mat, Row, Row) * SinSqr - 2 *
      GetMatrixElement(Mat, Row, Col) * SinCos + GetMatrixElement(Mat, Col, Col) * CosSqr;
    MatRowCol := (GetMatrixElement(Mat, Col, Col) -
      GetMatrixElement(Mat, Row, Row)) * SinCos + GetMatrixElement(Mat, Row, Col) * (CosSqr - SinSqr);
    FOR Index := 1 TO Dimen DO
      IF NOT (Index IN [Row, Col])
        THEN
          BEGIN
            MatRowIndex := GetMatrixElement(Mat, Row, Index) * CosTheta +
              GetMatrixElement(Mat, Col, Index) * SinTheta;
            MatColIndex := -GetMatrixElement(Mat, Row, Index) *
              SinTheta + GetMatrixElement(Mat, Col, Index) * CosTheta;
            SetMatrixElement(Mat, Row, Index, MatRowIndex);
            SetMatrixElement(Mat, Index, Row, MatRowIndex);
            SetMatrixElement(Mat, Col, Index, MatColIndex);
            SetMatrixElement(Mat, Index, Col, MatColIndex);
          END;
    SetMatrixElement(Mat, Row, Row, MatRowRow);
    SetMatrixElement(Mat, Col, Col, MatColCol);
    SetMatrixElement(Mat, Row, Col, MatRowCol);
    SetMatrixElement(Mat, Col, Row, MatRowCol);
  END; { procedure RotateMatrix }

  PROCEDURE RotateEigenvectors(SinTheta, CosTheta: double; Row, Col: WORD;
  VAR Eigenvectors: MatrixTyp);
{- This procedure rotates the Eigenvectors matrix through an angle Theta.  The
   rotation matrix is the identity matrix except for the Row,Row; Row,Col; Col,Col;
   and Col,Row elements.  The Eigenvectors matrix will be the product of all the
   rotation matrices which operate on Mat.  }

  VAR EigenvectorsRowIndex, EigenvectorsColIndex : double;
      Index                                      : WORD;

  BEGIN
    { Transform eigenvector matrix }
    FOR Index := 1 TO Dimen DO
      BEGIN
        EigenvectorsRowIndex :=  CosTheta * GetMatrixElement(Eigenvectors, Row, Index) +
          SinTheta * GetMatrixElement(Eigenvectors, Col, Index);
        EigenvectorsColIndex := -SinTheta * GetMatrixElement(Eigenvectors, Row, Index) +
          CosTheta * GetMatrixElement(Eigenvectors, Col, Index);
        SetMatrixElement(Eigenvectors, Row, Index, EigenvectorsRowIndex);
        SetMatrixElement(Eigenvectors, Col, Index, EigenvectorsColIndex);
      END;
  END; { procedure RotateEigenvectors }

  PROCEDURE NormalizeEigenvectors(VAR Eigenvectors : MatrixTyp);
  {- This procedure normalizes the eigenvectors so that the largest element in each vector is one. -}

  VAR Row      : WORD;
    Largest    : double;
    CurrentRow : VectorTyp;

  BEGIN { procedure NormalizeEigenvectors }
    FOR Row := 1 TO Dimen DO
      BEGIN
        GetRow(EigenVectors, Row, CurrentRow);
        Largest := FindLargest(CurrentRow);
        DivConstant(CurrentRow, Largest);
        SetRow(EigenVectors, CurrentRow, Row);
        DestroyVector(CurrentRow);     // was generated IN GetRow
     END;
  END; { procedure NormalizeEigenvectors }

BEGIN { procedure Jacobi }
  Dimen := MatrixRows(Mat);
  Result := TestData(Dimen, Mat, MaxIter);
  IF (Result <> 0) THEN EXIT;    // matrix NOT suitable FOR Jacobi
  Iter := 0;
  CreateIdentityMatrix(EigenVectors, Dimen);
  CreateVector(Eigenvalues, Dimen, 0);
  REPEAT
    Iter := Succ(Iter);
    SumSquareDiag := 0;
    FOR Diag := 1 TO Dimen DO
      SumSquareDiag := SumSquareDiag + Sqr(GetMatrixElement(Mat, Diag, Diag));
    Done := TRUE;
    FOR Row := 1 TO Pred(Dimen) DO
      FOR Column := Succ(Row) TO Dimen DO
        IF (Abs(GetMatrixElement(Mat, Row, Column)) > MaxError * SumSquareDiag)
          THEN
            BEGIN
              Done := FALSE;
              CalculateRotation(GetMatrixElement(Mat, Row, Row),
                GetMatrixElement(Mat, Row, Column),
                GetMatrixElement(Mat, Column, Column),
                SinTheta, CosTheta);
              RotateMatrix(SinTheta, CosTheta, Row, Column, Mat);
              RotateEigenvectors(SinTheta, CosTheta, Row, Column, Eigenvectors);
            END;
  UNTIL Done OR (Iter > MaxIter);
  FOR Diag := 1 TO Dimen DO
    SetVectorElement(Eigenvalues, Diag, GetMatrixElement(Mat, Diag, Diag));
  IF Iter > MaxIter THEN Result := 4;
END; { procedure Jacobi }


PROCEDURE SortEigenValues(VAR EigenValues: VectorTyp; VAR EigenVectors: MatrixTyp);

VAR n, nMax, i, j : WORD;
    Max           : double;
    ValuesInter   : VectorTyp;
    VectorsInter  : MatrixTyp;

BEGIN
  n := VectorLength(EigenValues);
  CopyVector(EigenValues, ValuesInter);
  CopyMatrix(EigenVectors, VectorsInter);
  FOR i := 1 TO n DO
    BEGIN
      Max := -MaxRealNumber;           // größten EigenValue identifizieren
      nMax := 0;
      FOR j := 1 TO n DO
        IF (GetVectorElement(ValuesInter, j) > Max)
          THEN
            BEGIN
              nMax := j;
              Max := GetVectorElement(ValuesInter, j);
            END;
      SetVectorElement(EigenValues, i, Max);     // sortierte Kopie erstellen
      FOR j := 1 TO n DO
        SetMatrixElement(EigenVectors, j, i, GetMatrixElement(VectorsInter, nMax, j));
      SetVectorElement(ValuesInter, nMax, -MaxRealNumber);
    // auf sehr kleinen Wert setzen
    END;
  DestroyVector(ValuesInter);
  DestroyMatrix(VectorsInter);
END;



END.    // Eigenvalues

