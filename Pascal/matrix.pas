UNIT Matrix;

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

USES Math, MathFunc, Vector;

CONST
  MatrixError: BOOLEAN = FALSE;   { toggle for error condition }

TYPE
  MatrixStruc = RECORD
    Columns, Rows: WORD;
    Data: ARRAY[1..MaxVector] OF float;
  END;
  MatrixTyp = ^MatrixStruc;

PROCEDURE CreateMatrix(VAR Mat: MatrixTyp; Rows, Columns: WORD; Value: float);

PROCEDURE DestroyMatrix(VAR Mat: MatrixTyp);

PROCEDURE ReadMatrix(MedStr: STRING; VAR A: MatrixTyp);

PROCEDURE WriteMatrix(MedStr: STRING; CONST A: MatrixTyp; ValidFigures: BYTE);

FUNCTION GetMatrixElement(CONST A: MatrixTyp; Row, Column: WORD): float;

PROCEDURE SetMatrixElement(VAR A: MatrixTyp; Row, Column: WORD; Value: float);

PROCEDURE CreateIdentityMatrix(VAR A: MatrixTyp; n: WORD);

PROCEDURE CreateNullMatrix(VAR A: MatrixTyp; n: WORD);

PROCEDURE CreateHilbertMatrix(VAR H: MatrixTyp; n: WORD);

FUNCTION MatrixRows(CONST A: MatrixTyp): WORD;

FUNCTION MatrixColumns(CONST A: MatrixTyp): WORD;

{ ************************************************************************* }

PROCEDURE GetRow(CONST A: MatrixTyp; Z: WORD; VAR Vek: VectorTyp);

PROCEDURE GetColumn(CONST A: MatrixTyp; S: WORD; VAR Vek: VectorTyp);

PROCEDURE ExchangeColumns(VAR A: MatrixTyp; Column1, Column2: WORD);

PROCEDURE ExchangeRows(VAR A: MatrixTyp; Row1, Row2: WORD);

PROCEDURE SetRow(VAR A: MatrixTyp; CONST Vek: VectorTyp; Z: WORD);

PROCEDURE SetColumn(VAR A: MatrixTyp; CONST Vek: VectorTyp; S: WORD);

PROCEDURE CopyMatrix(CONST Source: MatrixTyp; VAR Dest: MatrixTyp);

{ *************************** Matrixalgebra ******************************* }

PROCEDURE MatrixAdd(CONST A, B: MatrixTyp; VAR Res: MatrixTyp);

PROCEDURE SkalarMultiplikation(VAR A: MatrixTyp; x: float);

PROCEDURE MatrixInnerProduct(CONST A, B: MatrixTyp; VAR C: MatrixTyp);

PROCEDURE HadamardSchurProduct(CONST A, B: MatrixTyp; VAR C: MatrixTyp);

PROCEDURE MatrixDivision(CONST A, B: MatrixTyp; VAR C: MatrixTyp);

PROCEDURE CenterData(VAR A: MatrixTyp);

PROCEDURE ChangeMatrixNorm(VAR A: MatrixTyp; Norm: float);

FUNCTION Determinante(CONST A: MatrixTyp): float;

FUNCTION MatrixTrace(CONST A: MatrixTyp): float;

FUNCTION FrobeniusNorm(CONST A: MatrixTyp): float;

FUNCTION FrobeniusSkalarProduct(CONST A, B: MatrixTyp): float;

PROCEDURE InverseMatrix(VAR A: MatrixTyp);

PROCEDURE ERoMultAdd(VAR A: MatrixTyp; Faktor: float; Row1, Row2: WORD);

PROCEDURE MatrixTranspose(CONST A: MatrixTyp; VAR B: MatrixTyp);

PROCEDURE AntiSym(CONST A: MatrixTyp; VAR Symmetric, Antisymmetric: MatrixTyp);

PROCEDURE Diag(VAR Matrix: MatrixTyp);

PROCEDURE LeadingPrincipleMinors(VAR A: MatrixTyp; VAR V: VectorTyp);

PROCEDURE NegativeMatrix(VAR A: MatrixTyp);

{ ****************************** spezielle Matrizen *********************** }

FUNCTION MatrixSquare(CONST A: MatrixTyp): BOOLEAN;

FUNCTION MatrixSymmetric(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixLeftTrapezoid(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixRightTrapezoid(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixDiagonal(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixUpperTriangular(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixLowerTriangular(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixUpperHessenberg(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixLowerHessenberg(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixTridiagonal(VAR A: MatrixTyp): BOOLEAN;

FUNCTION MatrixPositivDefinite(CONST A: MatrixTyp): BOOLEAN;

FUNCTION NullMatrix(VAR A: MatrixTyp): BOOLEAN;

{ ***************************** Vector/Matrix ***************************** }

PROCEDURE DyadicVectorProduct(CONST A, B: VectorTyp; VAR C: MatrixTyp);

PROCEDURE MultMatrixVector(CONST Mat: MatrixTyp; CONST Vek: VectorTyp;
  VAR Result: VectorTyp);

PROCEDURE CangeVectorToMatrix(CONST Vek: VectorTyp; VAR Mat: MatrixTyp);

{ ************************* sorting ********************************** }

PROCEDURE ShellSortMatrix(VAR t: MatrixTyp; Column: WORD);
{ sorts the rows of a matrix by column "Column" }


IMPLEMENTATION

VAR
  CH: CHAR;

PROCEDURE CreateMatrix(VAR Mat: MatrixTyp; Rows, Columns: WORD; Value: float);

VAR
  i: WORD;
  x: longword;

BEGIN
  x := Rows * Columns * SizeOf(float) + 2 * SizeOf(WORD) + 4;
  TRY
    GetMem(Mat, x);
  except
    CH := WriteErrorMessage(' Not enough memory to create matrix');
    MatrixError := true;
    EXIT;
  END;
  Mat^.Columns := Columns;
  Mat^.Rows := Rows;
  FOR i := 1 TO (Rows * Columns) DO
    Mat^.Data[i] := Value;
END;


PROCEDURE DestroyMatrix(VAR Mat: MatrixTyp);

VAR
  x: longword;

BEGIN
  x := MatrixRows(Mat) * MatrixColumns(Mat) * SizeOf(float) + 2 * SizeOf(WORD) + 4;
  FreeMem(Mat, x);
END;

PROCEDURE ReadMatrix(MedStr: STRING; VAR A: MatrixTyp);

VAR
  i, j, n, p: WORD;
  x: float;
  Medium: TEXT;

BEGIN
  Assign(Medium, MedStr);
  Reset(Medium);
  IF IOResult <> 0
    THEN
      BEGIN
        CH := WriteErrorMessage(' Reading a matrix: File not found');
        MatrixError := true;
        EXIT;
      END;
  ReadLn(Medium, n);
  ReadLn(Medium, p);
  IF ((n * p > MaxVector))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Reading a matrix: Matrix too big');
        EXIT;
      END;
  CreateMatrix(A, n, p, 0.0);
  FOR i := 1 TO n DO
    BEGIN
      IF EoF(Medium)
        THEN
          BEGIN
            CH := WriteErrorMessage(' Reading a matrix: Unknown file format');
            MatrixError := true;
            EXIT;
          END;
      FOR j := 1 TO p DO
        BEGIN
          IF EoLn(Medium)
            THEN
              BEGIN
                CH := WriteErrorMessage(' Reading a matrix: Unknown file format');
                MatrixError := true;
                EXIT;
              END;
          Read(Medium, x);
          IF IOResult <> 0
            THEN
              BEGIN
                CH := WriteErrorMessage(' Reading a matrix: Unknown file format');
                MatrixError := true;
                EXIT;
              END;
          SetMatrixElement(A, i, j, x);
        END; { for j }
      ReadLn(Medium);
    END;  { for i }
  Close(Medium);
END;


PROCEDURE WriteMatrix(MedStr: STRING; CONST A: MatrixTyp; ValidFigures: BYTE);

VAR
  i, j: WORD;
  Medium: TEXT;

BEGIN
  Assign(Medium, MedStr);
  Rewrite(Medium);
  Writeln(Medium, A^.Rows);
  Writeln(Medium, A^.Columns);
  IF NOT (IOResult = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage(' Writing a matrix: Illegal File operation');
        MatrixError := true;
        EXIT;
      END;
  FOR i := 1 TO MatrixRows(A) DO
    BEGIN
      FOR j := 1 TO MatrixColumns(A) DO
        Write(Medium, FloatStr(GetMatrixElement(A, i, j), ValidFigures), ' ');
      Writeln(Medium);
    END;
  Close(Medium);
END;


PROCEDURE CopyMatrix(CONST Source: MatrixTyp; VAR Dest: MatrixTyp);

VAR
  i, j, p, n: WORD;

BEGIN
  n := MatrixRows(Source);
  p := MatrixColumns(Source);
  CreateMatrix(Dest, n, p, 0.0);
  IF MatrixError THEN EXIT;
  FOR i := 1 TO n DO
    FOR j := 1 TO p DO
      SetMatrixElement(Dest, i, j, GetMatrixElement(Source, i, j));
END;


FUNCTION MatrixSymmetric(VAR A: MatrixTyp): BOOLEAN;

VAR
  i, j, Dimen: WORD;
  x, y: float;

BEGIN
  Dimen := MatrixRows(A);
  IF NOT (MatrixSquare(A))
    THEN
      BEGIN
        MatrixSymmetric := FALSE;
          EXIT;
      END;
  FOR i := 1 TO Dimen DO
    FOR j := 1 TO Dimen DO
      BEGIN
        x := GetMatrixElement(A, i, j);
        y := GetMatrixElement(A, j, i);
        IF Abs(x - y) > Zero
          THEN
            BEGIN
              MatrixSymmetric := FALSE;
              EXIT;
            END
          ELSE   // make sure they are identical
            SetMatrixElement(A, i, j, GetMatrixElement(A, j, i));
      END;
  MatrixSymmetric := TRUE;
END;


PROCEDURE GetRow(CONST A: MatrixTyp; Z: WORD; VAR Vek: VectorTyp);

VAR
  j, m, n: WORD;

BEGIN
  m := MatrixRows(A);
  n := MatrixColumns(A);
  IF Z > m
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix-Error: accessing a non-existent row');
        MatrixError := TRUE;
        EXIT;
      END;
  CreateVector(Vek, n, 0.0);
  FOR j := 1 TO n DO
    SetVectorElement(Vek, j, GetMatrixElement(A, z, j));
END;


PROCEDURE SetRow(VAR A: MatrixTyp; CONST Vek: VectorTyp; Z: WORD);

VAR
  j, n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF (Z > n) OR (p <> VectorLength(Vek))
    THEN
      BEGIN
        CH := WriteErrorMessage('Matrix-error: set row with illegal parameter');
        MatrixError := TRUE;
        EXIT;
      END;
  FOR j := 1 TO p DO
    SetMatrixElement(A, z, j, GetVectorElement(Vek, j));
END;


PROCEDURE ExchangeRows(VAR A: MatrixTyp; Row1, Row2: WORD);

VAR
  Dummy: float;
  j, n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF ((Row1 > n) OR (Row2 > n))
    THEN
      BEGIN
        CH := WriteErrorMessage('Matrix-error: accessing non-existant row');
        MatrixError := TRUE;
        EXIT;
      END;
  FOR j := 1 TO p DO
    BEGIN
      Dummy := GetMatrixElement(A, Row1, j);
      SetMatrixElement(A, Row1, j, GetMatrixElement(A, Row2, j));
      SetMatrixElement(A, Row2, j, Dummy);
    END;
END; { procedure ExchangeRows }


PROCEDURE GetColumn(CONST A: MatrixTyp; S: WORD; VAR Vek: VectorTyp);

VAR
  i, n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF (S > p)
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix-error: accessing non-existant column');
        MatrixError := TRUE;
        EXIT;
      END;
  CreateVector(Vek, n, 0.0);
  FOR i := 1 TO n DO
    SetVectorElement(Vek, i, GetMatrixElement(A, i, S));
END;


PROCEDURE SetColumn(VAR A: MatrixTyp; CONST Vek: VectorTyp; S: WORD);

VAR
  i, n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF (S > p) OR (n <> VectorLength(Vek))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix-error: setting non-existant column');
        MatrixError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    SetMatrixElement(A, i, S, GetVectorElement(Vek, i));
END;


PROCEDURE ExchangeColumns(VAR A: MatrixTyp; Column1, Column2: WORD);
{ Columns n und m vertauschen }

VAR
  i, n, p: WORD;
  Dummy: float;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF ((Column1 > p) OR (Column2 > p))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix-errorr: accessing non-existant column');
        MatrixError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    BEGIN
      Dummy := GetMatrixElement(A, i, Column1);
      SetMatrixElement(A, i, Column1, GetMatrixElement(A, i, Column2));
      SetMatrixElement(A, i, Column2, Dummy);
    END;
END;


FUNCTION MatrixRows(CONST A: MatrixTyp): WORD;

BEGIN
  MatrixRows := A^.Rows;
END;


FUNCTION MatrixColumns(CONST A: MatrixTyp): WORD;

BEGIN
  MatrixColumns := A^.Columns;
END;


FUNCTION GetMatrixElement(CONST A: MatrixTyp; Row, Column: WORD): float;

VAR
  n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF (Row <= n) AND (Column <= p)
    THEN
      GetMatrixElement := A^.Data[Pred(Row) * p + Column]
    ELSE
      BEGIN
        MatrixError := true;
        CH := WriteErrorMessage(' Attempt to read a non-existent matrix element');
      END;
END;


PROCEDURE SetMatrixElement(VAR A: MatrixTyp; Row, Column: WORD; Value: float);

VAR
  n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF (Row <= n) AND (Column <= p)
    THEN
      A^.Data[Pred(Row) * p + Column] := Value
    ELSE
      BEGIN
        MatrixError := true;
        CH := WriteErrorMessage(' Attempt to write to a non-existent matrix element');
      END;
END;


PROCEDURE CreateIdentityMatrix(VAR A: MatrixTyp; n: WORD);
{ n*n identity matrix }

VAR
  i: WORD;

BEGIN
  CreateMatrix(A, n, n, 0.0);
  IF MatrixError THEN EXIT;
  FOR i := 1 TO n DO
    SetMatrixElement(A, i, i, 1.0);
END;


PROCEDURE CreateNullMatrix(VAR A: MatrixTyp; n: WORD);

BEGIN
  CreateMatrix(A, n, n, 0.0);
END;


PROCEDURE CreateHilbertMatrix(VAR H: MatrixTyp; n: WORD);

VAR
  i, j: WORD;

BEGIN
  CreateMatrix(H, n, n, 0.0);
  FOR i := 1 TO n DO
    FOR j := 1 TO n DO
      SetMatrixElement(H, i, j, 1 / Pred(i + j));
END;


FUNCTION MatrixTrace(CONST A: MatrixTyp): float;

VAR
  i, n, p: WORD;
  Sum: float;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF (n <> p)
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix-error: trace of a matrix that is not square');
        MatrixError := TRUE;
        EXIT;
      END;
  Sum := 0;
  FOR i := 1 TO n DO
    Sum := Sum + GetMatrixElement(A, i, i);
  MatrixTrace := Sum;
END;


FUNCTION FrobeniusNorm(CONST A: MatrixTyp): float;

VAR
  i, j, n, p: WORD;
  Sum: extended;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  Sum := 0.0;
  FOR i := 1 TO n DO
    FOR j := 1 TO p DO
      Sum := Sum + Sqr(GetMatrixElement(A, i, j));
  FrobeniusNorm := Sqrt(Sum);
END;


FUNCTION FrobeniusSkalarProduct(CONST A, B: MatrixTyp): float;

VAR
  CH: CHAR;
  C, D: MatrixTyp;
  n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  IF NOT ((n = MatrixRows(B)) AND (p = MatrixColumns(B)))
    THEN
      BEGIN
        CH := WriteErrorMessage(
          ' Matrix-error: Frobenius-Skalar product of incompatible matrices');
        MatrixError := TRUE;
        EXIT;
      END;
  MatrixTranspose(A, C);          // C = A^T
  MatrixInnerProduct(C, B, D);    // D = A^T B
  FrobeniusSkalarProduct := MatrixTrace(D);
  DestroyMatrix(C);
  DestroyMatrix(D);
END;


PROCEDURE ERoMultAdd(VAR A: MatrixTyp; Faktor: float; Row1, Row2: WORD);

VAR
  j: WORD;

BEGIN
  FOR j := 1 TO MatrixColumns(A) DO
    SetMatrixElement(A, Row2, j, GetmatrixElement(A, Row2, j) +
      GetMatrixElement(A, Row1, j) * Faktor);
END;


PROCEDURE MatrixTranspose(CONST A: MatrixTyp; VAR B: MatrixTyp);

VAR
  i, j, n, p: WORD;

BEGIN
  n := MatrixRows(A);
  p := MatrixColumns(A);
  CreateMatrix(B, p, n, 0.0);
  FOR i := 1 TO n DO
    FOR j := 1 TO p DO
      SetMatrixElement(B, j, i, GetMatrixElement(A, i, j));
END;


PROCEDURE Diag(VAR Matrix: MatrixTyp);

VAR
  n, p, i, j: WORD;

BEGIN
  n := MatrixRows(Matrix);
  p := MatrixColumns(Matrix);
  IF (p <> n)
    THEN
      BEGIN
        CH := WriteErrorMessage('Matrix-error: Diag of a non-square matrix');
        MatrixError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    FOR j := Succ(i) TO p DO
      BEGIN
        SetMatrixElement(Matrix, i, j, 0.0);
        SetMatrixElement(Matrix, j, i, 0.0);
      END;
END;


FUNCTION Determinante(CONST A: MatrixTyp): float;

VAR
  PartialDeter, Multiplier: float;
  Row, ReferenceRow: WORD;
  DetEqualsMaxError: BOOLEAN;
  Copy: MatrixTyp;

  PROCEDURE Pivot(ReferenceRow: WORD; VAR PartialDeter: float;
  VAR DetEqualsMaxError: BOOLEAN);
     {- This procedure searches the ReferenceRow column of the matrix Data for
        the first non-MaxError element below the diagonal. If it finds one, then
        the procedure switches rows so that the non-MaxError element is on the
        diagonal. Switching rows changes the determinant by a factor of -1;
        this change is returned in PartialDeter. If it doesn't find one, the
        matrix is singular and the Determinant is MaxError (DetEqualsMaxError = true
        is returned.  -}

  VAR
    NewRow: INTEGER;

  BEGIN
    DetEqualsMaxError := TRUE;
    NewRow := ReferenceRow;
    WHILE DetEqualsMaxError AND (NewRow < MatrixRows(Copy)) DO
      BEGIN  { Try to find a row with a non-MaxError    }
        NewRow := Succ(NewRow);
        IF Abs(GetMatrixElement(Copy, NewRow, ReferenceRow)) > MaxError
          THEN
            BEGIN
              ExchangeRows(Copy, NewRow, ReferenceRow); { Switch these two rows }
              DetEqualsMaxError := FALSE;
              PartialDeter := -PartialDeter;  { Switching rows changes }
            END; { the determinant by a factor of -1 }
      END;
  END; { procedure Pivot }

BEGIN  { Determinante }
  IF MatrixRows(A) = 1
    THEN
      BEGIN
        Result := GetMatrixElement(A, 1, 1);
        EXIT;
      END;
  CopyMatrix(A, Copy);
  IF MatrixError THEN EXIT;
  DetEqualsMaxError := FALSE;
  PartialDeter := 1;
  ReferenceRow := 0;
  { Make the matrix upper triangular }
  WHILE NOT (DetEqualsMaxError) AND (ReferenceRow < Pred(MatrixRows(Copy))) DO
    BEGIN
      INC(ReferenceRow);
      { If diagonal element is MaxError then switch rows }
      IF Abs(GetMatrixElement(Copy, ReferenceRow, ReferenceRow)) < MaxError
        THEN Pivot(ReferenceRow, PartialDeter, DetEqualsMaxError);
      IF NOT (DetEqualsMaxError)
        THEN
          FOR Row := Succ(ReferenceRow) TO MatrixRows(Copy) DO
          { Make the ReferenceRow element of this row MaxError }
            IF Abs(GetMatrixElement(Copy, Row, ReferenceRow)) > MaxError
              THEN
                BEGIN
                  Multiplier :=
                    -GetMatrixElement(Copy, Row, ReferenceRow) /
                    GetMatrixElement(Copy, ReferenceRow, ReferenceRow);
                  EROmultAdd(Copy, Multiplier, ReferenceRow, Row);
                END;
      { Multiply the diagonal Term into PartialDeter }
      PartialDeter := PartialDeter * GetMatrixElement(Copy, ReferenceRow, ReferenceRow);
    END; { while }
  IF DetEqualsMaxError
    THEN
      Result := 0
    ELSE
      Result := PartialDeter * GetMatrixElement(Copy, MatrixRows(Copy),
      MatrixColumns(Copy));
  DestroyMatrix(Copy);
END; { function Determinante }


PROCEDURE LeadingPrincipleMinors(VAR A: MatrixTyp; VAR V: VectorTyp);
// Implementierung nicht elegant, bei großen Matrizen Wiederholung vermeiden

VAR
  n, i, j, k: WORD;
  B: MatrixTyp;

BEGIN
  IF NOT MatrixSymmetric(A)
    THEN
      BEGIN
        CH := WriteErrorMessage( 'Matrix-error: leading principle minors of a non-square matrix');
        MatrixError := TRUE;
        EXIT;
      END;
  n := MatrixRows(A);
  CreateVector(V, n, 0.0);
  FOR i := 1 TO n DO
    BEGIN
      CreateMatrix(B, i, i, 0.0);
      FOR j := 1 TO i DO
        FOR k := 1 TO i DO
          SetMatrixElement(B, j, k, GetMatrixElement(A, j, k));
      SetVectorElement(V, i, Determinante(B));
      DestroyMatrix(B);
    END;
  FOR i := 1 TO n DO
    IF Abs(GetVectorElement(V, i)) < Zero THEN SetVectorElement(V, i, 0);
END;


PROCEDURE NegativeMatrix(VAR A: MatrixTyp);

VAR
  i, j: WORD;

BEGIN
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO MatrixColumns(A) DO
      SetMatrixElement(A, i, j, -GetMatrixElement(A, i, j));
END;


PROCEDURE InverseMatrix(VAR A: MatrixTyp);

VAR
  i, j, m, n, p: WORD;
  NoExch: 0..MaxVector;
  Exch: ARRAY [1..MaxVector, 1..2] OF INTEGER;


  PROCEDURE Transform(m: WORD);

  VAR
    B: MatrixTyp;
    i, j: WORD;
    Pivot: float;

  BEGIN
    CopyMatrix(A, B);
    Pivot := GetMatrixElement(B, m, m);
    IF Abs(Pivot) <= MaxError THEN
    BEGIN
      CH := WriteErrorMessage('Matrix error: inversion of singular matrix');
      MatrixError := true;
      EXIT;
    END;
    FOR i := 1 TO n DO
      FOR j := 1 TO n DO
        IF i <> m
          THEN
            IF j <> m
              THEN
                SetMatrixElement(B, i, j, GetMatrixElement(A, i, j) -
                   GetMatrixElement(A, i, m) * GetMatrixElement(A, m, j) / Pivot)
              ELSE
                SetMatrixElement(B, i, j, GetMatrixElement(A, i, j) / Pivot)
          ELSE IF j <> m
                 THEN SetMatrixElement(B, i, j, -GetMatrixElement(A, i, j) / Pivot)
                 ELSE SetMatrixElement(B, i, j, 1 / Pivot);
    DestroyMatrix(A);
    CopyMatrix(B, A);
    DestroyMatrix(B);
  END; { Transform }

BEGIN { InverseMatrix }
  NoExch := 0;
  m := MatrixRows(A);
  n := MatrixColumns(A);
  IF m <> n
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix error: inversion of non-quadratic matrix');
        MatrixError := true;
        EXIT;
      END;
  FOR i := 1 TO m DO
    BEGIN
      p := i;
      FOR j := Succ(i) TO n DO
        IF Abs(GetMatrixElement(A, j, i)) > Abs(p)
          THEN p := j;
      IF p <> i
        THEN
          BEGIN
            INC(NoExch);
            Exch[NoExch, 1] := i;
            Exch[NoExch, 2] := p;
            ExchangeRows(A, i, p);
          END;
      Transform(i);
      IF MatrixError THEN EXIT;
    END;
  FOR i := NoExch DOWNTO 1 DO
    ExchangeColumns(A, Exch[i, 2], Exch[i, 1]);
END; { InverseMatrix }


PROCEDURE MatrixAdd(CONST A, B: MatrixTyp; VAR Res: MatrixTyp);
{ Elementweise Addition zweier Matrizen. }

VAR
  i, j, Rows, Columns: WORD;

BEGIN
  Rows := MatrixRows(A);
  Columns := MatrixColumns(A);
  IF ((Rows <> MatrixRows(B)) OR (Columns <> MatrixColumns(B)))
    THEN
      BEGIN
        CH := WriteErrorMessage(
          ' Matrix error: addition of matrices not of the same size');
        MatrixError := true;
        EXIT;
      END;
  CreateMatrix(Res, Rows, Columns, 0);
  FOR i := 1 TO Rows DO
    FOR j := 1 TO Columns DO
      SetMatrixElement(Res, i, j, GetMatrixElement(A, i, j) + GetMatrixElement(B, i, j));
END;


PROCEDURE SkalarMultiplikation(VAR A: MatrixTyp; x: float);
{ Elementweise Multiplikation einer Matrix mit einem Skalar }

VAR
  i, j: WORD;

BEGIN
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO MatrixColumns(A) DO
      SetMatrixElement(A, i, j, GetMatrixElement(A, i, j) * x);
END;


PROCEDURE MatrixInnerProduct(CONST A, B: MatrixTyp; VAR C: MatrixTyp);

VAR
  i, j: WORD;
  Col, Row: VectorTyp;

BEGIN
  IF MatrixColumns(A) <> MatrixRows(B)
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix error: multiplication of A.Columns <> B.Rows');
        MatrixError := true;
        EXIT;
      END;
  CreateMatrix(C, MatrixRows(A), MatrixColumns(B), 0.0);
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO MatrixColumns(B) DO
    BEGIN
      GetRow(A, i, Row);
      GetColumn(B, j, Col);
      SetMatrixElement(C, i, j, VectorInnerProduct(Row, Col));
      DestroyVector(Row);
      DestroyVector(Col);
    END;
END;


PROCEDURE HadamardSchurProduct(CONST A, B: MatrixTyp; VAR C: MatrixTyp);

VAR
  i, j, n, p: WORD;

BEGIN
  n := MatrixRows(A);
  IF (n <> MatrixRows(B))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Hadamard-Schur multiplication: A.Rows <> B.Rows');
        MatrixError := true;
        EXIT;
      END;
  p := MatrixColumns(A);
  IF (p <> MatrixColumns(B))
    THEN
      BEGIN
        CH := WriteErrorMessage( 'Hadamard-Schur multiplication: A.Columns <> B.Columns');
        MatrixError := true;
        EXIT;
      END;
  CreateMatrix(C, n, p, 0.0);
  FOR i := 1 TO n DO
    FOR j := 1 TO p DO
      SetMatrixElement(C, i, j, GetMatrixElement(A, i, j) * GetMatrixElement(B, i, j));
END;


PROCEDURE MatrixDivision(CONST A, B: MatrixTyp; VAR C: MatrixTyp);

VAR
  Rows, Columns, i, j, k: WORD;
  M, N: MatrixTyp;

BEGIN
  Rows := MatrixRows(A);
  Columns := MatrixColumns(A);
  IF (Rows <> Columns) OR (Columns <> MatrixRows(B))
    THEN
      BEGIN
        CH := WriteErrorMessage(
          ' Matrix division: A.Columns <> B.Rows or A not quadratic');
        MatrixError := TRUE;
        EXIT;
      END;
  CopyMatrix(A, M);
  CopyMatrix(B, N);
  IF MatrixError THEN EXIT;
  FOR j := 1 TO Rows DO
    BEGIN
      IF GetMatrixElement(M, j, j) = 0
        THEN
          FOR k := j TO Rows DO
            IF GetMatrixElement(M, k, j) <> 0
              THEN
                BEGIN
                  ExchangeRows(M, j, k);
                  ExchangeRows(N, j, k);
                END
              ELSE IF k = Rows
                     THEN
                       BEGIN
                         CH := WriteErrorMessage(' Matrix division: no solution');
                         MatrixError := TRUE;
                         EXIT;
                       END;
      FOR i := Succ(j) TO Rows DO
        SetMatrixElement(M, i, j, GetMatrixElement(M, i, j) / GetMatrixElement(M, j, j));
      FOR i := Succ(j) TO Rows DO
        BEGIN
          FOR k := Succ(j) TO Rows DO
            SetMatrixElement(M, i, k, GetMatrixElement(M, i, k) -
              GetMatrixElement(M, j, k) * GetMatrixElement(M, i, j));
          FOR k := 1 TO MatrixColumns(N) DO
            SetMatrixElement(N, i, k, GetMatrixElement(N, i, k) -
              GetMatrixElement(N, j, k) * GetMatrixElement(M, i, j));
        END;
    END;
  FOR i := Rows DOWNTO 1 DO
    FOR k := 1 TO MatrixColumns(N) DO
      BEGIN
        FOR j := Succ(i) TO Columns DO
          SetMatrixElement(N, i, k, GetMatrixElement(N, i, k) -
            GetMatrixElement(N, j, k) * GetMatrixElement(M, i, j));
        SetMatrixElement(N, i, k, GetMatrixElement(N, i, k) /
          GetMatrixElement(M, i, i));
      END;
  CopyMatrix(N, C);
  DestroyMatrix(M);
  DestroyMatrix(N);
END;


PROCEDURE CenterData(VAR A: MatrixTyp);

VAR
  Data: VectorTyp;
  j, Columns: WORD;

BEGIN
  Columns := MatrixColumns(A);
  FOR j := 1 TO Columns DO
    BEGIN
      GetColumn(A, j, Data);
      Centre(Data);
      SetColumn(A, Data, j);
      DestroyVector(Data);
    END;
END;


PROCEDURE ChangeMatrixNorm(VAR A: MatrixTyp; Norm: float);

VAR
  SP: float;
  i, j: WORD;

BEGIN
  IF Abs(Norm) < MaxError
    THEN
      BEGIN
        CH := WriteErrorMessage(' Norm of a matrix must be greater than 0');
        MatrixError := true;
        EXIT;
      END;
  SP := 0;
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO MatrixColumns(A) DO
      SP := SP + Sqr(GetMatrixElement(A, i, j));
  IF Abs(SP) < MaxError
    THEN
      BEGIN
        CH := WriteErrorMessage(' Norm of a matrix: no solution');
        MatrixError := true;
        EXIT;
      END;
  Norm := Norm / Sqrt(SP);
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO MatrixColumns(A) DO
      SetMatrixElement(A, i, j, GetMatrixElement(A, i, j) * norm);
END;


PROCEDURE AntiSym(CONST A: MatrixTyp; VAR Symmetric, Antisymmetric: MatrixTyp);

VAR
  Dimen, i, j: WORD;

BEGIN
  Dimen := MatrixRows(A);
    IF Dimen <> MatrixColumns(A)
    THEN
      BEGIN
        CH := WriteErrorMessage('Antisymmetrical of an unsymmetric Matrix');
        MatrixError := true;
        EXIT;
      END;
  FOR i := 1 TO Dimen DO
    FOR j := 1 TO Dimen DO
      BEGIN
        SetMatrixElement(AntiSymmetric, i, j, 0.5 *
          (GetMatrixElement(A, i, j) + GetMatrixElement(A, j, i)));
        SetMatrixElement(Symmetric, i, j, GetMatrixElement(A, i, j) -
          GetMatrixElement(Antisymmetric, i, j));
      END;
END;

{ *********************** Matrizen und Vectoren *************************** }

PROCEDURE DyadicVectorProduct(CONST A, B: VectorTyp; VAR C: MatrixTyp);

VAR
  Row, Column: WORD;

BEGIN
  CreateMatrix(C, A^.Length, B^.Length, 0.0);
  FOR Row := 1 TO A^.Length DO
    FOR Column := 1 TO B^.Length DO
      SetMatrixElement(C, Row, Column, GetVectorElement(A, Row) *
        GetVectorElement(B, Column));
END;


PROCEDURE MultMatrixVector(CONST Mat: MatrixTyp; CONST Vek: VectorTyp;
  VAR Result: VectorTyp);

VAR
  Row, Column: WORD;
  Sum: float;

BEGIN
  IF NOT (MatrixColumns(Mat) = VectorLength(Vek))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Matrix and vector have different number of rows');
        MatrixError := true;
        EXIT;
      END;
  CreateVector(Result, MatrixRows(Mat), 0.0);
  FOR Row := 1 TO MatrixRows(Mat) DO
    BEGIN
      Sum := 0;
      FOR Column := 1 TO MatrixColumns(Mat) DO
        Sum := Sum + GetMatrixElement(Mat, Row, Column) + GetVectorElement(Vek, Column);
      SetVectorElement(Result, Row, Sum);
    END;
END;


PROCEDURE CangeVectorToMatrix(CONST Vek: VectorTyp; VAR Mat: MatrixTyp);

VAR
  j: WORD;

BEGIN
  CreateMatrix(Mat, VectorLength(Vek), 1, 0.0);
  FOR j := 1 TO VectorLength(Vek) DO
    SetMatrixElement(Mat, j, 1, GetVectorElement(Vek, j));
END;

FUNCTION MatrixSquare(CONST A: MatrixTyp): BOOLEAN;

BEGIN
  Result := MatrixRows(A) = MatrixColumns(A);
END;

FUNCTION MatrixLeftTrapezoid(VAR A: MatrixTyp): BOOLEAN;

VAR
  i, j: WORD;

BEGIN
  Result := FALSE;
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO Pred(i) DO
      IF Abs(GetMatrixElement(A, i, j)) > Zero
        THEN EXIT
        ELSE SetMatrixElement(A, i, j, 0.0);
  Result := TRUE;
END;

FUNCTION MatrixRightTrapezoid(VAR A: MatrixTyp): BOOLEAN;

VAR
  i, j: WORD;

BEGIN
  Result := FALSE;
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := Succ(i) TO MatrixColumns(A) DO
      IF Abs(GetMatrixElement(A, i, j)) > Zero
        THEN EXIT
        ELSE SetMatrixElement(A, i, j, 0.0); // make sure
  Result := TRUE;
END;


FUNCTION MatrixDiagonal(VAR A: MatrixTyp): BOOLEAN;

BEGIN
  Result := MatrixRightTrapezoid(A) AND MatrixLeftTrapezoid(A);
END;


FUNCTION MatrixUpperTriangular(VAR A: MatrixTyp): BOOLEAN;

BEGIN
  Result := MatrixLeftTrapezoid(A) AND MatrixSquare(A);
END;


FUNCTION MatrixLowerTriangular(VAR A: MatrixTyp): BOOLEAN;

BEGIN
  Result := MatrixRightTrapezoid(A) AND MatrixSquare(A);
END;


FUNCTION MatrixUpperHessenberg(VAR A: MatrixTyp): BOOLEAN;

VAR
  i, j: WORD;

BEGIN
  Result := FALSE;
  IF NOT (MatrixSquare(A)) THEN EXIT;
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO i - 2 DO
      IF Abs(GetMatrixElement(A, i, j)) > Null
        THEN EXIT
        ELSE SetMatrixElement(A, i, j, 0.0);
  Result := TRUE;
END;


FUNCTION MatrixLowerHessenberg(VAR A: MatrixTyp): BOOLEAN;

VAR
  i, j: WORD;

BEGIN
  Result := FALSE;
  IF NOT (MatrixSquare(A)) THEN EXIT;
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := i + 2 TO MatrixColumns(A) DO
      IF Abs(GetMatrixElement(A, i, j)) > Null
        THEN EXIT
        ELSE SetMatrixElement(A, i, j, 0.0);
  Result := TRUE;
END;



FUNCTION MatrixTridiagonal(VAR A: MatrixTyp): BOOLEAN;

BEGIN
  Result := MatrixUpperHessenberg(A) AND MatrixLowerHessenberg(A);
END;


FUNCTION MatrixPositivDefinite(CONST A: MatrixTyp): BOOLEAN;

VAR
  i, j: WORD;
  Akt, Max: float;

BEGIN
  Result := FALSE;
  Max := 0.0;
  IF NOT (MatrixSquare(A)) THEN EXIT;
  FOR i := 1 TO MatrixRows(A) DO
    BEGIN
      Akt := GetMatrixElement(A, i, i);
      IF Akt <= 0 THEN EXIT;
      IF Akt > Max THEN Max := Akt;
    END;
  FOR i := 1 TO MatrixRows(A) DO
    BEGIN
      FOR j := 1 TO Pred(i) DO
        BEGIN
          Akt := GetMatrixElement(A, i, j);
          IF Akt > Max
            THEN EXIT;
          IF Sqr(Akt) >= GetMatrixElement(A, i, i) * GetMatrixElement(A, j, j)
            THEN EXIT;
        END;
      FOR j := Succ(i) TO MatrixColumns(A) DO
      BEGIN
        Akt := GetMatrixElement(A, i, j);
        IF Akt > Max
          THEN EXIT;
        IF Sqr(Akt) >= GetMatrixElement(A, i, i) * GetMatrixElement(A, j, j)
          THEN EXIT;
      END;
    END;
  Result := TRUE;
END;


FUNCTION NullMatrix(VAR A: MatrixTyp): BOOLEAN;

VAR
  i, j: WORD;

BEGIN
  Result := TRUE;
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO MatrixColumns(A) DO
      IF Abs(GetMatrixElement(A, i, j)) > Zero
        THEN
          BEGIN
            Result := FALSE;
            EXIT;
          END
        ELSE
          SetMatrixElement(A, i, j, 0.0);
END;


PROCEDURE ShellSortMatrix(VAR t: MatrixTyp; Column: WORD);

LABEL
  10;

VAR
  i, j, k, l, m, nn, NaNs, n, LdN: INTEGER;
  s: float;
  tmp1, tmp2: VectorTyp;

BEGIN
  NaNs := 0;
  s := -MaxRealNumber;
  n := MatrixRows(t);
  FOR i := 1 TO n DO
    BEGIN
      IF IsNaN(GetMatrixElement(t, i, Column))    // check forR NaN data
        THEN INC(NaNs)
        ELSE IF (GetMatrixElement(t, i, Column)) > s
               THEN
                 s := GetMatrixElement(t, i, Column); // and find largest element of data vector
    END;
  s := 10 * s;
  IF (NaNs > 0)
    THEN
      FOR i := 1 TO n DO
        IF IsNaN(GetMatrixElement(t, i, Column))
          THEN SetMatrixElement(t, i, Column, s);
  // replace all NaN WITH very large number so they Move TO END OF vector
  LdN := Trunc(Ln(n) / Const_ln2);
  m := n;
  FOR nn := 1 TO LdN DO
    BEGIN
      m := m DIV 2;
      k := n - m;
      FOR j := 1 TO k DO
        BEGIN
          i := j;
10:       l := i + m;
          IF (GetMatrixElement(t, l, Column) < GetMatrixElement(t, i, Column))
            THEN
              BEGIN
                GetRow(t, i, tmp1);
                GetRow(t, l, tmp2);
                SetRow(t, tmp2, i);
                SetRow(t, tmp1, l);
                DestroyVector(tmp1);
                DestroyVector(tmp2);
                i := i - m;
                IF i >= 1 THEN GOTO 10;
              END;
        END;
    END;
  IF (NaNs > 0)
    THEN
      FOR i := Succ(n - NaNs) TO n DO  // change the top NaNs elements back TO NaN
        SetMatrixElement(t, i, Column, NaN);
END;

END. { unit Matrix }

