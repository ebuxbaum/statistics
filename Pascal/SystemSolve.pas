UNIT SystemSolve;
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

USES Math, MathFunc, Vector, Matrix;

TYPE
  StrArrayType = ARRAY[0..11] of STRING[30];

CONST
  SystemError: StrArrayType = ('ok', 'dimension error', 'matrix singular',
    'MaxIter < 1', 'no conversion', 'evaluation short-circuited',
    '', 'matrix not symmetric', 'not enough memory', '', '', '');

FUNCTION MakeUpperTriangular(VAR Koeff : MatrixTyp; VAR Right : VectorTyp) : BYTE;

FUNCTION HadamardConditionNumber (VAR Mat : MatrixTyp) : double;

FUNCTION GaussElimination(CONST Coefficients: MatrixTyp;
    VAR RightSide, Solution: VectorTyp): BYTE;

FUNCTION PartialPivoting(CONST Coefficients: MatrixTyp;
    VAR RightSide, Solution: VectorTyp): BYTE;

FUNCTION LU_Decompose(CONST Coefficients: MatrixTyp;
    VAR Decomp, Permute: MatrixTyp): BYTE;

FUNCTION LU_Solve(CONST Decomp, Permute: MatrixTyp; CONST RightSide: VectorTyp;
    VAR Solution: VectorTyp): BYTE;

FUNCTION Householder (CONST A : MatrixTyp; VAR Q, R : MatrixTyp) : BYTE;

IMPLEMENTATION

VAR CH  : CHAR;


FUNCTION MakeUpperTriangular(VAR Koeff : MatrixTyp; VAR Right : VectorTyp) : BYTE;

VAR Multiplier               : double;
    Row, ReferenceRow, Dimen : WORD;

  FUNCTION Pivot(ReferenceRow: INTEGER) : BYTE;
        { This procedure searches the ReferenceRow column of the Coefficients
          matrix for the first non-MaxError element below the diagonal. If it
          finds one, then the procedure switches rows so that the non-MaxError
          element is on the diagonal. It also switches the corresponding
          elements in the Constants vector. If it doesn't find one, the
          matrix is singular and no solution exists (Error = 2 is returned). }

  VAR NewRow : INTEGER;
      Dummy  : double;

  BEGIN
    Result := 2;          { No solution exists }
    NewRow := ReferenceRow;
    WHILE (Result > 0) AND (NewRow < Dimen) DO
      BEGIN    { Try to find a row with a non-MaxError diagonal element    }
        NewRow := Succ(NewRow);
        IF Abs(GetMatrixElement(Koeff, NewRow, ReferenceRow)) > MaxError
          THEN
            BEGIN
              ExchangeRows(Koeff, NewRow, ReferenceRow);
              Dummy := GetVectorElement(Right, NewRow);
              SetVectorElement(Right, NewRow, GetVectorElement(Right, ReferenceRow));
              SetVectorElement(Right, ReferenceRow, Dummy);
              Result := 0;    { Solution may exist }
            END;
      END;
  END; { procedure Pivot }

BEGIN { procedure UpperTriangular }
  Dimen := MatrixRows(Koeff);
  ReferenceRow := 0;
  Result := 0;
  WHILE ((Result = 0) AND (ReferenceRow < Pred(Dimen))) DO
    BEGIN
      INC(ReferenceRow);  { Check to see if the main diagonal element is MaxError }
      IF Abs(GetMatrixElement(Koeff, ReferenceRow, ReferenceRow)) < MaxError
        THEN Result := Pivot(ReferenceRow);
      IF Result = 0
        THEN
          FOR Row := Succ(ReferenceRow) TO Dimen DO
        { Make the ReferenceRow element of this row MaxError }
            IF Abs(GetMatrixElement(Koeff, Row, ReferenceRow)) > MaxError
              THEN
                BEGIN
                  Multiplier := -GetMatrixElement(Koeff, Row, ReferenceRow) / GetMatrixElement(Koeff, ReferenceRow, ReferenceRow);
                  EROmultAdd(Koeff, Multiplier, ReferenceRow, Row);
                  SetVectorElement(Right, Row, GetVectorElement(Right, Row) + Multiplier * GetVectorElement(Right, ReferenceRow));
                END;
    END; { while }
  IF Abs(GetMatrixElement(Koeff, MatrixRows(Koeff), MatrixColumns(Koeff))) < MaxError
    THEN Result := 2;    { No solution }
END; { procedure UpperTriangular }


FUNCTION HadamardConditionNumber (VAR Mat : MatrixTyp) : double;

VAR n, i, j           : WORD;
    res               : BYTE;
    Decomp, Permute   : MatrixTyp;
    temp, cond        : double;

BEGIN
  n := MatrixRows(Mat);
  IF NOT(MatrixSymmetric(Mat))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Hadamard condition number of a non-symmetric matrix');
        Result := NaN;
        EXIT;
      END;
  res := LU_Decompose(Mat, Decomp, Permute);
  IF MatrixError
    THEN
      BEGIN
        Result := NaN;
        EXIT;
      END;
  IF res = 2
    THEN
      Result := 0    // singular matrix
    ELSE
      BEGIN
        cond := 1.0;
        FOR i := 1 TO n DO
          BEGIN
            temp := 0.0;
            FOR j := 1 TO n DO
              temp := temp + Sqr(GetMatrixElement(Mat, i, j));
            cond := cond * GetMatrixElement(Decomp, i, i) / Sqrt(temp);
          END;
        Result := Abs(cond);
      END;
  DestroyMatrix(Decomp);
  DestroyMatrix(Permute);
END;

PROCEDURE BackwardSubst (VAR Koeff : MatrixTyp; VAR Right, Solution : VectorTyp);
{ This procedure applies backwards substitution to the upper triangular
  Coefficients matrix and Constants vector. The resulting vector is the
  solution to the set of equations and is stored in the vector Solution. }

VAR Term, Row, Dimen : WORD;
    Sum              : double;

BEGIN
  Dimen := MatrixRows(Koeff);
  FOR Term := Dimen DOWNTO 1 DO
    BEGIN
      Sum := 0;
      FOR Row := Succ(Term) TO Dimen DO
        Sum := Sum + GetMatrixElement(Koeff, Term, Row) * GetVectorElement(Solution, Row);
      SetVectorElement(Solution, Term, (GetVectorElement(Right, Term) - Sum) / GetMatrixElement(Koeff, Term, Term));
  END;
END; { procedure BackwardsSub }


FUNCTION Initial(CONST Coefficients : MatrixTyp; VAR RightSide, Solution : VectorTyp) : BYTE;

VAR Dimen : WORD;

BEGIN
  Result := 0;
  Dimen := MatrixRows(Coefficients);
  IF ((Dimen < 1) OR (Dimen <> MatrixColumns(Coefficients)) OR (Dimen <> VectorLength(RightSide)))
    THEN
      Result := 1
    ELSE
      IF Dimen = 1
        THEN
          IF Abs(GetMatrixElement(Coefficients, 1, 1)) < MaxError
            THEN Result := 2
            ELSE SetVectorElement(Solution, 1, GetVectorElement(RightSide, 1) / GetMatrixElement(Coefficients, 1, 1));
END; { procedure Initial }


FUNCTION GaussElimination(CONST Coefficients : MatrixTyp;
  VAR RightSide, Solution : VectorTyp) : BYTE;

VAR Error : BYTE;
    Dimen : WORD;
    Koeff : MatrixTyp;
    Right : VectorTyp;

BEGIN { procedure Gaussian_Elimination }
  Dimen := MatrixRows(Coefficients);
  Result := Initial(Coefficients, RightSide, Solution);
  IF result <> 0 THEN EXIT;
  CreateVector(Right, Dimen, 0.0);
  CreateVector(Solution, Dimen, 0.0);
  CopyVector(RightSide, Right);
  IF VectorError
    THEN
      BEGIN
        VectorError := FALSE;
        Result := 8;
        EXIT;
      END;
  CopyMatrix(Coefficients, Koeff);
  IF MatrixError
    THEN
      BEGIN
        MatrixError := FALSE;
        Result := 8;
        EXIT;
      END;
  IF MatrixRows(Koeff) > 1
    THEN
      BEGIN
        Error := MakeUpperTriangular(Koeff, Right);
        IF Error = 0 THEN BackwardSubst(Koeff, Right, Solution);
      END;
  Result := Error;
  DestroyVector(Right);
  DestroyMatrix(Koeff);
END; { procedure Gaussian_Elimination }

    { ************************************************************************** }

FUNCTION PartialPivoting(CONST Coefficients: MatrixTyp;
  VAR RightSide, Solution: VectorTyp): BYTE;

VAR Error : BYTE;
    Dimen : WORD;
    Koeff : MatrixTyp;
    Right : VectorTyp;

BEGIN  { procedure PartialPivoting }
  Dimen := MatrixRows(Coefficients);
  Result := Initial(Coefficients, RightSide, Solution);
  IF result <> 0 THEN EXIT;
  CopyMatrix(Coefficients, Koeff);
  IF MatrixError
    THEN
      BEGIN
        MatrixError := FALSE;
        Result := 8;
        EXIT;
      END;
  CreateVector(Solution, Dimen, 0.0);
  CopyVector(RightSide, Right);
  IF VectorError
    THEN
      BEGIN
        VectorError := FALSE;
        Result := 8;
        EXIT;
      END;
  IF Dimen > 1
    THEN
      BEGIN
        Error := MakeUpperTriangular(Koeff, Right);
        IF Error = 0 THEN BackwardSubst(Koeff, Right, Solution);
      END;
  Result := Error;
  DestroyMatrix(Koeff);
  DestroyVector(Right);
END; { Partial_Pivoting }

{ ************************************************************************** }

FUNCTION LU_Decompose(CONST Coefficients : MatrixTyp; VAR Decomp, Permute : MatrixTyp): BYTE;

VAR Error: BYTE;
    Koeff, Upper, Lower: MatrixTyp;
    Dimen: WORD;


  FUNCTION RowColumnMult(VAR Lower, Upper: MatrixTyp; Row, Column: WORD): double;
    { Function return: dot product of row Row of Lower and column Column of Upper }

  VAR Term : INTEGER;
      Sum  : double;

  BEGIN
    Sum := 0;
    FOR Term := 1 TO Pred(Row) DO
      Sum := Sum + GetMatrixElement(Lower, Row, Term) * GetMatrixElement(Upper, Term, Column);
    Result := Sum;
  END; { RowColumnMult }


  PROCEDURE Pivot(ReferenceRow: WORD; VAR Error: BYTE);
     { This procedure searches the ReferenceRow column of the Coefficients
       matrix for the element in the Row below the main diagonal which
       produces the largest value of Coefficients[Row, ReferenceRow] - Sum
       (for K=1 to pred(ReferenceRow) of Upper[Row, k] - Lower[k, ReferenceRow]
      If it finds one, then the procedure switches rows so that this element
      is on the main diagonal. The procedure also switches the corresponding
      elements in the Permute matrix and the Lower matrix. If the largest
      value of the above expression is MaxError, then the matrix is singular and
      no solution exists (Error = 2 is returned).   }

  VAR PivotRow, Row      : WORD;
      ColumnMax, TestMax : double;

  BEGIN { procedure Pivot }
    { First, find the row with the largest TestMax }
    PivotRow := ReferenceRow;
    ColumnMax := Abs(GetMatrixElement(Koeff, ReferenceRow, ReferenceRow) -
      RowColumnMult(Lower, Upper, ReferenceRow, ReferenceRow));
    FOR Row := Succ(ReferenceRow) TO Dimen DO
      BEGIN
        TestMax := Abs(GetMatrixElement(Koeff, Row, ReferenceRow) - RowColumnMult(Lower, Upper, Row, ReferenceRow));
        IF TestMax > ColumnMax
          THEN
            BEGIN
              PivotRow := Row;
              ColumnMax := TestMax;
            END;
      END;
    IF PivotRow <> ReferenceRow
      THEN   { Second, switch these two rows }
        BEGIN
          ExchangeRows(Koeff, PivotRow, ReferenceRow);
          ExchangeRows(Lower, PivotRow, ReferenceRow);
          ExchangeRows(Permute, PivotRow, ReferenceRow);
        END
      ELSE { If ColumnMax is MaxError, no solution exists }
        IF ColumnMax < MaxError
          THEN Error := 2;
  END; { procedure Pivot }


  PROCEDURE Decompose(VAR Error: BYTE);
     { This procedure decomposes the Coefficients matrix into two triangular
       matrices, a lower and an upper one.  The lower and upper matrices are
       combined into one matrix, Decomp.  The permutation matrix, Permute,
       records the effects of partial pivoting.  }

  VAR Term, Index: INTEGER;

    PROCEDURE Initialize;
          { This procedure initializes Lower and Upper to the MaxError matrix
            and Permute to the identity matrix. }

    BEGIN
      CreateMatrix(Upper, Dimen, Dimen, 0.0);
      CreateMatrix(Lower, Dimen, Dimen, 0.0);
      CreateIdentityMatrix(Permute, Dimen);
      IF MatrixError
        THEN
          BEGIN
            MatrixError := FALSE;
            Error := 8;
          END;
    END; { procedure Initialize }

  BEGIN { Decompose }
    Initialize;
    { partial pivoting on row 1 }
    Pivot(1, Error);
    IF Error = 0
      THEN
        BEGIN
          SetMatrixElement(Lower, 1, 1, 1);
          SetMatrixElement(Upper, 1, 1, GetMatrixElement(Koeff, 1, 1));
          FOR Term := 1 TO Dimen DO
            BEGIN
              SetMatrixElement(Lower, Term, 1, GetMatrixElement(Koeff, Term, 1) / GetMatrixElement(Upper, 1, 1));
              SetMatrixElement(Upper, 1, Term, GetMatrixElement(Koeff, 1, Term) / GetMatrixElement(Lower, 1, 1));
            END;
        END;
    Term := 1;
    WHILE (Error = 0) AND (Term < Pred(Dimen)) DO
      BEGIN
        Term := Succ(Term); { perform partial pivoting on row Term }
        Pivot(Term, Error);
        SetMatrixElement(Lower, Term, Term, 1);
        SetMatrixElement(Upper, Term, Term,
           GetMatrixElement(Koeff, Term, Term) - RowColumnMult(Lower, Upper, term, term));
        IF Abs(GetMatrixElement(Upper, Term, Term)) < MaxError
          THEN
            Error := 2   { no solutions }
          ELSE
            FOR Index := Succ(Term) TO Dimen DO
              BEGIN
                SetMatrixElement(Upper, Term, Index, GetMatrixElement(Koeff, Term, Index) - RowColumnMult(Lower, Upper, Term, Index));
                SetMatrixElement(Lower, Index, Term, (GetMatrixElement(Koeff, Index, Term) - RowColumnMult(Lower, Upper, Index, Term)) /
                    GetMatrixElement(Upper, Term, Term));
              END;
      END;
    SetMatrixElement(Lower, Dimen, Dimen, 1);
    SetMatrixElement(Upper, Dimen, Dimen, GetMatrixElement(Koeff, Dimen, Dimen) - RowColumnMult(Lower, Upper, Dimen, Dimen));
    IF Abs(GetMatrixElement(Upper, Dimen, Dimen)) < MaxError THEN Error := 2;
    { Combine the upper and lower triangular matrices into one }
    CopyMatrix(Upper, Decomp);
    FOR Term := 2 TO Dimen DO
      FOR Index := 1 TO Pred(Term) DO
        SetMatrixElement(Decomp, Term, Index, GetMatrixElement(Lower, Term, Index));
    DestroyMatrix(Upper);
    DestroyMatrix(Lower);
  END; { procedure Decompose }

BEGIN { LU_Decompose }
  Dimen := MatrixRows(Coefficients);
  CopyMatrix(Coefficients, Koeff);
  IF MatrixError
    THEN
      BEGIN
        MatrixError := FALSE;
        Result := 8;
        EXIT;
      END;
  IF Dimen < 1
    THEN
      Error := 1
    ELSE
      BEGIN
        Error := 0;
        IF Dimen = 1
          THEN
            BEGIN
              CopyMatrix(Koeff, Decomp);
              SetMatrixElement(Permute, 1, 1, 1);
            END
          ELSE
            Decompose(Error);
      END;
  Result := Error;
  DestroyMatrix(Koeff);
END; { LU_Decompose }

{ ************************************************************************** }

FUNCTION LU_Solve(CONST Decomp, Permute : MatrixTyp; CONST RightSide : VectorTyp;
  VAR Solution : VectorTyp): BYTE;

VAR Error    : BYTE;
    Dimen    : WORD;
    DEC, Per : MatrixTyp;
    Right    : VectorTyp;


  PROCEDURE FindSolution;
     { The Decom matrix contains a lower and upper triangular matrix. This
       procedure performs a two step backwards substitution to compute the
       solution to the system of equations.  First, forward substitution is
       applied to the lower triangular matrix and Constants vector yielding
       PartialSolution.  Then backwards substitution is applied to the Upper
       matrix and the PartialSolution vector yielding Solution.  }

  VAR PartialSolution : Vectortyp;
      Term, Index     : WORD;
      Sum             : double;

  BEGIN { FindSolution }
    { First solve the lower triangular matrix }
    CreateVector(PartialSolution, Dimen, 0.0);
    SetVectorElement(PartialSolution, 1, GetVectorElement(Right, 1));
    FOR Term := 2 TO Dimen DO
      BEGIN
        Sum := 0;
        FOR Index := 1 TO Pred(Term) DO
          IF Term = Index
          THEN
            Sum := Sum + GetVectorElement(PartialSolution, Index)
          ELSE
            Sum := Sum + GetMatrixElement(DEC, Term, Index) * GetVectorElement(PartialSolution, Index);
        SetVectorElement(PartialSolution, Term, GetVectorElement(Right, Term) - Sum);
      END;
    { Then solve the upper triangular matrix }
    SetVectorElement(Solution, Dimen, GetVectorElement(PartialSolution,Dimen) / GetMatrixElement(DEC, Dimen, Dimen));
    FOR Term := Pred(Dimen) DOWNTO 1 DO
      BEGIN
        Sum := 0;
        FOR Index := Succ(Term) TO Dimen DO
          Sum := Sum + GetMatrixElement(DEC, Term, Index) * GetVectorElement(Solution, Index);
        SetVectorElement(Solution, Term, (GetVectorElement(PartialSolution, Term) - Sum) / GetMatrixElement(DEC, Term, Term));
      END;
    DestroyVector(PartialSolution);
  END; { procedure FindSolution }


  PROCEDURE PermuteRightSide;

  VAR Row, Column   : WORD;
      Entry         : double;
      TempConstants : VectorTyp;

  BEGIN
    CreateVector(TempConstants, Dimen, 0.0);
    FOR Row := 1 TO Dimen DO
      BEGIN
        Entry := 0;
        FOR Column := 1 TO Dimen DO
          Entry := Entry + GetMatrixElement(Per, Row, Column) * GetVectorElement(Right, Column);
        SetVectorElement(TempConstants, Row, Entry);
      END;
    CopyVector(TempConstants, Right);
    DestroyVector(TempConstants);
  END; { PermuteRightSide }

BEGIN { Solve_LU_Decompostion }
  Dimen := MatrixRows(Decomp);
  CopyMatrix(Decomp, DEC);
  CopyMatrix(Permute, Per);
  CreateVector(Solution, Dimen, 0.0);
  IF MatrixError
    THEN
      BEGIN
        MatrixError := FALSE;
        Result := 8;
        EXIT;
      END;
  CopyVector(RightSide, Right);
  IF VectorError
    THEN
      BEGIN
        VectorError := FALSE;
        Result := 8;
        EXIT;
      END;
  IF Dimen < 1
    THEN
      Error := 1
    ELSE
      BEGIN
        Error := 0;
        PermuteRightSide;
        FindSolution;
      END;
  Result := Error;
  DestroyMatrix(DEC);
  DestroyMatrix(Per);
END; { procedure LU_Solve }

// ************* QR-Decomposition by Householder reflection *******************

FUNCTION Householder (CONST A : MatrixTyp; VAR Q, R : MatrixTyp) : BYTE;

VAR Rows, Columns,
    i, j, k          : WORD;
    c                : CHAR;
    alpha, beta      : double;
    x, u, e          : VectorTyp;
    Identity, Product,
    H, Hs            : MatrixTyp;

BEGIN
  Rows := MatrixRows(A);
  Columns := MatrixColumns(A);
  IF Rows < Columns
    THEN
      BEGIN
        c := WriteErrorMessage('QR decomposition of matrix with more columns than rows');
        MatrixError := TRUE;
        Householder := 1;                        // dimension error
        EXIT;
      END;
  CreateIdentityMatrix(Q, Rows);                 // neutral element OF matrix multiplication
  IF MatrixError THEN EXIT;
  CopyMatrix(A, R);                              // so that A IS unchanged
  IF MatrixError THEN EXIT;
  FOR j := 1 TO min(Pred(Rows), Columns) DO      // calculate j-th Householder Matrix
    BEGIN
      CreateVector(x, Succ(Rows-j), 0.0);
      FOR i := j TO Rows DO
        BEGIN
          k := Succ(i-j);
          SetVectorElement(x, k, GetMatrixElement(R, i, j));
        END;
      alpha := -Signum(GetVectorElement(x, 1)) * VectorEuklidianNorm(x);
      CreateVector(e, Succ(Rows-j), 0.0);
      SetVectorElement(e, 1, alpha);
      VectorAdd(x, e, u);                        // u = x + sgn(x_1) * ||x||_2 * (1, 0, ..., 0)
      DestroyVector(x);
      DestroyVector(e);
      DivConstant(u, Abs(GetVectorElement(u, 1))); // scale u by Abs(first element) -> v
      beta := 2 / VectorInnerProduct(u, u);      // \beta = 2 / (v \odot v)
      DyadicVectorProduct(u, u, Product);        // v \omul v
      DestroyVector(u);
      SkalarMultiplikation(Product, -beta);      //  beta * (v\omul v)
      CreateIdentityMatrix(Identity, Succ(Rows-j));
      MatrixAdd(Identity, Product, hs);          // Hs = I - beta * (v\omul v)
      DestroyMatrix(Identity);
      DestroyMatrix(Product);
      CreateIdentityMatrix(H, Rows);
      FOR i := j TO Rows DO                      // Copy Hs into lower right corner OF H
        FOR k := j TO Columns DO
          SetMatrixElement(H, i, k, GetMatrixElement(Hs, i-Pred(j), k-Pred(j)));
      MatrixInnerProduct(H, R, Product);         // R = HR
      DestroyMatrix(R);
      CopyMatrix(Product, R);
      DestroyMatrix(Product);
      MatrixInnerProduct(Q, H, Product);         // Q = QH
      DestroyMatrix(Q);
      CopyMatrix(Product, Q);
      DestroyMatrix(Product);
    END;
  Result := 0;                                   // everything ok
END;


END.   // SystemSolve

