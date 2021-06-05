UNIT CrossValidation;

INTERFACE

USES MathFunc, Vector, Matrix, Zufall;

CONST
  MaxK = 15;
  CrossValidationError: BOOLEAN = FALSE;

TYPE
  SplitDataTyp = ARRAY [1..MaxK] OF MatrixTyp;


PROCEDURE SplitDataMatrix(CONST Data: MatrixTyp;   // data matrix
  kValidate: WORD;                                 // no OF validation groups
  VAR SplitData: SplitDataTyp);                    // randomly distributed data
{ randomly splits the data matrix into kValidate sub-matrices }

PROCEDURE CreateTestData(CONST SplitData: SplitDataTyp;
  k, kValidate: WORD;
  VAR TestData: MatrixTyp);
{ combine all Groups except the out of box group k into test matrix }


IMPLEMENTATION

PROCEDURE SplitDataMatrix(CONST Data: MatrixTyp; kValidate: WORD;
                          VAR SplitData: SplitDataTyp);

VAR
  h, i, j, n, p : WORD;
  c             : CHAR;
  Available     : ARRAY [1..MaxK] OF WORD;
  CurrentRow    : VectorTyp;

BEGIN
  IF (kValidate > MaxK)
    THEN
      BEGIN
        CrossValidationError := TRUE;
        c := WriteErrorMessage('k-fold cross-validation: kValidate > maximum');
        EXIT;
      END;
  n := MatrixRows(Data);
  p := MatrixColumns(Data);
  j := n DIV kValidate;  // number OF elements OF all submatrices except last
  FOR i := 1 TO Pred(kValidate) DO
    BEGIN
      CreateMatrix(SplitData[i], j, p, 0.0);
      IF MatrixError
        THEN
          BEGIN
            CrossValidationError := TRUE;
            MatrixError := FALSE;
            c := WriteErrorMessage('k-fold cross-validation: not enough memory');
            EXIT;
          END;
      Available[i] := j;
    END;
  CreateMatrix(SplitData[kValidate], j + (n MOD kValidate), p, 0.0);
  // put left-overs into last group
  IF MatrixError
    THEN
      BEGIN
        CrossValidationError := TRUE;
        MatrixError := FALSE;
        c := WriteErrorMessage('k-fold cross-validation: not enough memory');
        EXIT;
      END;
  Available[kValidate] := j + (n MOD kValidate);
  FOR i := 1 TO n DO // randomly put each data row into one OF the kValidate submatrices
    BEGIN
      REPEAT
        j := Round(RandomLaplace(1, kValidate)); // select group
      UNTIL (Available[j] > 0);
      GetRow(Data, i, CurrentRow);
      h := Succ(MatrixRows(SplitData[j]) - Available[j]);
      SetRow(SplitData[j], CurrentRow, h);
      DEC(Available[j]);
      DestroyVector(CurrentRow);
    END;
END; { SplitDataMatrix }

PROCEDURE CreateTestData(CONST SplitData: SplitDataTyp; k, kValidate: WORD;
                         VAR TestData: MatrixTyp);

VAR
  j, n, p, Sum : WORD;
  i            : 1..MaxK;
  v            : VectorTyp;
  c            : CHAR;

BEGIN
  Sum := 0;
  FOR i := 1 TO kValidate DO  // calculate number OF rows OF TEST data
    IF i <> k
      THEN
        BEGIN
          n := MatrixRows(SplitData[i]);
          Sum := Sum + n;
        END;
  p := MatrixColumns(SplitData[1]);
  CreateMatrix(TestData, Sum, p, 0.0);
  IF MatrixError
    THEN
      BEGIN
        CrossValidationError := TRUE;
        MatrixError := FALSE;
        c := WriteErrorMessage('k-fold cross-validation: not enough memory');
        EXIT;
      END;
  Sum := 0;
  FOR i := 1 TO kValidate DO
    IF i <> k
      THEN
        BEGIN
          FOR j := 1 TO MatrixRows(SplitData[i]) DO
            BEGIN
              INC(Sum);
              GetRow(SplitData[i], j, v);
              SetRow(TestData, v, Sum);
              DestroyVector(v);
            END;
        END;
END; { CreateTestData }

END.

