UNIT PCA;

{ Performs routines in the context of principle component analysis }

INTERFACE

USES Math, Mathfunc, Stat, Correlations, Vector, Matrix;

CONST
  MaxVariables = 180;
  SepChar = ';';
  // IN Middle Europe variables IN CSV separated by ";" AS "," IS decimal separator
  SigVecs: WORD = 5;
  // number OF significant vectors, can be changed by calling PROGRAM
  ProblemName: STRING = ''; // first name OF all files produced;
  Border = 20;  // number OF ranges FOR statistical analysis OF correlations

TYPE
  DataTypes = (binary, nominal, ordinal, interval, rational);
  TypeVector = ARRAY[1..MaxVariables] OF DataTypes;
  FreqsType = ARRAY[DataTypes, DataTypes, -Border..Border] OF WORD;


PROCEDURE ReadCSV(VAR Types: TypeVector; VAR Data: MatrixTyp);
{ Read data from CSV file. }

PROCEDURE CalculateCorrelationMatrix(CONST Data: MatrixTyp;
  CONST Types: TypeVector; VAR Cor: MatrixTyp);
{ calculates correlation matrix for a mixed type data matrix. NaN-values
  are handled. }

PROCEDURE ReadCorrelations(VAR Cor: MatrixTyp);
{ If the correlation matrix has been calculated previously,
  read it from CSV file }

PROCEDURE WriteCorrelations(CONST Cor: MatrixTyp; Shrunk: BOOLEAN);
{ Writes the correlation matrix into a csv-file. File name will depend on
  whether the matrix has been shrunk.}

PROCEDURE AnalyseFrequencies(CONST Cor: MatrixTyp; CONST Types: TypeVector;
  VAR Freqs: FreqsType);
{ Determine distribution of correlation coefficients by variable type }

FUNCTION Bartlett(CONST Cor: MatrixTyp; Cases: WORD): double;
{ Bartlett's test for sphericity }

PROCEDURE ExplainedVariance(CONST EigenValues: VectorTyp;
  VAR Acc, Variance, CumVariance: VectorTyp);
{ calculate acceleration, explained variance and cumulative explained variance
  from eigenvalues and writes them into a csv-file }

PROCEDURE WriteEigenVector(CONST EigenVectors: MatrixTyp);
{ writes eigenvectors into a csv-file }

PROCEDURE MaximumLikelihood(CONST Data: MatrixTyp; CONST Types: TypeVector;
  VAR ImputVector: VectorTyp);
{ Determin average (cardinal) or most common (other) value of a column for
  imputation }

PROCEDURE RobustProduct(CONST A, B: MatrixTyp; CONST ImputVector: VectorTyp;
  VAR C: MatrixTyp);
{ robust matrix product, elements of A (original data) may be NaN and are replaced
  by the j-th element in ImputVector. For B (Eigenvectors) this precaution is not
  necessary. Only the first Max eigenvectors are used. }

PROCEDURE WriteComponentScores(CONST Scores: MatrixTyp);
{ writes component scores into a csv-file }

PROCEDURE CalculateLoadings(CONST Data, Scores: MatrixTyp;
  CONST Types: TypeVector; VAR Loadings: MatrixTyp);

PROCEDURE CalculateCommunalities(CONST Loadings: MatrixTyp;
  VAR Communalities, Uniqueness,
  VarianceAccounted: VectorTyp);
{ Row and column sums of squared loadings }

PROCEDURE WriteLoadings(CONST Loadings: MatrixTyp;
  CONST Communalities, Uniqueness,
  VarianceAccounted: VectorTyp;
  Rotated: BOOLEAN);
{ Writes loadings into a csv-file. File name will depend on whether the
  loadings have been rotated. }

IMPLEMENTATION

VAR
  Significance: SignificanceType;


PROCEDURE ReadCSV(VAR Types: TypeVector; VAR Data: MatrixTyp);

VAR
  InputFile: TEXT;
  c: CHAR;
  s: STRING;
  Variables, Cases, Error, i, j: WORD;
  Interim: MatrixTyp;
  IntTypes: TypeVector;
  x: double;

BEGIN
  Variables := 0;
  Cases := 0;
  Assign(InputFile, ProblemName + '.csv');
  Reset(InputFile);
  ReadLn(InputFile);            // ignore Line WITH variable names
  WHILE NOT EoLn(InputFile) DO  // Read data types from second Line, count Variables
  BEGIN
    INC(Variables);
    s := '';
    REPEAT
      Read(InputFile, c);
      IF c <> SepChar THEN
        S := s + c;
    UNTIL ((c = SepChar) OR EoLn(InputFile));
    IF s = 'binary'
      THEN IntTypes[Variables] := binary
      ELSE IF s = 'ordinal'
             THEN IntTypes[Variables] := ordinal
             ELSE IF s = 'interval'
                    THEN IntTypes[Variables] := interval
                    ELSE IF s = 'rational'
                           THEN IntTypes[Variables] := rational
                           ELSE IF s = 'nominal'
                                  THEN IntTypes[Variables] := nominal
                                  ELSE
                                    BEGIN              // abort PROGRAM WITH error message
                                      Write('unknown data type ', s, ' at n=', Variables: 3, '. Press <CR>:');
                                      ReadLn;
                                      HALT;
                                    END; { else }
  END; { while }
  ReadLn(InputFile);            // go TO next Line
  CreateMatrix(Interim, MaxVariables, Variables, 0.0);
  // Note: number OF persons NOT yet known
  WHILE NOT EoF(InputFile) DO   // now Read data from following lines AND count cases
    BEGIN
      INC(Cases);
      FOR i := 1 TO Variables DO
        BEGIN
          x := ReadFloat(InputFile);
          IF MathError
            THEN
              BEGIN
                MathError := FALSE;
                HALT;
              END;
          SetMatrixElement(Interim, Cases, i, x);

        END;
      ReadLn(InputFile);
    END; { while }
  Close(InputFile);
  CreateMatrix(Data, Cases, Pred(Variables), 0.0);
  // now generate correctly sized data matrix
  FOR i := 1 TO Cases DO
    FOR j := 2 TO Variables DO
      SetMatrixElement(Data, i, Pred(j), GetMatrixElement(Interim, i, j));
  DestroyMatrix(Interim);                         // destroy the interim matrix
  FOR j := 2 TO Variables DO
    // generate Types vector without 1ST column (CASE numbers)
    Types[Pred(j)] := IntTypes[j];
  Writeln(Cases: 4, ' cases with ', Variables: 3, ' variables read, ');
END;


PROCEDURE ReadCorrelations(VAR Cor: MatrixTyp);
// Read correlation matrix from CSV FILE

VAR
  InputFile: TEXT;
  c: CHAR;
  s: STRING;
  Variables, Error, i, j: WORD;
  x: double;

BEGIN
  Assign(InputFile, Problemname + '-NativeCorr.csv');
  Reset(InputFile);
  Variables := 0;
  WHILE NOT EoLn(InputFile) DO  // Read first Line WITH variable numbers, count Variables
    BEGIN   // nothing needs TO be done with these values
      INC(Variables);
      i := ReadInt(InputFile) ;
    END; { while }
  ReadLn(InputFile);            // go TO next Line
  DEC(Variables);               // Line number
  CreateMatrix(Cor, Variables, Variables, 0.0);
  FOR i := 1 TO Variables DO    // now Read data from following lines
    BEGIN
      j := ReadInt(InputFile);  // Line number
      FOR j := 1 TO Variables DO
        SetMatrixElement(Cor, i, j, ReadFloat(InputFile));   // values
      ReadLn(InputFile);
    END; { for i }
  Close(InputFile);
  Writeln(Variables: 3, ' variables read, ');
END;


PROCEDURE AnalyseFrequencies(CONST Cor: MatrixTyp; CONST Types: TypeVector;
  VAR Freqs: FreqsType);

CONST
  TypeText: ARRAY [DataTypes] OF STRING[8] =
    ('binary', 'nominal', 'ordinal', 'interval', 'rational');

VAR
  i, j: DataTypes;
  l, m, Variables: WORD;
  k: INTEGER;
  OutFile: TEXT;
  x: double;

BEGIN
  Assign(OutFile, ProblemName + '-Frequencies.csv');
  Rewrite(OutFile);
  Variables := MatrixRows(Cor);
  FOR i := LOW(DataTypes) TO High(DataTypes) DO
    FOR j := LOW(DataTypes) TO High(DataTypes) DO
      FOR k := -Border TO Border DO
        Freqs[i, j, k] := 0;
  FOR l := 1 TO Variables DO
    FOR m := 1 TO Variables DO
      BEGIN
        x := GetMatrixElement(Cor, l, m) * Border;
        INC(Freqs[Types[l], Types[m], Round(x)]);
      END;
  FOR i := LOW(DataTypes) TO High(DataTypes) DO
    BEGIN
      Writeln(OutFile, TypeText[i]: 8, SepChar);
      FOR j := LOW(DataTypes) TO High(DataTypes) DO
        BEGIN
          Write(OutFile, SepChar, TypeText[j]: 8, SepChar);
          FOR k := -Border TO Border DO
            Write(OutFile, Freqs[i, j, k]: 4, SepChar);
          Writeln(OutFile);
        END;
    END;
  Close(OutFile);
END;


PROCEDURE CalculateCorrelationMatrix(CONST Data: MatrixTyp;  CONST Types:
          TypeVector; VAR Cor: MatrixTyp);

VAR
  i, j, Variables: WORD;
  IVector, JVector: VectorTyp;
  Rxy: double;
  Cont: ContTable;

BEGIN
  Variables := MatrixColumns(Data);
  CreateIdentityMatrix(Cor, Variables);
  GetColumn(Data, i, IVector);
  for i := 1 to Variables do
    begin
      SetMatrixElement(Cor, i, i, 1.0);    // diagonal
      GetColumn(Data, i, IVector);
      FOR j := Succ(i) TO Variables DO
        BEGIN
          GetColumn(Data, j, JVector);
          CASE Types[i] OF
            nominal: CASE Types[j] OF
                      nominal:  BEGIN
                                  Contingency(IVector, JVector, Cont);
                                  Rxy := lambda(Cont);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                  DestroyContingency(Cont);
                                END;
                      ordinal:  BEGIN
                                  Contingency(IVector, JVector, Cont);
                                  Rxy := theta(Cont);
                                  Rxy := Sqrt(Rxy);
                                  // asymmetric: only positive values!
                                  DestroyContingency(Cont);
                                END;
                      binary:   BEGIN
                                  Contingency(IVector, JVector, Cont);
                                  Rxy := lambda(Cont);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                  DestroyContingency(Cont);
                                END;
                      interval,
                      rational: BEGIN
                                  Rxy := eta_sqr(IVector, JVector, Significance);
                                  Rxy := Sqrt(Rxy);
                                  // asymmetric: only positive values!
                                END;
                     END;
            ordinal: CASE Types[j] OF
                      nominal:  BEGIN
                                  Contingency(JVector, IVector, Cont);
                                  Rxy := theta(Cont);
                                  Rxy := Sqrt(Rxy);
                                  DestroyContingency(Cont);
                                END;
                      ordinal:  BEGIN
                                  Rxy :=
                                    OrdinalCorrelations(IVector, JVector, 'E', Significance);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                END;
                      binary:   BEGIN
                                  Rxy :=
                                    OrdinalCorrelations(IVector, JVector, 'E', Significance);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                END;
                      interval,
                      rational: Rxy := SpearmanRankCorrelation(IVector, JVector, Significance);
                    END;
            binary: CASE Types[j] OF
                      nominal:  BEGIN
                                  Contingency(IVector, JVector, Cont);
                                  Rxy := lambda(Cont);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                  DestroyContingency(Cont);
                                END;
                      ordinal:  BEGIN
                                  Rxy := OrdinalCorrelations(IVector, JVector, 'E', Significance);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                END;
                      binary:   BEGIN
                                  Contingency(IVector, JVector, Cont);
                                  Rxy := lambda(Cont);
                                  Rxy := sign(Rxy) * Sqrt(Abs(Rxy));
                                  DestroyContingency(Cont);
                                END;
                      interval,
                      rational:  Rxy := PointBiserialCorrelation(IVector, JVector, Significance);
              END;
          interval,
          rational: CASE Types[j] OF
                      nominal:  BEGIN
                                  Rxy := eta_sqr(JVector, IVector, Significance);
                                  Rxy := Sqrt(Rxy);
                                END;
                      ordinal: Rxy := SpearmanRankCorrelation(IVector, JVector, Significance);
                      binary:  Rxy := PointBiserialCorrelation(JVector, IVector, Significance);
                      interval,
                      rational: Rxy := PearsonProductMomentCorrelation(IVector, JVector, Significance);
                    END;
          END; { case Types[i] }
          SetMatrixElement(Cor, i, j, Rxy);    // upper half
          SetMatrixElement(Cor, j, i, Rxy);    // lower half
          Writeln(i:4, ' ', j:4, ' ', Rxy:3:3);
          DestroyVector(JVector);
        END;  { for j }
      DestroyVector(IVector);
    END; { for i }
  Write('Correlation matrix calculated, ');
END;


PROCEDURE WriteCorrelations(CONST Cor: MatrixTyp; Shrunk: BOOLEAN);

VAR
  i, j: WORD;
  OutFile: TEXT;

BEGIN
  IF Shrunk
    THEN Assign(OutFile, ProblemName + '-ShrunkCorr.csv')
    ELSE Assign(OutFile, Problemname + '-NativeCorr.csv');
  Rewrite(OutFile);
  Write(OutFile, SepChar);
  FOR j := 1 TO Cor^.Columns DO
    Write(OutFile, j: 4, SepChar); // LABEL columns
  Writeln(OutFile);
  FOR i := 1 TO Cor^.Rows DO // the matrix itself
    BEGIN
      Write(OutFile, i: 4, SepChar); // LABEL row
      FOR j := 1 TO Cor^.Columns DO
        Write(OutFile, GetMatrixElement(Cor, i, j):7:4, ';');
      Writeln(OutFile);
    END;  { for i }
  Close(OutFile);
  Writeln('Correlations written to file, ');
END; { WriteCorrelations}


FUNCTION Bartlett(CONST Cor: MatrixTyp; Cases: WORD): double;

VAR
  f, Variables: WORD;
  chi2, Det, P0: extended;

BEGIN
  Variables := MatrixRows(Cor);
  Det := Determinante(Cor);
  Write('Determinante = ', det: 10, ' ');
  IF (Det < Zero)                           // singular matrix
    THEN
      BEGIN
        Result := NaN;
        Writeln('Bartlet = NaN');
        EXIT;
      END;
  chi2 := -(Pred(Cases) - (2 * Variables + 5) / 6) * log(Det, 10);
  f := Round(Variables * Pred(Variables) / 2);
  Write('chi2 := ', chi2: 5: 3, ' ', 'f = ', f: 5, ' ');
  P0 := IntegralChi(chi2, f);
  Writeln('P0 = ', P0: 1: 4);
  Result := P0;
END;


PROCEDURE ExplainedVariance(CONST EigenValues: VectorTyp;
  VAR Acc, Variance, CumVariance: VectorTyp);

VAR
  i, Variables: WORD;
  x, Sum, Cummulative: double;
  OutFile: TEXT;

BEGIN
  Variables := VectorLength(EigenValues);
  CreateVector(Variance, Variables, 0.0);
  CreateVector(CumVariance, Variables, 0.0);
  CreateVector(Acc, Variables, 0.0);
  Sum := 0;
  i := 1;
  FOR i := 1 TO Variables DO
    Sum := Sum + GetVectorElement(EigenValues, i);
  Cummulative := 0;
  FOR i := 2 TO Pred(Variables) DO //  second derivative OF the scree-curve by finite differences
    SetVectorElement(Acc, i, GetVectorElement(EigenValues, Succ(i)) -
       2*GetVectorElement(EigenValues, i) +  GetVectorElement(EigenValues, Pred(i)));
  SetVectorElement(Acc, 1, NaN);
  SetVectorElement(Acc, Variables, NaN);
  Assign(OutFile, ProblemName + '-EigenValues.csv');
  Rewrite(OutFile);
  FOR i := 1 TO Variables DO
    BEGIN
      x := GetVectorElement(EigenValues, i) / Sum;
      Cummulative := Cummulative + x;
      SetVectorElement(Variance, i, x);
      SetVectorElement(CumVariance, i, Cummulative);
      Write(OutFile, i: 3, SepChar, GetVectorElement(EigenValues, i):3:3, SepChar,
        x:3:3, SepChar, Cummulative:3:3);
      Writeln(OutFile, SepChar, GetVectorElement(Acc, i):3:3);
    END;
  Close(OutFile);
END;


PROCEDURE WriteEigenVector(CONST EigenVectors: MatrixTyp);

VAR
  i, j: WORD;
  OutFile: TEXT;

BEGIN
  Assign(OutFile, ProblemName + '-EigenVectors.csv');
  Rewrite(OutFile);
  FOR i := 1 TO MatrixRows(EigenVectors) DO
    BEGIN
      Write(OutFile, i: 3, SepChar);
      FOR j := 1 TO SigVecs DO
        Write(OutFile, GetMatrixElement(EigenVectors, i, j): 3: 3, SepChar);
      Writeln(OutFile);
    END;
  Close(OutFile);
END;


PROCEDURE MaximumLikelihood(CONST Data: MatrixTyp; CONST Types: TypeVector;
  VAR ImputVector: VectorTyp);

VAR
  i, j, n, p, s, iMax: WORD;
  JVector: VectorTyp;
  x: double;
  Numbers: ARRAY [1..MaxSteps] OF WORD;

BEGIN
  n := MatrixRows(Data);
  p := MatrixColumns(Data);
  CreateVector(ImputVector, p, 0.0);
  FOR j := 1 TO p DO
    BEGIN
      GetColumn(Data, p, JVector);
      CASE Types[j] OF
        binary,
        nominal,
        ordinal:  BEGIN
                    FOR i := 1 TO MaxSteps DO Numbers[i] := 0;
                    FOR i := 1 TO n DO
                      BEGIN
                        x := GetVectorElement(JVector, i);
                        IF IsNaN(x)
                          THEN
                          ELSE INC(Numbers[Round(x)]);
                      END;
                    x := 0;
                    iMax := 0;
                    FOR i := 1 TO MaxSteps DO
                      IF (Numbers[i] > x)
                        THEN
                          BEGIN
                            x := Numbers[i];
                            iMax := i;
                          END;
                    SetVectorElement(ImputVector, j, iMax); // most common element
                  END;
        interval,
        rational: BEGIN
                    x := NeumaierSum(JVector)/ActualElements(JVector); //arithmetic mean
                    SetVectorElement(ImputVector, j, x);
                  END;
      END;  // case
      DestroyVector(JVector);
    END;
END;

PROCEDURE RobustProduct(CONST A, B: MatrixTyp; CONST ImputVector: VectorTyp;
  VAR C: MatrixTyp);

VAR
  i, j, k: WORD;
  Sum: double;

BEGIN
  IF MatrixColumns(A) <> MatrixRows(B)
    THEN
      BEGIN
        Write(' Matrix multiplication: A.Columns <> B.Rows');
        ReadLn;
        EXIT;
      END;
  IF MatrixColumns(B) < SigVecs
    THEN
      BEGIN
        Write(' Matrix multiplication: Number of Eigenvectors < SigVecs ');
        ReadLn;
        EXIT;
      END;
  CreateMatrix(C, MatrixRows(A), SigVecs, 0);
  FOR i := 1 TO MatrixRows(A) DO
    FOR j := 1 TO SigVecs DO
      BEGIN
        Sum := 0;
        FOR k := 1 TO MatrixColumns(A) DO
          IF IsNaN(GetMatrixElement(A, i, k))
            THEN Sum := Sum + GetVectorElement(ImputVector, k) *
                         GetMatrixElement(B, k, j)
          //  replace NaN WITH max likelyhood estimator
            ELSE Sum := Sum + GetMatrixElement(A, i, k) * GetMatrixElement(B, k, j);
        SetMatrixElement(C, i, j, Sum);
      END;
END;

PROCEDURE WriteComponentScores(CONST Scores: MatrixTyp);

VAR
  i, j: WORD;
  OutFile: TEXT;

BEGIN
  Assign(OutFile, ProblemName + '-Scores.csv');
  Rewrite(OutFile);
  FOR i := 1 TO MatrixRows(Scores) DO
    BEGIN
      FOR j := 1 TO SigVecs DO
        Write(OutFile, GetMatrixElement(Scores, i, j): 3: 3, SepChar);
      Writeln(OutFile);
    END;
  Close(OutFile);
END;

PROCEDURE CalculateLoadings(CONST Data, Scores: MatrixTyp;
  CONST Types: TypeVector; VAR Loadings: MatrixTyp);

VAR
  i, j, p : WORD;
  Rxy: double;
  DataVector, ScoreVector: VectorTyp;
  Significance: SignificanceType;

BEGIN
  p := MatrixColumns(Data);
  CreateMatrix(Loadings, p, SigVecs, 0.0);
  FOR i := 1 TO p DO                           // over all variables
    BEGIN
      GetColumn(Data, i, DataVector);
      FOR j := 1 TO SigVecs DO                  // over all Scores
        BEGIN
          GetColumn(Scores, j, ScoreVector);      // JVector IS always rational
          CASE Types[i] OF
            binary:   Rxy := PointBiserialCorrelation(DataVector,
                ScoreVector, Significance);
            nominal:  Rxy := Sqrt(eta_sqr(DataVector, ScoreVector, Significance));
            // by definition always positive
            ordinal:  Rxy := SpearmanRankCorrelation(DataVector,
                ScoreVector, Significance);
            interval,
            rational: Rxy := PearsonProductMomentCorrelation(DataVector, ScoreVector, Significance);
          END;
          SetMatrixElement(Loadings, i, j, Rxy);
          DestroyVector(ScoreVector);
        END;
      DestroyVector(DataVector);
    END;
END;


PROCEDURE CalculateCommunalities(CONST Loadings: MatrixTyp;
            VAR Communalities, Uniqueness, VarianceAccounted: VectorTyp);

VAR
  i, j, Rows, Columns: WORD;
  Sum1, Sum2, x, y: double;

BEGIN
  Rows := MatrixRows(Loadings);
  Columns := MatrixColumns(Loadings);
  CreateVector(Communalities, Rows, 0);
  FOR i := 1 TO Rows DO
    BEGIN
      Sum1 := 0;
      FOR j := 1 TO SigVecs DO
        Sum1 := Sum1 + Sqr(GetMatrixElement(Loadings, i, j));
      // row sum OF squared loadings
      SetVectorElement(Communalities, i, Sum1);
    END;
  CreateVector(VarianceAccounted, Columns + 2, 0.0);
  FOR j := 1 TO Columns DO
    BEGIN
      Sum1 := 0;
      FOR i := 1 TO Rows DO
        Sum1 := Sum1 + Sqr(GetMatrixElement(Loadings, i, j));
      // column sum OF squared loadings
      SetVectorElement(VarianceAccounted, j, Sum1);
    END;
  CreateVector(Uniqueness, Rows, 0);
  Sum1 := 0;
  Sum2 := 0;
  FOR i := 1 TO Rows DO
    BEGIN
      x := GetVectorElement(Communalities, i);
      y := 1 - x;
      SetVectorElement(Uniqueness, i, y);
      Sum1 := Sum1 + x;                    // explained
      Sum2 := Sum2 + y;                    // unexplained
    END;
  SetVectorElement(VarianceAccounted, Columns + 1, Sum1);
  SetVectorElement(VarianceAccounted, Columns + 2, Sum2);
END;


PROCEDURE WriteLoadings(CONST Loadings: MatrixTyp; CONST Communalities, Uniqueness,
                        VarianceAccounted: VectorTyp; Rotated: BOOLEAN);

VAR
  n, i, j: WORD;
  OutFile: TEXT;
  Variance: double;

BEGIN
  IF Rotated
    THEN Assign(OutFile, ProblemName + '-RotLoad.csv')
    ELSE Assign(OutFile, ProblemName + '-Loadings.csv');
  Rewrite(OutFile);
  n := MatrixRows(Loadings);
  FOR i := 1 TO n DO
    BEGIN
      FOR j := 1 TO SigVecs DO
        Write(OutFile, GetMatrixElement(Loadings, i, j):3:3, SepChar);
      Writeln(OutFile, SepChar, GetVectorElement(Communalities, i):3:3, SepChar,
        GetVectorElement(Uniqueness, i):3:3, SepChar);
    END;
  Writeln(OutFile);
  FOR j := 1 TO SigVecs DO
    Write(OutFile, GetVectorElement(VarianceAccounted, j):3:3, SepChar);
  Writeln(OutFile, GetVectorElement(VarianceAccounted, SigVecs + 1):3:3, SepChar,
    GetVectorElement(VarianceAccounted, SigVecs + 2):3:3, SepChar);
  Variance := GetVectorElement(VarianceAccounted, SigVecs + 1) +
    GetVectorElement(VarianceAccounted, SigVecs + 2);
  FOR j := 1 TO SigVecs DO
    Write(OutFile, GetVectorElement(VarianceAccounted, j) / Variance:3:3, SepChar);
  Writeln(OutFile, GetVectorElement(VarianceAccounted, SigVecs + 1) / Variance:3:3, SepChar,
    GetVectorElement(VarianceAccounted, SigVecs + 2) / Variance:3:3, SepChar);
  Close(OutFile);
END;


END.  // PCA

