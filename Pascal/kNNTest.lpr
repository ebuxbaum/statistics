PROGRAM kNNTest;

USES
  Math,            // free pascal standard math UNIT
  MathFunc,        // basic math routines
  Vector,          // vector algebra
  Matrix,          // matrix algebra
  Zufall,          // Random numbers
  Deskript,        // descriptive statistics
  kNN,             // k means clustering
  CrossValidation  // k-fold cross-validation
  ;

CONST
  n = 112;     // data sets
  p = 23;      // variables

VAR
  Data, ValidationResult, TestData, Centroids : MatrixTyp;
  SplitData                                   : SplitDataTyp;
  i, j, k, l, kValidate, kmeans, done, Count  : WORD;
  c                                           : CHAR;
  Distance                                    : float;
  Group                                       : GroupTyp;
  WSS                                         : ARRAY[1..kNNmax] OF float;

  PROCEDURE ReadCSV(n, p: WORD; FileName: STRING; VAR Data: MatrixTyp);
  // Read correlation matrix from CSV FILE

  VAR
    InputFile : TEXT;
    i, j      : WORD;
    x         : float;
    c         : CHAR;

  BEGIN
    Assign(InputFile, Filename);
    Reset(InputFile);
    IF IOResult <> 0
      THEN
        BEGIN
          c := WriteErrorMessage('Unable to open file, Press <CR>');
          HALT;
        END;
    ReadLn(InputFile);            // ignore first Line WITH headers
    ReadLn(InputFile);            // ignore second Line WITH types
    CreateMatrix(Data, n, p, 0.0);
    IF MatrixError
      THEN
        BEGIN
          c := WriteErrorMessage('program terminated');
          HALT;
        END;
    FOR i := 1 TO n DO            // now Read data from following lines
      BEGIN
        FOR j := 1 TO p DO
          BEGIN
            x := ReadFloat(InputFile);
            IF MathError
              THEN
                BEGIN
                  c := WriteErrorMessage('Unable to read datum, press <CR>');
                  HALT;
                END;
            SetMatrixElement(Data, i, j, x);
          END; { for j }
        ReadLn(InputFile);
      END; { for i }
    Close(InputFile);
  END; { ReadCSV }


BEGIN
  k := 2;
  ReadCSV(n, p, 'All.csv', Data);
  ShellSortMatrix(Data, 1);
  // ensure that data are sorted by study no
  kValidate := floor(Sqrt(n));                   // groups FOR k-fold cross-validation
  IF kValidate > MaxK THEN kValidate := MaxK;
  SplitDataMatrix(Data, kValidate, SplitData);
  IF CrossValidationError
    THEN
      BEGIN
        c := WriteErrorMessage('program terminated');
        HALT;
      END;
  CreateMatrix(ValidationResult, n, 3, 0.0);
  IF MatrixError
    THEN
      BEGIN
        c := WriteErrorMessage('program terminated');
        HALT;
      END;
  Count := 1;
  FOR i := 1 TO kvalidate DO
    BEGIN
      CreateTestData(SplitData, i, kValidate, TestData);  // create a data vector
      Writeln('Test data created ', i: 3);
      IF CrossValidationError
        THEN
          BEGIN
            c := WriteErrorMessage('program terminated');
            HALT;
          END;
      Distance := kMeansCluster(TestData, k, Centroids, Group);
      AssignTestData(SplitData[i], Centroids, Count, kValidate, ValidationResult);
      WSS[i] := WithinSumOfSquares(ValidationResult);
    END;
END.

