program Simplex;


USES
  MathFunc,            // basic math functions
  Vector,              // vector arithmetic
  Matrix,              // matrix arithmetic
  Zufall,              // random numbers
  SimplexFit           // curve fit by simplex
  ;

const MaxData = 100;

VAR Data                              : MatrixTyp;
    Calculated                        : VectorTyp;
    ProblemName, Formel, xName, yName : STRING;

{ *********************************************************************** }

PROCEDURE CreateDataSet (VAR Data : MatrixTyp);

VAR i : WORD;
    x, y : float;

BEGIN
  FOR i := 1 TO MatrixRows(Data) DO
    BEGIN
      x := i/10 ;
      SetMatrixElement(Data, i, 1, x);
      y := x / (1+x);                  // normalised Henri-Michaelis-Menten law
      y := y + RandomNormal(0, 0.1);   // add normal-distributed random noise
      SetMatrixElement(Data, i, 2, y);
    END;
END;


BEGIN {Hauptprogram}
  CreateMatrix(Data, MaxData, 2, 0.0);
  CreateVector(Calculated, MaxData, 0.0);
  CreateDataSet(Data);
  Approximation(Data, Calculated, ProblemName, Formel, xName, yName);
  DestroyMatrix(Data);
  DestroyVector(Calculated);
END.
