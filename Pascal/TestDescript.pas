program TestDescript;

uses math, MathFunc, Vector, Matrix, Zufall, Deskript;

const ProbSize = 50;
      Vars     = 10;

var Dm, Dr : VectorTyp;
    i, j : word;
    MData: MatrixTyp;

begin
  CreateMatrix(MData, ProbSize, Vars, 0.0);
  for i := 1 to ProbSize do
    begin
      if RandomLinear > 0.1  // 10% outliers
        then
          for j := 1 to Vars do
            SetMatrixElement(MData, i, j, 100 + RandomNormal(10-j, 2*j))
        else
          for j := 1 to Vars do
            SetMatrixElement(MData, i, j, 100 + RandomNormal(10+j, 3*j));
    end;
  writeln(1);
  RobustDistance (MData, Dr);
  writeln(2);
  MahalanobisDistance (MData, Dm);
  writeln(3);
  for i := 1 to ProbSize do
    writeln(i:3, ' ', FloatStr(GetVectorElement(Dm, i), 10), ' ', FloatStr(GetVectorElement(Dr, i), 10));

  readln;
  DestroyVector(Dr);
  DestroyVector(Dm);
  DestroyMatrix(MData);

end.

