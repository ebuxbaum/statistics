PROGRAM TestDescript;

USES math, MathFunc, Vector, Matrix, Zufall, Deskript;

CONST ProbSize = 100;
      Vars     =   2;

VAR Data, Weights, Dm, Dr       : VectorTyp;
    i, j                        : WORD;
    MinEle, MaxEle, sum, mean,
    std, q1, q2, q3, e2, e3, e4 : float;
    MData, DataWithWeights      : MatrixTyp;

BEGIN
  CreateVector(Data, ProbSize, 0.0);
  CreateVector(Weights, ProbSize, 1/ProbSize);
  FOR i := 1 TO ProbSize DO
    SetVectorElement(Data, i, 100 + RandomNormal(0, 10));

  Mean := ArithmeticMean(Data);
  Writeln('Arithmetic mean: ', FloatStr(Mean, 10));
  Writeln('Geometric mean: ', FloatStr(GeometricMean(Data), 10));
  Writeln('Harmonic mean: ', FloatStr(HarmonicMean(Data), 10));
  Writeln;
  Writeln('Power mean, k=1: ', FloatStr(GeneralMean(Data, Weights, 1), 10));
  Writeln('Power mean, k=0.001: ', FloatStr(GeneralMean(Data, Weights, 0.001), 10));
  Writeln('Power mean, k=-1: ', FloatStr(GeneralMean(Data, Weights, -1), 10));
  Writeln;
  std := StandardDeviation(Data);
  Writeln('Standard deviation', FloatStr(std, 10));
  Writeln('excess kurtosis', FloatStr(ExcessKurtosis(Data, mean, std), 10));
  Writeln('Skewness', FloatStr(Skewness(Data, mean, std), 10));
  Writeln;
  Writeln('Gini coefficient: ', FloatStr(Gini(Data, Mean), 10));
  Writeln('Herfindahl-Index: ', FloatStr(HerfindahlIndex(Data), 10));
  Writeln;
  q2 := Median(Data);
  q1 := Quantile(Data, 0.25);
  q3 := Quantile(Data, 0.75);
  Writeln('Median: ', FloatStr(q2, 10));
  Writeln('Trimedian: ', FloatStr(Trimedian(Data), 10)); // data have been sorted by Median
  Writeln('Inter-quartile distance: ', FloatStr(InterQuantilDistance(q1, q3), 10));
  Writeln('MAD: ', FloatStr(MAD(Data), 10));
  Writeln('std. error of median: ', FloatStr(StandardErrorOfMedian(Data), 10));
  Writeln;
  Writeln('Naive Hodges-Lehmann: ', FloatStr(NaiveHodgesLehmann(Data), 10),
  ' Hodges-Lehmann estimator: ', FloatStr(HodgesLehmann(Data), 10)); // data have been sorted by Median
  Writeln('Naive Sn: ', FloatStr(NaiveSn(Data), 10), ' Sn: ', FloatStr(Sn(Data), 10));
  Writeln('Naive Qn: ', FloatStr(NaiveQn(Data), 10));
//  Writeln('Qn: ', FloatStr(Qn(Data), 10));

  FOR i := 1 TO ProbSize DO SetVectorElement(Weights, i, RandomNormal(1/ProbSize, 0.01));
  MaxEle := FindLargest(Weights);
  MinEle := FindSmallest(Weights);
  Scale(Weights, MinEle, MaxEle); // weights IN 0..1
  Sum := TotalSum(Weights);
  FOR i := 1 TO ProbSize DO
    SetVectorElement(Weights, i, GetVectorElement(Weights, i)/Sum); // total weight = 1
  CreateMatrix(DataWithWeights, ProbSize, 2, 0.0);
  SetColumn(DataWithWeights, Data, 1);
  SetColumn(DataWithWeights, Weights, 2);
  Writeln;
  Write('Weighted median: ', FloatStr(WeightedMedian(DataWithWeights), 10));
  Writeln(', weighted LoMed:  ', FloatStr(WeightedLoMed(DataWithWeights), 10),
          ', weighted HiMed:  ', FloatStr(WeightedHiMed(DataWithWeights), 10));
  Writeln;
  e2 := ell2(Data);
  e3 := ell3(Data);
  e4 := ell4(Data);
  Writeln('ell2: ', FloatStr(e2, 10), ', ell3: ', FloatStr(e3, 10), ', ell4: ', FloatStr(e4, 10));
  Writeln('tau2: ', FloatStr(e2/mean, 10), ', tau3: ', FloatStr(e3/e2, 10), ', tau4: ', FloatStr(e4/e2, 10));
  Writeln('Ouartile coefficient of skewness: ', FloatStr(QuartileCoefficientOfSkewness(q1, q2, q3), 10));
  Writeln('Centile coefficient of kurtosis: ', FloatStr(CentilCoeffKurtosis(Data), 10));
  ReadLn;
  DestroyVector(Data);
  DestroyMatrix(DataWithWeights);

  CreateMatrix(MData, ProbSize, Vars, 0.0);
  FOR i := 1 TO (ProbSize DIV 10) DO              // 10% outliers
    BEGIN
      SetMatrixElement(MData, i, 1, 115 + RandomNormal(0, 3));
      SetMatrixElement(MData, i, 2, 110 + RandomNormal(0, 3));
    END;
  FOR i := Succ(ProbSize DIV 10) TO ProbSize DO   // 90% inliers
    BEGIN
      SetMatrixElement(MData, i, 1, 110 + RandomNormal(0, 3));
      SetMatrixElement(MData, i, 2, 115 + RandomNormal(0, 3));
    END;
  RobustDistance (MData, Dr);
  MahalanobisDistance (MData, Dm);
  FOR i := 1 TO ProbSize DO
    Writeln(i:3, ' ', FloatStr(GetMatrixElement(MData, i, 1), 11), ' ',
    FloatStr(GetMatrixElement(MData, i, 2), 11), ' ',
    FloatStr(GetVectorElement(Dm, i), 10), ' ', FloatStr(GetVectorElement(Dr, i), 10));

  ReadLn;
  DestroyVector(Dr);
  DestroyVector(Dm);
  DestroyMatrix(MData);
END.
