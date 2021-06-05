UNIT Circular;

INTERFACE

USES math,                  // Free Pascal math UNIT
     crt,                   // Free pascal UNIT
     MathFunc,              // basic math functions
     Dynam,                 // dynamic data structures
     Complex,               // Complex numbers
     Vector,                // vector arithmetic
     Matrix,                // matrix algebra
     Correlations,          // correlation coefficients
     Deskript,              // descriptive statistics
     Stat,                  // statistical distributions
     Zufall;                // pseudo-Random numbers

CONST CircleError : BOOLEAN = FALSE;
      Nintey = Const_pi / 2;          // 90° in rad

{ ****************** description of a single data vector ******************** }

PROCEDURE Transform (VAR Data : VectorTyp; FullCircle : float);

FUNCTION MedianDirection (TransformedData : VectorTyp) : float;

FUNCTION MeanVector (TransformedData : VectorTyp; p : word) : ComplexTyp;

FUNCTION CircularVariance (R : float) : float;

FUNCTION CircularStandardDeviation (R : float) : float;

FUNCTION CircularDispersion (TransformedData : VectorTyp; Mean : ComplexTyp) : float;

FUNCTION Kappa (R : float; n : WORD) : float;

FUNCTION TrigonometricMoment (Data : VectorTyp; MeanAngle : float; p : WORD) : ComplexTyp;

function CenteredMean(Moment : ComplexTyp) : ComplexTyp;

FUNCTION CircularSkew (TransformedData : VectorTyp; Mean : ComplexTyp) : float;

FUNCTION CenteredCircularSkew (Mean1, Mean2 : ComplexTyp) : float;

FUNCTION CircularKurtosis (TransformedData : VectorTyp; Mean : ComplexTyp) : float;

FUNCTION CenteredCircularKurtosis (Mean1, Mean2 : ComplexTyp) : float;

FUNCTION ConfidenceInterval (R, delta : float; n : WORD) : float;

{ ****************************** random numbers ***************************** }

FUNCTION RandomVonMieses(mu, kappa : float) : float;

FUNCTION RandomUniformCircular : float;

{ *********************** Tests for preferred direction ********************* }

FUNCTION Rayleigh (MeanVectorLength : float; n : WORD) : float;

FUNCTION HodgesAjne (TransformedData : VectorTyp; MeanAngle : float) : float;

PROCEDURE ChiSqrTest (Data : VectorTyp; Direction : float;
          Sectors : WORD; VAR ChiSqr : float; VAR dgf : WORD);

FUNCTION HomewardComponent (n : WORD; Mean : ComplexTyp;
         Expected : float) : float;

FUNCTION Rao (TransformedData : VectorTyp) : float;

PROCEDURE Kuipers (Data : MatrixTyp;   VAR K, U : float);

PROCEDURE CalculateCumulativeFrequencies (Data : VectorTyp;
          FullCircle : float; VAR Result : MatrixTyp);

FUNCTION OneSample (MeanAngle, R, delta, TestAngle : float; n : WORD) : BOOLEAN;

{ ******************************** grouped data ***************************** }

FUNCTION MeanVectorGrouped (Transformed : MatrixTyp;
         Difference : float) : ComplexTyp;


{ *********************** Compare two circular distributions **************** }

PROCEDURE Difference (Data1, Data2 : VectorTyp;
                      VAR KuipersV, WatsonsUSqr : float);

FUNCTION FTest (Data1, Data2 : VectorTyp) : float;

FUNCTION Wilcoxon (Data1, Data2 : VectorTyp;
                   Direction1, Direction2 : float) : float;

PROCEDURE WatsonWilliams (TransformedData1, TransformedData2 : VectorTyp);

FUNCTION KruskalWallis (TransformedData1, TransformedData2 : VectorTyp) : float;

{ ************************** linear-bivariate Data ************************** }

PROCEDURE DescribeBivariat (Data : MatrixTyp);

{ *********** linear dependent, circular independent Correlation ************ }

FUNCTION CAssociation (CONST theta, Y : VectorTyp; VAR Sig : SignificanceType) : float;

PROCEDURE CircularLinearCorrelation (const TransformedData, yData: VectorTyp;
                                     VAR Correlation : float; VAR Sig : SignificanceType);

PROCEDURE LinearPeriodicRankCorrelation (const TransformedData, yData: VectorTyp;
                                         VAR Correlation, U : float);

PROCEDURE TrigonometricPolynomial (const TransformedData, yData: VectorTyp;
                                   Periode : float);

{ ********* Correlation circular dependent and independent variable ********* }

PROCEDURE PeriodicRankCorrelation (CONST alpha, beta : VectorTyp;
                                     VAR rPlus, rMinus : float);

PROCEDURE CircularCircularCorrelation (Alpha, Beta : VectorTyp;
          VAR Correlation : float; var Sig : SignificanceType);

{ *************************************************************************** }

IMPLEMENTATION

VAR xMax       : WORD;
    ch         : CHAR;


FUNCTION ModuloTwoPi (Datum : float) : float;

VAR x : float;

BEGIN
  x := Datum;
  WHILE (x > Const_2pi) DO x := x - Const_2pi;
  WHILE (x < 0.0) DO x := x + Const_2pi;
  Result := x;
END;


PROCEDURE Transform (VAR Data : VectorTyp; FullCircle : float);

VAR n, i  : WORD;
    Datum : float;

BEGIN
  n := VectorLength(Data);
  FOR i := 1 TO n DO
     BEGIN
       Datum := GetVectorElement(Data, i) * Const_2pi / FullCircle;
       SetVectorElement(Data, i, ModuloTwoPi(Datum));
     END;
END;

{ ************************************************************************** }

FUNCTION MedianDirection (TransformedData : VectorTyp) : float;

BEGIN
  Result := Median(TransformedData);
END;


FUNCTION MeanVector (TransformedData : VectorTyp; p : word) : ComplexTyp;

VAR C, S, R, theta, Datum : float;
    n, i                  : WORD;

BEGIN
   n  := VectorLength(TransformedData);
   C := 0;
   S := 0;
   FOR i := 1 TO n DO
     BEGIN
       Datum := GetVectorElement(TransformedData, i);
       C := C + Cos(p * Datum);
       S := S + Sin(p * Datum);
     END;
  writeln('Sum       ', FloatStr(S, ValidFigures), '  ', FloatStr(C, ValidFigures));
  R := Sqrt(C*C + S*S) / n;
  CASE signum(C) OF
    0 : theta := Const_pi / 2;    // arccos(C/R) = arccos(0)
   -1 : theta := ArcTan(S/C) + Const_pi;
    1 : CASE signum(S) OF
          1 :  theta := ArcTan(S/C);
         -1 :  theta := ArcTan(S/C) + Const_2pi;
          0 :  theta := 0;
    END;
  END;
  Result := ComplexInit(R, theta); // convert to complex number
END;


FUNCTION CircularVariance (R : float) : float;

BEGIN
  Result := 1 - R;
END;


FUNCTION CircularStandardDeviation (R : float) : float;

BEGIN
  Result := sqrt(-2 * ln(R));
END;


FUNCTION Kappa (R : float; n : WORD) : float;

VAR k : float;

BEGIN
  IF R < 0.53
    THEN
      k := 2*R + pot(R,3) + 5*pot(R, 5)/6
    ELSE
      IF R < 0.85
        THEN
          k := -0.4 + 1.39*R + 0.43/(1-R)
        ELSE
          k := 1 / (pot(R,3) - 4*R*R + 3*R);
  IF n <= 15  // small sample correction
    THEN
      IF k < 2
        THEN k := max((k - 2/(n*k)), 0)
        ELSE k := pot(Pred(n), 3) * k / (n*n*n + n);
  Result := k;
END;


FUNCTION TrigonometricMoment (Data : VectorTyp; MeanAngle : float; p : WORD) : ComplexTyp;

VAR n, i : WORD;
    SumCos, SumSin, x : float;

BEGIN
  n := VectorLength(Data);
  SumCos := 0;
  SumSin := 0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(Data, i);
      SumSin := SumSin + Sin(p * (x - MeanAngle));
      SumCos := SumCos + Cos(p * (x - MeanAngle));
    END;
  Result := ComplexInit(SumCos/n, SumSin/n);
END;


function CenteredMean (Moment : ComplexTyp) : ComplexTyp;

var S, C, R, theta : float;

BEGIN
  C := Re(Moment);
  S := Im(Moment);
  R := sqrt(C*C + S*S);
  CASE signum(C) OF
    0 : theta := Const_pi / 2;    // arccos(C/R) = arccos(0)
   -1 : theta := ArcTan(S/C) + Const_pi;
    1 : CASE signum(S) OF
          1 :  theta := ArcTan(S/C);
         -1 :  theta := ArcTan(S/C) + Const_2pi;
          0 :  theta := 0;
    END;
  END;
  Result := ComplexInit(R, theta); // convert to complex number
END;


FUNCTION CircularDispersion (TransformedData : VectorTyp; Mean : ComplexTyp) : float;

var n, i          : word;
    theta, m1, m2 : float;

BEGIN
  n := VectorLength(TransformedData);
  theta := Im(Mean);
  m1 := Re(TrigonometricMoment(TransformedData, theta, 1));
  m2 := Re(TrigonometricMoment(TransformedData, theta, 2));
 Result := (1 - m2) / (2 * m1 * m1);
END;

FUNCTION CircularSkew (TransformedData : VectorTyp; Mean : ComplexTyp) : float;

VAR n, i          : WORD;
    theta, SumSin : float;

BEGIN
  n := VectorLength(TransformedData);
  theta := Im(Mean); // mean angle
  SumSin := 0.0;
  FOR i := 1 TO n DO
    SumSin := SumSin + sin(2 * (GetVectorElement(TransformedData, i) - theta));
  Result := SumSin / n;
END;


FUNCTION CenteredCircularSkew (Mean1, Mean2 : ComplexTyp) : float;

VAR theta1, theta2, R1, R2 : float;

BEGIN
  theta1 := Im(Mean1);
  R1 := Re(Mean1);
  theta2 := Im(Mean2);
  R2 := Re(Mean2);
  Result := pot(1-R1, 2/3);
  IF Result = 0
    THEN
      BEGIN
        CircleError := true;
        WriteLn('Error: Circular skew when resultant length is 1');
        EXIT;
      END;
  Result := R2 * sin(theta2 - 2*theta1) / Result;
END;


FUNCTION CircularKurtosis (TransformedData : VectorTyp; Mean : ComplexTyp) : float;

VAR n, i          : WORD;
    theta, SumCos : float;

BEGIN
  n := VectorLength(TransformedData);
  theta := Im(Mean); // mean angle
  SumCos := 0.0;
  FOR i := 1 TO n DO
    SumCos := SumCos + cos(2 * (GetVectorElement(TransformedData, i) - theta));
  Result := SumCos / n;
END;


FUNCTION CenteredCircularKurtosis (Mean1, Mean2 : ComplexTyp) : float;

VAR theta1, theta2, R1, R2 : float;

BEGIN
  theta1 := Im(Mean1);
  R1 := Re(Mean1);
  theta2 := Im(Mean2);
  R2 := Re(Mean2);
  Result := sqr(1-R1);
  IF Result = 0
    THEN
      BEGIN
        CircleError := true;
        WriteLn('Error: Circular kurtosis when resultant length is 1');
        EXIT;
      END;
  Result := (R2 * cos(theta2 - 2*theta1) - pot(R1, 4)) / Result;
END;


FUNCTION ConfidenceInterval (R, delta : float; n : WORD) : float;

VAR Rn, chi : float;
    c       : CHAR;

BEGIN
  Rn := R * n;
  chi := SignificanceLimit_Chi2(delta, 1);
  IF R > 0.9
    THEN
      Result := arccos(Sqrt(n - n*n - Rn*Rn * Exp(chi/n))/Rn)
    ELSE
      IF R > Sqrt(chi/(2*n))
        THEN
          Result := arccos(Sqrt((2 * n * (2*Rn*Rn - n*chi)) / (4*n - chi))/Rn)
        ELSE
          BEGIN
            c := WriteErrorMessage('Confidence interval of mean angle: vector length to small');
            CircleError := TRUE;
            EXIT;
          END;
END;

{ ************************************************************************** }

FUNCTION RandomVonMieses(mu, kappa : float) : float;

VAR s, u1, u2, theta : float;

BEGIN
  IF kappa > 1.3
    THEN s := 1/Sqrt(kappa)
    ELSE s := Const_pi * Exp(-kappa);
  REPEAT
    u1 := Random;
    u2 := Random;
    theta := s*(2*u2 - 1) / u1;    // generate prospective value
    IF Abs(theta) < Const_pi       // rejection pre-TEST
      THEN
        IF ((kappa*theta*theta) < (4 - 4*u1))  // acceptance pre-TEST
          THEN
            BEGIN
              Result := ModuloTwoPi(theta + mu);
              EXIT;
            END
          ELSE
            IF ((kappa*Cos(theta)) > (2*Ln(u1) + kappa))
              THEN
                BEGIN
                  Result := ModuloTwoPi(theta + mu);
                  EXIT;
                END
             ELSE  // IF < THEN TRY again
      ELSE // IF < THEN TRY again
  UNTIL FALSE;
END;

FUNCTION RandomUniformCircular : float;

BEGIN
  Result := Random * Const_2pi;
END;

{ ************************************************************************** }

FUNCTION Rayleigh (MeanVectorLength : float; n : WORD) : float;

VAR c : CHAR;
    x : float;

BEGIN
  IF (n < 10)
    THEN
      BEGIN
        c := WriteErrorMessage('Rayleigh test with less than 10 data points');
        CircleError := TRUE;
        EXIT;
      END;
  x := MeanVectorLength * n; // convert TO unscaled
  Result := Exp(Sqrt(1 + 4*n + 4*(n*n - x*x)) - (1 + 2*n));
END;


FUNCTION HodgesAjne (TransformedData : VectorTyp; MeanAngle : float) : float;

VAR StartAngle, StopAngle, x  : float;
    i, m1, m2, n              : WORD;

BEGIN
  StartAngle := ModuloTwoPi(MeanAngle - Nintey);
  StopAngle := ModuloTwoPi(MeanAngle + Nintey);
  IF StartAngle > StopAngle
    THEN
      BEGIN // exchange start- and stop angle
        x := StartAngle;
        StartAngle := StopAngle;
        StopAngle := x;
      END;
  n := VectorLength(TransformedData);
  m1 := 0;
  m2 := 0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(TransformedData, i);
      IF (x >= StartAngle) AND (x <= StopAngle)
        THEN INC(m1)  // count # OF data points on either side of the diameter
        ELSE INC(m2);
    END;
  IF (m1 > m2) THEN m1 := m2; // peak on the other side of diameter
  Result := 1/pot(2, Pred(n)) * (n - 2*m1) * BinomialCoef(n, m1);
END;

PROCEDURE ChiSqrTest (Data : VectorTyp; Direction : float; Sectors : WORD;
                      VAR ChiSqr : float; VAR dgf : WORD);

VAR i, j, Counter, n, max : WORD;
    Expected, Angle,
    StartAngle,
    Start, Finish, x      : float;
    ErrorStr, s           : STRING;

BEGIN
  n := VectorLength(Data);
  IF n < 8
    THEN
      BEGIN
        CircleError := TRUE;
        ch := WriteErrorMessage('Not enough data for ChiSqr-test (> 8 required)');
        EXIT;
      END;
  max := n DIV 4;
  IF max > 250 THEN max := 250;
  Expected := n / Sectors;
  IF Expected < 4
    THEN
      BEGIN
        CircleError := TRUE;
        Str(n:4, ErrorStr);
        Str(max:4, s);
        ErrorStr := ErrorStr + 'Data points allow at most ' + s + ' sectors';
        ch := WriteErrorMessage(ErrorStr);
        EXIT;
      END;
  ChiSqr := 0.0;
  Angle := Const_2pi / Sectors;
  StartAngle := Direction - 0.5 * Angle;
  dgf := Pred(Sectors);
  FOR i := 1 TO Sectors DO
    BEGIN
      Start := StartAngle + Pred(i) * Angle;
      Finish := Start + Angle;
      Counter := 0;
      FOR j := 1 TO n DO
        BEGIN
          x := GetVectorElement(Data, j);
          IF (Start <= Const_2pi) AND (Finish <= Const_2pi)
            THEN
              IF (x > Start) AND (x < Finish)
                THEN INC(Counter)
                ELSE // outside SECTOR
            ELSE
              IF (Start < Const_2pi) AND (Finish > Const_2pi)
                THEN
                  IF ((x > Start) AND (x < Const_2pi)) OR (x < ModuloTwoPi(Finish))
                    THEN INC(Counter)
                    ELSE
                ELSE
                  IF (x > ModuloTwoPi(Start)) AND (x < ModuloTwoPi(Finish))
                    THEN INC(Counter);
          END; // FOR j
       ChiSqr := ChiSqr + Sqr(Counter - Expected) / Expected;
    END;  // FOR i
END;


FUNCTION HomewardComponent (n : WORD; Mean : ComplexTyp; Expected : float) : float;

VAR Length, Phi, v : float;

BEGIN
  Polar(Mean, Length, Phi);
  Phi := ModuloTwoPi(Phi);
  v := Length * Cos(Phi - Expected);
  Result := Sqrt(2*n) * v;
END;


FUNCTION Rao (TransformedData : VectorTyp) : float;

VAR T, SumT, Expected : float;
    i, n              : WORD;

BEGIN
  n := VectorLength(TransformedData);
  Expected := Const_2pi / n;
  SumT := 0;
  FOR i := 2 TO n DO
    BEGIN
      T := GetVectorElement(TransformedData, i) - GetVectorElement(TransformedData, Pred(i));
      SumT := SumT + Abs(T - Expected);
    END;
  T := Const_2pi + GetVectorElement(TransformedData, 1) - GetVectorElement(TransformedData, n);
  SumT := SumT + Abs(T - Expected);
  Rao := 0.5 * SumT;
END;


PROCEDURE Kuipers (Data : MatrixTyp; VAR K, U : float);

VAR DPlus, DMinus, c, vSqr, cvn, vSum,
    Diff, Diff1, Diff2, Vn, vMean, x    : float;
    n, i                                : WORD;
    xVector, yVector,
    bVector                             : VectorTyp;

BEGIN
  n := MatrixRows(Data);
  GetColumn(Data, 1, xVector);
  GetColumn(Data, 2, yVector);
  GetColumn(Data, 3, bVector);
  vSum := NeumaierSum(bVector);
  vMean := vSum / n;
  DPlus := 0.0;
  DMinus := 0.0;
  vSqr := 0.0;
  cvn := 0.0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(bVector, i);
      Diff  := GetVectorElement(yVector, i) - x;
      Diff1 := GetVectorElement(yVector, Pred(i)) - x;
      Diff2 := GetVectorElement(yVector, Succ(i)) - x;
      IF Diff  > DPlus THEN DPlus := Diff;
      IF Diff1 > DPlus THEN DPlus := Diff1;
      IF Diff  < DMinus THEN DMinus := Diff;
      IF Diff1 < DMinus THEN DMinus := Diff1;
      c := Pred(2 * i);
      cvn := cvn + c * GetVectorElement(bVector, i) / n;
      vSqr := vSqr + Sqr(GetVectorElement(bVector, i));
    END;  // FOR
  Vn := (DPlus + Abs(DMinus)) * Sqrt(n);
  K := Vn;
  U := vSqr - cvn + n * (1/3 - Sqr(vMean - 0.5));
  DestroyVector(xVector);
  DestroyVector(yVector);
  DestroyVector(bVector);
END;


FUNCTION OneSample (MeanAngle, R, delta, TestAngle : float; n : WORD) : BOOLEAN;

VAR x : float;

BEGIN
  x := ConfidenceInterval(R, delta, n);
  Result := (TestAngle > (MeanAngle - delta)) AND (TestAngle < (MeanAngle + delta));
END;


PROCEDURE CalculateCumulativeFrequencies (Data : VectorTyp; FullCircle : float;
                                      VAR Result : MatrixTyp);

VAR i, n : WORD;
    Wert : float;
    D    : VectorTyp;

BEGIN
  n := VectorLength(Data);
  CopyVector(Data, D);
  ShellSort(D);
  FOR i := 1 TO n DO
   BEGIN
     Wert := GetVectorElement(D, i);
     SetMatrixElement(Result, i, 1, Wert);
     SetMatrixElement(Result, i, 2, i/n);
     SetMatrixElement(Result, i, 3, Wert/FullCircle);
   END;
  DestroyVector(D);
END;


FUNCTION MeanVectorGrouped (Transformed : MatrixTyp; Difference : float) : ComplexTyp;

VAR Total, n, i, j : WORD;
    x, y, z,
    Length, Angle  : float;
    Coord          : ComplexTyp;

BEGIN
    n := MatrixRows(Transformed);
    Total := 0;
    x := 0.0;
    y := 0.0;
    FOR i := 1 TO n DO
      BEGIN
        j := Round(GetMatrixElement(Transformed, i, 2));
        z := GetMatrixElement(Transformed, i, 1);
        Total := Total + j;
        x := x + j * Cos(z);
        y := y + j * Sin(z);
      END;
    x := x / Total;
    y := y / Total;
    Coord := ComplexInit(x, y);      // Re, Im
    Polar(Coord, Length, Angle);     // polar coordinates
    z := Difference / 2;
    z := z / Sin(z);
    Length := Length * z;            // correct quantisation error
    Result := Rect(Length, Angle);   // back TO cartesian coordinates }
END;


PROCEDURE WatsonWilliams (TransformedData1, TransformedData2 : VectorTyp);

VAR Transformed3                                   : VectorTyp;
    Mean1, Mean2, Mean3                            : ComplexTyp;
    Length1, Length2, Length3, Phi1, Phi2, Phi3,
    r, g, F, P0                                    : float;
    i, n1, n2, n3, nMean                           : WORD;

 BEGIN
  n1 := VectorLength(TransformedData1);
  n2 := VectorLength(TransformedData2);
  n3 := n1 + n2;
  CreateVector(Transformed3, n3, 0.0);
  FOR i := 1 TO n1 DO
    SetVectorElement(Transformed3, i, GetVectorElement(TransformedData1, i));
  FOR i := 1 TO n2 DO
    SetVectorElement(Transformed3, n1+i, GetVectorElement(TransformedData2, i));
  Mean1 := MeanVector(TransformedData1, 1);
  Polar(Mean1, Length1, Phi1);
  Phi1 := ModuloTwoPi(Phi1);
  Mean2 := MeanVector(TransformedData2, 1);
  Polar(Mean2, Length2, Phi2);
  Phi2 := ModuloTwoPi(Phi2);
  Mean3 := MeanVector(Transformed3, 1);
  Polar(Mean3, Length3, Phi3);
  Phi3 := ModuloTwoPi(Phi3);
  r := (Length1 + Length2) / 2;
  nMean := (n1 + n2) DIV 2;
  g := 1 + 3 / (8 * Kappa(r, nMean));
  F := g * (Length1 + Length2 - Length3) / (n3 - (Length1 + Length2));
  P0 := Integral_F(F, 1, n3-2);
  Writeln('F = ', FloatStr(F, ValidFigures), ', f1 = 1, f2 = ', n3-2:4,
          ' => P0 = ' + FloatStr(P0, ValidFigures));
  DestroyVector(Transformed3);
END;


FUNCTION FTest (Data1, Data2 : VectorTyp) : float;

VAR Mean1, Mean2                               : ComplexTyp;
    Angle1, Angle2, Length1, Length2, rMean    : float;
    n1, n2                                     : WORD;

BEGIN
  n1 := VectorLength(Data1);
  Mean1 := MeanVector(Data1, 1);
  Polar(Mean1, Length1, Angle1);
  n2 := VectorLength(Data2);
  Mean2 := MeanVector(Data2, 2);
  Polar(Mean2, Length2, Angle2);
  Length1 := Length1 * n1;
  Length2 := Length2 * n2;
  rMean := (Length1 + Length2) / (n1 + n2);
  IF rMean < 0.7
    THEN
     BEGIN
       CH := WriteErrorMessage('Mean vectors to short, F-Test not applicable');
       CircleError := TRUE;
       EXIT;
     END;
  Result := Pred(n2) * (n1 - Length1) / (Pred(n1) * (n2 - Length2));
END;


FUNCTION Wilcoxon (Data1, Data2 : VectorTyp; Direction1, Direction2 : float) : float;

VAR U1, U2, x                       : float;
    n1, n2, Rank, i, j, Sum1, Sum2  : WORD;
    Diff1, Diff2                    : VectorTyp;

BEGIN
  n1 := VectorLength(Data1);
  CreateVector(Diff1, n1, 0.0);
  FOR i := 1 TO n1 DO
   BEGIN
     x := Abs(GetVectorElement(Data1, i) - Direction1);
     IF x > Pi
       THEN SetVectorElement(Diff1, i, Const_2pi-x)
       ELSE SetVectorElement(Diff1, i, x);
   END;
  ShellSort(Diff1);
  n2 := VectorLength(Data2);
  CreateVector(Diff2, n2, 0.0);
  FOR i := 1 TO n2 DO
   BEGIN
     x := Abs(GetVectorElement(Data2, i) - Direction2);
     IF x > Pi
       THEN SetVectorElement(Diff2, i, Const_2pi-x)
       ELSE SetVectorElement(Diff2, i, x);
   END;
  ShellSort(Diff2);
  Rank := 1;
  j := 1;
  Sum1 := 0;
  Sum2 := 0;
  IF n1 < n2
    THEN
     BEGIN
       FOR i := 1 TO n1 DO
         BEGIN
           WHILE (GetVectorElement(Diff2, j) < GetVectorElement(Diff1, i)) DO
             BEGIN
               INC(Sum2, Rank);
               INC(Rank);
               INC(j);
             END;
           INC(Sum1, Rank);
           INC(Rank);
         END;
       FOR i := j TO n2 DO
        BEGIN
          INC(Sum2, Rank);
          INC(Rank);
        END;
     END
    ELSE
     BEGIN
       FOR i := 1 TO n2 DO
         BEGIN
           WHILE (GetVectorElement(Diff1, j) < GetVectorElement(Diff2, i)) DO
             BEGIN
               INC(Sum1, Rank);
               INC(Rank);
               INC(j);
             END;
           INC(Sum2, Rank);
           INC(Rank);
         END;
       FOR i := j TO n1 DO
        BEGIN
          INC(Sum1, Rank);
          INC(Rank);
        END;
     END;
  DestroyVector(Diff1);
  DestroyVector(Diff2);
  U1 := Sum1 - n1 * Succ(n1) / 2;
  U2 := Sum2 - n2 * Succ(n2) / 2;
  IF U1 < U2
    THEN Result := U1
    ELSE Result := U2;
END;


FUNCTION KruskalWallis (TransformedData1, TransformedData2 : VectorTyp) : float;

VAR Data3                            : VectorTyp;
    Median3, x, xmin, xmax, P           : float;
    i, n1, n2, n3, m1, m2, m                          : WORD;

BEGIN
  n1 := VectorLength(TransformedData1);
  n2 := VectorLength(TransformedData2);
  n3 := n1 + n2;
  IF (n1 <= 10) OR (n2 <= 10)
    THEN
      BEGIN
        ch := WriteErrorMessage('Circular Kruskal-Wallis test: n > 10 required for all groups');
        CircleError := TRUE;
        EXIT;
      END;
  CreateVector(Data3, n3, 0.0);
  FOR i := 1 TO n1 DO
    SetVectorElement(Data3, i, GetVectorElement(TransformedData1, i));
  FOR i := 1 TO n2 DO
    SetVectorElement(Data3, n1+i, GetVectorElement(TransformedData2, i));
  ShellSort(Data3);
  Median3 := MedianDirection(Data3);
  xmin := ModuloTwoPi(Median3 + Nintey);
  xmax := ModuloTwoPi(Median3 - Nintey);
  m1 := 0;
  FOR i := 1 TO n1 DO
    BEGIN
      x := GetVectorElement(TransformedData1, i);
      IF (x > xmin) AND (x < xmax) THEN INC(m1);
    END;
  m2 := 0;
  FOR i := 1 TO n2 DO
    BEGIN
      x := GetVectorElement(TransformedData2, i);
      IF (x > xmin) AND (x < xmax) THEN INC(m2);
    END;
  m := m1 + m2;
  P := n3*n3 / (m * (n3 - m)) * (m1*m1/n1 + m2*m2/n2) - n3 * m / (n3 - m);
  Result := IntegralChi(P, 1);
  Writeln('Chi^2 = ', FloatStr(P, ValidFigures), ' with 1 degrees of freedom, P0 = ',
          FloatStr(x, ValidFigures));
  DestroyVector(Data3);
END;


PROCEDURE Difference (Data1, Data2 : VectorTyp; VAR KuipersV, WatsonsUSqr : float);

TYPE ListenRecTyp = RECORD
                      Angle, k1, k2 : float;
                    END;

VAR n1, n2, Cumulative1, Cumulative2                  : WORD;
    Diff, Dif, DifSqr, DifPos, DifNeg, Angle1, Angle2 : float;
    ListenRec                                         : ListenRecTyp;
    Liste                                             : FiFo;

BEGIN
  n1 := VectorLength(Data1);
  n2 := VectorLength(Data2);
  ShellSort(Data1);
  ShellSort(Data2);
  Cumulative1 := 0;
  Cumulative2 := 0;
  InitFiFo(Liste);
  WHILE (Cumulative1 < n1) AND (Cumulative2 < n2) DO          { Data durchlaufen, }
    BEGIN                                                     { cumulative Frequenzen }
      Angle1 := GetVectorElement(Data1, Succ(Cumulative1));   { erzeugen und in FIFO  }
      Angle2 := GetVectorElement(Data2, Succ(Cumulative2));   { Liste speichern }
      IF Angle1 < Angle2
        THEN
         BEGIN
             INC(Cumulative1);
             ListenRec.Angle := Angle1;
         END
        ELSE
         BEGIN
             INC(Cumulative2);
             ListenRec.Angle := Angle2;
         END;
      ListenRec.k1 := 1.0 * Cumulative1 / n1;
      ListenRec.k2 := 1.0 * Cumulative2 / n2;
      Put(Liste, ListenRec, SizeOf(ListenRec));
    END; { while }
  WHILE (Cumulative1 < n1) DO                               { eventuell übrig gebliebene }
   BEGIN                                                    { Elemente aus 1. Datasatz  }
     INC(Cumulative1);                                      { bearbeiten }
     ListenRec.k1 := 1.0 * Cumulative1 / n1;
     ListenRec.k2 := 1.0 * Cumulative2 / n2;
     ListenRec.Angle := GetVectorElement(Data1, Cumulative1);
     Put(Liste, ListenRec, SizeOf(ListenRec));
   END;
  WHILE (Cumulative2 < n2) DO                               { dito für 2. Datasatz }
   BEGIN
     INC(Cumulative2);
     ListenRec.k1 := 1.0 * Cumulative1 / n1;
     ListenRec.k2 := 1.0 * Cumulative2 / n2;
     ListenRec.Angle := GetVectorElement(Data2, Cumulative2);
     Put(Liste, ListenRec, SizeOf(ListenRec));
   END;
  DifPos := 0.0;
  DifNeg := 0.0;
  Dif    := 0.0;
  DifSqr := 0.0;
  WHILE NOT EmptyFiFo(Liste) DO                             { Liste durchlaufen und }
    BEGIN                                                   { größte positive und   }
      Get(Liste, ListenRec, SizeOf(ListenRec));             { negative Difference suchen }
      Diff := ListenRec.k1 - ListenRec.K2;
      IF Diff > DifPos THEN DifPos := Diff;
      IF Diff < DifNeg THEN DifNeg := Diff;
      Dif := Dif + Diff;                                    { aufsummieren für Watson's Test }
      DifSqr := DifSqr + Sqr(Diff);
    END;
  KuipersV := n1 * n2 * (DifPos - DifNeg);
  WatsonsUSqr := n1 * n2 / Sqr(n1+n2) * (DifSqr - Sqr(Dif) / (n1+n2));
END;


PROCEDURE DescribeBivariat (Data : MatrixTyp);

VAR i, n1                                                     : WORD;
    SumX, SumY, a, b, c, d, r, SumXSqr, SumYSqr, SumXDevYDev,
    xSta, ySta, xMean, yMean, CoVariance, Correlation, Phi,
    Major, Minor, xMin, xMax, yMin, yMax                      : float;

BEGIN
  SumX := 0.0;
  SumY := 0.0;
  SumXSqr := 0.0;
  SumYSqr := 0.0;
  SumXDevYDev := 0.0;
  n1 := MatrixRows(Data);
  xMin := MaxRealNumber;
  xMax := MinRealNumber;
  yMin := MaxRealNumber;
  yMax := MinRealNumber;
  FOR i := 1 TO n1 DO
    BEGIN
      a := GetMatrixElement(Data, i, 1);
      IF a > xMax THEN xMax := a;
      IF a < xMin THEN xMin := a;
      SumX := SumX + a;
      SumXSqr := SumXSqr + Sqr(a);
      a := GetMatrixElement(Data, i, 2);
      IF a > yMax THEN yMax := a;
      IF a < yMin THEN yMin := a;
      SumY := SumY + a;
      SumYSqr := SumYSqr + Sqr(a);
    END;
  xMean := SumX / n1;
  yMean := SumY / n1;
  xSta := Sqrt((SumXSqr - Sqr(SumX) / n1) / Pred(n1));
  ySta := Sqrt((SumYSqr - Sqr(SumY) / n1) / Pred(n1));
  Write('Mean x = ', FloatStr(xMean, ValidFigures), ' ± ', FloatStr(xSta, ValidFigures));
  Writeln('Mean y = ', FloatStr(yMean, ValidFigures), ' ± ', FloatStr(ySta, ValidFigures));
  FOR i := 1 TO n1 DO
    BEGIN
      a := GetMatrixElement(Data, i, 1) - xMean;
      b := GetMatrixElement(Data, i, 2) - yMean;
      SumXDevYDev := SumXDevYDev + a * b;
    END;
  CoVariance := SumXDevYDev / Pred(n1);
  Correlation := CoVariance / (xSta * ySta);
  A := Sqr(ySta);
  B := - CoVariance;
  C := Sqr(xSta);
  D := (1 - Sqr(Correlation)) * Sqr(xSta) * Sqr(ySta);
  R := Sqrt(Sqr(A - C) + 4 * Sqr(B));
  Phi := ArcTan(2 * B / (A - C - R));
  Major := Sqrt(2 * D / (A + C - R));
  Minor := Sqrt(2 * D / (A + C + R));
  Writeln('r = ', FloatStr(Correlation, ValidFigures), ' Φ = ', FloatStr(Phi, ValidFigures));
  Writeln('Major = ', FloatStr(Major, ValidFigures), ' Minor = ', FloatStr(Minor, ValidFigures));
  Write('Press any key: ');
END;

// regression and correlation betweeen circular independent and linear dependent variable

FUNCTION CAssociation (CONST theta, Y : VectorTyp; VAR Sig : SignificanceType) : float;

VAR Tc, Ts, an, Un, P, xi, ti : float;
    i, n                      : WORD;

  FUNCTION InterpolateCritical (Un : float; n : WORD) : float;
  // interpolates critical values in table A10 from Fisher (1993)

  const a10 = 4.608; a05 = 5.909; a01 = 8.912;
        b10 = 0.668; b05 = 0.773; b01 = 0.865;

  BEGIN
    if Un < (a10 * (1 - pot(b10, n)))                // critical value for 10%
      THEN Result := 1.0                             // to give it some value
      ELSE if Un < (a05 * (1 - pot(b05, n)))         // critical value for 5%
             THEN Result := 0.1
             ELSE if Un < (a01 * (1 - pot(b01, n)))  // critical value for 1%
                    THEN Result := 0.05
                    ELSE Result := 0.01;
  END;

BEGIN
  n := VectorLength(theta);
  IF VectorLength(Y) <> n
    THEN
     BEGIN
       ch := WriteErrorMessage('Linear-circular association: length of theta and x not equal');
       CircleError := true;
       EXIT;
     END;
  Tc := 0;
  Ts := 0;
  FOR i := 1 TO n DO
    BEGIN
      xi := GetVectorElement(Y, i);
      ti := GetVectorElement(theta, i);
      IF IsNaN(xi) OR IsNaN(ti)
        THEN // ignore data pair if either is NaN
        ELSE
          BEGIN
            Tc := Tc + xi*cos(ti);
            Ts := Ts + xi*sin(ti);
          END;
    END;
  IF ODD(n)
    THEN an := 2*pot(sin(Const_pi/n), 4) / pot(1+cos(Const_pi/n), 3)
    ELSE an := 1 / (1 + 5*pot(cot(Const_pi/n), 2) + 4*pot(cot(Const_pi/n), 4));
  Result := an * (Tc*Tc * Ts*Ts); // correlation coefficient
  WITH Sig DO
    BEGIN
      TestValue := 24 * (Tc*Tc * Ts*Ts) / (pot(n, 3) + pot(n, 2)); // test statistics Un
      Freedom := n;
      IF n > 100
        THEN P0 := exp(-pot(Un, 2)/2)
        ELSE P0 := InterpolateCritical(Un, n);
    END;
END;

PROCEDURE RankXY (const TransformedData, yData : VectorTyp; var Ranks : VectorTyp);

VAR n, i, Rank, MaxPos : WORD;
    Kopie              : MatrixTyp;
    Maximum, x         : float;

BEGIN
  n := VectorLength(TransformedData);
  CreateMatrix(Kopie, n, 2, 0.0);
  SetColumn(Kopie, TransformedData, 1);
  SetColumn(Kopie, yData, 2);
  Rank := n;                       { Feld nach aufsteigenden x sortieren }
  REPEAT
    Maximum := MinRealNumber;
    FOR i := 1 TO Rank DO
      BEGIN
        x := GetMatrixElement(Kopie, i, 1);
        IF x > Maximum
          THEN
            BEGIN
              Maximum := x;
              MaxPos := i;
            END;
      END;
    x := GetMatrixElement(Kopie, MaxPos, 2);
    SetMatrixElement(Kopie, MaxPos, 1, GetMatrixElement(Kopie, Rank, 1));
    SetMatrixElement(Kopie, MaxPos, 2, GetMatrixElement(Kopie, Rank, 2));
    SetMatrixElement(Kopie, Rank, 1, Rank);
    SetMatrixElement(Kopie, Rank, 2, x);
    DEC(Rank);
  UNTIL (Rank = 0);
  Rank := n;                       { jetzt die Ränge der y-Elemente bestimmen }
  REPEAT
    Maximum := MinRealNumber;
    FOR i := 1 TO Rank DO
      BEGIN
        x := GetMatrixElement(Kopie, i, 2);
        IF x > Maximum
          THEN
            BEGIN
              Maximum := x;
              MaxPos := i;
            END;
      END;
    x := GetMatrixElement(Kopie, MaxPos, 1);
    SetMatrixElement(Kopie, MaxPos, 1, GetMatrixElement(Kopie, Rank, 1));
    SetMatrixElement(Kopie, MaxPos, 2, GetMatrixElement(Kopie, Rank, 2));
    SetMatrixElement(Kopie, Rank, 2, Maximum);
    SetMatrixElement(Kopie, Rank, 1, x);
    SetVectorElement(Ranks, Round(x), Rank);
    DEC(Rank);
  UNTIL (Rank = 0);
  DestroyMatrix(Kopie);
END;


PROCEDURE LinearPeriodicRankCorrelation (const TransformedData, yData: VectorTyp;
                                         VAR Correlation, U : float);

VAR Ranks               : VectorTyp;
    n, i                : WORD;
    epsilon, x, y, c, s : float;

BEGIN
  n := VectorLength(TransformedData);
  IF n <> VectorLength(yData)
    THEN
      BEGIN
        ch := WriteErrorMessage('Linear-periodic correlation: Vectors of different lengths');
        CircleError := TRUE;
        EXIT;
      END;
  epsilon := Const_2pi / n;
  CreateVector(Ranks, n, 0.0);
  RankXY(TransformedData, yData, Ranks);
  c := 0.0;
  s := 0.0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(Ranks, i);
      y := epsilon * i;
      c := c + x * Cos(y);
      s := s + x * Sin(y);
    END;
  IF Odd(n)
    THEN x := 2 * pot(Sin(Pi/n),4) / pot(1 + Cos(Pi/n), 3)
    ELSE x := 1 / (1 + 5 * pot(cot(Pi/n),2) + 4 * pot(cot(Pi/n),4));
  y := 24 / (x * Sqr(n) * (n+1));
  Correlation := x * (Sqr(c) + Sqr(s));
  U := y * Correlation;
  DestroyVector(Ranks);
END;


PROCEDURE CircularLinearCorrelation (const TransformedData, yData: VectorTyp;
                                     VAR Correlation : float; VAR Sig : SignificanceType);

VAR SinPhi, CosPhi, y : VectorTyp;
    n, i              : WORD;
    x, ryc, rys, rcs  : float;

BEGIN
  n := VectorLength(TransformedData);
  IF n <> VectorLength(yData)
    THEN
      BEGIN
        ch := WriteErrorMessage('Linear-periodic correlation: Vectors of different lengths');
        CircleError := TRUE;
        EXIT;
      END;
  CreateVector(SinPhi, n, 0.0);
  CreateVector(CosPhi, n, 0.0);
  CreateVector(y, n, 0.0);
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(TransformedData, i);
      SetVectorElement(SinPhi, i, Sin(x));
      SetVectorElement(CosPhi, i, Cos(x));
    END;
  ryc := PearsonProductMomentCorrelation(yData, CosPhi, sig);
  rys := PearsonProductMomentCorrelation(yData, SinPhi, sig);
  rcs := PearsonProductMomentCorrelation(CosPhi, SinPhi, sig);
  DestroyVector(SinPhi);
  DestroyVector(CosPhi);
  Correlation := (Sqr(ryc) + Sqr(rys) - 2 * ryc * rys * rcs) / (1 - Sqr(rcs)); // r^2
  WITH Sig DO
    BEGIN
      TestValue := n * Correlation;
      Freedom := 2;
      P0 := IntegralChi(TestValue, 2);
    END;
  Correlation := Sqrt(Correlation);
END;


PROCEDURE TrigonometricPolynomial (const TransformedData, yData: VectorTyp;
                                   Periode : float);

VAR i, n                                                     : WORD;
    Omega, t, s, c, y, tMax, tMin, yMax, yMin,
    SumY, SumC, SumS, SumCSqr, SumSSqr, SumYC,
    SumYS, SumCS, M, X, Z, Hilfs1, Hilfs2, Hilfs3, Hilfs4,
    Correlation, U, A, Phi, r, TEST, P0                       : float;
    Sig                                                      : SignificanceType;
    hstr                                                     : STRING;

BEGIN
  n := VectorLength(TransformedData);
  IF n <> VectorLength(yData)
    THEN
      BEGIN
        ch := WriteErrorMessage('Linear-periodic correlation: Vectors of different lengths');
        CircleError := TRUE;
        EXIT;
      END;
  SumY := 0.0;
  SumC := 0.0;
  SumS := 0.0;
  SumCSqr := 0.0;
  SumSSqr:= 0.0;
  SumYC := 0.0;
  SumYS := 0.0;
  SumCS := 0.0;
  tMax := MinRealNumber;
  tMin := MaxRealNumber;
  yMax := MinRealNumber;
  yMin := MaxRealNumber;
  FOR i := 1 TO n DO
    BEGIN
      t := GetVectorElement(TransformedData, i);  // circular datum
      IF t > tMax THEN tMax := t;
      IF t < tMin THEN tMin := t;
      c := Cos(Periode*t);
      s := Sin(Periode*t);
      SumC := SumC + c;
      SumS := SumS + s;
      y := GetVectorElement(yData, i);            // linear datum
      IF y > yMax THEN yMax := y;
      IF y < yMin THEN yMin := y;
      SumY := SumY + y;
      SumCSqr := SumCSqr + Sqr(c);
      SumSSqr:= SumSSqr + Sqr(s);
      SumYC := SumYC + y * c;
      SumYS := SumYS + y * s;
      SumCS := SumCS + c * s;
    END;
  Hilfs1 := Sqr(SumC) + n * SumSSqr;
  Hilfs2 := (SumS * SumC * SumCS - Sqr(SumS) - n * Sqr(SumCS))
            / Hilfs1 - Sqr(SumS) / n + SumSSqr;
  Hilfs3 := (Sqr(SumS) * SumYC - n * SumYC * SumCS +
            SumCS * SumC * SumY - SumCS * SumC * SumS) / Hilfs1;
  Hilfs4 := (Sqr(SumS) * Sqr(SumC) - SumS * Sqr(SumC) * SumY)
            / (n * Hilfs1);
  Z := (SumYS + Hilfs3 + Hilfs4 - SumS * SumY / n) / Hilfs2;
  X := (n * SumYC - SumC * SumY + SumC * SumS - n * SumCS * Z)
       / Hilfs1;
  M := (SumY - SumC * X - SumS * Z) / n;
  A := Sqrt(Sqr(X) + Sqr(Z));
  CASE signum(X) OF
   1 : Phi := ArcTan(Z/X);
  -1 : IF x < 0 THEN Phi := Const_pi + ArcTan(Z/X);
   0 : CASE signum(y) OF
        1 : Phi := Nintey;                       {  90° }
       -1 : IF Y < 0 THEN Phi := Const_pi*3/2;   { 270° }
        0 :  BEGIN
               CH := WriteErrorMessage('unknown akrophase angle ');
               CircleError := TRUE;
               EXIT;
             END;
       END;
  END;
  IF tMin > 0 THEN tmin := 0;
  IF yMin > 0 THEN yMin := 0;
  IF A < 0
    THEN hstr := '-'
    ELSE hstr := '+';
  Writeln('y = ', FloatStr(M, ValidFigures), ' ', HStr, ' ', FloatStr(A, ValidFigures),
          ' * cos(', FloatStr(Periode, ValidFigures), ' * t - ', FloatStr(Phi, ValidFigures), ')');
  CircularLinearCorrelation(TransformedData, yData, r, Sig);
  WITH Sig DO
    Writeln('Parametric: r = ', FloatStr(r, ValidFigures), ' chi2 = ',
     FloatStr(TestValue, ValidFigures), ' with 2 dgf, P(H0: r = 0) = ', FloatStr(P0, ValidFigures));
  LinearPeriodicRankCorrelation(TransformedData, yData, Correlation, u);
  Writeln('Non-Parametric r = ', FloatStr(Correlation, ValidFigures), ' U = ', FloatStr(U, ValidFigures));
END;

// ********** correlation between two circular variables **********************

PROCEDURE PeriodicRankCorrelation (CONST alpha, beta : VectorTyp;
                                     VAR rPlus, rMinus : float);

VAR Ranks                                             : VectorTyp;
    n, i                                              : WORD;
    Epsilon, DeltaPlus, DeltaMinus, x, y,
    SumSinPlus, SumSinMinus, SumCosPlus, SumCosMinus  : float;

BEGIN
  n := VectorLength(alpha);
  IF (n <> VectorLength(beta))
    THEN
      BEGIN
        ch := WriteErrorMessage('Circular-circular correlation: Vectors of unequal length');
        CircleError := TRUE;
        EXIT;
      END;
  CreateVector(Ranks, n, 0.0);
  RankXY(alpha, beta, Ranks);
  SumCosPlus := 0.0;
  SumCosMinus := 0.0;
  SumSinPlus := 0.0;
  SumSinMinus := 0.0;
  Epsilon := Const_2pi / n;
  FOR i := 1 TO n DO
    BEGIN
      DeltaPlus := epsilon * (i + GetVectorElement(Ranks, i));
      SumCosPlus := SumCosPlus + Cos(DeltaPlus);
      SumSinPlus := SumSinPlus + Sin(DeltaPlus);
      DeltaMinus := epsilon * (i - GetVectorElement(Ranks, i));
      SumCosMinus := SumCosMinus + Cos(DeltaMinus);
      SumSinMinus := SumSinMinus + Sin(DeltaMinus);
    END;
  DestroyVector(Ranks);
  rPlus := Sqrt(Sqr(SumCosMinus) + Sqr(SumSinMinus)) / n;
  rMinus := Sqrt(Sqr(SumCosPlus) + Sqr(SumSinPlus)) / n;
END;


FUNCTION JuppMardia (CONST alpha, beta : VectorTyp) : float;

VAR rCC, rCs, rSC, rSS, r1, r2, x, y : float;
    cosX, sinX, cosY, sinY           : VectorTyp;
    n, i                             : WORD;
    Sig                              : SignificanceType;

BEGIN
  n := VectorLength(alpha);
  IF (n <> VectorLength(beta))
    THEN
      BEGIN
        ch := WriteErrorMessage('Circular-circular correlation: Vectors of unequal length');
        CircleError := TRUE;
        EXIT;
      END;
  CreateVector(cosX, n, 0.0);
  CreateVector(sinX, n, 0.0);
  CreateVector(cosY, n, 0.0);
  CreateVector(sinY, n, 0.0);
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(alpha, i);
      y := GetVectorElement(beta, i);
      SetVectorElement(cosX, i, Cos(x));
      SetVectorElement(sinX, i, Sin(x));
      SetVectorElement(cosY, i, Cos(Y));
      SetVectorElement(sinY, i, Sin(y));
    END;
  rCC := PearsonProductMomentCorrelation(cosX, cosY, Sig);
  rCS := PearsonProductMomentCorrelation(cosX, sinY, Sig);
  rSC := PearsonProductMomentCorrelation(sinX, cosY, Sig);
  rSS := PearsonProductMomentCorrelation(sinX, sinY, Sig);
  r1  := PearsonProductMomentCorrelation(cosX, sinX, Sig);
  r2  := PearsonProductMomentCorrelation(cosY, sinY, Sig);
  DestroyVector(cosX);
  DestroyVector(sinX);
  DestroyVector(cosY);
  DestroyVector(sinY);
  Result := (Sqr(rCC) + Sqr(rCS) + Sqr(rSC) + Sqr(rSS) + 2 * (rCC * rSS + rCS * rSC) * r1 * r2
          - 2 * (rCC * rCS + rSC * rSS) * r2 - 2 * (rCC * rSC + rCS * rSS) * r1)
          / ((1 - Sqr(r1)) * (1 - Sqr(r2)));
END;


PROCEDURE CircularCircularCorrelation (Alpha, Beta : VectorTyp;
          VAR Correlation : float; var Sig : SignificanceType);

VAR SumSinSin, SumSinSinSqr, SumSinA, SumSinB,
    MeanA, MeanB, f, x, y                     : float;
    n, na, i                                  : WORD;

BEGIN
  n := VectorLength(alpha);
  IF (n <> VectorLength(beta))
    THEN
      BEGIN
        ch := WriteErrorMessage('Circular-circular correlation: Vectors of unequal length');
        CircleError := TRUE;
        EXIT;
      END;
  MeanA := MedianDirection(Alpha);
  MeanB := MedianDirection(Beta);
  SumSinSin := 0;
  SumSinSinSqr := 0;
  SumSinA := 0;
  SumSinB := 0;
  na := 0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(Alpha, i);
      y := GetVectorElement(Beta, i);
      IF IsNaN(x) OR IsNaN(y)
        THEN // DO nothing, allows data containing NaNs
        ELSE
          BEGIN
            INC(na);  // actual n, excluding NaNs
            x := Sin(x - MeanA);
            y := Sin(y - MeanB);
            SumSinSin := SumSinSin + x*y;
            SumSinSinSqr := SumSinSinSqr + x*x * y*y;
            SumSinA := SumSinA + x*x;
            SumSinB := SumSinB + y*y;
          END;
    END;
  Correlation :=  SumSinSin / Sqrt(SumSinSinSqr);
  f := na * SumSinA * SumSinB / SumSinSinSqr;
  WITH Sig DO
    BEGIN
      TestValue := Sqrt(f) * Correlation;
      P0 := IntegralGauss(TestValue);
    END;
END;


END.

