UNIT LinearRegression;

INTERFACE

USES math, MathFunc, Vector, Stat, Correlations;

CONST RegressionError : BOOLEAN = FALSE;

TYPE CurveTyp  = (Origin, Linear, Exponential, Power, Hyperbola, Inverse, Maximum,
                  Sigmoidal, ExpSigmoidal, ModPower, Hill);
     ResultTyp = RECORD
                   a, sa,                            // intercept
                   b, sb,                            // slope
                   m, sm,                            // exponent
                   r, t, P0 : double;                // TEST statistics
                 END;

PROCEDURE Approximation(CONST x, y, weight : VectorTyp; ct : CurveTyp;
          VAR yCalc : VectorTyp; VAR Res : ResultTyp);
// linear regression FOR least sum OF squares

PROCEDURE Transformation(CONST xOrig, yOrig : VectorTyp; VAR xTrans, yTrans : VectorTyp;
                         ct : CurveTyp; Schaetzwert : double);
// transformation OF x, y TO linearise

PROCEDURE Retransform (CONST x, yTrans : Vectortyp; VAR TRes, Res : ResultTyp;
                       ct : CurveTyp; Schaetzwert : double;
                       VAR yCalc : VectorTyp);
// re-transformation OF linear results TO curve

PROCEDURE LinFit (Data : MatrixTyp;  mwt: BOOLEAN; VAR a, b, siga, sigb, chi2, q: real);
{ Robuste Daten-Modellierung nach Press et al.: Numerical Recepies in Pascal,
  Cambridge 1989, pp 590-8
  Anpassung von Daten an eine Ausgleichsgerade. Dabei wird jedoch die uebliche
  Annahme fallengelassen, alle Datenpunkte haetten die gleiche Standardabweichung.
  Eingabe: Data enthaelt die Daten, dabei stehen in jeder Zeile x, y und die
           Standardabweichung fuer y (falls letztere nicht zur Verfuegung
           stehen, muá mwt = false gesetzt werden).
  Ausgabe: a und b sind die Parameter des Modells, siga und sigb die Fehler-
           grenzen. chi2 ist ein Mass fuer die Guete des Modells, q die War-
           scheinlichkeit gegen dieses Modell. }

PROCEDURE Deming (CONST x, y : VectorTyp; delta : double;
                  VAR xCalc, yCalc : VectorTyp; VAR Res : ResultTyp);
// Deming regression FOR data WITH errors IN x AND y

PROCEDURE TheilSenKendall (CONST x, y : VectorTyp; VAR yCalc : VectorTyp;
                           VAR Res : ResultTyp);
// non-parametric regression

PROCEDURE Eisenthal (CONST Substrate, Velocity : VectorTyp; VAR yCalc : VectorTyp;
                           VAR Res : ResultTyp);
// Eisenthal & Cornish-Bowden 1974



IMPLEMENTATION

PROCEDURE Approximation(CONST x, y, weight : VectorTyp; ct : CurveTyp;
          VAR yCalc : VectorTyp; VAR Res : ResultTyp);

VAR Sx, Sy, Sw, Sxx, Syy, Sxy, SDxDy, SDx2, SDy2,xMean, yMean : double;
    c                                                         : CHAR;
    n                                                         : WORD;

  PROCEDURE DataSums;

  VAR
    i: INTEGER;
    DX, dy, xi, yi, wi : double;

  BEGIN
    Sx := 0;
    Sy := 0;
    Sw := 0;
    Sxy := 0;
    Sxx := 0;
    Syy := 0;
    SDxDy := 0;
    SDx2 := 0;
    SDy2 := 0;
    FOR i := 1 TO VectorLength(x) DO
      IF IsNaN(GetVectorElement(x, i)) OR IsNaN(GetVectorElement(y, i))
        THEN
        ELSE
          BEGIN
            xi := GetVectorElement(x, i);
            yi := GetVectorElement(y, i);
            wi := GetVectorElement(weight, i);
            Sx := Sx + xi * wi;
            Sy := Sy + yi * wi;
            Sw := Sw + wi;
            Sxy := Sxy + xi * yi * wi;
            Sxx := Sxx + xi * xi * wi;
            Syy := Syy + yi * yi * wi;
          END;
    xMean := Sx / Sw;
    yMean := Sy / Sw;
    FOR i := 1 TO VectorLength(x) DO
      IF IsNaN(GetVectorElement(x, i)) OR IsNaN(GetVectorElement(y, i))
        THEN
        ELSE
          BEGIN
            DX := GetVectorElement(x, i) - xMean;
            dy := GetVectorElement(y, i) - yMean;
            wi := GetVectorElement(Weight, i);
            SDxDy  := SDxDy + DX * dy * wi;
            SDx2 := SDx2  + DX * DX * wi;
            SDy2 := SDy2  + dy * dy * wi;
          END;
    SDxDy := SDxDy / Sw;
    SDx2  := SDx2  / Sw;
    SDy2  := SDy2  / Sw;
  END;


  PROCEDURE Parameters;

  VAR i: INTEGER;
      a1, a2, a3, a4 : double;

  BEGIN
    IF (ct = Origin)
      THEN
        BEGIN   // Line through origin
          Res.b := Sxy / Sxx;
          Res.a := 0;
        END
      ELSE
        BEGIN   // all others
          a1 := Sw * Sxy;
          a2 := Sx * Sy;
          a3 := Sw * Sxx;
          a4 := Sx * Sx;
          a4 := (a1-a2) / (a3-a4);
          Res.b := (Sw * Sxy - Sx * Sy) / (Sw * Sxx - Sx * Sx);
          Res.a := yMean - Res.b * xMean;
        END;
    FOR i := 1 TO VectorLength(x) DO
      IF IsNaN(GetVectorElement(x, i))
        THEN SetVectorElement(yCalc, i, NaN)
        ELSE SetVectorElement(yCalc, i, Res.a + Res.b * GetVectorElement(x, i));
  END;


  PROCEDURE ErrorEstimates;

  VAR i, f : INTEGER;
      Se : double;

  BEGIN
    CASE ct OF
      Origin            : f := Round(n - 1);    {degrees of freedom = data points - parameters}
      Linear..Maximum   : f := Round(n - 2);
      Sigmoidal..Hill   : f := Round(n - 3);
    END;
    WITH Res DO
      BEGIN
        se := 1/(Sw * f) * (Sw*Syy - Sy*Sy - Res.b*Res.b*(Sw*Sxx - Sx*Sx));
        sb := Sw*se / (Sw*Sxx - Sx*Sx);
        IF (ct = Origin)
          THEN sa := 0
          ELSE sa := (1/Sw) * sb *Sxx;
        se := Sqrt(se);
        sa := Sqrt(Res.sa);
        sb := Sqrt(Res.sb);
        r := SDxDy / Sqrt(SDx2 * SDy2);
        IF (r < 1)  // prevent division by 0
          THEN
            BEGIN
              t := r * Sqrt(f/(1-r*r));
              P0 := Integral_t(t, f);
            END
          ELSE
            BEGIN
              t := NaN;
              P0 := 0;
            END;
      END;
  END;

BEGIN
  IF VectorLength(x) <> VectorLength(y)
    THEN
      BEGIN
        c := WriteErrorMessage('Linear regression: unequal length of dependent and independent data vector');
        RegressionError := TRUE;
        EXIT;
      END;
  n := VectorLength(x);
  CreateVector(yCalc, VectorLength(x), 0.0);
  DataSums;
  Parameters;
  ErrorEstimates;
END;


PROCEDURE Transformation(CONST xOrig, yOrig : VectorTyp; VAR xTrans, yTrans : VectorTyp;
                         ct : CurveTyp; Schaetzwert : double);
// transformation OF x, y TO linearise

VAR i, n : INTEGER;

BEGIN
    n := VectorLength(xOrig);
    CreateVector(xTrans, n, 0.0);
    CreateVector(yTrans, n, 0.0);
    FOR i := 1 TO n DO
      IF IsNaN(GetVectorelement(xOrig, i))
        THEN
          SetVectorElement(xTrans, i, NaN)
        ELSE
          CASE ct OF                        {Transformation of x-values}
            Origin, Linear, Exponential, Inverse, Maximum, ExpSigmoidal :
               SetVectorElement(xTrans, i, GetVectorElement(xOrig, i));
            Power, Sigmoidal, ModPower :
               SetVectorElement(xTrans, i, Ln(GetVectorElement(xOrig, i)));
            Hyperbola :
               SetVectorElement(xTrans, i, 1 / GetVectorElement(xOrig, i));
            Hill :
               SetVectorElement(xTrans, i, log(GetVectorElement(xOrig, i), 10));
          END;
    FOR i := 1 TO n DO
      IF IsNaN(GetVectorelement(yOrig, i))
        THEN
          SetVectorElement(yTrans, i, NaN)
        ELSE
          CASE ct OF                 {Transformation of y-values}
            Origin, Linear          :  SetVectorElement(yTrans, i, GetVectorElement(yOrig, i));
            Exponential, Power      :  SetVectorElement(yTrans, i, Ln(GetVectorElement(yOrig, i)));
            Hyperbola, Inverse      :  SetVectorElement(yTrans, i, 1 / GetVectorElement(yOrig, i));
            Maximum                 :  IF IsNaN(GetVectorelement(xOrig, i))
                                         THEN SetVectorElement(yTrans, i, NaN)
                                         ELSE SetVectorElement(yTrans, i, Ln(GetVectorElement(xOrig, i) / GetVectorElement(yOrig, i)));
            Sigmoidal, ExpSigmoidal :  SetVectorElement(yTrans, i, Ln(Schaetzwert / GetVectorElement(yOrig, i) - 1));
            ModPower                :  SetVectorElement(yTrans, i, Ln(GetVectorElement(yOrig, i) - Schaetzwert));
            Hill                    :  SetVectorElement(yTrans, i, log(GetVectorElement(yOrig, i) / (Schaetzwert - GetVectorElement(yOrig, i)), 10));
          END;
END;


PROCEDURE Retransform (CONST x, yTrans : Vectortyp; VAR TRes, Res : ResultTyp;
                       ct : CurveTyp; Schaetzwert : double;
                       VAR yCalc : VectorTyp);

VAR i, n : WORD;

BEGIN
  n := VectorLength(x);
  CreateVector(yCalc, n, 0.0);
  CASE ct OF
    Origin, Linear : BEGIN
                       CopyVector(yTrans, yCalc);
                       Res := TRes;
                       Res.m := 0;
                       Res.sm := 0;
                     END;
    Exponential    : BEGIN
                       Res.a  := Exp(TRes.a);
                       Res.sa := TRes.sa/Tres.a * Res.a;
                       Res.b  := 0;
                       Res.sb := 0;
                       Res.m  := TRes.b;
                       Res.sm := TRes.sb;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a *
                                    Exp(Res.m * GetVectorElement(x, i)));
                         END;
                     END;
    Power          : BEGIN
                       Res.a  := Exp(TRes.a);
                       Res.sa := TRes.sa/Tres.a * Res.a;
                       Res.b  := 0;
                       Res.sb := 0;
                       Res.m  := TRes.b;
                       Res.sm := TRes.sb;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a *
                                    Pot(GetVectorElement(x, i), Res.m)) ;
                         END;
                     END;
    Hyperbola      : BEGIN
                       Res.a  := 1/TRes.a;
                       Res.sa := TRes.sa/Tres.a * Res.a;
                       Res.b  := Tres.b * Res.a;
                       Res.sb := TRes.sb/Tres.b * Res.b;
                       Res.m  := 0;
                       Res.sm := 0;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a * GetVectorElement(x, i)
                                     / (Res.b + GetVectorElement(x, i)));
                         END;
                     END;
    Inverse        : BEGIN
                       Res.a  := 1/TRes.b;
                       Res.sa := TRes.sb/Tres.b * Res.a;
                       Res.b  := Tres.a * Res.a;
                       Res.sb := TRes.sa/Tres.a * Res.b;
                       Res.m  := 0;
                       Res.sm := 0;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a /
                                    (Res.b + GetVectorElement(x, i)));
                         END;
                     END;
    Maximum        : BEGIN
                       Res.a  := Exp(-TRes.a);
                       Res.sa := Abs(TRes.sa/Tres.a * Res.a);
                       Res.b  := 0;
                       Res.sb := 0;
                       Res.m  := -TRes.b;
                       Res.sm := Abs(TRes.sb/Tres.b * Res.m);
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a * GetVectorElement(x, i)
                                     * Exp(Res.m * GetVectorElement(x, i)));
                         END;
                     END;
    Sigmoidal      : BEGIN
                       Res.a  := Schaetzwert;
                       Res.sa := NaN;
                       Res.b  := Exp(TRes.a);
                       Res.sb := TRes.sa/Tres.a * Res.a;
                       Res.m  := TRes.b;
                       Res.sm := TRes.sb/Tres.b * Res.m;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a /
                                    (1 + Res.b * Pot(GetVectorElement(x, i), Res.m)));
                         END;
                     END;
    ExpSigmoidal   : BEGIN
                       Res.a  := Schaetzwert;
                       Res.sa := NaN;
                       Res.b  := Exp(TRes.a);
                       Res.sb := TRes.sa/Tres.a * Res.a;
                       Res.m  := TRes.b;
                       Res.sm := TRes.sb/Tres.b * Res.m;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a /
                                    (1 + Res.b * Exp(GetVectorElement(x, i) * Res.m)));
                         END;
                     END;
    ModPower       : BEGIN
                       Res.a  := Schaetzwert;
                       Res.sa := NaN;
                       Res.b  := Exp(TRes.a);
                       Res.sb := TRes.sa/Tres.a * Res.a;
                       Res.m  := TRes.b;
                       Res.sm := TRes.sb/Tres.b * Res.m;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a + Res.b *
                                    Pot(GetVectorElement(x, i), Res.m));
                         END;
                     END;
    Hill           : BEGIN
                       Res.a  := Schaetzwert;
                       Res.sa := NaN;
                       Res.b  := pot(10, -TRes.a);
                       Res.sb := Abs(TRes.sa/Tres.a * Res.a);
                       Res.m  := TRes.b;
                       Res.sm := TRes.sb/Tres.b * Res.m;
                       CreateVector(yCalc, VectorLength(yTrans), 0.0);
                       FOR i := 1 TO n DO
                         BEGIN
                           IF IsNaN(GetVectorElement(yTrans, i))
                             THEN SetVectorElement(yCalc, i, NaN)
                             ELSE SetVectorElement(yCalc, i, Res.a * pot(GetVectorElement(x, i), Res.m)
                                    / (Res.b + Pot(GetVectorElement(x, i), Res.m)));
                         END;
                     END;
  END; // CASE
  Res.r := TRes.r;
  Res.t := TRes.t;
  Res.P0 := TRes.P0;
END;


PROCEDURE LinFit (Data : MatrixTyp;  mwt: BOOLEAN; VAR a, b, siga, sigb, chi2, q: real);

VAR i, nData            : INTEGER;
    wt, t, sy, sxoss,
    sx, st2, ss, sigdat : float;

BEGIN
   sx := 0.0;
   sy := 0.0;
   st2 := 0.0;
   b := 0.0;
   nData := MatrixRows(Data);
   IF mwt
     THEN
       BEGIN
         ss := 0.0;
         FOR i := 1 TO ndata DO
           BEGIN
             wt := 1.0 / Sqr(GetMatrixElement(Data, i, 3));
             ss := ss + wt;
             sx := sx + GetMatrixElement(Data, i, 1) * wt;
             sy := sy + GetMatrixElement(Data, i, 2) * wt
           END
      END
     ELSE
       BEGIN
         FOR i := 1 TO ndata DO
            BEGIN
                sx := sx + GetMatrixElement(Data, i, 1);
                sy := sy + GetMatrixElement(Data, i, 2);
            END;
         ss := ndata
       END;
   sxoss := sx/ss;
   IF mwt
     THEN
       BEGIN
         FOR i := 1 TO ndata DO
           BEGIN
             t := (GetMatrixElement(Data, i, 1) -sxoss) / GetMatrixElement(Data, i, 3);
             st2 := st2+t*t;
             b := b + t * GetMatrixElement(Data, i, 2) / GetMatrixElement(Data, i, 3);
           END
       END
     ELSE
       BEGIN
         FOR i := 1 TO ndata DO
           BEGIN
             t := GetMatrixElement(Data, i, 1) - sxoss;
             st2 := st2 + t * t;
             b := b + t * GetMatrixElement(Data, i, 2);
           END
       END;
   b := b / st2;
   a := (sy - sx * b) / ss;
   siga := Sqrt((1.0 + sx * sx / (ss * st2)) / ss);
   sigb := Sqrt(1.0 / st2);
   chi2 := 0.0;
   IF NOT mwt
     THEN
       BEGIN
         FOR i := 1 TO ndata DO
           chi2 := chi2 + Sqr(GetMatrixElement(Data, i, 2) - a - b * GetMatrixElement(Data, i, 1));
         q := 1.0;
         sigdat := Sqrt(chi2 / (ndata - 2));
         siga := siga * sigdat;
         sigb := sigb * sigdat;
       END
     ELSE
       BEGIN
         FOR i := 1 TO ndata DO
           chi2 := chi2 + Sqr((GetMatrixElement(Data, i, 2) - a - b *
                   GetMatrixElement(Data, i, 1)) / GetMatrixElement(Data, i, 3));
         q := IntegralChi(chi2, ndata-2);
       END;
END;


PROCEDURE Deming (CONST x, y : VectorTyp; delta : double;
                  VAR xCalc, yCalc : VectorTyp; VAR Res : ResultTyp);

VAR i, n, s                                : WORD;
    Sx, Sy, sxx, syy, sxy, xMean, yMean,b_ : double;
    Significance                           : SignificanceType;
    c                                      : CHAR;

BEGIN
  IF VectorLength(x) <> VectorLength(y)
    THEN
      BEGIN
        RegressionError := TRUE;
        c := WriteErrorMessage('Deming regression: unequal length of data vectors');
        EXIT;
      END;
  n := VectorLength(x);
  s := 0;                                                   // calculate means
  Sx := 0;
  Sy := 0;
  FOR i := 1 TO n DO
    BEGIN
      IF IsNaN(GetVectorElement(x, i)) OR IsNaN(GetVectorElement(y, i))
        THEN
        ELSE
          BEGIN
            INC(s);                                        // valid data pairs
            Sx := Sx + GetVectorElement(x, i);
            Sy := Sy + GetVectorElement(y, i);
          END;
    END;
  xMean := Sx / s;                                    // calculate 2nd moments
  yMean := Sy / s;
  sxx := 0;
  syy := 0;
  sxy := 0;
  FOR i := 1 TO n DO
    BEGIN
      IF IsNaN(GetVectorElement(x, i)) OR IsNaN(GetVectorElement(y, i))
        THEN
        ELSE
          BEGIN
            sxx := sxx + Sqr(GetVectorElement(x, i) - xMean);
            syy := syy + Sqr(GetVectorElement(y, i) - yMean);
            sxy := sxy + (GetVectorElement(x, i) - xMean) * (GetVectorElement(y, i) - yMean)
          END;
    END;
   sxx := sxx/s; // sample variance x
   syy := syy/s; // sample variance y
   sxy := sxy/s; // sample covariance x,y                  // calculate params
   Res.b := (syy - delta*sxx + Sqrt((syy-delta*sxx)*(syy-delta*sxx) + 4*delta*sxy*sxy))
             / (2*sxy);
   Res.sb := NaN; // Deming regression doesn't provide error estimates
   Res.a := yMean - Res.b*xMean;
   Res.sa := NaN;
   Res.r := QuadrantCorrelation(x, y, Significance);
   Res.P0 := Significance.P0;
   Res.t := Significance.Testvalue; // actually chi^2, not t
   CreateVector(xCalc, n, 0.0);                    // calculate xCalc and yCalc
   CreateVector(yCalc, n, 0.0);
   b_ := Res.b/(Res.b*Res.b + delta);
   FOR i := 1 TO n DO
     BEGIN
       IF IsNaN(GetVectorElement(x, i)) OR IsNaN(GetVectorElement(y, i))
         THEN
         ELSE
           BEGIN
             SetVectorElement(xCalc, i, GetVectorElement(x, i) + (GetVectorElement(y, i) -
                              Res.a - Res.b * GetVectorElement(x, i)));
             SetVectorElement(yCalc, i, Res.a + Res.b * GetVectorElement(xCalc, i));
           END;
     END;
END;


PROCEDURE TheilSenKendall (CONST x, y : VectorTyp; VAR yCalc : VectorTyp;
                           VAR Res : ResultTyp);

VAR i, j, n, s                                                       : WORD;
    Slopes, Slopes_big, Intercepts, Intercepts_big, xSorted, ySorted : VectorTyp;
    x1, x2, y1, y2, xmin, xmax, xdiff, Sx, Sy, xMed, yMed            : double;
    Significance                                                     : SignificanceType;
    c                                                                : CHAR;
    unknown                                                          : BOOLEAN;

BEGIN
  IF VectorLength(x) <> VectorLength(y)
    THEN
      BEGIN
        RegressionError := TRUE;
        c := WriteErrorMessage('Linear regression: unequal Length OF dependent AND independent data vector');
        EXIT;
      END;
  n := VectorLength(x);
  s := 0;
  CreateVector(Slopes_big, Round(n*(n-1)/2 + 0.5), 0.0);   // maximal possible number OF slopes
  xmax := FindLargest(x);
  xmin := FindSmallest(x);
  xdiff := (xmax - xmin) / (5 * n);  // factor 5 IS arbitrary
  FOR i := 1 TO n DO
     FOR j := Succ(i) TO n DO
       BEGIN
         x1 := GetVectorElement(x, i);
         x2 := GetVectorElement(x, j);
         y1 := GetVectorElement(y, i);
         y2 := GetVectorElement(y, j);
         unknown :=  IsNaN(x1) OR IsNaN(x2) OR  IsNaN(y1) OR IsNaN(y2);
         IF (Abs(x1 - x2) < xdiff) OR unknown
           THEN
           ELSE
             BEGIN
               INC(s);
               SetVectorElement(Slopes_big, s, (y1-y2)/(x1-x2));
             END;
       END;
   CreateVector(slopes, s, 0.0);
   FOR i := 1 TO s DO  // remove any empty values from slope vector
     SetVectorElement(slopes, i, GetVectorElement(slopes_big, i));
   DestroyVector(Slopes_big);
   Res.b := Median(slopes);
   Res.sb := (Quantile(slopes, 0.75) - Quantile(slopes, 0.25)) / 2;
   CopyVector(x, xSorted);
   xMed := Median(xSorted);
   CopyVector(y, ySorted);
   yMed := Median(ySorted);
   Res.a := yMed - Res.b * xMed; // most stable estimator FOR intercept
   CreateVector(Intercepts_big, n, 0.0);
   s := 0;
   FOR i := 1 TO n DO  // calculate intercepts FOR all x/y pairs
     IF IsNan(GetVectorElement(x, i)) OR IsNaN(GetVectorElement(y, i))
       THEN
       ELSE
         BEGIN
           SetVectorElement(Intercepts_big, i, GetVectorElement(y, i) - Res.b * GetVectorElement(x, i));
           INC(s);
         END;
   CreateVector(Intercepts, s, 0.0);
   FOR i := 1 TO s DO  // remove any empty values from intercept vector
     SetVectorElement(Intercepts, i, GetVectorElement(Intercepts_big, i));
   DestroyVector(Intercepts_big);
   Res.sa := (Quantile(intercepts, 0.75) - Quantile(intercepts, 0.25)) / 2;
   Res.r := QuadrantCorrelation(x, y, Significance);
   Res.P0 := Significance.P0;
   Res.t := Significance.Testvalue; // actually chi^2, NOT t
   CreateVector(yCalc, n, 0.0);
   FOR i := 1 TO n DO
     SetVectorElement(yCalc, i, Res.a + Res.b * GetVectorElement(x, i));
   DestroyVector(xSorted);
   DestroyVector(ySorted);
   DestroyVector(Slopes);
   DestroyVector(Intercepts);
END;

PROCEDURE Eisenthal (CONST Substrate, Velocity : VectorTyp; VAR yCalc : VectorTyp;
                           VAR Res : ResultTyp);

VAR Km, Vmax, Km_big, Vmax_big : VectorTyp;
    i, j, n, s                 : WORD;
    c                          : CHAR;
    SI, vi, Sj, vj             : double;
    Significance               : SignificanceType;

BEGIN
  IF VectorLength(Substrate) <> VectorLength(Velocity)
    THEN
      BEGIN
        RegressionError := TRUE;
        c := WriteErrorMessage('Direct plot: unequal Length OF dependent AND independent data vector');
        EXIT;
      END;
  n := VectorLength(Substrate);
  s := 0;
  CreateVector(Km_big, Round(n*(n-1)/2 + 0.5), 0.0);   // maximal possible number OF pairs
  CreateVector(Vmax_big, Round(n*(n-1)/2 + 0.5), 0.0);
  FOR i := 1 TO n DO
     FOR j := Succ(i) TO n DO
       BEGIN
         SI := GetVectorElement(Substrate, i);
         Vi := GetVectorElement(Velocity, i);
         Sj := GetVectorElement(Substrate, j);
         Vj := GetVectorElement(Velocity, j);
         IF (IsNaN(SI) OR IsNaN(vi) OR  IsNaN(Sj) OR IsNaN(vj))
           THEN
             // ignore
           ELSE
             BEGIN
               INC(s);
               SetVectorElement(Km_big,   s, (vj - vi) / (vi/SI - vj/Sj));
               SetVectorElement(Vmax_big, s, (SI - Sj) / (SI/vi - Sj/vj));
             END;
       END;
   CreateVector(Km, s, 0.0);
   CreateVector(Vmax, s, 0.0);
   FOR i := 1 TO s DO  // remove any empty values from result vectors
     BEGIN
       SetVectorElement(Km, i, GetVectorElement(Km_big, i));
       SetVectorElement(Vmax, i, GetVectorElement(Vmax_big, i));
     END;
   DestroyVector(Km_big);
   DestroyVector(Vmax_big);
   Res.a := Median(Km);
   Res.sa := (Quantile(Km, 0.75) - Quantile(Km, 0.25)) / 2;
   Res.b := Median(Vmax);
   Res.sb := (Quantile(Vmax, 0.75) - Quantile(Vmax, 0.25)) / 2;
   FOR i := 1 TO s DO
     Writeln(i:3, ' ', GetVectorElement(Km, i):3:3, ' ', GetVectorElement(Vmax, i):3:3);
   Res.r := QuadrantCorrelation(Substrate, Velocity, Significance);
   Res.P0 := Significance.P0;
   Res.t := Significance.Testvalue; // actually chi^2, NOT t
   CreateVector(yCalc, n, 0.0);
   FOR i := 1 TO n DO
     BEGIN
       SI := GetVectorElement(Substrate, i);
       IF IsNaN(SI)
         THEN
           SetVectorElement(yCalc, i, NaN)
         ELSE
           SetVectorElement(yCalc, i, Res.b * SI / (Res.a + SI));
     END;
  DestroyVector(Km);
  DestroyVector(Vmax);
END;

end.







