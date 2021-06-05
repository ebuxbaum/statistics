PROGRAM TestCircular;

USES math,         // standard math LIBRARY
     MathFunc,     // mathematical functions
     Complex,      // complex numbers
     Vector,       // vector arithmetic
     Matrix,       // matrix arithmetic
     Zufall,       // Random numbers
     Stat,         // statistical significance
     Deskript,     // descriptive statistics
     Nonparam,     // non-parametric tests
     Circular;     // circular statistics

CONST Length =   67;
      Steps  =   30;

TYPE Counts = ARRAY [0..Steps] OF WORD;

VAR xR, yR, zR, xTheta, yTheta, zTheta, a : float;
    xMean, yMean, zMean, xMean2, yMean2,
    zMean2                                : ComplexTyp;
    xVec, yVec, zVec, aVec                : VectorTyp;
    i                                     : WORD;
    xCounts, yCounts, zCounts             : Counts;

  PROCEDURE CreateRandomVectors (VAR  xVec, yVec, zVec : VectorTyp;
                                 VAR  xCounts, yCounts, zCounts : Counts);
  VAR i       : WORD;
      x, y, z : float;

  BEGIN
    CreateVector(xVec, Length, 0.0);
    CreateVector(yVec, Length, 0.0);
    CreateVector(zVec, Length, 0.0);
    FOR i := 1 TO Length DO
      BEGIN
        x := RandomVonMieses(0.0, 1.0);
        SetVectorElement(xVec, i, x);
        INC(xCounts[Round(x/a)]);
        y := RandomVonMieses(3.0, 10.0);
        SetVectorElement(yVec, i, y);
        INC(yCounts[Round(y/a)]);
        z := RandomUniformCircular;
        SetVectorElement(zVec, i, z);
        INC(zCounts[Round(z/a)]);
      END;
  END;


  PROCEDURE ReadRandomVectors (VAR  xVec, yVec, zVec : VectorTyp;
                               VAR  xCounts, yCounts, zCounts : Counts);

  VAR i       : WORD;
      x, y, z : float;
      InFile  : TEXT;

  BEGIN
    CreateVector(xVec, Length, 0.0);
    CreateVector(yVec, Length, 0.0);
    CreateVector(zVec, Length, 0.0);
    ASSIGN(InFile, 'Kreis.csv');
    TRY
      RESET(InFile);
    EXCEPT
      Write('could not open file "Kreis.csv"');
      ReadLn;
      HALT;
    END;
    ReadLn(InFile); // headline
    FOR i := 1 TO Length DO
      BEGIN
        x := ReadFloat(Infile);
        SetVectorElement(xVec, i, x);
        INC(xCounts[Round(x/a)]);
        y := ReadFloat(Infile);
        SetVectorElement(yVec, i, y);
        INC(yCounts[Round(y/a)]);
        z := ReadFloat(Infile);
        SetVectorElement(zVec, i, z);
        INC(zCounts[Round(z/a)]);
        ReadLn(InFile);
      END;
    CLOSE(InFile);
  END;

PROCEDURE ReadBP (VAR  Time, Sys, Dia, Pulse, PP : VectorTyp);

VAR i       : WORD;
    InFile  : TEXT;

BEGIN
  CreateVector(Time, Length, 0.0);
  CreateVector(Sys, Length, 0.0);
  CreateVector(Dia, Length, 0.0);
  CreateVector(Pulse, Length, 0.0);
  CreateVector(PP, Length, 0.0);
  ASSIGN(InFile, 'Blutdruck.csv');
  TRY
    RESET(InFile);
  EXCEPT
    Write('could not open file "Blutdruck.csv"');
    ReadLn;
    HALT;
  END;
  ReadLn(InFile); // headline
  FOR i := 1 TO Length DO
    BEGIN
      SetVectorElement(Time, i, ReadFloat(Infile));
      SetVectorElement(Sys, i, ReadFloat(Infile));
      SetVectorElement(Dia, i, ReadFloat(Infile));
      SetVectorElement(Pulse, i, ReadFloat(Infile));
      SetVectorElement(PP, i, ReadFloat(Infile));
      ReadLn(InFile);
    END;
  CLOSE(InFile);
END;


BEGIN
  inc(ValidFigures);
  a := Const_2pi/Steps;
  CreateVector(aVec, Length, 0.0);
  FOR i := 0 TO Steps DO
    BEGIN
      xCounts[i] := 0;
      yCounts[i] := 0;
      zCounts[i] := 0;
    END;
//  CreateRandomVectors(xVec, yVec, zVec, xCounts, yCounts, zCounts);
  ReadRandomVectors(xVec, yVec, zVec, xCounts, yCounts, zCounts);
  Writeln('        x           y           z');
  xCounts[Steps] := xCounts[Steps] + xCounts[0];  // deal WITH cross-over
  yCounts[Steps] := yCounts[Steps] + yCounts[0];
  zCounts[Steps] := zCounts[Steps] + zCounts[0];
  FOR i := 1 TO Steps DO
    Writeln(a*(Pred(i)+i)/2:1:3, ' ', xCounts[i]:4, '        ', yCounts[i]:4, '        ', zCounts[i]:4);
  Writeln;
  Writeln('Median:  ', FloatStr(MedianDirection(xVec), ValidFigures), '  ',
                       FloatStr(MedianDirection(yVec), ValidFigures), '  ',
                       FloatStr(MedianDirection(zVec), ValidFigures));
  xMean := MeanVector(xVec,1);
  yMean := MeanVector(yVec,1);
  zMean := MeanVector(zVec,1);
  xR    := Re(xMean);
  xTheta  := Im(xMean);
  yR    := Re(yMean);
  yTheta  := Im(yMean);
  zR    := Re(zMean);
  zTheta  := Im(zMean);
  WriteLn('R:       ', FloatStr(xR, ValidFigures), '  ',
                       FloatStr(yR, ValidFigures), '  ',
                       FloatStr(zR, ValidFigures));
  WriteLn('Theta:   ', FloatStr(xTheta, ValidFigures), '  ',
                       FloatStr(yTheta, ValidFigures), '  ',
                       FloatStr(zTheta, ValidFigures));
  writeLn('Variance ', FloatStr(CircularVariance(xR), ValidFigures), '  ',
                       FloatStr(CircularVariance(yR), ValidFigures), '  ',
                       FloatStr(CircularVariance(zR), ValidFigures));
  writeLn('std.dev. ', FloatStr(CircularStandardDeviation(xR), ValidFigures), '  ',
                       FloatStr(CircularStandardDeviation(yR), ValidFigures), '  ',
                       FloatStr(CircularStandardDeviation(zR), ValidFigures));
  WriteLn('kappa:   ', FloatStr(Kappa(xR, Length), ValidFigures), '  ',
                       FloatStr(Kappa(yR, Length), ValidFigures), '  ',
                       FloatStr(Kappa(zR, Length), ValidFigures));
  xMean2 := CenteredMean(TrigonometricMoment(xVec, xTheta, 2));
  yMean2 := CenteredMean(TrigonometricMoment(yVec, yTheta, 2));
  zMean2 := CenteredMean(TrigonometricMoment(zVec, zTheta, 2));
  WriteLn('skew:    ', FloatStr(CenteredCircularSkew(xMean, xMean2), ValidFigures), '  ',
                       FloatStr(CenteredCircularSkew(yMean, yMean2), ValidFigures), '  ',
                       FloatStr(CenteredCircularSkew(zMean, zMean2), ValidFigures));
  WriteLn('kurtosis:', FloatStr(CenteredCircularKurtosis(xMean, xMean2), ValidFigures), '  ',
                       FloatStr(CenteredCircularKurtosis(yMean, yMean2), ValidFigures), '  ',
                       FloatStr(CenteredCircularKurtosis(zMean, zMean2), ValidFigures));
  WriteLn('delta   :', FloatStr(CircularDispersion(xVec, xMean), ValidFigures), '  ',
                       FloatStr(CircularDispersion(yVec, yMean), ValidFigures), '  ',
                       FloatStr(CircularDispersion(zVec, zMean), ValidFigures));

  WriteLn('Raileigh  ', 100*Rayleigh(Re(xMean), Length):3:3, '%      ', 100*Rayleigh(Re(yMean), Length):3:3,
          '%      ', 100*Rayleigh(Re(zMean), Length):3:3, '%');
  WriteLn('Hodges:   ', 100*HodgesAjne(xVec, 0.0):3:3, '%      ', 100*HodgesAjne(yVec, 3.0):1:3,
          '%      ', 100*HodgesAjne(zVec, 0.0):1:3, '%');
  Write('Press return:');
  ReadLn;
  DestroyVector(xVec);
  DestroyVector(yVec);
  DestroyVector(zVec);
END.

