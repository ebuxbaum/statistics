UNIT Correlations;

{ Various types of correlation coefficients. Missing data should be encoded
  NaN. Nominal, binary and ordinal values are given as integer numbers,
  but in vectors and fields of type double (alternatively, the type definitions
  in units Vector and Matrix would have to be overlayed).

  The record Significance contains the appropriate Testvariable, degrees of
  freedom and probability for H_0: r = 0 vs H_1: r <> 0 }

INTERFACE

USES Math, MathFunc, Vector, Stat, Deskript;

CONST
  MaxCases = 1360;
  MaxSteps = 100; // maximal # OF different values FOR ordinal AND nominal data
  CorrelationsError: BOOLEAN = FALSE;

TYPE
  ContStruc = RECORD
    Table: ARRAY [0..MaxSteps, 0..MaxSteps] OF WORD;
    RowSums, ColumnSums: ARRAY [0..MaxSteps] OF WORD;
    r, c, n: WORD;
    // number OF rows, columns AND legal (no NaN) comparisons
  END;
  ContTable = ^ContStruc;

  Confusion = RECORD
    a, b, c, d: WORD;
  END;

  SignificanceType = RECORD
    TestValue, Freedom, P0: float;
  END;

PROCEDURE CorrelationSignificance(CONST Vector1, Vector2: VectorTyp;
  r: float; VAR Significance: SignificanceType);
{ t-test for significance of r_p, r_s, r_j  }

FUNCTION PearsonProductMomentCorrelation(CONST Vector1, Vector2: VectorTyp;
  VAR Significance: SignificanceType): float;
{ Pearson's product moment correlation for two data vectors that are both at
  least interval scale. n is number of cases.
  Pearson: Proc. Royal Soc. London 58 (1895) 240-2}

FUNCTION WeightedPearson(CONST Vector1, Vector2, weight: VectorTyp;
  VAR Significance: SignificanceType): float;
{ use a weight vector to emphasise or deemphasise certain data pairs }

PROCEDURE Rank(VAR Data: VectorTyp);
{ converts an interval or rational data vector into ranks, taking care of ties }

FUNCTION SpearmanRankCorrelation(CONST Vector1, Vector2: VectorTyp;
  VAR Significance: SignificanceType): float;
{ Spearman's rank correlation for two interval scaled vectors with no ties.
  Spearman: Am. J. Psychol. 1 (1904) 72-101
  Internal copies of Vectors are ranked, then PearsonProductMomentCorrelation is called.
  Ties are handled. }

FUNCTION QuadrantCorrelation(CONST Vector1, Vector2: VectorTyp;
  VAR Significance: SignificanceType): float;
{ Nonparametric correlation between two cardinal vectors, does not require linearity.
  A chi^2 test is performed for the significance of r_q }

FUNCTION PointBiserialCorrelation(CONST BinaryVector, CardinalVector: VectorTyp;
  VAR Significance: SignificanceType): float;
{ Correlation in dichotomous exams. BinaryVector is a vector of 0/1 of one question,
  CardinalVector contains the overall grades.}

FUNCTION ConfusionTable(CONST Vector1, Vector2: VectorTyp): Confusion;
{ a : joint presence, d: joint absence, b or c : one case only }

FUNCTION McConnaugheyCorrelation(c: confusion): float;
{ Correlation between two 0/1 vectors according to McConnaughey: Penelitian
  Laut di Indonesia [Marine Res. Indonesia] 7 (1964) 1-40. According to Hubalek:
  Biol. Rev. 57:4 (1982) 669–89 this is the only coefficient of 42 suggested in
  the literature that meets desirable statistical properties and maps to -1..1 }

FUNCTION MatthewsCorrelation(c: Confusion): float;
{ Characterises confusion tables. Matthews: Biochim. Biophys. Acta 405:2 (1975) 442 }

FUNCTION Rk(Cont: ContTable): float;
{ Expansion of Matthews correlation to k times k tables,
  J. Gorodkin: Comput. Biol. Chem. 28:5 (2004) 367–374. }

PROCEDURE Contingency(Data1, Data2: VectorTyp; VAR Cont: ContTable);
{ Generates a contingency table for two data vectors of length Length. }

PROCEDURE WriteContingency(CONST Cont: ContTable; MedStr: STRING; ValidFigures: BYTE);
{ Writes the contingency table Cont onto the file named MedStr (CON: for screen).
  ValidFigures is the number of places for each figure. }

PROCEDURE DestroyContingency(VAR Cont: ContTable);
{ releases the memory occupied by Cont, its value becomes nil. }

FUNCTION Chi2(CONST Contingency: ContTable;
  VAR Significance: SignificanceType): float;
{ Calculates chi^2 for a contingency table
  Pearson, K.: Philosophical Magazine 50:302 (1900) 157-175}

FUNCTION Phi2(CONST Contingency: ContTable;
  VAR Significance: SignificanceType): float;
{ phi^2 = chi^2 / n }

FUNCTION CramersVTilde(chi2: float; r, c, n: WORD): float;
{ Cramér's V~ with bias correction, given  chi^2, r rows, c columns and n data.
  Cramér, H.: Mathematical Methods of Statistics. Princeton (Princeton University Press)
  1946, ISBN 0-691-08004-6, pp. 282 }

FUNCTION lambda(CONST Contingency: ContTable): float;
{  Symmetric coefficient of relative predictability for nominal/nominal association
   Goodman & Kruskal: J. Am. Statistical Assoc. 49:268 (1954) 732-764}

FUNCTION lambda_r(CONST Contingency: ContTable): float;
{ Asymmetric coefficient of relative predictability for nominal/nominal association
  for row case }

FUNCTION lambda_c(CONST Contingency: ContTable): float;
{ Asymmetric coefficient of relative predictability for nominal/nominal association
  for column case }

FUNCTION OrdinalCorrelations(CONST Data1, Data2: VectorTyp; Formula: CHAR;
  VAR Significance: SignificanceType): float;
{ calculates correlation between ordinal/ordinal and ordinal/binary data vectors.
  case Formula = T return Kendall's tau
                 G return Goodman & Kruskal's gamma
                 D return Kim's d_(xy) (which is <> d_(yx), if need be, rerun with exchanged vectors)
                 E return Wilson's e
  These measures are distinguished by how ties in Data1 and Data2 are handled.
  Literatur: Freeman, L.C: J. Math. Sociol. 12:1 (1986) 49–69.}

FUNCTION theta(CONST Contingency: ContTable): float;
{ calculates theta for nominal/ordinal association.
  Literatur: Freeman, L.C.: Elementary applied Statistics. New York (Wiley) 1965 }

FUNCTION eta_sqr(CONST NominalVector, CardinalVector: VectorTyp;
  VAR Significance: SignificanceType): float;
{ nominal/cardinal associations. Nominal data  als continuous numbers 1..Classes
  Literatur https://de.wikipedia.org/wiki/Fehlerreduktionsmaße 2017-07-29}

{ ************************* calculations with r ***************************** }

FUNCTION AverageCorrelations(CONST Correlations: VectorTyp): float;
{ Average of k correlation coefficients
  Bushman & Wang: Psychol. Bulletin 117:3 (1995) 530,
  Corey et al.: J. Gen. Psychol. 125:3 (1998) 245-61 }

FUNCTION Angle(r: float): float;
{ convert a correlation coefficient to the angle between vectors (in degrees) }


IMPLEMENTATION

FUNCTION TestDataVectorLength(Data1, Data2: VectorTyp): BOOLEAN;

VAR c : char;

BEGIN
  IF (VectorLength(Data1) <> VectorLength(Data2))
    THEN
      BEGIN
        c := WriteErrorMessage(' Correlation: vectors of unequal length');
        Result := FALSE;
      END
    ELSE
      Result := TRUE;
END;

PROCEDURE CorrelationSignificance(CONST Vector1, Vector2: VectorTyp;
  r: float; VAR Significance: SignificanceType);

VAR
  i, n, m: WORD;

BEGIN
  n := VectorLength(Vector1);
  IF NOT (TestDataVectorLength(Vector1, Vector2))
    THEN EXIT;
  m := 0;
  FOR i := 1 TO n DO
    IF (IsNaN(GetVectorElement(Vector1, i)) OR IsNaN(GetVectorElement(Vector1, i)))
      THEN // ignore
      ELSE INC(m);
  Significance.TestValue := (m - 2) / (1 - r * r);
  IF signum(Significance.TestValue) >= 0
    THEN Significance.TestValue := r * Sqrt(Significance.TestValue)
    ELSE Significance.TestValue := 0;
  Significance.Freedom := m - 2;
  Significance.P0 := Integral_t(Significance.TestValue, Round(Significance.Freedom));
END;


FUNCTION PearsonProductMomentCorrelation(CONST Vector1, Vector2: VectorTyp;
  VAR Significance: SignificanceType): float;

VAR
  i, j: WORD;
  SumXY, SumX, SumY, SumX2, SumY2, varX, varY, covXY, x, y, r: float;

BEGIN
  IF NOT (TestDataVectorLength(Vector1, Vector2)) THEN EXIT;
  SumXY := 0;
  SumX := 0;
  SumY := 0;
  SumX2 := 0;
  SumY2 := 0;
  j := 0;
  FOR i := 1 TO VectorLength(Vector1) DO
    BEGIN
      x := GetVectorElement(Vector1, i);
      y := GetVectorElement(Vector2, i);
      IF (IsNaN(x) OR IsNaN(y))
        THEN // ignore
        ELSE
          BEGIN
            SumXY := SumXY + x * y;
            SumX := SumX + x;
            SumY := SumY + y;
            SumX2 := SumX2 + x * x;
            SumY2 := SumY2 + y * y;
            INC(j); // actual number OF valid comparisons
          END; { else }
    END; { for }
  varX := j * SumX2 - SumX * SumX;
  varY := j * SumY2 - SumY * SumY;
  covXY := j * SumXY - SumX * SumY;
  IF (varX * varY) = 0
    THEN r := 0     // ???
    ELSE r := covXY / Sqrt(varX * varY);
  CorrelationSignificance(Vector1, Vector2, r, Significance);
  Result := r;
END; { PearsonProductMoment }

FUNCTION WeightedPearson(CONST Vector1, Vector2, weight: VectorTyp;
  VAR Significance: SignificanceType): float;

VAR
  i, n: WORD;
  BarX, BarY, SumWX, SumWY, SumW, r, SumWXd, SumWYd, SumWXYd, CovXY,
  VarX, VarY, x, y, w: float;

BEGIN
  IF NOT (TestDataVectorLength(Vector1, Vector2)) THEN EXIT;
  n := VectorLength(Vector1);
  SumWX := 0;
  SumWY := 0;
  SumW := 0;
  SumWXd := 0;
  SumWYd := 0;
  SumWXYd := 0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(Vector1, i);
      y := GetVectorElement(Vector2, i);
      w := GetVectorElement(Weight, i);
      IF (IsNaN(x) OR IsNaN(x) OR IsNaN(w))
        THEN // ignore
        ELSE
          BEGIN
            SumWX := SumWX + w * x;
            SumWY := SumWY + w * y;
            SumW := SumW + w;
          END; { else }
    END; { for }
  BarX := SumWX / SumW;     // calculate averages
  BarY := SumWY / SumW;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(Vector1, i);
      y := GetVectorElement(Vector2, i);
      w := GetVectorElement(Weight, i);
      IF (IsNaN(x) OR IsNaN(y) OR IsNaN(w))
        THEN  // ignore
        ELSE
          BEGIN
            SumWXd := SumWXd + w * (x - BarX) * (x - BarX);
            SumWYd := SumWYd + w * (y - BarY) * (y - BarY);
            SumWXYd := SumWXYd + w * (x - BarX) * (y - BarY);
          END; { else }
    END; { for }
  CovXY := SumWXYd / SumW;
  VarX := SumWXd / SumW;
  VarY := SumWYd / SumW;
  IF (varX * varY) = 0
    THEN r := 0   // ???
    ELSE r := CovXY / Sqrt(VarX * VarY);
  CorrelationSignificance(Vector1, Vector2, r, Significance);
  Result := r;
END; { WeightedPearson }

PROCEDURE Rank(VAR Data: VectorTyp);

LABEL
  done;

VAR
  i, j, k, Valid, n: WORD;
  r, x, y: float;
  Sorted: VectorTyp;

BEGIN
  n := VectorLength(Data);
  CopyVector(Data, Sorted);
  ShellSort(Sorted);            // create a sorted Copy OF data  (any NaNs TO highest #)
  Valid := n;
  FOR i := 1 TO n DO
    IF IsNaN(GetVectorElement(Sorted, i))
      THEN DEC(Valid);         // determine number OF non-NaN data points
  i := 1;
  REPEAT
    x := GetVectorElement(Sorted, i);
    IF IsNaN(x)
      THEN
      ELSE
        BEGIN
          j := i;
          REPEAT             // find number OF Elements that are equal TO x
            y := GetVectorElement(Sorted, Succ(j));
            IF IsNaN(y)
              THEN GOTO Done // leave counting LOOP, no more valid data
              ELSE IF (x = y)
                     THEN INC(j);
          UNTIL ((x < y) OR (j = n));
Done:     r := (1.0 * i + j) / 2;     // rank IS average OF lowest + highest element
          FOR k := 1 TO Valid DO      // change all elements OF Data equal TO x -> rank
            BEGIN
              y := GetVectorElement(Data, k);
              IF IsNaN(y)
                THEN
                ELSE IF (y = x)
                       THEN SetVectorElement(Data, k, r);
            END;
          i := Succ(j);
        END;
  UNTIL (i >= Valid) OR IsNaN(x);
  DestroyVector(Sorted);
END;


FUNCTION SpearmanRankCorrelation(CONST Vector1, Vector2: VectorTyp;
  VAR Significance: SignificanceType): float;

VAR
  Ranked1, Ranked2: VectorTyp;

BEGIN
  IF NOT (TestDataVectorLength(Vector1, Vector2)) THEN EXIT;
  CopyVector(Vector1, Ranked1);
  CopyVector(Vector2, Ranked2);
  Rank(Ranked1);
  Rank(Ranked2);
  Result := PearsonProductMomentCorrelation(Ranked1,
    Ranked2, Significance);
  DestroyVector(Ranked1);
  DestroyVector(Ranked2);
END;

FUNCTION QuadrantCorrelation(CONST Vector1, Vector2: VectorTyp;
  VAR Significance: SignificanceType): float;

VAR
  Sorted1, Sorted2            : VectorTyp;
  n, i, n_plus, n_minus       : WORD;
  Median1, Median2, x, y, n_e : float;
  c                           : char;

BEGIN
  IF NOT (TestDataVectorLength(Vector1, Vector2)) THEN EXIT;
  n := VectorLength(Vector1);
  CopyVector(Vector1, Sorted1);
  CopyVector(Vector2, Sorted2);
  Median1 := Median(Sorted1);
  Median2 := Median(Sorted2);
  DestroyVector(Sorted1);
  DestroyVector(Sorted2);
  n_plus := 0;
  n_minus := 0;
  FOR i := 1 TO n DO
    BEGIN
      x := GetVectorElement(Vector1, i);
      y := GetVectorElement(Vector2, i);
      IF IsNaN(x) OR IsNaN(y)
        THEN // ignore NaNs
        ELSE IF ((x = Median1) OR (y = Median2))
               THEN // ignore data ON coordinate axis OF system WITH origin Median1/Median2
               ELSE IF (((x > Median1) AND (y > Median2)) OR ((x < Median1) AND (y < Median2)))
                      THEN INC(n_plus)   // concordance: quadrant I OR III
                      ELSE INC(n_minus); // discordance: quadrant II OR IV
    END;
  Result := (n_plus - n_minus) / (n_plus + n_minus);
  n_e := 1.0 * (n_plus + n_minus) / 2;
  IF (n_e < 5)
    THEN
      BEGIN
        c := WriteErrorMessage(' Quadrant correlation: cannot calculate significance because n_e < 5');
        CorrelationsError := TRUE;
        EXIT;
      END;
  Significance.TestValue := ((n_plus - n_e) * (n_plus - n_e) +
    (n_minus - n_e) * (n_minus - n_e)) / n_e;
  Significance.Freedom := 1;
  Significance.P0 := IntegralChi(Significance.TestValue, Round(Significance.Freedom));
END;


FUNCTION PointBiserialCorrelation(CONST BinaryVector, CardinalVector: VectorTyp;
  VAR Significance: SignificanceType): float;

VAR
  n, i, TotalI, GoodI, BadI                 : WORD;
  SumX, SumXGood, SumXBad, s, barX, x, y, r : float;
  c                                         : CHAR;

BEGIN
  IF NOT (TestDataVectorLength(BinaryVector, CardinalVector)) THEN EXIT;
  n := VectorLength(CardinalVector);
  TotalI := 0;
  SumX := 0;
  FOR i := 1 TO n DO       // calculate average x
    IF NOT (IsNaN(GetVectorElement(CardinalVector, i)) OR
            IsNaN(GetVectorElement(BinaryVector, i)))
      THEN
        BEGIN
          INC(TotalI);
          SumX := SumX + GetVectorElement(CardinalVector, i);
        END;
  BarX := SumX / TotalI;
  SumX := 0;
  FOR i := 1 TO n DO       // calculate standard deviation OF x
    IF NOT (IsNaN(GetVectorElement(CardinalVector, i)) OR
            IsNaN(GetVectorElement(BinaryVector, i)))
      THEN
        SumX := SumX + (GetVectorElement(CardinalVector, i) - barX) *
               (GetVectorElement(CardinalVector, i) - barX);
  s := Sqrt(SumX / TotalI);
  GoodI := 0;
  BadI := 0;
  SumXGood := 0;
  SumXBad := 0;
  FOR i := 1 TO n DO
    IF (IsNaN(GetVectorElement(CardinalVector, i)) OR
        IsNaN(GetVectorElement(BinaryVector, i)))
      THEN  // ignore
      ELSE IF (GetVectorElement(BinaryVector, i) = 1)
             THEN
               BEGIN
                 INC(GoodI);
                 SumXGood := SumXGood + GetVectorElement(CardinalVector, i);
               END
             ELSE IF (GetVectorElement(BinaryVector, i) = 0)
                    THEN
                      BEGIN
                        INC(BadI);
                        SumXBad := SumXBad + GetVectorElement(CardinalVector, i);
                      END
                    ELSE
                      BEGIN
                        c := WriteErrorMessage('Point biserial correlation: illegal value in exam table');
                        CorrelationsError := TRUE;
                        EXIT;
                      END;
  IF ((GoodI = 0) OR (BadI = 0))
    THEN
      r := 0
    ELSE
      BEGIN
        SumXGood := SumXGood / GoodI;
        SumXBad := SumXBad / BadI;
        x := Sqrt(GoodI * BadI / TotalI / TotalI);
        y := (SumXGood - SumXBad) / s;
        r := x * y;
      END;
  CorrelationSignificance(BinaryVector, CardinalVector, r, Significance);
  Result := r;
END;

FUNCTION ConfusionTable(CONST Vector1, Vector2: VectorTyp): Confusion;

VAR
  i, n: WORD;      // contingency table values

BEGIN
  IF NOT (TestDataVectorLength(Vector1, Vector2)) THEN EXIT;
  n := VectorLength(Vector1);
  WITH ConfusionTable DO
    BEGIN
      a := 0;
      b := 0;
      c := 0;
      d := 0;
      FOR i := 1 TO n DO
        IF (IsNaN(GetVectorElement(Vector1, i)) OR IsNaN(GetVectorElement(Vector2, i)))
          THEN // ignore
          ELSE IF (GetVectorElement(Vector1, i) = 1)   // count contingency table elements
                 THEN IF (GetVectorElement(Vector2, i) = 1)
                        THEN INC(a)  // joint presence
                        ELSE INC(b)  // discordance
                 ELSE IF (GetVectorElement(Vector2, i) = 1)
                        THEN INC(c)  // discordance
                        ELSE INC(d); // joint absence
    END;

END;

FUNCTION McConnaugheyCorrelation(c: confusion): float;

BEGIN
  WITH c DO
    BEGIN
      Result := (1.0 * (a + b) * (a + c));
      IF (Result = 0)
        THEN // shouldn't happen, but result will be 0
        ELSE Result := (1.0 * a * a - 1.0 * b * c) / Result;
    END;
END;

FUNCTION MatthewsCorrelation(c: Confusion): float;

BEGIN
  WITH c DO
    BEGIN
      Result := 1.0 * (a + c) * (a + b) * (d + c) * (d + b);
      // calculate denominator first
      IF Result <= 0
        THEN // shouldn't happen, but result will be 0, the correct limiting result
        ELSE
        Result := (a * d - c * b) / Sqrt(Result);
    END;
END;

FUNCTION Rk(Cont: ContTable): float;

VAR
  k, i                   : WORD;
  Sum1, Sum2, Sum3, Sum4 : float;
  c                      : char;

BEGIN
  k := Cont^.r;
  IF k <> Cont^.c
    THEN
      BEGIN
        c := WriteErrorMessage('Gorodkin''s Rk: Contingency table not square ');
        CorrelationsError := TRUE;
        EXIT;
      END;
  Sum1 := 0;
  Sum2 := 0;
  Sum3 := 0;
  Sum4 := 0;
  FOR i := 1 TO k DO
    BEGIN
      Sum1 := Sum1 + (Cont^.RowSums[i] * Cont^.ColumnSums[i]);
      Sum2 := Sum2 + (Cont^.ColumnSums[i] * Cont^.ColumnSums[i]);
      Sum3 := Sum3 + (Cont^.RowSums[i] * Cont^.RowSums[i]);
      Sum4 := Sum4 + Cont^.Table[i, i];
    END;
  Result := (Sum4 * Cont^.n - Sum1) / Sqrt((Cont^.n * Cont^.n - Sum2) *
    (Cont^.n * Cont^.n - Sum3));
END;


PROCEDURE Contingency(Data1, Data2: VectorTyp; VAR Cont: ContTable);

VAR
  i, j, l, m, Length : WORD;
  c                  : char;

BEGIN
  IF NOT (TestDataVectorLength(Data1, Data2)) THEN EXIT;
  Length := VectorLength(Data1);
  TRY
    GetMem(Cont, SizeOf(ContStruc));
  except
    c := WriteErrorMessage(' Contingency table: Not enough memory to create table');
    EXIT;
  END;
  FOR i := 0 TO MaxSteps DO
    // we don't know number of different values yet, so use maximum
    begin
      for j := 0 to MaxSteps do
        Cont^.Table[i, j] := 0;
      Cont^.RowSums[i] := 0;
      Cont^.ColumnSums[i] := 0;
    end;
  Cont^.n := 0;
  Cont^.r := 0;
  Cont^.c := 0;
  for i := 1 to Length do
    begin
      if (IsNaN(GetVectorElement(Data1, i)) or IsNaN(GetVectorElement(Data2, i)))
        then  // ignore
        else
          begin
            l := trunc(GetVectorElement(Data1, i));
            m := trunc(GetVectorElement(Data2, i));
            Inc(Cont^.Table[l, m]);
            if (l > Cont^.r)
              then Cont^.r := l; // set r to highest value of Data1
            if (m > Cont^.c)
              then Cont^.c := m; // set c to highest value of Data2
            Inc(Cont^.n);
          end;
    end;
  for i := 0 to Cont^.r do
    for j := 0 to Cont^.c do
      begin
        Cont^.RowSums[i] := Cont^.RowSums[i] + Cont^.Table[i, j];
        Cont^.ColumnSums[j] := Cont^.ColumnSums[j] + Cont^.Table[i, j];
      end;
end;


procedure DestroyContingency(var Cont: ContTable);

begin
  FreeMem(Cont, SizeOf(ContStruc));
end;


PROCEDURE WriteContingency(CONST Cont : ContTable; MedStr : String; ValidFigures : byte);

VAR i, j    : WORD;
    Medium  : TEXT;
    c       : CHAR;

BEGIN
  IF (MedStr = 'CON')
    THEN
      BEGIN
        FOR i := 0 TO Cont^.r DO
          BEGIN
            FOR j := 0 TO Cont^.c Do
              write(Cont^.Table[i,j]:ValidFigures, ' ');
            writeln('| ', Cont^.RowSums[i]:ValidFigures);
          END;
        FOR j := 0 TO Cont^.c DO
          FOR i := 1 TO succ(ValidFigures) DO write('_');
        write('|');
        FOR i := 1 TO succ(ValidFigures) DO write('_');
        writeln;
        FOR j := 0 TO Cont^.c DO write(Cont^.ColumnSums[j]:ValidFigures, ' ');
        writeln('| ', Cont^.n:ValidFigures);
      END
    ELSE
      BEGIN
        TRY
          assign(Medium, MedStr);
          rewrite(Medium);
        EXCEPT
          c := WriteErrorMessage(' Writing contingency table: could not open file ');
          CorrelationsError := true;
          EXIT;
        END;
        FOR i := 0 TO Cont^.r DO
          BEGIN
            FOR j := 0 TO Cont^.c DO
              write(Medium, Cont^.Table[i,j]:ValidFigures, ' ');
            writeln(Medium, '| ', Cont^.RowSums[i]:ValidFigures);
          END;
        FOR j := 0 TO Cont^.c DO
          FOR i := 1 TO succ(ValidFigures) DO write(Medium, '_');
        write(Medium, '|');
        FOR i := 1 TO succ(ValidFigures) DO write(Medium, '_');
        writeln(Medium);
        FOR j := 0 TO Cont^.c DO write(Medium, Cont^.ColumnSums[j]:ValidFigures, ' ');
        writeln(Medium, '| ', Cont^.n:ValidFigures);
      END;
END;


FUNCTION Chi2(const Contingency: ContTable;
  var Significance: SignificanceType): float;

VAR
  H: array[0..MaxSteps] of float;
  i, j: word;
  A: float;

BEGIN
  FOR i := 1 TO MaxSteps DO
    H[i] := 0.0;
  FOR i := 0 TO Contingency^.r DO
    FOR j := 0 TO Contingency^.c DO
      BEGIN
        A := Contingency^.Table[i, j];
        IF Contingency^.ColumnSums[j] = 0
          THEN A := 0
          ELSE A := A * A / Contingency^.ColumnSums[j];
        H[i] := H[i] + A;
      END;
  A := 0;
  FOR i := 0 TO Contingency^.r DO  // calculate first phi^2 = chi^2 / n
    BEGIN
      IF (Contingency^.RowSums[i] = 0)
        THEN
        ELSE A := A + H[i] / Contingency^.RowSums[i];
    END;
  A := (A - 1) * Contingency^.n;
  Significance.TestValue := A;
  Significance.Freedom := pred(Contingency^.r) * pred(Contingency^.c);
  Significance.P0 := IntegralChi(Significance.TestValue, round(Significance.Freedom));
  Result := A;
END;

FUNCTION Phi2(CONST Contingency: ContTable;
  VAR Significance: SignificanceType): float;

VAR
  x: float;

BEGIN
  x := chi2(Contingency, Significance);
  Result := x / Contingency^.n;
END;

FUNCTION CramersVTilde(chi2: float; r, c, n: word): float;

VAR
  phi2, hilfs, ct, rt: float;

BEGIN
  phi2 := chi2 / n;
  hilfs := phi2 - pred(c) * pred(r) / pred(n);     // now calculate tilde-version
  IF hilfs > 0
    THEN phi2 := hilfs
    ELSE phi2 := 0;
  ct := c - PRED(c) * PRED(c) / PRED(n);
  rt := r - PRED(r) * PRED(r) / PRED(n);
  IF rt > ct
    THEN rt := ct;                        // determine min(r~,c~)
  Result := sqrt(phi2 / (rt - 1));
end;

FUNCTION lambda(CONST Contingency: ContTable): float;

VAR
  i, j: word;
  MaxR, MaxC, SumMaxTi, SumMaxTj: float;
  MaxTi, MaxTj: ARRAY[0..MaxSteps] OF WORD;
  // can't be vectorType AS 0 must be legal index

BEGIN
  FOR i := 0 TO MaxSteps DO
    BEGIN
      MaxTi[i] := 0;
      MaxTj[i] := 0;
    END;
  SumMaxTi := 0;
  FOR i := 0 TO Contingency^.r DO
    BEGIN
      FOR j := 0 TO Contingency^.c DO
        IF MaxTi[i] < Contingency^.Table[i, j] THEN
          MaxTi[i] := Contingency^.Table[i, j];
      SumMaxTi := SumMaxTi + MaxTi[i];
    END;
  SumMaxTj := 0;
  FOR j := 0 TO Contingency^.c DO
    BEGIN
      FOR i := 0 TO Contingency^.r DO
        IF MaxTj[j] < Contingency^.Table[i, j]
          THEN MaxTj[j] := Contingency^.Table[i, j];
      SumMaxTj := SumMaxTj + MaxTj[j];
    END;
  MaxR := Contingency^.RowSums[0];
  FOR i := 1 TO Contingency^.r DO
    IF (Contingency^.RowSums[i] > MaxR)
      THEN MaxR := Contingency^.RowSums[i];
  MaxC := Contingency^.ColumnSums[0];
  FOR j := 1 TO Contingency^.c DO
    IF (Contingency^.ColumnSums[j] > MaxC)
      THEN MaxC := Contingency^.ColumnSums[j];
  Result := (2 * Contingency^.n - (MaxR + MaxC));
  IF Result = 0
    THEN           // ???
    ELSE lambda := (SumMaxTi + SumMaxTj - MaxR - MaxC) / lambda;
END;


FUNCTION lambda_r(CONST Contingency: ContTable): float;

VAR
  i, j: WORD;
  MaxR, MaxT, SumCT: float;

BEGIN
  MaxR := 0;
  FOR i := 0 TO Contingency^.r DO
    IF (Contingency^.RowSums[i] > MaxR)
      THEN MaxR := Contingency^.RowSums[i];
  SumCT := 0;
  FOR j := 0 TO Contingency^.c DO
    BEGIN
      MaxT := Contingency^.Table[0, j];
      FOR i := 1 TO Contingency^.r DO
        IF (MaxT < Contingency^.Table[i, j])
          THEN MaxT := Contingency^.Table[i, j];
      SumCT := SumCT + Contingency^.ColumnSums[j] - MaxT;
    END;
  Result := Contingency^.n - MaxR;
  IF (Result = 0)
    THEN          // ???
    ELSE lambda_r := (Contingency^.n - MaxR - SumCT) / lambda_r;
END;


FUNCTION lambda_c(CONST Contingency: ContTable): float;

VAR
  i, j: WORD;
  MaxC, MaxT, SumRT: float;

BEGIN
  MaxC := 0;
  FOR j := 0 TO Contingency^.c DO
    IF (Contingency^.ColumnSums[j] > MaxC)
      THEN MaxC := Contingency^.ColumnSums[j];
  SumRT := 0;
  FOR i := 0 TO Contingency^.r DO
    BEGIN
      MaxT := Contingency^.Table[i, 0];
      FOR j := 1 TO Contingency^.c DO
        IF (MaxT < Contingency^.Table[i, j])
          THEN MaxT := Contingency^.Table[i, j];
      SumRT := SumRT + Contingency^.RowSums[i] - MaxT;
    END;
  lambda_c := (Contingency^.n - MaxC);
  IF lambda_c = 0
    THEN          // ???
    ELSE lambda_c := (Contingency^.n - MaxC - SumRT) / lambda_c;
END;


FUNCTION OrdinalCorrelations(CONST Data1, Data2: VectorTyp; Formula: CHAR;
  VAR Significance: SignificanceType): float;

VAR
  i, j, k, C, D, X0, Y0, Z, Length: WORD;
  x1, x2, y1, y2, U2, U3, V2, V3: float;
  Cont: ContTable;

BEGIN
  IF NOT (TestDataVectorLength(Data1, Data2)) THEN EXIT;
  Length := VectorLength(Data1);
  Contingency(Data1, Data2, Cont);
  C := 0;
  D := 0;
  X0 := 0;
  Y0 := 0;
  Z := 0;
  FOR i := 1 TO Length DO
    FOR j := Succ(i) TO Length DO
      BEGIN
        x1 := GetVectorElement(Data1, i);
        x2 := GetVectorElement(Data1, j);
        y1 := GetVectorElement(Data2, i);
        y2 := GetVectorElement(Data2, j);
        IF (IsNaN(x1) OR IsNaN(x2) OR IsNaN(y1) OR IsNaN(y2))
          THEN // ignore NaNs
          ELSE IF (Trunc(x1) = Trunc(x2))
                 THEN IF (Trunc(y1) = Trunc(y2))
                        THEN INC(Z)                   // tie IN X AND Y
                        ELSE INC(X0)                  // tie IN X only
                 ELSE IF (Trunc(y1) = Trunc(y2))
                        THEN INC(Y0)                  // tie IN Y only
                        ELSE IF (Trunc(x1) > Trunc(x2))
                              THEN IF (Trunc(y1) > Trunc(y2))
                                  THEN INC(C)
                                  ELSE INC(D)
                              ELSE IF (Trunc(y1) > Trunc(y2))
                                     THEN INC(D)
                                     ELSE INC(C);
      END; // FOR j
  CASE UpCase(Formula) OF
    'T', 'G': BEGIN  // calculate Kendal's tau or Goodman & Kruskal's gamma
                Result := (1.0 * C + D);
                IF (Result = 0)
                  THEN      // ???
                  ELSE Result := (1.0 * C - D) / Result;
              END;
    'D':      BEGIN  // calculate Kim's d_{xy} (note that d_{xy} <> d_{yx})
                Result := C + D + X0;
                if (Result = 0)
                  then      // ???
                  else Result := (1.0 * C - D) / Result;
              end;
    'E':      begin  // calculate Wilson's e
                Result := C + D + X0 + Y0;
                IF (Result = 0)
                  THEN      // ???
                  ELSE Result := (1.0 * C - D) / Result;
              END;
  END; // case
  U2 := 0; // now start WITH significance
  FOR i := 0 TO Cont^.r DO
    FOR j := Succ(i) TO Cont^.r DO
      U2 := U2 + Cont^.RowSums[i] * Cont^.RowSums[j];
  // sum OF product OF row sums 2 at a time
  V2 := 0;
  FOR i := 0 TO Cont^.c DO
    FOR j := Succ(i) TO Cont^.c DO
      U2 := U2 + Cont^.ColumnSums[i] * Cont^.ColumnSums[j];
  // sum OF product OF column sums 2 at a time
  U3 := 0;
  FOR i := 0 TO Cont^.c DO
    FOR j := Succ(i) TO Cont^.c DO
      FOR k := Succ(j) TO Cont^.c DO
        U3 := U3 + Cont^.RowSums[i] * Cont^.RowSums[j] * Cont^.RowSums[k];
  // sum OF product OF row sums 3 at a time
  V3 := 0;
  FOR i := 0 TO Cont^.c DO
    FOR j := Succ(i) TO Cont^.c DO
      FOR k := Succ(j) TO Cont^.c DO
        U3 := U3 + Cont^.ColumnSums[i] * Cont^.ColumnSums[j] * Cont^.ColumnSums[k];
  // sum OF product OF column sums 3 at a time
  Significance.Freedom := Sqrt(U2 * V2 / Pred(Cont^.n) - (U2 * V3 + V2 * U3) /
    (Cont^.n * Pred(Cont^.n)) + U3 * V3 / (Cont^.n * Pred(Cont^.n) * (Cont^.n - 2)));
  // standard deviation
  Significance.P0 := 0;
  // IF it can't be calculated
  Significance.TestValue := (2 * pred(Cont^.r) * pred(Cont^.c));
  if (Significance.TestValue <> 0) then
    Significance.TestValue := abs(C - D) - Cont^.n / Significance.TestValue; // average
  if (Significance.Freedom <> 0) then
    Significance.P0 := IntegralGauss(Significance.TestValue / Significance.Freedom);
  // error corresponding to z
end;


FUNCTION theta(CONST Contingency: ContTable): float;

VAR
  i, j, k, l: word;
  SumT, fa, fb, SumDi, nr: float;

BEGIN
  SumT := 0;
  FOR i := 0 TO Contingency^.r DO
    FOR j := succ(i) TO Contingency^.r DO
      SumT := SumT + Contingency^.RowSums[i] * Contingency^.RowSums[j];
  SumDi := 0;
  FOR i := 0 TO Contingency^.r DO
    BEGIN
      FOR k := succ(i) TO Contingency^.r DO
        BEGIN
          fa := 0;
          fb := 0;
          FOR j := 0 TO Contingency^.c DO
            BEGIN
              nr := Contingency^.Table[i, j];
              IF (j < Contingency^.c)
                THEN
                  FOR l := succ(j) TO Contingency^.c DO
                    fa := fa + Contingency^.Table[k, l] * nr;
              IF (j > 0)
                THEN // this check actually is necessary
                  FOR l := pred(j) DOWNTO 0 DO
                    fb := fb + Contingency^.Table[k, l] * nr;
            END;
          SumDi := SumDi + abs(fb - fa);
        END;
    END;
  IF (SumT = 0)
    THEN Result := 0    // ???
    ELSE Result := SumDi / SumT;
END;


FUNCTION eta_sqr(CONST NominalVector, CardinalVector: VectorTyp;
  VAR Significance: SignificanceType): float;

VAR
  i, k, n, Classes, Length, Value1: WORD;
  E1, E2, total, r: float;
  SumOfK: array [0..MaxSteps] of float;
  NofK: ARRAY [0..MaxSteps] OF WORD;

BEGIN
  IF NOT (TestDataVectorLength(NominalVector, CardinalVector)) THEN EXIT;
  Length := VectorLength(NominalVector);
  Classes := 0;
  FOR i := 1 TO Length DO
    BEGIN
      IF IsNaN(GetVectorElement(NominalVector, i))
        THEN
        ELSE IF (round(GetVectorElement(NominalVector, i)) > Classes)
               THEN Classes := ROUND(GetVectorElement(NominalVector, i));
      // find largest value for nominal vector
    END;
  FOR i := 0 TO Classes DO
    BEGIN
      NofK[i] := 0;
      SumOfK[i] := 0;
    END;
  Total := 0;
  n := 0;
  FOR i := 1 TO Length DO
    IF ((IsNaN(GetVectorElement(NominalVector, i))) OR
      (IsNaN(GetVectorElement(CardinalVector, i))))
      THEN   // ignore NaNs
      ELSE
        BEGIN
          Value1 := trunc(GetVectorElement(NominalVector, i));
          Inc(NofK[Value1]);
          // how many x are there for any of the 0..Classes values of x?
          SumOfK[Value1] := SumOfK[Value1] + GetVectorElement(CardinalVector, i);
        END;
  FOR k := 0 TO Classes DO
    BEGIN
      n := n + NofK[k];                         // number of data pairs that are not NaN
      Total := Total + SumOfK[k];
      IF (NofK[k] <> 0)
        THEN SumOfK[k] := SumOfK[k] / NofK[k]; // group-averages
    END;
  Total := Total / n;                          // average over all data
  E1 := 0;
  E2 := 0;
  FOR i := 1 TO Length DO
    BEGIN
      IF ((IsNaN(GetVectorElement(NominalVector, i))) OR
        (IsNaN(GetVectorElement(CardinalVector, i))))
        THEN
        ELSE E1 := E1 + sqr(GetVectorElement(CardinalVector, i) - Total);
    END;
  FOR i := 1 TO Length DO
    BEGIN
      IF ((IsNaN(GetVectorElement(NominalVector, i))) OR
          (IsNaN(GetVectorElement(CardinalVector, i))))
        THEN
        ELSE E2 := E2 + sqr(GetVectorElement(CardinalVector, i) -
                   SumOfK[trunc(GetVectorElement(NominalVector, i))]);
    END;
  IF (E1 = 0)
    THEN r := 0.0    // ???
    ELSE r := 1 - E2 / E1;
  IF (r = 1)
    THEN
      BEGIN
        Significance.P0 := NaN;
        Significance.Freedom := n - r;
        Significance.TestValue := NaN;
      END
    ELSE
      BEGIN
        Significance.TestValue := r / (1 - r) * (n - r) / (r - 1);
        // calculate F
        Significance.Freedom := n - r;
        // other is r-1 and can be calculated in calling program
        Significance.P0 := Integral_F(Significance.TestValue, round(n - r), round(r - 1));
        // calculate probability for 0-hypotheses
        Result := r;
      END;
END;


FUNCTION AverageCorrelations(CONST Correlations: VectorTyp): float;

VAR
  i, j, k: WORD;
  SumR: float;

BEGIN
  k := VectorLength(Correlations);
  SumR := 0;
  j := 0;
  FOR i := 1 TO k DO
    IF IsNaN(GetVectorElement(Correlations, i))
      THEN // ignore
      ELSE
        BEGIN
          SumR := SumR + tanh(GetVectorElement(Correlations, i));
          Inc(j);
        END;
  Result := arctanh(SumR / j);
END;

FUNCTION Angle (r : float) : float;

BEGIN
   Result := ArcCos(r) * 180 / Const_pi;
END;


end.  // correlations

