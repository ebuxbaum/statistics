PROGRAM Similarity;

USES math, MathFunc, Vector, Matrix, BigSet, Correlations;

CONST MaxVariables =  180;
      MaxCases     = 1400;
      SepChar      =  ';';    // separates data IN csv-fles

TYPE DataTypes  = (binary, nominal, ordinal, interval, rational);
     TypeVector = ARRAY[1..MaxVariables] OF DataTypes;

VAR Types      : TypeVector;
    Data, Dist1, Dist2 : MatrixTyp;
    Variables,
    Cases,i,j  : WORD;
    r          : double;

TYPE ClusterFeld = ARRAY [1..MaxCases] OF SetType;     // declarations should be local TO ClusterAnalysis

VAR HilfsMatrix  : MatrixTyp;    { akt. Distanzen der Cluster oben, coph. Matrix unten und Clusterzustand diagonal }
    Cluster      : ClusterFeld;
    v            : VectorTyp;

PROCEDURE ReadCSV(VAR Cases, Variables : WORD; VAR Types: TypeVector;
                  VAR Data: MatrixTyp);

CONST MaxCases = 200;

VAR
  InputFile: TEXT;
  c: CHAR;
  s: STRING;
  Error, i, j: WORD;
  Interim: MatrixTyp;
  x: double;

BEGIN
  Variables := 0;
  Cases := 0;
  Assign(InputFile, 'All.csv');
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
    CASE s OF
      'binary'   : Types[Variables] := binary;
      'ordinal'  : Types[Variables] := ordinal;
      'nominal'  : Types[Variables] := nominal;
      'interval' : Types[Variables] := interval;
      'rational' : Types[Variables] := rational
      ELSE
        BEGIN              // Format error: abort PROGRAM WITH error message
          Write('unknown data type ', s, ' at n=', Variables: 3, '. Press <CR>:');
          ReadLn;
          HALT;
       END; { else }
    END; { case }
  END; { while }
  ReadLn(InputFile);            // go TO next Line
  CreateMatrix(Interim, MaxCases, Variables, 0.0);  // Note: number OF persons NOT yet known
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
  CreateMatrix(Data, Cases, Variables, 0.0);  // now generate correctly sized data matrix
  FOR i := 1 TO Cases DO
    FOR j := 1 TO Variables DO
      SetMatrixElement(Data, i, j, GetMatrixElement(Interim, i, j));
  DestroyMatrix(Interim);                         // destroy the interim matrix
  Writeln(Cases: 4, ' cases with ', Variables: 3, ' variables read, ');
END;

PROCEDURE Gowers(CONST Types : TypeVector; CONST Data : MatrixTyp;
                  VAR Dist : MatrixTyp);
// calculate Gower's modified general similarity coefficients (see Sneath & Sokal 1973, section 4.4)
// for ordinal data see: Padani, Taxon 48 (1999) 331-40

var SumW, i, j, k, Cases, Variables : word;
    S, SumWS, x, y          : double;
    Maximum, Minimum, Range : array [2..MaxVariables] of double; // as first column is case number

begin
  Variables := MatrixColumns(Data);
  Cases     := MatrixRows(Data);
  for k := 2 to Variables do  // calculate range of intervall and rational data, ignore StudyNamber
    begin
      Maximum[k] := -MaxRealNumber;
      Minimum[k] := MaxRealNumber;
      case Types[k] of
        ordinal, interval, rational : begin
              for i := 1 to Cases do
                  begin
                    x := GetMatrixElement(Data, i, k);
                    if not isNaN(x)
                      then
                        begin
                          if (x > Maximum[k]) then Maximum[k] := x;
                          if (x < Minimum[k]) then Minimum[k] := x;
                        end;
                  end; { for i }
                Range[k] := Maximum[k] - Minimum[k];
              end
        else Range[k] := 0.0;  // binary, nominal: just give it a defined value, won't be used
      END;  { case }
    END; { for k }
  Writeln('ranges calculated');
  Write('Distances: ');
  CreateIdentityMatrix(Dist, Cases);
  FOR i := 1 TO Cases DO
    BEGIN
      FOR j := Succ(i) TO Cases DO // calculate upper half OF similarity matrix
        BEGIN
           SumW := 0;
           SumWS := 0;
           FOR k := 2 TO Variables DO  // ignore CASE numbes
             BEGIN
               x := GetMatrixElement(Data, i, k);
               y := GetMatrixElement(Data, j, k);
//               Write(OutFile, i:4, ';', j:4, ';', k:4, ';', FloatStr(x, 8), ';', FloatStr(y, 8), ';');
               IF (isNaN(x) OR isNaN(y))
                 THEN // SumW AND SumWS both increase by 0
                 ELSE
                   CASE Types[k] OF
                     binary, nominal : BEGIN
                                         IF x = y
                                           THEN S := 0
                                           ELSE S := 1;
                                         INC(SumW);
                                         SumWS := SumWS + S; // AS W = 1 -> WS = S
                                       END;
                     ordinal, interval, rational :
                                       BEGIN
                                         S := Abs(x - y) / Range[k];
                                         INC(SumW);
                                         SumWS := SumWS + S;
                                       END;
                   END; { case }
             END; { for k }
           IF SumW = 0
             THEN x := NaN                     // so that such cases can be identified IN the distance matrix
             ELSE x := SumWS / SumW;
           SetMatrixElement(Dist, i, j, x);
           SetMatrixElement(Dist, j, i, SumW); // put number OF valid compares into lower half OF similarity matrix
        END;  { for j }
      IF (i MOD 20 = 0) THEN Write('.');
    END; { for i }
  Writeln(' calculated, ');
END;  { Gowers }

PROCEDURE CalcCaseCorrelations (CONST Data : MatrixTyp; VAR Corr : MatrixTyp);
{ Calculate row-wise correlation as similarity measure. Requires cardinal data. }

VAR i, j, k, n, Cases, Variables : WORD;
    Vec1, Vec2                : VectorTyp;
    Significance              : SignificanceType;
    r                         : float;

BEGIN
  Variables := MatrixColumns(Data);
  Cases := MatrixRows(Data);
  CreateVector(Vec1, Pred(Variables), 0.0);
  CreateVector(Vec2, Pred(Variables), 0.0);
  CreateIdentityMatrix(Corr, Cases);
  FOR i := 1 TO Cases DO
    BEGIN
      FOR k := 2 TO Variables DO  // ignore the study # IN column 1
        SetVectorElement(Vec1, Pred(k), GetMatrixElement(Data, i, k));
      FOR j := Succ(i) TO cases DO
        BEGIN
          n := 0;
          FOR k := 2 TO Variables DO
            BEGIN
              SetVectorElement(Vec2, Pred(k), GetMatrixElement(Data, j, k));
              IF IsNaN(GetVectorElement(Vec1, Pred(k))) OR IsNaN(GetVectorElement(Vec1, Pred(k)))
                THEN
                ELSE INC(n);
            END;
          r := PearsonProductMomentCorrelation(Vec1, Vec2, Significance); // Q-mode
          r := 0.5 * (r + 1);  // convert range from -1..1 TO 0..1
          r := 1 - r;          // similarity TO distance
          SetMatrixElement(Corr, i, j, r);
          SetMatrixElement(Corr, j, i, n);
        END;
      IF (i MOD 20 = 0) THEN Write('.');
    END; { for i }
  Writeln(' calculated, ');
  DestroyVector(Vec1);
  DestroyVector(Vec2);
END;

PROCEDURE WriteDistances (CONST Data, Dist : MatrixTyp);

VAR i, j, cases  : WORD;
    OutFile : TEXT;
    MaxS, MinS, MaxC, MinC, x : double;

BEGIN
  Cases := MatrixRows(Dist);
  Assign(OutFile, 'Distances.csv');
  Rewrite(OutFile);
  MaxS := MinRealNumber;
  MaxC := MinRealNumber;
  MinS := MaxRealNumber;
  MinC := MaxRealNumber;
  Write(OutFile, 'StudyNum;');  // LABEL columns
  FOR i := 1 TO Cases DO Write(OutFile, Round(GetMatrixElement(Data, i, 1)):4, ';');
  Writeln(OutFile);
  FOR i := 1 TO Cases DO // the matrix itself
    BEGIN
      Write(OutFile, Round(GetMatrixElement(Data, i, 1)):4, ';');
      FOR j := 1 TO Cases DO
        IF (j>i)
          THEN
            BEGIN
              x := GetMatrixElement(Dist, i, j);
              Write(OutFile, x:6:4, ';');   // upper half WITH similarity
              IF IsNaN(x)
                THEN
                ELSE
                  BEGIN
                    IF (x > MaxS) THEN MaxS := x;
                    IF (x < MinS) THEN MinS := x;
                  END;
            END
          ELSE
             BEGIN
              x := GetMatrixElement(Dist, i, j);
              Write(OutFile, Round(x):4, ';');  // lower half WITH # OF comparisons
              IF (j < i) // ignore i=j
                THEN
                  BEGIN
                    IF (x > MaxC) THEN MaxC := x;
                    IF (x < MinC) THEN MinC := x;
                  END;
            END;
       Writeln(OutFile)
    END;  { for i }
  Close(OutFile);
  Writeln;
  Writeln('Distances written to file, ');
  Writeln('Range ', MinS:5:3, '-', MaxS:5:3, ' with ', MinC:3:0, '-', MaxC:3:0, ' comparisons');
END; { WriteDistances }


PROCEDURE AnalyseFrequencies(CONST Dist: MatrixTyp);

CONST Border = 50;    // number OF ranges FOR statistical analysis OF correlations

TYPE FreqsType  = ARRAY[-Border..Border] OF WORD;

VAR
  i, j, Cases     : WORD;
  k               : INTEGER;
  OutFile         : TEXT;
  x               : double;
  Freqs           : FreqsType;

BEGIN
  Assign(OutFile, 'DistFreq.csv');
  Rewrite(OutFile);
  Cases := MatrixRows(Dist);
  FOR k := -Border TO Border DO
    Freqs[k] := 0;
  FOR i := 1 TO Cases DO
    FOR j := Succ(i) TO Cases DO
    BEGIN
      x := GetMatrixElement(Dist, i, j) * Border;
      INC(Freqs[Round(x)]);
    END;
  FOR k := -Border TO Border DO
     Writeln(OutFile, k/Border:1:4, SepChar, Freqs[k]:6, SepChar);
  Writeln(OutFile);
  Close(OutFile);
END;


FUNCTION ClusterAnalysis (CONST Dist : MatrixTyp) : double;

VAR i, j,
    AnzCluster   : WORD;
    MinAbst      : double;
    ClusterDatei : TEXT;

     PROCEDURE Minimum (CONST HilfsMatrix : MatrixTyp; VAR MinAbst : double);

     VAR i, j     : WORD;

     BEGIN
       FOR i := 1 TO Pred(Cases) DO
         FOR j := Succ(i) TO Cases DO
           IF GetMatrixElement(HilfsMatrix, i, j) < MinAbst
             THEN MinAbst := GetMatrixElement(HilfsMatrix, i, j);
   END;


     PROCEDURE UniteCluster (VAR Cluster : ClusterFeld; i, j : WORD; MinAbst : double);

     VAR k : WORD;

     BEGIN
       SetUnion(Cluster[i], Cluster[i], Cluster[j]);
       ClearAllBits(Cluster[j]);
       Writeln(ClusterDatei);
       Writeln;
       Write(ClusterDatei, MinAbst:6:4, ' ');
       Write(MinAbst:6:4, ' ');
       FOR k := 1 TO Cases DO
         IF InSet(Cluster[i], k)
           THEN
             BEGIN
               Write(ClusterDatei, Round(GetMatrixElement(Data, k, 1)):4, ' ');
               Write(Round(GetMatrixElement(Data, k, 1)):4, ' ');
             END;
       Writeln(ClusterDatei);
       Writeln;
     END;


     PROCEDURE Cophen (VAR HilfsMatrix : MatrixTyp; CONST Cluster : ClusterFeld;
                       MinAbst : double);

     VAR i, j, l : WORD;

     BEGIN
       FOR i := 1 TO Cases DO
         IF NOT EmptySet(Cluster[i])
           THEN
             FOR j := 1 TO Pred(Cases) DO
                IF InSet(Cluster[i], j)
                  THEN
                    FOR l := Succ(i) TO Cases DO
                       IF InSet(Cluster[i], l)
                         THEN
                           IF (GetMatrixElement(HilfsMatrix, l, j) = 0)
                             THEN SetMatrixElement(HilfsMatrix, l, j, MinAbst);
     END; { Cophen }


     PROCEDURE NewDistances (CONST Cluster : ClusterFeld; CONST Dist : MatrixTyp;
                                VAR HilfsMatrix : MatrixTyp);

     VAR i, j, k, l, b : WORD;
         a             : double;

     BEGIN
       FOR i := 1 TO Pred(Cases) DO
         FOR j := Succ(i) TO Cases DO
           BEGIN
             IF ((NOT EmptySet(Cluster[i])) AND (NOT EmptySet(Cluster[j])))
               THEN
                 BEGIN
                   a := 0.0;
                   b := 0;
                   FOR k := 1 TO Cases DO
                     IF InSet(Cluster[i], k)
                       THEN
                         FOR l := 1 TO Cases DO
                           IF InSet(Cluster[j], l)
                             THEN
                               BEGIN
                                 a := a + GetMatrixElement(Dist, k, l);
                                 INC(b);
                               END;
                   SetMatrixElement(HilfsMatrix, i, j, a/b);
                 END
               ELSE
                 SetMatrixElement(HilfsMatrix, i, j, MaxInt);
           END; { for }
     END; {NeueSimularitaeten}


     FUNCTION NewDiagonal (CONST Cluster: ClusterFeld; VAR HilfsMatrix : MatrixTyp) : WORD;

     VAR i, j : WORD;

     BEGIN
       j := 0;
       FOR i := 1 TO Cases DO
         IF NOT EmptySet(Cluster[i])
           THEN
             BEGIN
               SetMatrixElement(HilfsMatrix, i, i, 1);
               INC(j);
             END
           ELSE
             SetMatrixElement(HilfsMatrix, i, i, 0);
       NewDiagonal := j;
     END;


     FUNCTION Correlation (CONST Dist, HilfsMatrix : MatrixTyp) : double;

     VAR DistMittel, CoMittel,
         DistKor, CoKor,
         SumDistKorSqr, SumCoKorSqr,
         SumDistKorCoKor             : double;


         PROCEDURE Mittelwerte (CONST Dist, HilfsMatrix : MatrixTyp);

         VAR SumDist, SumCo : double;
             i, j, z       : WORD;

         BEGIN
           SumDist := 0;
           SumCo := 0;
           z := 0;
           FOR i := 1 TO Pred(Cases) DO
             FOR j := Succ(i) TO Cases DO
               BEGIN
                 SumDist := SumDist + GetMatrixElement(Dist, i, j);
                 SumCo := SumCo + GetMatrixElement(HilfsMatrix, j, i);
                 INC(z);
               END;
           DistMittel := SumDist / z;
           CoMittel := SumCo / z;
         END;


         FUNCTION Summen (CONST Dist, HilfsMatrix : MatrixTyp) : double;

         VAR i, j : WORD;

         BEGIN
           SumDistKorSqr := 0;
           SumCoKorSqr := 0;
           SumDistKorCoKor := 0;
           FOR i := 1 TO Pred(Cases) DO
             FOR j := Succ(i) TO Cases DO
               BEGIN
                 DistKor := GetMatrixElement(Dist, i, j) - DistMittel;
                 CoKor := GetMatrixElement(HilfsMatrix, j, i) - CoMittel;
                 SumDistKorSqr := SumDistKorSqr + Sqr(DistKor);
                 SumCoKorSqr := SumCoKorSqr + Sqr(CoKor);
                 SumDistKorCoKor := SumDistKorCoKor + DistKor * CoKor;
               END;
           Writeln(SumDistKorCoKor:4:1, ' ', SumDistKorSqr:4:1, ' ',  SumCoKorSqr:4:1);
           Summen := SumDistKorCoKor / (Sqrt(SumDistKorSqr) * Sqrt(SumCoKorSqr));
         END;

     BEGIN
       Mittelwerte(Dist, HilfsMatrix);
       Correlation := Summen(Dist, HilfsMatrix);
     END;

BEGIN
  ClusterAnalysis := 0.0;
  Assign(ClusterDatei, 'Similarities.CLU');
  Rewrite(ClusterDatei);
  CopyMatrix(Dist, HilfsMatrix);
  FOR i := 1 TO Cases DO         // jedem OTU sein eigenes Cluster
    BEGIN
      ClearAllBits(Cluster[i]);
      SetBit(Cluster[i], i);
    END;
  REPEAT                                        { Beginn der Analysenschleife }
    MinAbst := MaxRealNumber;
    Minimum(HilfsMatrix, MinAbst);
    FOR i := 1 TO Pred(Cases) DO               { neue Cluster bilden }
      FOR j := Succ(i) TO Cases DO
        IF (GetMatrixElement(HilfsMatrix, i, j) = MinAbst)
          THEN UniteCluster(Cluster, i, j, MinAbst);
    Cophen(HilfsMatrix, Cluster, MinAbst);
    NewDistances(Cluster, Dist, HilfsMatrix);
    AnzCluster := NewDiagonal(Cluster, HilfsMatrix);
  UNTIL AnzCluster = 1;
  Close(ClusterDatei);
  Result := Correlation(Dist, Hilfsmatrix);
END;

FUNCTION MatrixCorrelation (CONST Dist1, Dist2 : MatrixTyp) : float;

VAR i, j, k : WORD;
  SumXY, SumX, SumY, SumX2, SumY2, varX, varY, covXY, x, y, r : float;
  OutFile : TEXT;

BEGIN
  Assign(OutFile, 'DistCorr.csv');
  Rewrite(OutFile);
  Cases := MatrixRows(Dist1);
  SumXY := 0;
  SumX := 0;
  SumY := 0;
  SumX2 := 0;
  SumY2 := 0;
  k := 0;
  FOR i := 1 TO Cases DO
    FOR j := Succ(i) TO Cases DO
      BEGIN
        x := GetMatrixElement(Dist1, i, j);
        y := GetMatrixElement(Dist2, i, j);
        Writeln(OutFile, i:3, ';', j:3, ';', FloatStr(x, 10), ';', FloatStr(y, 10), ';');
        SumXY := SumXY + x * y;
        SumX := SumX + x;
        SumY := SumY + y;
        SumX2 := SumX2 + x * x;
        SumY2 := SumY2 + y * y;
        INC(k); // actual number OF valid comparisons
       END; { for }
  varX := k * SumX2 - SumX * SumX;
  varY := k * SumY2 - SumY * SumY;
  covXY := k * SumXY - SumX * SumY;
  IF Abs(varX * varY) < Zero
    THEN Result := 0     // ???
    ELSE Result := covXY / Sqrt(varX * varY);
  Writeln(OutFile, FloatStr(Result, 10));
  Close(OutFile);
END; { MatrixCorrelation }

BEGIN
  ReadCSV (Cases, Variables, Types, Data);
  Gowers(Types, Data, Dist1);
  CalcCaseCorrelations(Data, Dist2);
  WriteDistances(Data, Dist2); // Data IS source OF study-number
  AnalyseFrequencies(Dist2);
  FOR i := 1 TO Cases DO
     FOR j := Succ(i) TO Cases DO
       BEGIN
         SetMatrixElement(Dist1, j, i, GetMatrixElement(Dist1, i, j));       // symmetrieren
         SetMatrixElement(Dist2, j, i, GetMatrixElement(Dist2, i, j));
       END;
  LeadingPrincipleMinors(Dist2, V);
  FOR i := 1 TO Cases DO
    Writeln(FloatStr(GetVectorElement(V, i), 10));
  DestroyVector(V);
  r := MatrixCorrelation (Dist1, Dist2);
  r := ClusterAnalysis(Dist2);
  Write('Cophenetic correlation', r:1:3, ' Press <CR> to finish:');
  ReadLn;
END.

