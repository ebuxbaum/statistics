unit Nonparam;
{ Non-Parametrische Tests. Fuer diese Tests ist es nicht erforderlich, daß
  die Daten normalverteilt sind. Sie koennen daher angewendet werden, wo
  t-Tests nicht zulaessig sind }

interface

uses windows, mathfunc, Vector, stat;

const NonParamError : boolean = false;

procedure WilcoxonSignedRank (x, y : VectorTyp);
{ nicht-parametrischer Test bei gepaarten Daten }

procedure MannWitneyU (x, y : VectorTyp);
{ nicht-parametrischer Test bei nicht gepaarten Daten }

procedure VorzeichenTest (x, y : VectorTyp);
{ testet die Haeufigkeit, daá x[i] - y[i] > 0 gegen die Binominal-Verteilung }

implementation


procedure FindRanks (x, Rank : VectorTyp);

var i, j, n, NrGleich, NrKleiner : word;
    xWert                        : float;

begin
  n := VectorLength(x);
  for i := 1 to n do
    begin
      NrKleiner := 0;
      NrGleich := 0;
      xWert := GetVectorElement(x, i);
      for j := 1 to n do
        if GetVectorElement(x, j) = xWert
          then
            inc(NrGleich)
          else
            if GetVectorElement(x, j) < xWert
              then
                inc(NrKleiner);
      SetVectorElement(Rank, i, NrKleiner + succ(NrGleich) / 2);
    end;
end;



procedure WilcoxonSignedRank (x, y : VectorTyp);

var i, n, MinusCount, ZeroCount,
    NonZeroCount, PlusCount     : word;
    AbsDiff, Diff, Rank         : VectorTyp;
    Difference, TMean, TVar,
    TMinus, TPlus,
    NormalDerivative            : float;
    ch				: char;

begin
  if VectorLength(x) <> VectorLength(y)    { keine gepaarten Daten }
    then
      begin
        ch := WriteErrorMessage('Wilcoxon-Test: Vectors have different length! ');
        NonparamError := true;
        exit;
      end;
  n := VectorLength(x);
  ZeroCount := 0;
  MinusCount := 0;
  CreateVector(AbsDiff, n, 0.0);
  CreateVector(Diff, n, 0.0);
  CreateVector(Rank, n,0.0);
  for i := 1 to n do
    begin
      Difference := GetVectorElement(x, i) - GetVectorElement(y, i);
      if Difference = 0.0
        then
          begin
            inc(ZeroCount);
            SetVectorElement(AbsDiff, i, MaxRealNumber);
	  end
	else
          SetVectorElement(AbsDiff, i, abs(Difference));
      SetVectorElement(Diff, i, Difference);
    end;
  FindRanks(AbsDiff, Rank);
  TMinus := 0.0;
  for i := 1 to n do
    if GetVectorElement(Diff, i) < 0
      then
        begin
          inc(MinusCount);
          TMinus := TMinus + GetVectorElement(Rank, i);
        end
      else
        if GetVectorElement(Diff, i) > 0
          then
            begin
              inc(PlusCount);
              TPlus := TPlus + GetVectorElement(Rank, i);
            end;
  NonZeroCount := n - ZeroCount;
  TMean := NonZeroCount * succ(NonZeroCount) / 4;
  TVar := TMean * succ(2 * (NonZeroCount)) / 6;
  NormalDerivative := (TMinus - TMean) / sqrt(TVar);
  DestroyVector(AbsDiff);
  DestroyVector(Diff);
  DestroyVector(Rank);
  writeln('Difference         Number       Sum of ranks ');
  writeln('1 < 2              ', MinusCount:5, '             ', TMinus:ValidFigures);
  writeln('1 > 2              ', PlusCount:5,  '             ', TPlus:ValidFigures);
  writeln('1 = 2              ', ZeroCount:5,  '             ');
  writeln('                   _____');
  writeln('                   ', n:5);
  writeln;
  writeln('Standard normal derivative: ', NormalDerivative:ValidFigures);
end;


procedure MannWitneyU (x, y : VectorTyp);

var i, n1, n2, nTotal  : word;
    Rank, AlleDaten    : VectorTyp;
    UMean, UVariance,
    RankSum1, Ranksum2,
    U1, U2, Normal     : float;

begin
    n1 := VectorLength(x);
    n2 := VectorLength(y);
    nTotal := n1 + n2;
    CreateVector(AlleDaten, nTotal, 0.0);
    CreateVector(Rank, nTotal, 0.0);
    for i := 1 to n1 do
      SetVectorElement(Alledaten, i, GetVectorElement(x, i));
    for i := 1 to n2 do
      SetVectorElement(Alledaten, i+n1, GetVectorElement(y, i));
    FindRanks(AlleDaten, Rank);
    RankSum1 := 0.0;
    for i := 1 to n1 do
      RankSum1 := RankSum1 + GetVectorElement(Rank, i);
    RankSum2 := nTotal * succ(nTotal)/2 - RankSum1;
    U1 := RankSum1 - n1*succ(n1) / 2;
    U2 := n1 * n2 - U1;
    UMean := n1 * n2 / 2;
    UVariance := n1 * n2 * succ(n1) / 12;
    Normal := (U1 - UMean) / sqrt(UVariance);
    writeln('               Vector 1     Vector 2', Ranksum1, '     ', RankSum2:ValidFigures);
    writeln('Sum of ranks: ', n1:5, '        ', n2:5);
    writeln('Sample size:  ', U1:ValidFigures, '     ', U2:ValidFigures);
    writeln;
    writeln('U-statistics: ', Normal:ValidFigures);
end;


procedure VorzeichenTest (x, y : VectorTyp);

var i, nTotal, nPositiv, nZero : word;
    p0, Diff                   : float;
    ch			       : char;

begin
    if VectorLength(x) <> VectorLength(y)    { keine gepaarten Daten }
     then
       begin
         ch := WriteErrorMessage('Sign-Test: Data vectors have unequal lengths!');
         NonparamError := true;
         exit;
       end;
    nTotal := VectorLength(x);
    nPositiv := 0;
    nZero := 0;
    for i := 1 to nTotal do
      begin
        Diff := GetVectorElement(x, i) - GetVectorElement(y, i);
        if Diff = 0.0
          then
            inc(nZero)
          else
            if Diff > 0
              then
                inc(nPositiv);
      end;
    P0 := BinominalIntegral(nTotal-nZero, nPositiv, 0.5);
    writeln('Total number of data pairs:         ', nTotal:5);
    writeln('Number of positiv differences:      ', nPositiv:5);
    writeln('Number of negativ differences:      ', (nTotal-nPositiv-nZero):5);
    writeln('Number of differences equal to 0:   ', nZero:5);
    writeln('Probability for equal means:        ', P0:ValidFigures, '(', (P0*100):6:2, '%)');
end;


end.

