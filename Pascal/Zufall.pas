UNIT Zufall;

{ Berechnung von Zufallszahlen unterschiedlicher Verteilung
  Literatur: W.H.Press et al., Numerical Recepies in Pascal, Cambridge 1989 }

INTERFACE

USES MathFunc, dos;

FUNCTION RandomLaplace(LowerLimit, UpperLimit: float): float;
{ gleichverteilte ZufallsZahlen aus [LowerLimit..UpperLimit] }

FUNCTION RandomLaplace(LowerLimit, UpperLimit: LONGINT): LONGINT;

FUNCTION RandomBinominal(p: double; n: WORD): WORD;
{ Number der Treffer bei n Versuchen mit Trefferwarscheinlichkeit p }

FUNCTION RandomNormal(Average, SigmaSqr: double): double;
{ Normalverteilte Zufallszahlen.  }

FUNCTION RandomExponential(Average: double): double;
{ Exponential-verteilte Zufallszahlen. Simuliert z. B. die Zeit zwischen
  Ereignissen, die mit constanter Rate, aber zufaellig verteilt eintreten
  (in dieser Situation waehre die Number der Ereignisse pro Zeiteinheit
  Poisson-verteilt!) }

FUNCTION RandomPoisson(Average: double): WORD;
{ Poisson-verteilte Zufallszahlen. Simuliert das Eintreten seltener Ereignisse,
  pro Zeiteinheit (der Abstand zwischen den Ereignissen ist dann Exponential-
  verteil!) }

FUNCTION RandomGamma(Form: WORD; Average: double): double;
{ Verallgemeinerung der Exponential-Verteilten Zufallszahlen. Erlaubt die
  Simulation des Eintreffens von Ereignissen, deren Modus > 0 ist, z. B.
  die Lebensdauer von Produkten bis zum Versagen }

FUNCTION RandomBernoulli(P: double): BYTE;
{ Bernoulli-verteilte Zufallszahlen, immer 1 oder 0. Simuliert binaere Ereignisse
  wie z. B. Muenzwurf. }

FUNCTION RandomGeometric(P: double): WORD;
{ Number der Versuche, die erforderlich sind, bis eine Bernoulli-Verteilung 1
  wird }

FUNCTION RandomPareto(Minimum, Shape: double): double;
{ Schadenshoehe, die mit (Zufall * 100) %iger Warscheinlichkeit nicht
  ueberschritten wird.
  Dabei ist Minimum der Minimalschaden pro Claim (>= 0) und Shape
  der Exponent der Pareto-Verteilung (> 1) }

FUNCTION RandomChiSqr(DegFreedom: WORD): double;
{ Chi-Quadrat verteilte Zufallszahlen }

FUNCTION RandomT(DegFreedom: WORD): double;
{ t-Verteilte Zufallszahlen }

FUNCTION RandomF(f1, f2: WORD): double;
{ F-Verteilte Zufallszahlen }

IMPLEMENTATION

FUNCTION RandomLaplace(LowerLimit, UpperLimit: float): float;

VAR
  dummy: float;

BEGIN
  IF LowerLimit > UpperLimit
    THEN
      BEGIN
        Dummy := LowerLimit;
        LowerLimit := UpperLimit;
        UpperLimit := Dummy;
      END;
  Result := LowerLimit + (UpperLimit - LowerLimit) * Random;
END;


FUNCTION RandomLaplace (LowerLimit, UpperLimit : LONGINT) : LONGINT;

VAR dummy : LONGINT;

BEGIN
  IF LowerLimit > UpperLimit
    THEN
      BEGIN
        Dummy := LowerLimit;
        LowerLimit := UpperLimit;
        UpperLimit := Dummy;
      END;
  Result := LowerLimit + Trunc((UpperLimit-Pred(LowerLimit)) * Random);
END;


FUNCTION RandomBinominal(p: double; n: WORD): WORD;

VAR
  i, Sum: WORD;

BEGIN
  Sum := 0;
  FOR i := 1 TO n DO
    IF Random <= p THEN INC(Sum);
  Result := Sum;
END;


FUNCTION RandomNormal(Average, SigmaSqr: double): double;

BEGIN
  Result := Sqrt(-2 * Ln(Random)) * Sin(2 * Pi * Random) *
    SigmaSqr + Average;
END;


FUNCTION RandomExponential(Average: double): double;

BEGIN
  Result := -Ln(Random) * Average;
END;


FUNCTION RandomPoisson(Average: double): WORD;

VAR
  Count: WORD;
  ProbZero, Product: double;

BEGIN
  Count := 0;
  Product := Random;
  ProbZero := Exp(-Average);
  WHILE Product > ProbZero DO
    BEGIN
      INC(Count);
      Product := Product * Random;
    END;
  Result := Count;
END;


FUNCTION RandomGamma(Form: WORD; Average: double): double;

VAR
  Zaehler: WORD;
  Product: double;

BEGIN
  Product := 1;
  FOR Zaehler := 1 TO Form DO
    Product := Product * Random;
  Result := -Average * Ln(Product);
END;

FUNCTION RandomBernoulli(P: double): BYTE;

BEGIN
  Result := Ord(Random < P);
END;


FUNCTION RandomGeometric(P: double): WORD;

VAR
  Count: WORD;

BEGIN
  Count := 1;
  WHILE Random > P DO
    INC(Count);
  Result := Count;
END;


FUNCTION RandomChiSqr(DegFreedom: WORD): double;

VAR
  Zaehler: WORD;
  Sum: double;

BEGIN
  Sum := 0;
  FOR Zaehler := 1 TO DegFreedom DO
    Sum := Sum + Sqr(RandomNormal(0, 1));
  Result := Sum;
END;


FUNCTION RandomT(DegFreedom: WORD): double;

VAR
  Zaehler: WORD;
  Sum: double;

BEGIN
  Sum := 0;
  FOR Zaehler := 1 TO DegFreedom DO
    Sum := Sum + Sqr(RandomNormal(0, 1));
  Result := RandomNormal(0, 1) * Sqrt(DegFreedom / Sum);
END;


FUNCTION RandomF(f1, f2: WORD): double;

BEGIN
  Result := (RandomChiSqr(f1) / f1) / (RandomChiSqr(f2) / f2);
END;

FUNCTION RandomPareto(Minimum, Shape: double): double;

BEGIN
  Result := Minimum / Pot((1 - Random), 1 / Shape);
END;

END.    // Zufall

