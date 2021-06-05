UNIT SimplexFit;

{Regressionsanalyse nach M.Caceci, W.P. Cacheris: Fitting Curves to Data,
      Byte Magazine May 1984, S. 340-350
 Unter Benutzung des Formel-Compilers aus Pascal International 8/1987, S. 52-60
 Die Bestimmung der Fehlergrenzen durch Monte-Carlo-Simulation ist beschrieben
      von Straume & Johnson, Meth. Enzymol. 210 (1992) 117-129
 Die Idee der Minimierung von Chi2 für Daten mit bekannten Fehlergrenzen stammt
      aus Press et al.: Numerical Recepies in Pascal, Cambridge 1989
 Gesamtprogramm copyright 1988-1994 by Dr. Engelbert Buxbaum }

INTERFACE

USES MathFunc,            // basic maths
     crt,                 // low level system calls
     Calc,                // formula compiler
     Vector,              // vector arithmetic
     Matrix,              // matrix arithmetic
     Deskript,            // descriptive statistics
     Zufall               // random numbers
     ;

CONST
  SimplexError: BOOLEAN = FALSE;

PROCEDURE Approximation(Data: MatrixTyp; VAR yBerech: VectorTyp;
  VAR ProblemName, Formel, xName, yName: STRING);


IMPLEMENTATION

CONST
  alfa                  = 1.0;                 // Reflektion coefficient
  beta                  = 0.5;                 // Kontraktion coefficient
  gamma                 = 2.0;                 // Expansion coefficient
  MaxParameter          = 10;
  MaxVariablen          = 10;
  MaxN                  = MaxParameter + 1;    // dimension of simplex
  lw                    = 5;                   // linewidth of data field + 1
  Page                  = 12;
  Ja                    = 'Y';                 // key for yes
  Nein                  = 'N';                 // key for no

var ch : char;                                     // for error handling

PROCEDURE Approximation(Data: MatrixTyp; VAR yBerech: VectorTyp;
  VAR ProblemName, Formel, xName, yName: STRING);

TYPE
  AVektor = ARRAY[1..MaxN] OF float;
  DataRow = ARRAY [1..MaxVariablen] OF float;
  simpl   = ARRAY [1..MaxN] OF AVektor;            // Simplex
  ParFeld = ARRAY [1..MaxParameter] OF STRING[10];
  VarFeld = ARRAY [1..MaxVariablen] OF STRING[10];

VAR
  FehlerGrenzen,                    // Letzte Spalte der Daten-Matrix Fehlergrenzen?
  done: BOOLEAN;                    // Konvergenz

  i, j, n: BYTE;

  h, l: ARRAY [1..MaxN] OF BYTE;    // Zahl Hoch/Niedrig Parameter

  Daten,                            // Zahl der Datenpunkte
  MaxIter,                          // max.Zahl der Iterationen
  NIter: WORD;                      // Zahl der Iterationen

  Next,                             // neuer Vertex zu testen
  center,                           // Minimum aller Vertexe}
  mean, error, MaxErr,              // Maximal zulaessiger Fehler}
  p, q,                             // um ersten Simplex zu berechnen
  step: AVektor;                    // Eingabe Startschritte

  simp: simpl;                      // Simplex

  NameA,                            // Name des Eingabefiles
  NameB: STRING[64];                // Name des Ausgabefiles

  Eindat,                           // Eingabefile
  Ausdat: TEXT;                     // Ausgabefile}

  x: Calc_VarTab;                   // FOR formula compiler
  dummy, Sigma: float;              // y-Standardabweichung
  Formula: Calc_String;
  FormProg: Calc_Prog;

  Variablen, Parameter: BYTE;

  ErrorStat,                        // TRUE wenn Fehlerstatistik gewuenscht
  Erster: BOOLEAN;                  // Startsimplex nur einmal ausgeben

  Antwort: CHAR;

  Methode: (Summe, Median, ChiSqr);  // was soll minimiert werden?


  {****************************************************************************}

  PROCEDURE LiesFunction(VAR x: Calc_VarTab; VAR Formula: Calc_String;
  VAR FormProg: Calc_Prog);

  VAR
    i: BYTE;
    dummy: float;
    Name: Calc_IdStr;
    Input : TEXT;


    PROCEDURE Hilfe;

    BEGIN
      Writeln('The function used for fitting to the data should be entered in the ');
      Writeln('following manner: ');
      Writeln('1) Number and names of the variables (i. e. measured data)');
      Writeln('2) Number and name(s) of the parameters');
      Writeln('3) The formula itself. All names must be entered exactly as defined.');
      Writeln('   No undefined names are allowed. The formula must end with a ');
      Writeln('   semicolon.');
      Writeln;
      Writeln('The compiler ''knows'' the following constants and and functions, which ');
      Writeln('can not be redefined:');
      Writeln('Constants: e, pi                        Basic operators: +, -, *, /, ^');
      Writeln('Integer: div, mod, ggt, kgv             Logarithms: ln, lg, ld, exp');
      Writeln('sin, cos, tan, cot and the equivalent hyperbolic and arcus functions');
      Writeln('Various Functions: abs, deg, rad, fak, sgn');
      Writeln;
    END;


  BEGIN
    Hilfe;
    Assign(Input, 'CON');
    Reset(Input);
    IF SimplexError THEN EXIT;
    CalcDecMod := TRUE;                     {nur definierte Vars und Parms}
    x := NewVarTab;
    Writeln('Data file has ', MatrixColumns(Data), 'columns');
    IF MatrixColumns(Data) > 2
      THEN
        BEGIN
          Write('Last column independend variable or error margin (v/e):');
          REPEAT
            Readln(Antwort);
            Antwort := UpCase(Antwort);
          UNTIL (Antwort IN ['V', 'F', 'E', #27]);
          IF Antwort = #27
            THEN
              BEGIN
                SimplexError := TRUE;
                EXIT;
              END;
          FehlerGrenzen := (Antwort = 'F') OR (Antwort = 'E');
        END
      ELSE
        FehlerGrenzen := FALSE;
    IF FehlerGrenzen
      THEN Variablen := Pred(MatrixColumns(Data))
      ELSE Variablen := MatrixColumns(Data);
    Daten := MatrixRows(Data);
    FOR i := 1 TO Pred(Variablen) DO
      BEGIN
        Write('Name of the ', i, '. (independent) variable: ');
        Readln(Name);
        dummy := AddToVarTab(x, Name);
      END;
    Write('Name of the ', Variablen, '. (dependent) variable: ');
    Readln(Name);
    dummy := AddToVarTab(x, Name);
    Write('How many parameters do you want to use [1..', MaxParameter, ']: ');
    ReadLn(Parameter);
    N := Parameter + 1;               {Dimensionen des Simplex}
    FOR i := 1 TO Parameter DO
      BEGIN
        Write('Name of the ', i, '. parameter: ');
        ReadLn(Name);
        dummy := AddToVarTab(x, Name);
      END;
    REPEAT
      Writeln('Please enter the equation: ');
      Writeln;
      Write(x^[Variablen].VarId, ' = ');
      ReadLn(Formula);
      CompileExpression(Formula, x, FormProg);
      IF NOT CalcResult
        THEN Writeln('Unable to compile the equation, please try again: ');
    UNTIL CalcResult;
    Formel := x^[Variablen].VarID + ' = ' + Formula;
    xName := x^[1].VarID;
    yName := x^[Variablen].VarID;
  END;


  FUNCTION f(p: AVektor; d: MatrixTyp; Zeile: WORD): float;

  VAR
    i: BYTE;

  BEGIN
    FOR i := 1 TO Variablen DO
      AssignVar(x, x^[i].VarId, GetMatrixElement(d, Zeile, i));
    FOR i := 1 TO Parameter DO
      AssignVar(x, x^[Variablen + i].VarId, p[i]);
    Result := CalcExpression(FormProg, x);
  END;


  PROCEDURE Inparam(VAR MaxIter: WORD; VAR simp: Simpl; VAR Step, MaxErr: AVektor);
  {Einlesen aller benutzerdefinierten Parameter}

  VAR
    FalscheEingabe: BOOLEAN;
    Quelle: CHAR;
    i: WORD;
    c: CHAR;

  BEGIN
    Writeln('This routine calculates curve fits by the simplex algorithm');
    Writeln;
    LiesFunction(x, Formula, FormProg);
    IF SimplexError THEN EXIT;
    REPEAT
      FalscheEingabe := FALSE;
      Writeln;
      Write('Please enter the name of the output file: ');
      ReadLn(NameB);
      Assign(Ausdat, NameB);
      Rewrite(Ausdat);
      IF IOResult <> 2
        THEN   { d. h., Datei existiest schon }
          BEGIN
            REPEAT
              Write(NameB, ' already exists. Overwrite (y/n): ');
              ReadLn(c);
              c := UpCase(c);
            UNTIL (c = JA) OR (c = Nein) OR (c = #27);
          IF (c = #27)
            THEN
              BEGIN
                SimplexError := TRUE;
                EXIT;
              END;
          FalscheEingabe := (c = Nein);
        END;
    UNTIL NOT FalscheEingabe;
    Rewrite(Ausdat);
    Write(Ausdat, 'best fit for equation: ');
    Writeln(Ausdat, x^[Variablen].VarId, ' = ', Formula);
    Writeln(Ausdat);
    REPEAT
      FalscheEingabe := FALSE;
      REPEAT
        IF FehlerGrenzen
          THEN Write('Minimise sum of squares, median of squares or Chi2 (S/M/X): ')
          ELSE Write('Minimise sum or median of squares (S/M): ');
        ReadLn(c);
        c := UpCase(c);
        IF NOT ((c = 'S') OR (c = 'M') OR (c = 'X') OR (c = #27))
          THEN
            BEGIN
              Sound(400);
              Delay(50);
              NoSound;
            END;
      UNTIL (c = 'S') OR (c = 'M') OR (c = 'X') OR (c = #27);
      CASE c OF
        'S': BEGIN
               Methode := Summe;
               Writeln(Ausdat, 'Minimising sum of squared residuals ');
             END;
        'M': BEGIN
               Methode := Median;
               Writeln(Ausdat, 'Minimising median of squared residuals ');
             END;
        'X': BEGIN
               IF FehlerGrenzen
                 THEN
                   BEGIN
                     Methode := ChiSqr;
                     Writeln(Ausdat, 'Minimizing Chi2');
                   END
                 ELSE
                   BEGIN
                     ch := WriteErrorMessage('Chi2 erfordert Fehlergrenzen in der letzten Daten-Spalte');
                     SimplexError := TRUE;
                     EXIT;
                   END;
             END;
        #27: BEGIN
               SimplexError := TRUE;
               EXIT;
             END;
      END; { case }
      Writeln(Ausdat);
      Write(Ausdat, '                            ');
      FOR i := 1 TO Parameter DO
        Write(Ausdat, x^[i + Variablen].VarId,
          ' ': (ValidFigures + 2 - Length(x^[i + Variablen].VarId)));
      CASE Methode OF
        Summe: Writeln(Ausdat, 'sum of squares ');
        Median: Writeln(Ausdat, 'median of squares');
        ChiSqr: Writeln(Ausdat, 'Chi2');
      END;
    UNTIL NOT FalscheEingabe;
    REPEAT
      FalscheEingabe := FALSE;
      Write('Please enter maximal number of iterations: ');
      ReadLn(MaxIter);
    UNTIL NOT FalscheEingabe;
    REPEAT
      FalscheEingabe := FALSE;
      Writeln('Please enter initial guesses for all parameters: ');
      Write(Ausdat, 'Start coordinates:      ');
      FOR i := 1 TO Parameter DO
        BEGIN
          Write(x^[i + Variablen].VarId, ' = ');
          ReadLn(simp[1, i]);
          IF (i MOD lw) = 0 THEN Writeln(Ausdat);
          Write(Ausdat, FloatStr(simp[1, i], ValidFigures), '  ');
        END;
      Writeln(Ausdat);
      Writeln(Ausdat);
    UNTIL NOT FalscheEingabe;
    REPEAT
      FalscheEingabe := FALSE;
      Writeln('Please enter starting step width for all parameters ');
      Writeln('(ca. 1/10 to 1/2 of initial value)');
      Write(Ausdat, 'Starting step width:    ');
      FOR i := 1 TO Parameter DO
        BEGIN
          Write('Step width:     ', x^[i + Variablen].VarId, ' = ');
          ReadLn(step[i]);
          IF (i MOD lw) = 0 THEN Writeln(Ausdat);
          Write(Ausdat, FloatStr(step[i], ValidFigures), '  ');
        END;
      Writeln(Ausdat);
      Writeln(Ausdat);
      Write(Ausdat, 'max. alowable error:    ');
      FOR i := 1 TO n DO
        BEGIN
          MaxErr[i] := MaxError;
          IF (i MOD lw) = 0 THEN Writeln(Ausdat);
          Write(Ausdat, FloatStr(MaxErr[i], ValidFigures), '  ');
        END;
      Writeln(Ausdat);
      Writeln(Ausdat);
    UNTIL NOT FalscheEingabe;
    REPEAT
      Write('Do you want error margins for the parameters (y/n): ');
      ReadLn(c);
      c := UpCase(c);
    UNTIL (c = JA) OR (c = Nein) OR (c = #27);
    IF (c = #27)
      THEN
        BEGIN
          SimplexError := TRUE;
          EXIT;
        END;
    Writeln(UpCase(c));
    ErrorStat := (c = JA);
  END;


  PROCEDURE sum_of_residuals(VAR z: AVektor; Data: MatrixTyp);
  {Berechnet die Summe der Fehlerquadrate}

  VAR
    i: WORD;

  BEGIN
    z[n] := 0.0;
    FOR i := 1 TO Daten DO
      z[n] := z[n] + Sqr(f(z, Data, i) - GetMatrixElement(Data, i, Variablen));
  END;


  PROCEDURE Median_Of_Squares(VAR z: AVektor; Data: MatrixTyp);
  { berechnet den Median der Fehlerquadrate }

  VAR
    i: WORD;
    Residuals: VectorTyp;

  BEGIN
    CreateVector(Residuals, Daten, 0.0);
    FOR i := 1 TO Daten DO
      SetVectorElement(Residuals, i, Sqr(f(z, Data, i) -
        GetMatrixElement(Data, i, Variablen)));
    ShellSort(Residuals);
    z[n] := Quantile(Residuals, 0.5);
    DestroyVector(Residuals);
  END;

  PROCEDURE Chi(VAR z: AVektor; Data: MatrixTyp);
  {Berechnet chi2}

  VAR
    i: WORD;

  BEGIN
    z[n] := 0.0;
    FOR i := 1 TO Daten DO
      z[n] := z[n] + Sqr(f(z, Data, i) - GetMatrixElement(Data, i, Variablen) /
        GetMatrixElement(Data, i, Succ(Variablen)));
  END;


  PROCEDURE report;
  {berichtet Programstatistik}

   VAR
    dy, h, Zaehler, Nenner, Fehler, Mittel, rSqr: float;
    d1, d2: TEXT;
    HilfsStr: STRING[14];
    i, j: WORD;

  BEGIN
    Writeln(Ausdat);
    Writeln(Ausdat, 'Routine was left after ', NIter: 5, ' iterations');
    Writeln(Ausdat);
    Writeln(Ausdat);
    Writeln(Ausdat, 'The final simplex is: ');
    Write(Ausdat, '                        ');
    FOR j := 1 TO n DO
      BEGIN
        FOR i := 1 TO n DO
          BEGIN
            IF (i MOD lw) = 0 THEN
              Writeln(Ausdat);
            Write(Ausdat, FloatStr(simp[j, i], ValidFigures), '  ');
          END;
        Writeln(Ausdat);
        Write(Ausdat, '                        ');
      END;
    Writeln(Ausdat);
    Write(Ausdat, 'The mean is:            ');
    FOR i := 1 TO n DO
      BEGIN
        IF (i MOD lw) = 0 THEN
          Writeln(Ausdat);
        Write(Ausdat, FloatStr(mean[i], ValidFigures), '  ');
      END;
    Writeln(Ausdat);
    Writeln(Ausdat);
    Write(Ausdat, 'error:                  ');
    FOR i := 1 TO n DO
      BEGIN
        IF (i MOD lw) = 0 THEN
          Writeln(Ausdat);
        Write(Ausdat, FloatStr(error[i], ValidFigures), '  ');
      END;
    Writeln(Ausdat);
    Writeln(Ausdat);
    Write(Ausdat, '  #    ');
    FOR i := 1 TO Variablen DO
      Write(Ausdat, x^[i].VarId, ' ': (ValidFigures + 2 - Length(x^[i].VarId)));
    IF Fehlergrenzen THEN
      Write(Ausdat, 'Delta ', x^[Variablen].VarId, ' ': (ValidFigures -
        2 - Length(x^[i].VarId)));
    Write(Ausdat, x^[Variablen].VarId, ' ber.',
      ' ': (ValidFigures - 2 - Length(x^[Variablen].VarId)));
    Writeln(Ausdat, 'error:                  ');
    Zaehler := 0;
    sigma := 0.0;
    FOR i := 1 TO Daten DO
      BEGIN
        h := f(mean, Data, i);
        SetVectorElement(yBerech, i, h);
        dy := GetMatrixElement(Data, i, Variablen) - h;
        sigma := sigma + Sqr(dy);
        Write(Ausdat, i: 4, '  ');
        FOR j := 1 TO MatrixColumns(Data) DO
          Write(Ausdat, FloatStr(GetMatrixElement(Data, i, j), ValidFigures), '  ');
        Writeln(Ausdat, FloatStr(h, ValidFigures), '  ', FloatStr(dy, ValidFigures));
        Zaehler := Zaehler + GetMatrixElement(Data, i, Variablen);
      END;
    Writeln(Ausdat);
    sigma := Sqrt(sigma / Daten);
    Writeln(Ausdat, 'The standard deviation is:    ', FloatStr(sigma, ValidFigures));
    Fehler := sigma / Sqrt(Daten - Parameter);
    Writeln(Ausdat, 'The error of the function is: ', FloatStr(Fehler, ValidFigures));
    Mittel := Zaehler / Daten;
    Zaehler := 0;
    Nenner := 0;
    FOR i := 1 TO Daten DO
      BEGIN
        Zaehler := Zaehler + Sqr(GetVectorElement(yBerech, i) - Mittel);
        Nenner := Nenner + Sqr(GetMatrixElement(Data, i, Variablen) - Mittel);
      END;
    rSqr := Zaehler / Nenner;
    Writeln(Ausdat, 'r2:                           ', FloatStr(rSqr, ValidFigures));
  END;


  PROCEDURE First;

  VAR
    i, j: WORD;

  BEGIN
    Write(Ausdat, 'Start simplex ');
    FOR j := 1 TO n DO
      BEGIN
        Write(Ausdat, ' simp[', j: 3, ']');
        FOR i := 1 TO n DO
          BEGIN
            IF (i MOD lw) = 0 THEN Writeln(Ausdat);
            Write(Ausdat, FloatStr(simp[j, i], ValidFigures), '  ');
          END;
        Writeln(Ausdat);
        Write(Ausdat, '              ');
      END;
    Writeln(Ausdat);
  END;


  PROCEDURE new_vertex;
  {ersetzt worst durch next}

  VAR
    i, j: WORD;

  BEGIN
    IF erster THEN Write(' --- ', NIter: 4);
    FOR i := 1 TO n DO
      BEGIN
        simp[h[n], i] := Next[i];
        IF erster THEN Write(FloatStr(Next[i], ValidFigures));
      END;
    IF erster THEN Writeln;
  END;


  PROCEDURE order;
  {Highs und Lows fuer jeden Parameter}

  VAR
    i, j: BYTE;

  BEGIN
    FOR j := 1 TO n DO
      BEGIN
        FOR i := 1 TO n DO
          BEGIN
            IF simp[i, j] < simp[l[j], j] THEN l[j] := i;
            IF simp[i, j] > simp[h[j], j] THEN h[j] := i;
          END;
      END;
  END;


  PROCEDURE Iteration(Data: MatrixTyp);

  VAR
    i, j: WORD;

  BEGIN
    NIter := 0;
    REPEAT
      done := TRUE;
      NIter := Succ(NIter);
      FOR i := 1 TO n DO center[i] := 0.0;
      FOR i := 1 TO n DO
        IF i <> h[n]
          THEN
            FOR j := 1 TO Parameter DO
              center[j] := center[j] + simp[i, j];
      FOR i := 1 TO n DO
        BEGIN
          center[i] := center[i] / Parameter;
          Next[i] := (1.0 + alfa) * center[i] - alfa * simp[h[n], i];
        END;
      CASE Methode OF
        Summe: sum_of_residuals(Next, Data);
        Median: Median_Of_Squares(Next, Data);
        ChiSqr: Chi(Next, Data);
      END;
      IF Next[n] <= simp[l[n], n]
        THEN
          BEGIN
            new_vertex;
            FOR i := 1 TO n DO
              Next[i] := gamma * simp[h[n], i] + (1.0 - gamma) * center[i];
            CASE Methode OF
              Summe: sum_of_residuals(Next, Data);
              Median: Median_Of_Squares(Next, Data);
              ChiSqr: Chi(Next, Data);
            END;
            IF Next[n] <= simp[l[n], n] THEN new_vertex;
          END
        ELSE
          BEGIN
            IF Next[n] <= simp[h[n], n]
              THEN
                new_vertex
              ELSE
                BEGIN
                  FOR i := 1 TO Parameter DO
                    Next[i] := beta * simp[h[n], i] + (1.0 - beta) * center[i];
                  CASE Methode OF
                    Summe: sum_of_residuals(Next, Data);
                    Median: Median_Of_Squares(Next, Data);
                    ChiSqr: Chi(Next, Data);
                  END;
                  IF Next[n] <= simp[h[n], n]
                    THEN
                      new_vertex
                    ELSE
                      BEGIN
                        FOR j := 1 TO Parameter DO
                          BEGIN
                            simp[i, j] := (simp[i, j] + simp[l[n], j]) * beta;
                            CASE Methode OF
                              Summe: sum_of_residuals(simp[i], Data);
                              Median: Median_Of_Squares(simp[i], Data);
                              ChiSqr: Chi(simp[i], Data);
                            END; // case
                          END; // for j
                      END; // else
                END; // else
          END; //else
      order;
      FOR j := 1 TO n DO
        BEGIN
          error[j] := (simp[h[j], j] - simp[l[j], j]) / simp[h[j], j];
          IF done
            THEN IF error[j] > MaxErr[j]
                   THEN done := FALSE;
        END;
    UNTIL (done OR (NIter = MaxIter));
  END;


  PROCEDURE DoIteration(VAR Simp: Simpl; Data: MatrixTyp; VAR Mean: AVektor);

  VAR
    i, j: WORD;

  BEGIN
    CASE Methode OF
      Summe: sum_of_residuals(simp[1], Data);
      Median: Median_Of_Squares(simp[1], Data);
      ChiSqr: Chi(simp[1], Data);
    END;
    FOR i := 1 TO Parameter DO
      BEGIN
        p[i] := step[i] * (Sqrt(n) + Parameter - 1) / (Parameter * Sqrt(2));
        q[i] := step[i] * (Sqrt(n) - 1) / (Parameter * Sqrt(2));
      END;
    FOR i := 2 TO n DO
      BEGIN
        FOR j := 1 TO Parameter DO simp[i, j] := simp[1, j] + q[j];
        simp[i, i - 1] := simp[1, i - 1] + p[i - 1];
        CASE Methode OF
          Summe: sum_of_residuals(simp[i], Data);
          Median: Median_Of_Squares(simp[i], Data);
          ChiSqr: Chi(simp[i], Data);
        END;
      END;
    FOR i := 1 TO n DO
      BEGIN
        l[i] := 1;
        h[i] := 1;
      END;
    order;
    IF Erster THEN First;
    Iteration(Data);
    FOR i := 1 TO n DO
      BEGIN
        mean[i] := 0.0;
        FOR j := 1 TO n DO
          mean[i] := mean[i] + simp[j, i];
        mean[i] := mean[i] / n;
      END;
  END;


  PROCEDURE CalculateErrorMargins;

  CONST
    Anzahl = 100;

  VAR
    Counter, Differenz, Ergebniss: AVektor;
    SimulatedData, Ergebnisse: MatrixTyp;
    i, j, k: WORD;


    PROCEDURE Simulate(VAR SimulatedData: MatrixTyp; yBerech: VectorTyp;
      Sigma: float);

    VAR
      i: WORD;
      Wert: float;

    BEGIN
      FOR i := 1 TO MatrixColumns(SimulatedData) DO
        BEGIN
          Wert := RandomLaplace(GetVectorElement(yBerech, i), Sigma);
          SetMatrixElement(SimulatedData, i, Variablen, Wert);
        END;
    END;

  BEGIN
    Sigma := Sqr(Sigma);
    CopyMatrix(Data, SimulatedData);
    CreateMatrix(Ergebnisse, Anzahl, N, 0.0);
    FOR j := 1 TO N DO
      Counter[j] := 0.0;
    FOR i := 1 TO 100 DO
      BEGIN
        Write(i: 3, '  ');
        Simulate(SimulatedData, yBerech, Sigma);
        DoIteration(Simp, SimulatedData, Ergebniss);
        FOR j := 1 TO N DO
          BEGIN
            SetMatrixElement(Ergebnisse, i, j, Ergebniss[j]);
            Counter[j] := Counter[j] + Ergebniss[j];
            Write(FloatStr(Ergebniss[j], ValidFigures), '  ');
          END;
        Writeln;
      END;
    DestroyMatrix(SimulatedData);
    FOR j := 1 TO N DO
      BEGIN
        Counter[j] := Counter[j] / Anzahl;       { Mittelwerte der Parameter }
        Differenz[j] := 0.0;
      END;
    FOR i := 1 TO Anzahl DO
      FOR j := 1 TO N DO
        Differenz[j] := Sqr(Counter[j] - GetMatrixElement(Ergebnisse, i, j));
    Writeln(Ausdat);
    Writeln(Ausdat, 'Standard deviation of the parameters: ');
    Write(Ausdat, '                        ');
    FOR j := 1 TO N DO
      BEGIN
        Differenz[j] := Sqrt(Differenz[j] / Pred(Anzahl));
        IF (j MOD lw) = 0 THEN Writeln(Ausdat);
        Write(Ausdat, FloatStr(Sqrt(Differenz[j]), ValidFigures), '  ');
      END;
    Writeln(Ausdat);
  END;

BEGIN {Approximation}
  Erster := TRUE;
  Inparam(MaxIter, simp, step, MaxErr);
  IF SimplexError THEN EXIT;
  DoIteration(Simp, Data, Mean);
  CreateVector(yBerech, Daten, 0.0);
  report;
  Erster := FALSE;
  IF ErrorStat THEN CalculateErrorMargins;
  Close(Ausdat);
END;

END.

