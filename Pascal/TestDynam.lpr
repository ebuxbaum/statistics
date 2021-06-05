PROGRAM TestDynam;

USES dynam, mathfunc, crt;

{$F+} FUNCTION Kleiner(VAR a, b: STRING): BOOLEAN; {$F+}
{diese Funktion darf nicht lokal zu BinaerBaumTesten deklariert werden }
BEGIN
  Kleiner := a < b
END;

PROCEDURE BinaerBaumTesten;

TYPE str20 = STRING[20];

VAR Anzahl, i   : LONGINT;
    TempStr     : str20;
    Baum, temp  : bbtree;
    OK          : BOOLEAN;

  PROCEDURE ZufallsString(VAR ts: str20);

  VAR i: WORD;

  BEGIN
    ts[0] := Chr(20);
    FOR i := 1 TO 20 DO
      ts[i] := Chr(Random(27) +64);
  END;

BEGIN
  CheckBreak := TRUE;
  NewBBTree(Baum);
  Write('Bitte Anzahl der Zufallsstrings angeben: ');
  ReadLn(Anzahl);
  Writeln(Anzahl, ' Zufallsstrings werden erzeugt und im Baum abgelegt.');
  FOR i := 1 TO Anzahl DO
    BEGIN
      ZufallsString(TempStr);
      InsVarBBTree(Baum, TempStr, SizeOf(TempStr), ComparisonTyp(@Kleiner), OK)
    END;
  FirstVarBBTree(Baum, TempStr, OK);
  WHILE OK DO
    BEGIN
      Writeln(TempStr);
      NextVarBBTree(Baum, TempStr, ComparisonTyp(@Kleiner), OK);
    END;
  Writeln;
  Writeln('Ende des Tests Binär-Baum, weiter mit RETURN:');
  ReadLn;
END;

PROCEDURE LifoFifoTesten;

VAR s          : STRING[20];
    LifoListe  : Lifo;
    FifoListe  : Fifo;
    i          : WORD;

BEGIN
  ClrScr;
  Writeln('Demoprogramm fr die Arbeit mit Fifo und Lifo-Strukturen');
  Writeln;
  InitLifo(LifoListe);
  InitFifo(FifoListe);
  Writeln('Bitte geben Sie jetzt 10 kurze Strings ein (z.B. Namen)');
  Writeln;
  Writeln('EINGABE', 'AUSGABE von LIFO': 30, 'AUSGABE von FIFO': 30);
  Writeln;
  FOR i := 1 TO 10 DO
    BEGIN
      ReadLn(s);
      PUSH(LifoListe, s, SizeOf(s));
      Put(FifoListe, s, SizeOf(s))
    END;
  GotoXY(1, 7);
  WHILE NOT EmptyLifo(LifoListe) DO
    BEGIN
      GotoXY(30, WhereY);
      POP(LifoListe, s, SizeOf(s));
      Writeln(s)
    END;
  GotoXY(1, 7);
  WHILE NOT EmptyFifo(FifoListe) DO
    BEGIN
      GotoXY(60, WhereY);
      Get(FifoListe, s, SizeOf(s));
      Writeln(s)
    END;
  Writeln;
  Writeln('weiter mit RETURN:');
  ReadLn;
END; {LifoFifoTesten}

TYPE KeyType = WORD;

{$F+} FUNCTION less (VAR X, Y: KeyType): BOOLEAN;{$F-}
{auch diese Funktion darf nicht lokal zu BayerBaumTesten deklariert sein}
BEGIN
  Result := X < Y
END; (* less *)

PROCEDURE BayerBaumTesten;

VAR  B: BTree;
     k, t: KeyType;
     i, IO: INTEGER;
     Okay: BOOLEAN;

BEGIN
  NewTree(B, 'temp', 16, SizeOf(KeyType), ComparisonTyp(@less), IO);
  IF IO <> 0
    THEN
      BEGIN
        Writeln('Fehler bei der Baum-Erstellung!');
        HALT
      END;
  FOR i := 1 TO 500 DO
    BEGIN
      REPEAT
        k := 1 + Random(1000);
        t := k;
        FindKey(B, t, Okay)
      UNTIL NOT Okay;
      InsertKey(B, k, Okay);
      Writeln(i: 3, k: 5, Okay: 10)
    END;
  CloseTree(B, IO);
  IF IO <> 0
    THEN
      BEGIN
        Writeln('Baum kann nicht abgeschlossen werden!');
        HALT
      END;
  Write('Weiter mit RETURN');
  ReadLn;
  OpenTree(B, 'temp', ComparisonTyp(@less), IO);
  Writeln('Höhe des Baumes: ', Height(B): 1: 0);
  Writeln('Kapazität      : ', MaxEntries(B): 1: 0);
  Writeln('Mittlerer Zugr.: ', Access(B): 1: 2);
  Writeln('Ordnung        : ', B.n);
  Writeln('Seiten         : ', B.Pages);
  Write('Weiter mit RETURN');
  ReadLn;
  IF IO <> 0
    THEN
      BEGIN
        Writeln('Fehler bei der Baum-Öffnung!');
        HALT
      END;
  FirstKey(B, k, Okay);
  WHILE Okay DO
    BEGIN
      Write(k: 4);
      NextKey(B, k, Okay)
    END;
  Write('Weiter mit RETURN');
  ReadLn;
  FOR i := 1 TO 300 DO
    BEGIN
      REPEAT
        k := 1 + Random(1000);
        FindKey(B, k, Okay)
      UNTIL Okay;
      DeleteKey(B, k, Okay);
      Writeln(i: 3, k: 5, Okay: 10)
    END;
  Write('Weiter mit RETURN');
  ReadLn;
  LastKey(B, k, Okay);
  WHILE Okay DO
    BEGIN
      Write(k: 4);
      PrevKey(B, k, Okay)
    END;
  CloseTree(B, IO)
END;

BEGIN  // TEST PROGRAM
  LifoFifoTesten;
  BinaerBaumTesten;
  BayerBaumTesten;
END.
