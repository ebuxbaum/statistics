UNIT Dynam;

INTERFACE

USES MathFunc; // here used only for error handling

CONST DynamError : BOOLEAN = FALSE;
      MaxElem    = 2000;

TYPE
  Space = ARRAY[0..MaxElem] OF BYTE;

{ ***************************** LiFo ***************************** }

  LList = ^LListNode;

  LListNode = RECORD
    Size: WORD;
    Info: ^Space;
    Next: LList;
  END;

  LIFO = LList;

PROCEDURE InitLIFO(VAR L: Lifo);
{einen leeren Stapel initialisieren}

FUNCTION EmptyLIFO(L: LIFO): BOOLEAN;
{prüfen, ob der Stapel L leer ist}

PROCEDURE PUSH(VAR L: LIFO; VAR I; S: WORD);
{einen Wert auf den Stapel legen}

PROCEDURE POP(VAR L: Lifo; VAR I; S: WORD);
{einen Wert aus dem Stapel nehmen}

{ ***************************** FiFo ***************************** }

TYPE
  FList = ^FListNode;

  FListNode = RECORD
    Size: WORD;
    Info: ^Space;
    Next: FList;
  END;

  FIFO = RECORD
    First, last: FList;
  END;

PROCEDURE InitFIFO(VAR F: FIFO);
{legt eine leere Schlange an}

FUNCTION EmptyFIFO(F: FIFO): BOOLEAN;
{true, wenn die Schlange leer ist}

PROCEDURE Put(VAR F: FIFO; VAR I; S: WORD);
{einen Wert I mit Größe S in die Schlange schieben}

PROCEDURE Get(VAR F: FIFO; VAR I; S: WORD);
{einen Wert I der Größe S aus der Schlange holen}


IMPLEMENTATION

VAR ch : char; // for error handling

{*********************************  LiFo  ***********************************}

PROCEDURE InitLIFO(VAR L: LIFO);


BEGIN
  L := NIL;
END; (* InitLIFO *)


FUNCTION EmptyLIFO(L: LIFO): BOOLEAN;

BEGIN
  Result := L = NIL;
END; (* EmptyLIFO *)


PROCEDURE PUSH(VAR L: LIFO; VAR I; S: WORD);

VAR
  p: LList;

BEGIN
  NEW(p);
  WITH p^ DO
    BEGIN
      Size := S;
      Next := L;
      TRY
        GetMem(Info, Size);
      EXCEPT
        ch := WriteErrorMessage('not enough memory');
        DynamError := TRUE;
        EXIT;
      END;
      Move(I, Info^, Size);
    END; (* WITH *)
  L := p;
END; (* Push *)


PROCEDURE POP(VAR L: LIFO; VAR I; S: WORD);

VAR
  p: LList;

BEGIN
  p := L;
  WITH p^ DO
    BEGIN
      IF Size <> S
        THEN
          BEGIN
            Writeln('TYPE-MISMATCH-ERROR');
            HALT;
          END;
      Move(Info^, I, Size);
      FreeMem(Info, Size);
      L := Next;
    END; (* WITH *)
  DISPOSE(p);
END; (* Pop *)

{********************************* FiFo *************************************}

PROCEDURE InitFiFo(VAR F: FiFo);

BEGIN
  WITH F DO
    BEGIN
      First := NIL;
      last := NIL;
    END;
END;


FUNCTION EmptyFiFo(F: FiFo): BOOLEAN;

BEGIN
  Result := F.First = NIL;
END;


PROCEDURE Put(VAR F: FiFo; VAR i; s: WORD);

VAR p: FList;

BEGIN
  NEW(p);
  WITH p^ DO
    BEGIN
      Size := S;
      Next := NIL;
      TRY
        GetMem(Info, Size);
      EXCEPT
        ch := WriteErrorMessage('not enough memory');
        DynamError := TRUE;
        EXIT;
      END;
      Move(I, Info^, Size);
    END;
  IF EmptyFiFo(F)
    THEN
      BEGIN
        F.First := p;
        F.last := p;
      END
    ELSE
      BEGIN
        F.last^.Next := p;
        F.last := p;
      END;
END;


PROCEDURE Get(VAR F: FiFo; VAR i; s: WORD);

VAR
  p: FList;

BEGIN
  p := F.First;
  WITH p^ DO
    BEGIN
      IF Size <> s
        THEN
          BEGIN
            Writeln('Type-Mismatch-Error');
            HALT;
          END;
      Move(Info^, I, Size);
      FreeMem(Info, Size);
      F.First := Next;
    END;
  DISPOSE(p);
END;


END.  { Dynam }

