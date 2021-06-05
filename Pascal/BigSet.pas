UNIT BigSet;
{ Sets with 0..pred(MaxCardinality) members. MaxCardinality should be a
  multiple of the WordSize (16 bit on Windows).
  The calling program needs to keep tab on the actual size of the set used. }

INTERFACE

USES Classes;

CONST MaxCardinality = 2048;                                       // change AS needed
      WordSize       = SizeOf(WORD) * 8;                           // IN bits
      MaxMembers     = Pred(MaxCardinality DIV (WordSize));        // number OF words needed TO fit elements

TYPE SetType = ARRAY [0..MaxMembers] OF WORD;

PROCEDURE SetBit (VAR S : SetType; Bit : WORD);
{ set respective bit to true (= member of set) }

PROCEDURE ClearBit (VAR S : SetType; Bit : WORD);
{ set respective bit to false (= not member of set) }

PROCEDURE ToggleBit (VAR S : SetType; Bit : WORD);
{ change value of respective bit from true to false or vice versa }

PROCEDURE ClearAllBits (VAR S : SetType);
{ initialize all elements to 0 (false, non-member) }

PROCEDURE SetAllBits (VAR S : SetType);
{ initialize all elements to 1 (true, member) }

FUNCTION InSet (CONST S : SetType; Bit : WORD) : BOOLEAN;
{ return status of respective bit }

PROCEDURE SetUnion (VAR S: SetType; CONST T, U : SetType);
{ calculates union S of sets T and U ( S = T or U)  }

PROCEDURE SetIntersection (VAR S: SetType; CONST T, U : SetType);
{ calculates intersection of T and U (S = T and U) }

PROCEDURE SetComplement (VAR S: SetType; CONST T, U : SetType);
{ calculates complement of T and U (all elements in T that are not in U, S = T - U)  }

PROCEDURE SetSymDifference (VAR S: SetType; CONST T, U : SetType);
{ calculates  ((T-U) or (U-T)) }

FUNCTION SubSet (CONST S, T : SetType) : BOOLEAN;
{ true if S is a subset of T (all members of S are also members of T, includes S = T) }

FUNCTION TrueSubSet (CONST S, T : SetType) : BOOLEAN;
{ true if S is a true subset of T (all members of S are also members of T, but S < T) }

FUNCTION EqualSet (CONST S, T : SetType) : BOOLEAN;
{ true if S and T are equal }

FUNCTION EmptySet (CONST S : SetType) : BOOLEAN;
{ true if all elements of S are 0 }

FUNCTION Cardinality (CONST S : SetType) : WORD;
{ Number of set members }

IMPLEMENTATION

CONST Bits : ARRAY [0..15] OF WORD = (   1,    2,    4,    8,   16,   32,    64,   128,
                                       256,  512, 1024, 2048, 4096, 8192, 16384, 32768);

PROCEDURE SetBit(VAR S : SetType; Bit : WORD);

VAR n : WORD;

BEGIN
  n := Bit DIV WordSize;         // identify the WORD that stores the bit
  S[n] := S[n] OR Bits[Bit MOD 16];
END;


PROCEDURE ClearBit (VAR S : SetType; Bit : WORD);

VAR n : WORD;

BEGIN
  n := Bit DIV WordSize;
  S[n] := S[n] AND NOT Bits[Bit MOD 16];
END;


PROCEDURE ToggleBit (VAR S : SetType; Bit : WORD);

VAR n : WORD;

BEGIN
  n := Bit DIV WordSize;
  S[n] := S[n] XOR Bits[Bit MOD 16];
END;


FUNCTION InSet (CONST S : SetType; Bit : WORD) : BOOLEAN;

VAR n : WORD;

BEGIN
  n := Bit DIV WordSize;
  Result := S[n] AND Bits[Bit MOD 16] <> 0;
END;

PROCEDURE ClearAllBits (VAR S : SetType);

VAR n : WORD;

BEGIN
  FOR n := 0 TO MaxMembers DO
    S[n] := 0;
END;

PROCEDURE SetAllBits (VAR S : SetType);

VAR n : WORD;

BEGIN
  FOR n := 0 TO MaxMembers DO
    S[n] := Pred(MaxCardinality);
END;

PROCEDURE SetUnion (VAR S: SetType; CONST T, U : SetType);

VAR n : WORD;
BEGIN
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF (InSet(T, n) OR InSet(U, n))
      THEN SetBit(S, n)
      ELSE ClearBit(S, n);
END;

PROCEDURE SetIntersection (VAR S: SetType; CONST T, U : SetType);

VAR n : WORD;
BEGIN
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF (InSet(T, n) AND InSet(U, n))
      THEN SetBit(S, n)
      ELSE ClearBit(S, n);
END;

PROCEDURE SetComplement (VAR S: SetType; CONST T, U : SetType);

VAR n : WORD;
BEGIN
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF (InSet(T, n) AND NOT(InSet(U, n)))
      THEN SetBit(S, n)
      ELSE ClearBit(S, n);
END;

PROCEDURE SetSymDifference (VAR S: SetType; CONST T, U : SetType);

VAR H1, H2 : SetType;

BEGIN
  SetComplement(H1, T, U);
  SetComplement(H2, U, T);
  SetUnion(S, H1, H2)
END;

FUNCTION SubSet (CONST S, T : SetType) : BOOLEAN;

VAR n : WORD;
BEGIN
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF (InSet(S, n) AND NOT InSet(T, n))  // found one member OF S that IS NOT IN T?
      THEN
        BEGIN
          SubSet := FALSE;
          EXIT;
        END;
  Result := TRUE;  // IF none found
END;

FUNCTION EqualSet (CONST S, T : SetType) : BOOLEAN;

VAR n : WORD;
BEGIN
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF (InSet(S, n) <> InSet(T, n))  // found one inequality?
      THEN
        BEGIN
          Result := FALSE;
          EXIT;
        END;
  Result := TRUE;
END;

FUNCTION TrueSubSet (CONST S, T : SetType) : BOOLEAN;

BEGIN
  Result := Subset(S, T) AND NOT EqualSet(S, T);
END;

FUNCTION EmptySet (CONST S : SetType) : BOOLEAN;

VAR n : WORD;
BEGIN
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF InSet(S, n)  // found one member?
      THEN
        BEGIN
          Result := FALSE;
          EXIT;
        END;
  Result := TRUE;
END;

FUNCTION Cardinality (CONST S : SetType) : WORD;

VAR n, c : WORD;
BEGIN
  c := 0;
  FOR n := 0 TO Pred(MaxCardinality) DO
    IF InSet(S, n) THEN INC(c);   // count members
  Result := c;
END;

END.

