UNIT Vector;

INTERFACE

USES Math, mathfunc;

CONST
  MaxVector = 8000;             // PowerOfTwo[16] DIV SizeOf(float);
  VectorError: BOOLEAN = FALSE; // toggle IF error occurs.


TYPE
  VectorStruc = RECORD
    Length: WORD;
    Data: ARRAY [1..MaxVector] OF float;
  END;
  VectorTyp = ^VectorStruc;

PROCEDURE CreateVector(VAR Vec: VectorTyp; Length: WORD; Value: float);
{ Generates a vector of length Length and sets all elements to Value }

PROCEDURE DestroyVector(VAR Vec: VectorTyp);
{ release memory occupied by vector }

PROCEDURE ReadVector(VAR Vec: VectorTyp; MedStr: STRING);
{ read vector from the file indicated by MedStr }

PROCEDURE WriteVector(MedStr: STRING; CONST Vec: VectorTyp; ValidFigures: BYTE);
{ write Vector to file indicated by MedStr }

FUNCTION VectorLength(CONST Vec: VectorTyp): WORD;

FUNCTION GetVectorElement(CONST Vec: VectorTyp; n: WORD): float;

PROCEDURE SetVectorElement(VAR Vec: VectorTyp; n: WORD; c: float);

PROCEDURE CopyVector(CONST Source: VectorTyp; VAR Dest: VectorTyp);
{ use instead of Dest := Source! }

PROCEDURE LoadConstant(VAR Vec: VectorTyp; C: float);
{ set all elements of vector to C }

PROCEDURE AddConstant(VAR A: VectorTyp; C: float);
{ Add C to all elements of A and save result in B }

PROCEDURE SubConstant(VAR A: VectorTyp; C: float);

PROCEDURE MulConstant(VAR A: VectorTyp; C: float);

PROCEDURE DivConstant(VAR A: VectorTyp; C: float);

PROCEDURE VectorAdd(CONST A, B: VectorTyp; VAR C: VectorTyp);
{ add to all i-th element of A the i-th element of B and store in the i-th element of C }

PROCEDURE VectorSub(CONST A, B: VectorTyp; VAR C: VectorTyp);
{ subtract from i-th element of A the i-th element of B and store in the i-th element of C }

FUNCTION VectorInnerProduct(CONST A, B: VectorTyp): float;
{ inner product (standard skalar product) of two vectors }

FUNCTION TotalSum(CONST Vec: VectorTyp): float;
{ TotalSum := a[1] + a[2] + ... + a[n] }

FUNCTION NeumaierSum(CONST A: VectorTyp): float;
{ same as TotalSum, but reduces accumulation error }

FUNCTION TotalProduct(CONST Vec: VectorTyp): float;
{ TotalProduct := a[1] * a[2] * ... * a[n] }

FUNCTION ActualElements(CONST A: VectorTyp): WORD;
{ count number of non-NaN elements }

FUNCTION VectorSignum(Vec: VectorTyp): INTEGER;
{ returns 1 if all elements > 0 (vector strictly positive),
          0 if at least one element is zero (vector positive)
         -1 if at least one element is negative }

FUNCTION FindLargest(Vec: VectorTyp): float;
{ returns largest element of vector }

FUNCTION FindSmallest(Vec: VectorTyp): float;
{ returns smallest element of vector }

FUNCTION FindLargestAbsolute(Vec: VectorTyp): float;
{ returns absolute largest element of vector }

FUNCTION FindSmallestAbsolute(Vec: VectorTyp): float;
{ returns absolute smallest element of vector }

PROCEDURE Scale(VAR Vec: VectorTyp; min, max: float);
{ calculate (x - min) / (max - min) to scale the data. If min and max are the
  actual minimum and maximum of data vector, the vector is scaled to [0..1].
  If data have outlayers, use, say, Q95 and Q5 instead. }

PROCEDURE Centre(VAR Vec: VectorTyp);
{ subtract arithmetic mean from all data }

FUNCTION VectorAngle(CONST A, B: VectorTyp): float;
{ angle between two vectors }

{ ------------ Vector norms, normalisation  ------------------- }

FUNCTION VectorEuklidianNorm(CONST Vec: VectorTyp): float;
{ Euklidian length of vector }

FUNCTION VectorMaximumNorm(CONST Vec: VectorTyp): float;
{ maximum absolute value }

FUNCTION VectorAbsoluteSumNorm(CONST Vec: VectorTyp): float;
{ Sum of the absolute values }

PROCEDURE NormaliseVector(VAR Vec: VectorTyp; Norm: float);
{ Divide all elements of the vector by its norm }

{ ---------------- Elemente sortieren ---------------- }

PROCEDURE ShellSort(VAR t: VectorTyp);
{ sorts data in vector in increasing order. NaN will be sorted to
  the highest positions. Modified from
  http://jean-pierre.moreau.pagesperso-orange.fr/Pascal/sort2_pas.txt }

{ ------------------ Distance ------------------------ }

FUNCTION SquaredEuklidianDistance(CONST A, B: VectorTyp;
  IgnoreFirst: BOOLEAN): float;
{ square of Euklidian distance between two vectors. If IgnoreFirst is true, than
  case numbers in the first column are ignored. }


{ ---------------------------------------------------- }

IMPLEMENTATION

VAR
  CH: CHAR;

PROCEDURE CreateVector(VAR Vec: VectorTyp; Length: WORD; Value: float);

VAR
  i: WORD;

BEGIN
  TRY
    GetMem(Vec, Length * SizeOf(float) + SizeOf(WORD) + 6);
  EXCEPT
    CH := WriteErrorMessage(' Not enough memory to create vector');
    VectorError := true;
    EXIT;
  END;
  Vec^.Length := Length;
  FOR i := 1 TO Length DO
    Vec^.Data[i] := Value;
END;

PROCEDURE DestroyVector(VAR Vec: VectorTyp);

VAR
  x: WORD;

BEGIN
  x := Vec^.Length * SizeOf(float) + SizeOf(WORD) + 6;
  FreeMem(Vec, x);
END;


PROCEDURE ReadVector(VAR Vec: VectorTyp; MedStr: STRING);

VAR
  j, l: WORD;
  Medium: TEXT;

BEGIN
  Assign(Medium, MedStr);
  Reset(Medium);
  IF IOResult <> 0 THEN
  BEGIN
    CH := WriteErrorMessage(' File with vector not found');
    VectorError := true;
    EXIT;
  END;
  ReadLn(Medium, l);
  FOR j := 1 TO l DO
  BEGIN
    IF EoF(Medium) THEN
    BEGIN
      CH := WriteErrorMessage(' Unknown file format');
      VectorError := true;
      EXIT;
    END;
    ReadLn(Medium, Vec^.Data[j]);
    IF IOResult <> 0 THEN
    BEGIN
      CH := WriteErrorMessage(' Unknown file format');
      VectorError := true;
      EXIT;
    END;
  END;
  Close(Medium);
END;


PROCEDURE WriteVector(MedStr: STRING; CONST Vec: VectorTyp; ValidFigures: BYTE);

VAR
  j: WORD;
  Medium: TEXT;

BEGIN
  Assign(Medium, MedStr);
  Rewrite(Medium);
  Writeln(Medium, Vec^.Length);
  FOR j := 1 TO Vec^.Length DO
    Write(Medium, FloatStr(Vec^.Data[j], ValidFigures), ' ');
  Close(Medium);
END;


FUNCTION VectorLength(CONST Vec: VectorTyp): WORD;

BEGIN
  Result := Vec^.Length;
END;


FUNCTION GetVectorElement(CONST Vec: VectorTyp; n: WORD): float;

VAR
  s1, s2: STRING;

BEGIN
  IF ((n <= VectorLength(Vec)) AND (n > 0))
    THEN
      Result := Vec^.Data[n]
    ELSE
      BEGIN
        Str(n: 4, s1);
        Str(VectorLength(Vec): 4, s2);
        CH := WriteErrorMessage(' Attempt to read non-existend vector element No ' +
               s1 + ' of ' + s2);
      END;
END;


PROCEDURE SetVectorElement(VAR Vec: VectorTyp; n: WORD; C: float);

BEGIN
  IF ((n <= VectorLength(Vec)) AND (n > 0)) THEN
    Vec^.Data[n] := C
  ELSE
    BEGIN
      CH := WriteErrorMessage(' Attempt to write to non-existend vector element');
      VectorError := true;
    END;
END;


PROCEDURE CopyVector(CONST Source: VectorTyp; VAR Dest: VectorTyp);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(Source);
  TRY
    GetMem(Dest, n * SizeOf(float) + SizeOf(WORD) + 6);
  EXCEPT
    CH := WriteErrorMessage(' Copy vector: Not enough memory to create copy');
    VectorError := true;
    EXIT;
  END;
  CreateVector(Dest, n, 0.0);
  FOR i := 1 TO n DO
    SetVectorElement(Dest, i, GetVectorElement(Source, i));
END;


PROCEDURE LoadConstant(VAR Vec: VectorTyp; C: float);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(Vec);
  IF (n = 0) THEN EXIT;
  FOR i := 1 TO n DO
    SetVectorElement(Vec, i, C);
END;


PROCEDURE AddConstant(VAR A: VectorTyp; C: float);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(A);
  IF (n = 0) THEN EXIT;
  FOR i := 1 TO n DO
    SetVectorElement(A, i, GetVectorElement(A, i) + C);
END;


PROCEDURE SubConstant(VAR A: VectorTyp; C: float);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(A);
  IF (n = 0) THEN EXIT;
  FOR i := 1 TO n DO
    SetVectorElement(A, i, GetVectorElement(A, i) - C);
END;



PROCEDURE MulConstant(VAR A: VectorTyp; C: float);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(A);
  IF (n = 0) THEN EXIT;
  FOR i := 1 TO n DO
    SetVectorElement(A, i, GetVectorElement(A, i) * C);

END;


PROCEDURE DivConstant(VAR A: VectorTyp; C: float);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(A);
  IF (n = 0) THEN EXIT;
  FOR i := 1 TO n DO
    SetVectorElement(A, i, GetVectorElement(A, i) / C);
END;

PROCEDURE VectorAdd(CONST A, B: VectorTyp; VAR C: VectorTyp);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(A);
  IF (n = 0) OR (n <> VectorLength(B))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Addition of vectors: vectors of different lengths');
        VectorError := true;
        EXIT;
      END;
  CreateVector(C, n, 0.0);
  FOR i := 1 TO n DO
    SetVectorElement(C, i, GetVectorElement(A, i) + GetVectorElement(B, i));
END;


PROCEDURE VectorSub(CONST A, B: VectorTyp; VAR C: VectorTyp);

VAR
  i, n: WORD;

BEGIN
  n := VectorLength(A);
  IF (n = 0) OR (n <> VectorLength(B))
    THEN
      BEGIN
        CH := WriteErrorMessage(' Subtraction of vectors: vectors of different lengths');
        VectorError := true;
        EXIT;
      END;
  CreateVector(C, n, 0.0);
  FOR i := 1 TO n DO
    SetVectorElement(C, i, GetVectorElement(A, i) - GetVectorElement(B, i));
END;


FUNCTION TotalSum(CONST Vec: VectorTyp): float;

VAR
  i, n: WORD;
  x: float;

BEGIN
  n := VectorLength(Vec);
  IF (n = 0)
    THEN
      BEGIN
        TotalSum := 0;
        EXIT;
      END;
  x := 0;
  FOR i := 1 TO n DO
    IF IsNaN(GetVectorElement(Vec, i))
      THEN
      ELSE x := x + GetVectorElement(Vec, i);
  Result := x;
END;


FUNCTION NeumaierSum(CONST A: VectorTyp): float;

VAR
  Sum, c, t, x: float;
  i, j: WORD;

BEGIN
  IF (VectorLength(A) = 0)
    THEN
      BEGIN
        NeumaierSum := 0;
        EXIT;
      END;
  i := 1;
  REPEAT     // search FOR the first non-NaN element
    Sum := GetVectorElement(A, i);
    INC(i)
  UNTIL NOT (IsNaN(Sum));
  c := 0.0;
  FOR j := i TO VectorLength(A) DO
    IF IsNaN(GetVectorElement(A, j))
      THEN
      ELSE
        BEGIN
          x := GetVectorElement(A, j);
          t := Sum + x;
          IF Abs(Sum) >= Abs(x)
            THEN c := c + ((Sum - t) + x)   // LOW-order digits OF A(j) are lost
            ELSE c := c + ((x - t) + Sum);  // LOW-order digits OF Sum are lost
          Sum := t;
        END;
  Result := Sum + c;
END;


FUNCTION TotalProduct(CONST Vec: VectorTyp): float;

VAR
  i, n: WORD;
  x: float;

BEGIN
  n := VectorLength(Vec);
  IF (n = 0)
    THEN
      BEGIN
        TotalProduct := 0;
        EXIT;
      END;
  x := 0;
  FOR i := 1 TO n DO
    IF IsNaN(Vec^.Data[i])
      THEN
      ELSE x := x + Ln(GetVectorElement(Vec, i));
  TotalProduct := Exp(x);
END;

FUNCTION VectorSignum(Vec: VectorTyp): INTEGER;

VAR
  x: float;
  i: WORD;

BEGIN
  VectorSignum := 1;
  FOR i := 1 TO VectorLength(Vec) DO
    BEGIN
      x := GetVectorElement(Vec, i);
      IF IsNaN(x)
        THEN
        ELSE
          CASE signum(x) OF
            -1: BEGIN
                  VectorSignum := -1; // one negative element IS enough
                  EXIT;
                END;
            0: BEGIN
                 Result := 0; // continue IN CASE a negative element comes later
               END;
            1: BEGIN
                // do nothing
               END;
          END; { case }
    END; { for }
END;


FUNCTION VectorEuklidianNorm(CONST Vec: VectorTyp): float;

VAR
  i, n: WORD;
  x: VectorTyp;

BEGIN
  n := VectorLength(Vec);
  IF (n = 0)
    THEN
      BEGIN
        Result := 0;
        EXIT;
      END;
  CreateVector(x, n, 0.0);
  IF VectorError THEN EXIT;
  FOR i := 1 TO Vec^.Length DO
    SetVectorElement(x, i, Sqr(GetVectorElement(Vec, i)));
  VectorEuklidianNorm := Sqrt(NeumaierSum(x)); // avoid rounding error during summation
  DestroyVector(x);
END;


FUNCTION VectorMaximumNorm(CONST Vec: VectorTyp): float;

VAR
  Max, Element: float;
  i: WORD;

BEGIN
  IF (VectorLength(Vec) = 0)
    THEN
      BEGIN
        Result := 0;
        EXIT;
      END;
  Max := -1e300;
  FOR i := 1 TO VectorLength(Vec) DO
    BEGIN
      Element := Abs(GetVectorElement(Vec, i));
      IF Element > Max THEN Max := Element;
    END;
  Result := Max;
END;


FUNCTION VectorAbsoluteSumNorm(CONST Vec: VectorTyp): float;

VAR
  i: WORD;
  Sum: float;

BEGIN
  IF (VectorLength(Vec) = 0)
    THEN
      BEGIN
        Result := 0;
        EXIT;
      END;
  Sum := 0;
  FOR i := 1 TO VectorLength(Vec) DO
    Sum := Sum + Abs(GetVectorElement(Vec, i));
  Result := Sum;
END;


PROCEDURE NormaliseVector(VAR Vec: VectorTyp; Norm: float);

BEGIN
  DivConstant(Vec, Norm);
END;


FUNCTION VectorInnerProduct(CONST A, B: VectorTyp): float;

VAR
  i, n1, n2: WORD;
  x: VectorTyp;

BEGIN
  n1 := VectorLength(A);
  n2 := VectorLength(B);
  IF (n1 <> n2)
    THEN
      BEGIN
        CH := WriteErrorMessage('Vector error: Inner Product of vectors of unequal length');
        VectorError := true;
        EXIT;
      END;
  IF (n1 = 0)
    THEN
      BEGIN
        Result := 0;
        EXIT;
      END;
  CreateVector(x, n1, 0.0);
  FOR i := 1 TO n1 DO
    SetVectorElement(x, i, GetVectorElement(A, i) * GetVectorElement(B, i));
  Result := NeumaierSum(x);
  DestroyVector(x);
END;


FUNCTION FindLargestAbsolute(Vec: VectorTyp): float;

VAR
  i: WORD;
  a: float;

BEGIN
  a := GetVectorElement(Vec, 1);
  FOR i := 2 TO VectorLength(Vec) DO
    IF Abs(a) < Abs(GetVectorElement(Vec, i))
      THEN a := GetVectorElement(Vec, i);
  Result := a;
END;


FUNCTION FindSmallestAbsolute(Vec: VectorTyp): float;

VAR
  i: WORD;
  a: float;

BEGIN
  a := GetVectorElement(Vec, 1);
  FOR i := 2 TO VectorLength(Vec) DO
    IF Abs(a) > Abs(GetVectorElement(Vec, i))
      THEN a := GetVectorElement(Vec, i);
  Result := a;
END;


FUNCTION FindLargest(Vec: VectorTyp): float;

VAR
  i: WORD;
  a: float;

BEGIN
  a := GetVectorElement(Vec, 1);
  FOR i := 2 TO VectorLength(Vec) DO
    IF a < GetVectorElement(Vec, i)
      THEN a := GetVectorElement(Vec, i);
  Result := a;
END;


FUNCTION FindSmallest(Vec: VectorTyp): float;

VAR
  i: WORD;
  a: float;

BEGIN
  a := GetVectorElement(Vec, 1);
  FOR i := 2 TO VectorLength(Vec) DO
    IF a > GetVectorElement(Vec, i)
      THEN a := GetVectorElement(Vec, i);
  Result := a;
END;


FUNCTION VectorAngle(CONST A, B: VectorTyp): float;

VAR
  n1, n2: WORD;

BEGIN
  n1 := VectorLength(A);
  n2 := VectorLength(B);
  IF (n1 <> n2)
    THEN
      BEGIN
        CH := WriteErrorMessage('Vector-Error: Angle between vectors of unequal dimension');
        VectorError := TRUE;
        EXIT;
      END;
  IF (n1 = 0)
    THEN Result := 0
    ELSE Result := ArcCos(VectorInnerProduct(A, B) /
      (VectorEuklidianNorm(A) * VectorEuklidianNorm(B)));
END;


FUNCTION ActualElements(CONST A: VectorTyp): WORD;

VAR
  i, j: WORD;

BEGIN
  j := 0;
  FOR i := 1 TO VectorLength(A) DO
    IF IsNaN(GetVectorElement(A, i))
      THEN
      ELSE INC(j);
  Result := j;
END;


PROCEDURE Scale(VAR Vec: VectorTyp; min, max: float);

VAR
  i: WORD;
  Range: float;

BEGIN
  Range := max - min;
  FOR i := 1 TO VectorLength(Vec) DO
    SetVectorElement(Vec, i, (GetVectorElement(Vec, i) - min) / Range);
END;


PROCEDURE Centre(VAR Vec: VectorTyp);

VAR
  Mean: float;
  i: WORD;

BEGIN
  Mean := NeumaierSum(Vec) / ActualElements(Vec);
  FOR i := 1 TO VectorLength(Vec) DO
    SetVectorElement(Vec, i, GetVectorElement(Vec, i) - Mean);
END;


FUNCTION SquaredEuklidianDistance(CONST A, B: VectorTyp;
  IgnoreFirst: BOOLEAN): float;

VAR
  p, i, start: WORD;
  c: CHAR;
  ai, bi, Sum: float;

BEGIN
  p := VectorLength(A);
  IF VectorLength(B) <> p
    THEN
      BEGIN
        VectorError := TRUE;
        c := WriteErrorMessage('Euklidian distance of two vectors: vectors have unequal length');
        EXIT;
      END;
  Sum := 0;
  IF IgnoreFirst
    THEN Start := 2   // first column study number
    ELSE Start := 1;  // first column data
  FOR i := Start TO p DO
    BEGIN
      ai := GetVectorElement(A, i);
      bi := GetVectorElement(B, i);
      IF IsNaN(ai) OR IsNaN(Bi)
        THEN
          BEGIN
            VectorError := TRUE;
            c := WriteErrorMessage(
              'Euklidian distance of two vectors: vectors contain NaN');
            EXIT;
          END;
      Sum := Sum + Sqr(ai - bi);
    END;
  Result := Sum;
END;


PROCEDURE ShellSort(VAR t: VectorTyp);

LABEL
  10;

VAR
  i, j, k, l, m, nn, NaNs, n, LdN: INTEGER;
  tmp, s: float;

BEGIN
  NaNs := 0;
  s := -1e300;
  n := VectorLength(t);
  FOR i := 1 TO n DO
    BEGIN
      IF IsNaN(GetVectorElement(t, i))    // check FOR NaN data
        THEN INC(NaNs)
        ELSE IF (GetVectorElement(t, i)) > s
               THEN s := GetVectorElement(t, i); // AND find largest element OF data vector
    END;
  s := 10 * s;
  IF (NaNs > 0)
    THEN
      FOR i := 1 TO n DO
        IF IsNaN(GetVectorElement(t, i))
          THEN SetVectorElement(t, i, s);
  // replace all NaN WITH very large number so they Move TO END OF vector
  LdN := Trunc(Ln(n) / Const_ln2);
  m := n;
  FOR nn := 1 TO LdN DO
    BEGIN
      m := m DIV 2;
      k := n - m;
      FOR j := 1 TO k DO
        BEGIN
          i := j;
10:       l := i + m;
          IF (GetVectorElement(t, l) < GetVectorElement(t, i))
            THEN
              BEGIN
                tmp := GetVectorElement(t, i);
                SetVectorElement(t, i, GetVectorElement(t, l));
                SetVectorElement(t, l, tmp);
                i := i - m;
                IF i >= 1 THEN GOTO 10;
              END;
        END;
    END;
  IF (NaNs > 0) THEN
    FOR i := Succ(n - NaNs) TO n DO  // change the top NaNs elements back TO NaN
      SetVectorElement(t, i, NaN);
END;


END.

