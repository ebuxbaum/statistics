PROGRAM TestComplex;

USES Mathfunc, Complex;

VAR a, b, c, d : ComplexTyp;
    x, y, z    : double;


PROCEDURE EinOperand(a, r : ComplexTyp; Operation : STRING);

BEGIN
  Write(Operation,'(');
  Write(ComplexToStr(a, 6, 3));
  Write(') = ');
  Write(ComplexToStr(r, 6, 3));
  Writeln;
END;


PROCEDURE ZweiOperanden (a1, a2, r : ComplexTyp; Operation : STRING);

BEGIN
  Write('(');
  Write(ComplexToStr(a1, 6, 3));
  Write(') ', Operation, ' (');
  Write(ComplexToStr(a2, 6, 3));
  Write(') = (');
  Write(ComplexToStr(r, 6, 3));
  Writeln(')');
END;


PROCEDURE ComplexMitReel (a1, r : ComplexTyp; a2 : double; Operation : STRING);

BEGIN
  Write('(');
  Write(ComplexToStr(a1, 6, 3));
  Write(') ', Operation, ' ', a2:6:3);
  Write(' = (');
  Write(ComplexToStr(r, 6, 3));
  Writeln(')');
END;


PROCEDURE TestAddSub;

BEGIN
  c := a + b;
  ZweiOperanden(a, b, c, '+');
  d := c - b;
  ZweiOperanden(c, b, d, '-');
END;


PROCEDURE TestMulDiv;

BEGIN
  c := a * b;
  ZweiOperanden(a, b, c, '*');
  d := c / b;
  IF ComplexError
    THEN ComplexError := FALSE
    ELSE ZweiOperanden(c, b, d, '/');
END;


PROCEDURE TestMulDivMitReel;

BEGIN
  b := a * x;
  ComplexMitReel(a, b, x, '*');
  c := b / x;
  IF ComplexError
    THEN ComplexError := FALSE
    ELSE ComplexMitReel(b, c, x, '/');
END;


PROCEDURE TestPolarRect;

BEGIN
  Polar(a, y, z);
  c := Rect(y, z);
  Write('Polar(');
  Write(ComplexToStr(a, 6, 3));
  Writeln(') = (', y:6:3, ', ', z:6:3, ')');
  Write('Rect(', y:6:3, ', ', z:6:3, ') = (');
  Write(ComplexToStr(c, 6, 3));
  Writeln(')');
END;


PROCEDURE TestExpLn;

BEGIN
  b := ComplexExp(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'exp');
  c := ComplexLn(b);
  IF ComplexError
    THEN ComplexError := FALSE
    ELSE EinOperand(b, c, 'ln');
END;


PROCEDURE TestPotReelComplex;

BEGIN
  b := ComplexPower(x, a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  Write(x:6:3, '^(');
  Write(ComplexToStr(a, 6, 3));
  Write(') = ');
  Write(ComplexToStr(b, 6, 3));
  Writeln;
  c := ComplexRoot(b, a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  Write('(');
  Write(ComplexToStr(b, 6, 3));
  Write(')^1/(');
  Write(ComplexToStr(a, 6, 3));
  Write(') = ');
  Write(ComplexToStr(c, 6, 3));
  Writeln;
END;


PROCEDURE TestPotComplexReel;

BEGIN
  b := ComplexPower(a, x);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  ComplexMitReel(a, b, x, '^');
  c := ComplexRoot(b, x);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  ComplexMitReel(b, c, x, '^1/');
END;


PROCEDURE TestPotComplexComplex;

BEGIN
  c := ComplexPower(a, b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  ZweiOperanden(a, b, c, '^');
  d := ComplexRoot(c, b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  ZweiOperanden(c, b, d, '^1/');
END;


PROCEDURE TestSin;

BEGIN
  b := ComplexSin(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'sin');
  c := ComplexArcSin(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arcsin');
END;


PROCEDURE TestCos;

BEGIN
  b := ComplexCos(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'cos');
  c := ComplexArcCos(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arccos');
END;


PROCEDURE TestTan;

BEGIN
  b := ComplexTan(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'tan');
  c := ComplexArcTan(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arctan');
END;


PROCEDURE TestCot;

BEGIN
  b := ComplexCot(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'cot');
  c := ComplexArcCot(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arctan');
END;


PROCEDURE TestSinh;

BEGIN
  b := ComplexSinh(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'sinh');
  c := ComplexArSinh(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arsinh');
END;


PROCEDURE TestCosh;

BEGIN
  b := ComplexCosh(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'cosh');
  c := ComplexArCosh(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arcosh');
END;


PROCEDURE TestTanh;

BEGIN
  b := ComplexTanh(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'tanh');
  c := ComplexArTanh(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
       END;
  EinOperand(b, c, 'artanh');
END;


PROCEDURE TestCoth;

BEGIN
  b := ComplexCoth(a);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(a, b, 'coth');
  c := ComplexArCoth(b);
  IF ComplexError
    THEN
      BEGIN
        ComplexError := FALSE;
        EXIT;
      END;
  EinOperand(b, c, 'arcoth');
END;


PROCEDURE TestAll;

VAR c : ComplexTyp;

BEGIN
  TestAddSub; Writeln;
  TestMulDiv; Writeln;
  TestMulDivMitReel; Writeln;
  TestPolarRect; Writeln;
  TestExpLn; Writeln;
  TestPotReelComplex; Writeln;
  TestPotComplexReel; Writeln;
  TestPotComplexComplex; Writeln;
  TestSin; Writeln;
  TestCos; Writeln;
  TestTan; Writeln;
  TestCot; Writeln;
  TestSinh; Writeln;
  TestCosh; Writeln;
  TestTanh; Writeln;
  TestCoth; Writeln;
  ReadLn;
END;


BEGIN
  a := ComplexInit(0.621, 0.567);
  b := ComplexInit(0.5, 0.4);
  x := 1.5;
  TestAll;
  a := ComplexInit(1.0, -1.0);
  x := 0.2;
  TestAll;
  a := ComplexInit(-0.5, 0.2);
  x := 0.2;
  TestAll;
  a := ComplexInit(-1.0, -1.5);
  x := 1.5;
  TestAll;
  a := ComplexInit(1.0, 0.0);
  x := 1.5;
  TestAll;
  a := ComplexInit(0.0, 1.0);
  x := 1.5;
  TestAll;
  a := ComplexInit(0.0, 0.0);
  x := 1.5;
  TestAll;
END.

