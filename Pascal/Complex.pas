UNIT Complex;

INTERFACE

USES Math, mathfunc;

TYPE
  ComplexTyp = RECORD
    RealPart,
    ImagPart: float;
  END;

CONST Const_i    : ComplexTyp = (RealPart : 0.0; ImagPart : 1.0);
      Const_Null : ComplexTyp = (RealPart : 0.0; ImagPart : 0.0);

VAR
  ComplexError: BOOLEAN = FALSE;     // error-Flag

{ ************************** Type conversion ******************************** }

FUNCTION ComplexInit(RePart, ImPart: float): ComplexTyp;

FUNCTION Re(z: ComplexTyp): float;
{ real part of a ComplexTyp number }

FUNCTION Im(z: ComplexTyp): float;
{ imaginary part of a ComplexTyp number }

{ *************************** E/A- Routines ***************************** }

FUNCTION ComplexToStr(z: ComplexTyp; m, n: BYTE): STRING;

{ ***************************** Conversions *************************** }

FUNCTION ComplexIsNull(z: ComplexTyp): BOOLEAN;

FUNCTION ComplexConj(z: ComplexTyp): ComplexTyp;

FUNCTION ComplexAbs(z: ComplexTyp): float;

FUNCTION ComplexArg(z: ComplexTyp): float;

FUNCTION ComplexInv(z: ComplexTyp): ComplexTyp;

FUNCTION ComplexNeg(z: ComplexTyp): ComplexTyp;

PROCEDURE Polar(z: ComplexTyp; VAR r, phi: float);

FUNCTION Rect(r, phi: float): ComplexTyp;

{ ***************** Calculation with ComplexTyp ************** }

 OPERATOR = (z, b: ComplexTyp) : BOOLEAN;

 OPERATOR + (z, b: ComplexTyp): ComplexTyp;

 OPERATOR - (z, b: ComplexTyp): ComplexTyp;

 OPERATOR - (z: ComplexTyp): ComplexTyp;   // unary -

 OPERATOR * (z, b: ComplexTyp): ComplexTyp;

 OPERATOR / (z, b: ComplexTyp): ComplexTyp;

 OPERATOR ** (Basis, Exponent: ComplexTyp): ComplexTyp;

 { *********** Calculation with ComplexTyp and Real *********** }

 OPERATOR = (z: ComplexTyp; x: float) : BOOLEAN;

 OPERATOR = (x: float; z: ComplexTyp) : BOOLEAN;

 OPERATOR := (x: real) : ComplexTyp;

 OPERATOR := (z: ComplexTyp) : real;

 OPERATOR + (z: ComplexTyp; x: float) : ComplexTyp;

 OPERATOR + (x: float; z: ComplexTyp) : ComplexTyp;

 OPERATOR - (z : ComplexTyp; r : real) : ComplexTyp;

 OPERATOR - (r : real; z : ComplexTyp) : ComplexTyp;

 OPERATOR * (z: ComplexTyp; x: float): ComplexTyp;

 OPERATOR * (x: float; z: ComplexTyp) : ComplexTyp;

 OPERATOR / (z: ComplexTyp; x: float): ComplexTyp;

 OPERATOR / (x: float; z: ComplexTyp) : ComplexTyp;

 OPERATOR ** (Basis : ComplexTyp; Exponent : float) : ComplexTyp;

 OPERATOR ** (Basis : float; Exponent : ComplexTyp) : ComplexTyp;

{ ********************** Powers, roots ****************************** }

FUNCTION ComplexLn(z: ComplexTyp): ComplexTyp;
{ ln(z) = ln(abs(z)) + i * arg(z) }

FUNCTION ComplexExp(z: ComplexTyp): ComplexTyp;
{ exp(x + iy) = exp(x)*cos(y) + i*exp(x)*sin(y) }

FUNCTION ComplexPower(Basis, Exponent: ComplexTyp): ComplexTyp;
{ b^e = exp(e * ln(b)) }

FUNCTION ComplexPower(Basis: ComplexTyp; Exponent: float): ComplexTyp;
{ z^x = exp(x * ln(z)) }

FUNCTION ComplexPower(Basis: float; Exponent: ComplexTyp): ComplexTyp;
{ x^z = exp(z * ln(x)) }

FUNCTION ComplexRoot(z, b: ComplexTyp): ComplexTyp;
{ b-te Wurzel aus z = z^(1/b) }

FUNCTION ComplexRoot(z: ComplexTyp; x: float): ComplexTyp;
{ x-te Wurzel aus z = z^(1/x) }

FUNCTION ComplexRoot(x: float; z: ComplexTyp): ComplexTyp;
{ z-te Wurzel aus x = x^(1/z) }

{ *********************** Cyclometric functions ************************** }

FUNCTION ComplexSin(z: ComplexTyp): ComplexTyp;
{ sin(x +- iy) = sin(x) * cosh(y) +- i * cos(x) * sinh(y) }

FUNCTION ComplexCos(z: ComplexTyp): ComplexTyp;
{ cos(x +- iy) = cos(x) * cosh(y) +- i * sin(x) * sinh(y) }

FUNCTION ComplexTan(z: ComplexTyp): ComplexTyp;
{ tan(x +- iy) = (sin(2x) / (cos(2x) + cosh(2x)) +- i * (sinh(2x) / (cos(2x) + cosh(2x)) }

FUNCTION ComplexCot(z: ComplexTyp): ComplexTyp;
{ cot(z) = cos(z) / sin(z) }

FUNCTION ComplexArcSin(z: ComplexTyp): ComplexTyp;
{ arcsin(z) = 1/i * ln(i*z + sqrt(1 - sqr(z))) }

FUNCTION ComplexArcCos(z: ComplexTyp): ComplexTyp;
{ arccos(z) = 1/i * ln(z + sqrt(sqr(z) - 1)) }

FUNCTION ComplexArcTan(z: ComplexTyp): ComplexTyp;
{ arctan(z) = 1/2i * ln((1+iz) / (1-iz)) }

FUNCTION ComplexArcCot(z: ComplexTyp): ComplexTyp;
{ arccot(z) = -1/2i * ln((1+iz) / (iz - 1)) }

{ *********************** Hyperbolic funktions ********************** }

FUNCTION ComplexSinh(z: ComplexTyp): ComplexTyp;
{ sinh(x + iy) = sinh(x)*cos(y) +- i*cosh(x)*sin(y)
               = (exp(z) - exp(-z)) / 2 }

FUNCTION ComplexCosh(z: ComplexTyp): ComplexTyp;
{ cosh(x + iy) = cosh(x)*cos(y) +- i*sinh(x)*sin(y)
               = (exp(z) + exp(-z)) / 2 }

FUNCTION ComplexTanh(z: ComplexTyp): ComplexTyp;
{ tanh(x + iy) = (sinh(2x)/(cosh(2x)+cos(2y)) + i*(sin(2y)/(cosh(2x)+cos(2y)))
               = (exp(z) - exp(-z)) / (exp(z) + exp(-z)) }

FUNCTION ComplexCoth(z: ComplexTyp): ComplexTyp;
{ coth(z) = cosh(z) / sinh(z) }

FUNCTION ComplexArSinh(z: ComplexTyp): ComplexTyp;
{ arsinh(z) = ln(z + sqrt(sqr(z) + 1)) }

FUNCTION ComplexArCosh(z: ComplexTyp): ComplexTyp;
{ arcosh(z) = ln(z + sqrt(sqr(z) - 1)) }

FUNCTION ComplexArTanh(z: ComplexTyp): ComplexTyp;
{ artanh(z) := -1/2 * ln((1 + z) / (1 - z)) }

FUNCTION ComplexArCoth(z: ComplexTyp): ComplexTyp;
{ arcoth(z) := 1/2 * ln((1 + z) / (z - 1)) }

{ *************************** Implementation ****************************** }

IMPLEMENTATION

VAR ch : CHAR;  // used FOR error handling

{ **************************** Typ-Wandlungen ****************************** }

FUNCTION Re(z: ComplexTyp): float;

BEGIN
  Result := z.RealPart;
END;


FUNCTION Im(z: ComplexTyp): float;

BEGIN
  Result := z.ImagPart;
END;


FUNCTION ComplexIsNull(z: ComplexTyp): BOOLEAN;

BEGIN
  Result := (Abs(Re(z)) < Zero) AND (Abs(Im(z)) < Zero);
END;


FUNCTION ComplexInit(RePart, ImPart: float): ComplexTyp;

BEGIN
  Result.RealPart := RePart;
  Result.ImagPart := ImPart;
END;


{ ***************************** I/O-Routinen ******************************* }

FUNCTION ComplexToStr(z: ComplexTyp; m, n: BYTE): STRING;
  {Komplexe Zahl ausgeben}

VAR
  OutStr, ImagStr: STRING;
  k: INTEGER;

BEGIN
  Str(Re(z): m: n, OutStr);
  WHILE OutStr[1] = ' ' DO
    Delete(OutStr, 1, 1);
  IF Im(z) >= 0
    THEN OutStr := OutStr + '+'
    ELSE OutStr := OutStr + '-';
  Str(Abs(Im(z)): m: n, ImagStr);
  WHILE ImagStr[1] = ' ' DO
    Delete(ImagStr, 1, 1);
  OutStr := OutStr + ImagStr + 'i';
  FOR k := 1 TO (m - Length(OutStr)) DO
    OutStr := ' ' + OutStr;
  Result := OutStr;
END;

{ ************************** Hilfsfunktionen ******************************* }

FUNCTION ComplexConj(z: ComplexTyp): ComplexTyp;

BEGIN
  Result := ComplexInit(Re(z), -Im(z));
END;


FUNCTION ComplexAbs(z: ComplexTyp): float;

VAR
  x: float;

BEGIN
  x := (Sqr(Re(z)) + Sqr(Im(z)));
  IF Abs(x) < Zero
    THEN Result := 0
    ELSE Result := Sqrt(x);
END;


FUNCTION ComplexArg(z: ComplexTyp): float;

VAR
  i, r: float;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage(' Complex numbers: Argument of (0 + 0i)');
        ComplexError := TRUE;
        EXIT;
      END;
  i := Im(z);
  r := Re(z);
  IF r = 0
    THEN Result := Signum(i) * 0.5 * Pi
    ELSE
      IF r > 0
        THEN Result := ArcTan(i / r)
        ELSE Result := ArcTan(i / r) + Signum(i) * Pi;
END;


FUNCTION ComplexInv (z : ComplexTyp) : ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex Numbers: Attempt to invert (0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE Result := ComplexConj(z) / Sqr(ComplexAbs(z));
END;


FUNCTION ComplexNeg(z: ComplexTyp): ComplexTyp;

BEGIN
  Result := ComplexInit(-Re(z), -Im(z));
END;


PROCEDURE Polar(z: ComplexTyp; VAR r, phi: float);

BEGIN
  r := ComplexAbs(z);
  IF r <> 0
    THEN IF Abs(Im(z) / r) = 1
           THEN phi := Pi / 2 * (Im(z) / r) / Abs(Im(z) / r)
           ELSE phi := ArcTan(Im(z) / Re(z))
    ELSE
      BEGIN
        phi := 0;
        IF Re(z) < 0
          THEN
            IF Im(z) <> 0
              THEN phi := phi + Pi * Abs(Im(z)) / Im(z)
              ELSE phi := Pi;
      END;
END;


FUNCTION Rect(r, phi: float): ComplexTyp;

BEGIN
  Result := ComplexInit(r * Cos(phi), r * Sin(phi));
END;

{ *********************** Operators for complex ************************ }

OPERATOR = (z, b: ComplexTyp) : BOOLEAN;

BEGIN
  Result := (Re(z) = Re(b)) AND (Im(z) = Im(b));
END;


OPERATOR + (z, b: ComplexTyp): ComplexTyp;

BEGIN
  Result := ComplexInit(Re(z) + Re(b), Im(z) + Im(b));
END;


OPERATOR - (z, b: ComplexTyp): ComplexTyp;

BEGIN
  Result := ComplexInit(Re(z) - Re(b), Im(z) - Im(b));
END;


OPERATOR - (z: ComplexTyp): ComplexTyp;   // unary -

BEGIN
  Result := ComplexInit(-Re(z), -Im(z));
END;


OPERATOR * (z, b: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := Re(z) * Re(b) - Im(z) * Im(b);
  i := Re(z) * Im(b) + Im(z) * Re(b);
  Result := ComplexInit(r, i);
END;


OPERATOR / (z, b: ComplexTyp): ComplexTyp;

VAR r, i : float;

BEGIN
  IF ComplexIsNull(b)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Division by (0 + 0i) ');
        ComplexError := TRUE;
        EXIT;
      END;
  r := (Re(z) * Re(b) + Im(z) * Im(b)) / (Sqr(Re(b)) + Sqr(Im(b)));
  i := (Im(z) * Re(b) - Re(z) * Im(b)) / (Sqr(Re(b)) + Sqr(Im(b)));
  Result := ComplexInit(r, i);
END;


OPERATOR ** (Basis, Exponent: ComplexTyp): ComplexTyp;

BEGIN
  IF ComplexIsNull(Basis)
    THEN Result := ComplexInit(0, 0)
    ELSE Result := ComplexExp(Exponent * ComplexLn(Basis));
END;

{ ***************** Mixed Operators for complex and float ******************* }

OPERATOR = (z: ComplexTyp; x: float) : BOOLEAN;

BEGIN
  Result := (Abs(Im(z)) < Zero) AND (Re(z) = x);
END;

OPERATOR = (x: float; z: ComplexTyp) : BOOLEAN;

BEGIN
  Result := (Abs(Im(z)) < Zero) AND (Re(z) = x);
END;


OPERATOR := (x: real) : ComplexTyp;

BEGIN
  Result := ComplexInit(x, 0);
END;


OPERATOR := (z: ComplexTyp) : real;

BEGIN
  IF Abs(Im(z)) > Zero
    THEN
      BEGIN
        ch := WriteErrorMessage('Assignment of complex to real variable: Im(z) <> 0');
        ComplexError := TRUE;
      END
    ELSE
      Result := Re(z);
END;

OPERATOR + (z: ComplexTyp; x: float) : ComplexTyp;

BEGIN
  Result := ComplexInit(x + Re(z), Im(z));
END;

OPERATOR + (x: float; z: ComplexTyp) : ComplexTyp;

BEGIN
  Result := ComplexInit(x + Re(z), Im(z));
END;


OPERATOR - (z : ComplexTyp; r : real) : ComplexTyp;

BEGIN
  Result := ComplexInit(Re(z) - r, Im(z));
END;


OPERATOR - (r : real; z : ComplexTyp) : ComplexTyp;

BEGIN
  Result := ComplexInit(r - Re(z), -Im(z));
END;


OPERATOR * (z: ComplexTyp; x: float) : ComplexTyp;

BEGIN
  Result := ComplexInit(Re(z) * x, Im(z) * x);
END;


OPERATOR * (x: float; z: ComplexTyp) : ComplexTyp;

BEGIN
  Result := ComplexInit(Re(z) * x, Im(z) * x);
END;


OPERATOR / (z: ComplexTyp; x: float): ComplexTyp;

BEGIN
  IF x = 0
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to divide by zero');
        ComplexError := TRUE;
      END
    ELSE
      Result := ComplexInit(Re(z) / x, Im(z) / x);
END;


OPERATOR / (x: float; z: ComplexTyp) : ComplexTyp;

VAR d : float;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Division by (0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE
      BEGIN
        d := Sqr(Re(z)) + Sqr(Im(z));
        Result := ComplexInit(x*Re(z)/d, -x*Im(z)/d);
      END;
END;

OPERATOR ** (Basis : ComplexTyp; Exponent : float) : ComplexTyp;

BEGIN
  IF ComplexIsNull(Basis)
    THEN Result := ComplexInit(0, 0)
    ELSE Result := ComplexExp(ComplexLn(Basis) * Exponent);
END;


OPERATOR ** (Basis : float; Exponent : ComplexTyp) : ComplexTyp;

BEGIN
  IF Abs(Basis) < Zero
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Complex power of 0');
        ComplexError := TRUE;
      END
    ELSE
      Result := ComplexExp(Exponent * Ln(Basis));
END;

{ ************************************************************************* }

FUNCTION ComplexExp(z: ComplexTyp): ComplexTyp;

VAR
  i, r: float;

BEGIN
  i := Im(z);
  r := Exp(Re(z));
  Result := ComplexInit(r * Cos(i), r * Sin(i));
END;


FUNCTION ComplexLn(z: ComplexTyp): ComplexTyp;

VAR r, i: float;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate ln(0+0i) ');
        ComplexError := TRUE;
      END
    ELSE
      BEGIN
        r := Ln(ComplexAbs(z));
        i := ComplexArg(z);
        Result := ComplexInit(r, i);
      END;
END;


FUNCTION ComplexPower(Basis, Exponent: ComplexTyp): ComplexTyp;

BEGIN
  IF ComplexIsNull(Basis)
    THEN Result := ComplexInit(0, 0)
    ELSE Result := ComplexExp(Exponent * ComplexLn(Basis));
END;


FUNCTION ComplexRoot(z, b: ComplexTyp): ComplexTyp;

BEGIN
  IF ComplexIsNull(b)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate the (0 + 0i)-th root');
        ComplexError := TRUE;
        EXIT;
      END
    ELSE
      Result := ComplexPower(z, ComplexInv(b));
END;


FUNCTION ComplexPower(Basis: ComplexTyp; Exponent: float): ComplexTyp;

BEGIN
  IF ComplexIsNull(Basis)
    THEN Result := ComplexInit(0, 0)
    ELSE Result := ComplexExp(ComplexLn(Basis) * Exponent);
END;


FUNCTION ComplexRoot(z: ComplexTyp; x: float): ComplexTyp;

BEGIN
  IF x = 0
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: 0-th root');
        ComplexError := TRUE;
      END
    ELSE
      IF ComplexIsNull(z)
        THEN Result := ComplexInit(0, 0)
        ELSE Result := ComplexPower(z, 1 / x);
END;


FUNCTION ComplexPower(Basis: float; Exponent: ComplexTyp): ComplexTyp;

BEGIN
  IF Abs(Basis) < Zero
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Complex power of 0');
        ComplexError := TRUE;
      END
    ELSE
      Result := ComplexExp(Exponent * Ln(Basis));
END;


FUNCTION ComplexRoot(x: float; z: ComplexTyp): ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: (0 + 0i)-th root');
        ComplexError := TRUE;
      END
    ELSE IF x = 0
           THEN Result := ComplexInit(0, 0)
           ELSE Result := ComplexPower(x, ComplexInv(z));
END;

{ ***************************** Kreis-Funktionen *************************** }

FUNCTION ComplexSin(z: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := Sin(Re(z)) * cosh(Im(z));
  i := Cos(Re(z)) * sinh(Im(z));
  Result := ComplexInit(r, i);
END;


FUNCTION ComplexCos(z: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := Cos(Re(z)) * cosh(Im(z));
  i := -Sin(Re(z)) * sinh(Im(z));
  Result := ComplexInit(r, i);
END;


FUNCTION ComplexTan(z: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := Sin(2 * Re(z)) / (Cos(2 * Re(z)) + cosh(2 * Im(z)));
  i := sinh(2 * Im(z)) / (Cos(2 * Re(z)) + cosh(2 * Im(z)));
  Result := ComplexInit(r, i);
END;


FUNCTION ComplexCot(z: ComplexTyp): ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate cot(0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE
      Result := ComplexCos(z) / ComplexSin(z);
END;

{ *************************** Arcus-Funktionen ***************************** }

FUNCTION ComplexArcSin(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1, c2, c3: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex Numbers: Attempt to calculate arcsin(0 + 0i) ');
        ComplexError := TRUE;
        EXIT;
      END;
  c1 := ComplexInit(1, 0);
  c2 := ComplexInit(0, 1);
  c3 := ComplexInit(0, -1);
  c0 := ComplexRoot(c1 - ComplexPower(z, 2), 2);
  c0 := ComplexLn(c2 * z + c0);
  Result := c3 * c0;
END;


FUNCTION ComplexArcCos(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1, c2: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate arccos(0 + 0i) ');
        ComplexError := TRUE;
        EXIT;
      END;
  c1 := ComplexInit(1.0, 0.0);
  c2 := ComplexInit(0.0, -1.0);
  c0 := ComplexRoot(ComplexPower(z, 2) - c1, 2);
  Result := ComplexLn(z + c0) * c2;
END;


FUNCTION ComplexArcTan(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1, c2, c3, a1, a2: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate arctan(0 + 0i) ');
        ComplexError := TRUE;
        EXIT;
      END;
  c1 := ComplexInit(0, 1);
  c2 := ComplexInit(1, 0);
  c3 := ComplexInit(0, 2);
  c0 := c1 * z;
  a1 := c2 + c0;
  a2 := c2 - c0;
  c0 := ComplexLn(a1 / a2);
  a1 := c2 / c3;
  Result := a1 * c0;
END;


FUNCTION ComplexArcCot(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1, c2, c3, a1, a2: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate arccot(0 + 0i) ');
        ComplexError := TRUE;
        EXIT;
      END;
  c1 := ComplexInit(0, 1);
  c2 := ComplexInit(1, 0);
  c3 := ComplexInit(0, -2);
  c0 := c1 * z;
  a1 := c0 + c2;
  a2 := c0 - c2;
  c0 := ComplexLn(a1 / a2);
  a1 := c2 / c3;
  Result := a1 * c0;
END;


{ **************************** Hyperbolic funktions ************************* }

FUNCTION ComplexSinh(z: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := sinh(Re(z)) * Cos(Im(z));
  i := cosh(Re(z)) * Sin(Im(z));
  Result := ComplexInit(r, i);
END;


FUNCTION ComplexCosh(z: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := cosh(Re(z)) * Cos(Im(z));
  i := sinh(Re(z)) * Sin(Im(z));
  Result := ComplexInit(r, i);
END;


FUNCTION ComplexTanh(z: ComplexTyp): ComplexTyp;

VAR
  r, i: float;

BEGIN
  r := sinh(2 * Re(z)) / (cosh(2 * Re(z)) + Cos(2 * Im(z)));
  i := Sin(2 * Im(z)) / (cosh(2 * Re(z)) + Cos(2 * Im(z)));
  Result := ComplexInit(r, i);
END;


FUNCTION ComplexCoth (z : ComplexTyp) : ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate coth(0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE
      Result := ComplexCosh(z) / ComplexSinh(z);
END;

{ ***************************** Area-Funktionen **************************** }

FUNCTION ComplexArSinh(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate arsinh(0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE
      BEGIN
        c1 := ComplexInit(1, 0);
        c0 := ComplexRoot(ComplexPower(z, 2) + c1, 2);
        Result := ComplexLn(z + c0);
      END;
END;


FUNCTION ComplexArCosh(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate arcosh(0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE
      BEGIN
        c1 := ComplexInit(1, 0);
        c0 := ComplexRoot(ComplexPower(z, 2) - c1, 2);
        Result := ComplexLn(z + c0);
      END;
END;


FUNCTION ComplexArTanh(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1, c2, a1, a2: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
        ch := WriteErrorMessage('Complex numbers: Attempt to calculate artanh(0 + 0i) ');
        ComplexError := TRUE;
      END
    ELSE
      BEGIN
        c1 := ComplexInit(1, 0);
        c2 := ComplexInit(2, 0);
        a1 := c1 + z;
        a2 := c1 - z;
        c0 := ComplexLn(a1 / a2);
        Result := c0 / c2;
      END;
END;


FUNCTION ComplexArCoth(z: ComplexTyp): ComplexTyp;

VAR
  c0, c1, c2, a1, a2: ComplexTyp;

BEGIN
  IF ComplexIsNull(z)
    THEN
      BEGIN
       ch := WriteErrorMessage('Complex numbers: Attempt to calculate arcoth(0 + 0i) ');
       ComplexError := TRUE;
      END
    ELSE
      BEGIN
        c1 := ComplexInit(1, 0);
        c2 := ComplexInit(2, 0);
        a1 := z + c1;
        a2 := z - c1;
        c0 := ComplexLn(a1 / a2);
        Result := c0 / c2;
      END;
END;


END.
