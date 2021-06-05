UNIT SingularValue;
{ Literatur:
    W.H. Press et. al.: Mumerical recipes in Pascal: The art of scientific computation
    Cambridge (Cambridge University Press) 1989
}

INTERFACE

USES MathFunc, Vector, Matrix;

PROCEDURE svdcmp(VAR a: MatrixTyp; VAR w: VectorTyp; VAR v: MatrixTyp);

PROCEDURE svdSynth(VAR A: MatrixTyp; CONST W: VectorTyp; CONST U, V: MatrixTyp);
{ Berechnet A aus U, W, und V (nicht V^T). }


IMPLEMENTATION

VAR CH: CHAR;

PROCEDURE svdcmp(VAR a: MatrixTyp; VAR w: VectorTyp; VAR v: MatrixTyp);
{ Berechnet die Singular Value Decomposition A = U * W * V^T zur
  m*n Matrix A. Dabei überschreibt U die Matrix A, welche bei Bedarf
  vor Aufruf der Prozedur kopiet werden muß. Die n*n Diagonalmatrix W
  wird zurckgegeben als Vector der Diagonal-Elemente. Wichtig: In
  V steht die n*n Matrix V, nicht V^T. V und W werden in der Routine initialisiert.
  Literatur: Press et al., Numerical Methods in Pascal, Cambridge 1989  }

LABEL 1, 2, 3;

CONST MaxIter = 50;

VAR m, n, nm, mnmin, l, k, j, jj, its, i       : INTEGER;
  z, y, x, scale, s, h, g, f, c, anorm, Number : double;
  rv1                                          : VectorTyp;


  FUNCTION sign(a, b: double): double;

  BEGIN
    IF (b >= 0.0)
      THEN Result := Abs(a)
      ELSE Result := -Abs(a);
  END;

BEGIN
  { Householder reduction to bidiagonal form }
  m := MatrixRows(A);
  n := MatrixColumns(A);
  CreateMatrix(V, n, n, 0.0);
  CreateVector(W, n, 0.0);
  CreateVector(rv1, n, 0.0);
  g := 0.0;
  scale := 0.0;
  anorm := 0.0;
  FOR i := 1 TO n DO
    BEGIN
      l := i + 1;
      SetVectorElement(rv1, i, scale * g);
      g := 0.0;
      s := 0.0;
      scale := 0.0;
      IF (i <= m)
        THEN
          BEGIN
            FOR k := i TO m DO
              scale := scale + Abs(GetMatrixElement(A, k, i));
            IF (scale <> 0.0)
              THEN
                BEGIN
                  FOR k := i TO m DO
                    BEGIN
                      Number := GetMatrixElement(A, k, i) / scale;
                      SetMatrixElement(A, k, i, Number);
                      s := s + Sqr(Number);
                    END;
                  f := GetMatrixElement(A, i, i);
                  g := -sign(Sqrt(s), f);
                  h := f * g - s;
                  SetMatrixElement(A, i, i, f - g);
                  IF (i <> n)
                    THEN
                      FOR j := l TO n DO
                        BEGIN
                          s := 0.0;
                          FOR k := i TO m DO
                          s := s + GetMatrixElement(A, k, i) * GetMatrixElement(A, k, j);
                          f := s / h;
                          FOR k := i TO m DO
                            SetMatrixElement(A, k, j, GetMatrixElement(A, k, j) +
                               f * GetMatrixElement(A, k, i));
                        END;
                  FOR k := i TO m DO
                    SetMatrixElement(A, k, i, scale * GetMatrixElement(A, k, i));
                END;
          END;
      SetVectorElement(W, i, scale * g);
      g := 0.0;
      s := 0.0;
      scale := 0.0;
      IF ((i <= m) AND (i <> n))
        THEN
          BEGIN
            FOR k := l TO n DO
              scale := scale + Abs(GetMatrixElement(A, i, k));
            IF (scale <> 0.0)
              THEN
                BEGIN
                  FOR k := l TO n DO
                    BEGIN
                      Number := GetMatrixElement(A, i, k) / scale;
                      SetMatrixElement(A, i, k, Number);
                      s := s + Sqr(Number);
                    END;
                  f := GetMatrixElement(A, i, l);
                  g := -sign(Sqrt(s), f);
                  h := f * g - s;
                  SetMatrixElement(A, i, l, f - g);
                  FOR k := l TO n DO
                    SetVectorElement(rv1, k, GetMatrixElement(A, i, k) / h);
                  IF (i <> m)
                    THEN
                      FOR j := l TO m DO
                        BEGIN
                          s := 0.0;
                          FOR k := l TO n DO
                            s := s + GetMatrixElement(A, j, k) * GetMatrixElement(A, i, k);
                          FOR k := l TO n DO
                            SetMatrixElement(A, j, k, GetMatrixElement(A, j, k) +
                              s * GetVectorElement(rv1, k));
                       END;
                  FOR k := l TO n DO
                    SetMatrixElement(A, i, k, scale * GetMatrixElement(A, i, k));
                END;
          END;
      anorm := max(anorm, (Abs(GetVectorElement(W, i)) + Abs(GetVectorElement(rv1, i))));
    END;
  { Accumulation of right-hand transform }
  FOR i := n DOWNTO 1 DO
    BEGIN
      IF (i < n)
        THEN
          BEGIN
            IF (g <> 0.0)
              THEN
                BEGIN
                  FOR j := l TO n DO
                    SetMatrixElement(V, j, i, GetMatrixElement(A, i, j) /
                      GetMatrixElement(A, i, l) / g);
                  FOR j := l TO n DO
                    BEGIN
                      s := 0.0;
                      FOR k := l TO n DO
                        s := s + GetMatrixElement(A, i, k) * GetMatrixElement(V, k, j);
                      FOR k := l TO n DO
                        SetMatrixElement(V, k, j, GetMatrixElement(V, k, j) +
                           s * GetMatrixElement(V, k, i));
                    END;
                END;
            FOR j := l TO n DO
              BEGIN
                SetMatrixElement(V, i, j, 0.0);
                SetMatrixElement(V, j, i, 0.0);
              END;
          END;
      SetMatrixElement(V, i, i, 1.0);
      g := GetVectorElement(rv1, i);
      l := i;
    END;
  { Accumulation of left-hand transformation }
  IF m < n
    THEN mnmin := m
    ELSE mnmin := n;
  FOR i := mnmin DOWNTO 1 DO
    BEGIN
      l := i + 1;
      g := GetVectorElement(W, i);
      IF (i < n)
        THEN FOR j := l TO n DO SetMatrixElement(A, i, j, 0.0);
      IF (g <> 0.0)
        THEN
          BEGIN
            g := 1.0 / g;
            IF (i <> n)
              THEN
                FOR j := l TO n DO
                  BEGIN
                    s := 0.0;
                    FOR k := l TO m DO
                      s := s + GetMatrixElement(A, k, i) * GetMatrixElement(A, k, j);
                    f := (s / GetMatrixElement(A, i, i)) * g;
                    FOR k := i TO m DO
                      SetMatrixElement(A, k, j, GetMatrixElement(A, k, j) + f *
                         GetMatrixElement(A, k, i));
                  END;
            FOR j := i TO m DO
              SetMatrixElement(A, j, i, GetMatrixElement(A, j, i) * g);
          END
        ELSE
          FOR j := i TO m DO SetMatrixElement(A, j, i, 0.0);
    SetMatrixElement(A, i, i, GetMatrixElement(A, i, i) + 1.0);
    END;
  { Diagonalisation of bidiagonal form }
  FOR k := n DOWNTO 1 DO                { loop over singular values }
    BEGIN
      FOR its := 1 TO MaxIter DO        { loop over allowed iterations }
        BEGIN
          FOR l := k DOWNTO 1 DO       { test for splitting }
            BEGIN
              nm := l - 1;
              IF ((Abs(GetVectorElement(rv1, l)) + anorm) = anorm)
                THEN GOTO 2;
              IF ((Abs(GetVectorElement(W, nm)) + anorm) = anorm)
                THEN GOTO 1;
            END;
1:        c := 0.0;               { Cancellation of rv1[l] if l > 1  }
          s := 1.0;
          FOR i := l TO k DO
            BEGIN
              f := s * GetVectorElement(rv1, i);
              IF ((Abs(f) + anorm) <> anorm)
                THEN
                  BEGIN
                    g := GetVectorElement(W, i);
                    h := Pythag(f, g);
                    SetVectorElement(W, i, h);
                    h := 1.0 / h;
                    c := (g * h);
                    s := -(f * h);
                    FOR j := 1 TO m DO
                      BEGIN
                        y := GetMatrixElement(A, j, nm);
                        z := GetMatrixElement(A, j, i);
                        SetMatrixElement(A, j, nm, (y * c) + (z * s));
                        SetmatrixElement(A, j, i, -(y * s) + (z * c));
                      END;
                  END;
            END;
2:        z := GetVectorElement(W, k);
          IF (l = k)
            THEN
              BEGIN               { Convergence }
                IF (z < 0.0)    { Singular Value is made non-negative }
                  THEN
                    BEGIN
                      SetVectorElement(W, k, -z);
                      FOR j := 1 TO n DO
                        SetMatrixElement(V, j, k, -GetMatrixElement(V, j, k));
                    END;
                GOTO 3;
              END;
          IF (its = MaxIter)
            THEN
              BEGIN
                Writeln('no convergence in ', Maxiter: 2, ' SVDCMP iterations');
                ReadLn;
              END;
          x := GetVectorElement(W, l);     { Shift from bottom 2 by 2 minor }
          nm := k - 1;
          y := GetVectorElement(W, nm);
          g := GetVectorElement(rv1, nm);
          h := GetVectorElement(rv1, k);
          f := ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
          g := Pythag(f, 1.0);
          f := ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
          { Next QR-transformation }
          c := 1.0;
          s := 1.0;
          FOR j := l TO nm DO
            BEGIN
              i := j + 1;
              g := GetVectorElement(rv1, i);
              y := GetVectorElement(W, i);
              h := s * g;
              g := c * g;
              z := Pythag(f, h);
              SetVectorElement(rv1, j, z);
              c := f / z;
              s := h / z;
              f := (x * c) + (g * s);
              g := -(x * s) + (g * c);
              h := y * s;
              y := y * c;
              FOR jj := 1 TO n DO
                BEGIN
                  x := GetMatrixElement(V, jj, j);
                  z := GetMatrixElement(V, jj, i);
                  SetMatrixElement(V, jj, j, (x * c) + (z * s));
                  SetMatrixElement(V, jj, i, -(x * s) + (z * c));
                END;
              z := Pythag(f, h);
              SetVectorElement(W, j, z);   { Rotation can be arbitrary if z=0 }
              IF (z <> 0.0)
                THEN
                  BEGIN
                    z := 1.0 / z;
                    c := f * z;
                    s := h * z;
                  END;
              f := (c * g) + (s * y);
              x := -(s * g) + (c * y);
              FOR jj := 1 TO m DO
                BEGIN
                  y := GetMatrixElement(A, jj, j);
                  z := GetMatrixElement(A, jj, i);
                  SetMatrixElement(A, jj, j, (y * c) + (z * s));
                  SetMatrixElement(A, jj, i, -(y * s) + (z * c));
                END;
            END;
          SetVectorElement(rv1, l, 0.0);
          SetVectorElement(rv1, k, f);
          SetVectorElement(W, k, x);
        END;
3:     ;
    END;
  DestroyVector(rv1);
END;

PROCEDURE svdSynth(VAR A : MatrixTyp; CONST W : VectorTyp; CONST U, V : MatrixTyp);

VAR j, k, l, n, p: INTEGER;

BEGIN
  n := MatrixRows(U);
  p := MatrixColumns(U);
  IF VectorLength(W) <> p
    THEN
      BEGIN
        CH := WriteErrorMessage('Singular value synthesis: Dimension error');
        EXIT;
      END;
  IF ((MatrixRows(V) <> p) OR (MatrixColumns(V) <> p))
    THEN
      BEGIN
        CH := WriteErrorMessage('Singular value synthesis: Dimension error');
        EXIT;
      END;
  CreateMatrix(A, n, p, 0.0);
  FOR k := 1 TO n DO
    FOR l := 1 TO p DO
      BEGIN
        FOR j := 1 TO p DO
          SetMatrixElement(A, k, l, GetMatrixElement(A, k, l) +
              GetMatrixElement(U, k, j) * GetVectorElement(W, j) * GetMatrixElement(V, l, j));
     END;
END;


END. // SingularValue

