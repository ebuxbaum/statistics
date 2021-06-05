unit StatFunc;

{ Statistische Funktionen. Signifikanzschranken und Integrale für Beta-,
  Gamma-, Normal-, F-, t- und X2-Verteilung.
  Nach einer Vorlage von PC-SIG auf eine Turbo-Pascal Unit umgeschrieben
  von E. Buxbaum 1993 }

interface

uses mathfunc, stat;

function CDBeta(X, Alpha, Beta: float; var Cprec: float): float;
(*-------------------------------------------------------------------------*)
(*       Purpose:   Evaluates CPDF of Incomplete Beta Function             *)
(*                                                                         *)
(*       Calling Sequence:                                                 *)
(*                 X      --- Upper percentage point of PDF                *)
(*                 Alpha  --- First shape parameter                        *)
(*                 Beta   --- Second shape parameter                       *)
(*                 Cprec  --- Actual resulting precision                   *)
(*       Result:   Resultant probability                                   *)
(*                                                                         *)
(*       Calls: AlGama                                                     *)
(*                                                                         *)
(*       Method:                                                           *)
(*            The continued fraction expansion as given by                 *)
(*            Abramowitz and Stegun (1964) is used.  This                  *)
(*            method works well unless the minimum of (Alpha, Beta)        *)
(*            exceeds about 70000.                                         *)
(*                                                                         *)
(*            An error in the input arguments results in a returned        *)
(*            probability of -1.                                           *)
(*-------------------------------------------------------------------------*)


function BetaInv(P, Alpha, Beta: float; var Cprec: float): float;
(*--------------------------------------------------------------------------*)
(*       Purpose:  Calculates inverse Beta distribution                     *)
(*                                                                          *)
(*       Calling Sequence:                                                  *)
(*                 P       --- Cumulative probability for which percentage  *)
(*                             point is to be found.                        *)
(*                 Alpha   --- First parameter of Beta distribution         *)
(*                 Beta    --- Second parameter of Beta distribution        *)
(*                 Cprec   --- no. of digits of precision actually obtained *)
(*                                                                          *)
(*       Result:  percentage point of Beta dist.                            *)
(*                                                                          *)
(*       Calls:  CDBeta  (cumulative Beta distribution)                     *)
(*                                                                          *)
(*       Remarks:                                                           *)
(*            Because of the relationship between the Beta                  *)
(*            distribution and the F distribution, this                     *)
(*            routine can be used to find the inverse of                    *)
(*            the F distribution.                                           *)
(*                                                                          *)
(*            The percentage point is returned as -1.0 if any error occurs  *)
(*                                                                          *)
(*            Newtons' method is used to search for the percentage point.   *)
(*            The algorithm is based upon AS110 from Applied Statistics.    *)
(*                                                                          *)
(*--------------------------------------------------------------------------*)


function ALGama(Arg: float): float;
(*-------------------------------------------------------------------------*)
(*  ALGama  -- Logarithm of Gamma Distribution                             *)
(*            Minimax polynomial approximations are used over the          *)
(*            intervals [-inf,0], [0,.5], [.5,1.5], [1.5,4.0],             *)
(*            [4.0,12.0], [12.0,+inf].                                     *)
(*                                                                         *)
(*            See Hart et al, "Computer Approximations",                   *)
(*            Wiley(1968), p. 130F, and also                               *)
(*                                                                         *)
(*            Cody and Hillstrom, "Chebyshev approximations for            *)
(*            the natural logarithm of the Gamma function",                *)
(*            Mathematics of Computation, 21, April, 1967, P. 198F.        *)
(*                                                                         *)
(*            There are some machine-dependent constants used --           *)
(*               MaxRealNumber   --- Largest value for which ALGama                 *)
(*                          can be safely computed.                        *)
(*-------------------------------------------------------------------------*)


function GammaIn(Y, P: float; var Cprec: float): float;
(*-------------------------------------------------------------------------*)
(*       Purpose:   Evaluates Incomplete Gamma integral                    *)
(*       Calling Sequence:                                                 *)
(*                 Y      --- Gamma distrib. value                         *)
(*                 P      --- Degrees of freedom                           *)
(*                 Ifault --- error indicator                              *)
(*       Result: probability                                               *)
(*       Remarks:                                                          *)
(*            Either an infinite series summation or a continued fraction  *)
(*            approximation is used, depending upon the argument range.    *)
(*-------------------------------------------------------------------------*)


function Erf(Z: float): float;
(*--------------------------------------------------------------------------*)
(*       Purpose:   Calculates the error function.                          *)
(*                                                                          *)
(*       Calls:  None                                                       *)
(*                                                                          *)
(*       Method:                                                            *)
(*                                                                          *)
(*          A rational function approximation is used, adjusted for the     *)
(*          argument size.  The approximation gives 13-14 digits of         *)
(*          accuracy.                                                       *)
(*                                                                          *)
(*       Remarks:                                                           *)
(*                                                                          *)
(*          The error function can be used to find normal integral          *)
(*          probabilities because of the simple relationship between        *)
(*          the error function and the normal distribution:                 *)
(*                                                                          *)
(*              Norm(X) := (1 + Erf( X / Sqrt(2))) / 2, X >= 0;             *)
(*              Norm(X) := (1 - Erf(-X / Sqrt(2))) / 2, X < 0.              *)
(*--------------------------------------------------------------------------*)


implementation


const
  Xln2sp = 9.18938533204673E-01;    // LogE(Sqrt(2 * Const_pi))
  MaxPrec = 16;                     // Max. precision
  LnTenInv = 0.4342944819032520;    // 1 / LN(10)

function ALGama(Arg: float): float;

var Rarg, Alinc, Scale, Top, Bot, Frac, Algval : float;
    I, Iapprox, Iof, Ilo, Ihi                  : integer;
    Qminus, Qdoit                              : boolean;

const
  P: array [1..29] of float =
    ( 4.12084318584770E+00,  8.56898206283132E+01,  2.43175243524421E+02,
     -2.61721858385614E+02, -9.22261372880152E+02, -5.17638349802321E+02,
     -7.74106407133295E+01, -2.20884399721618E+00,  5.15505761764082E+00,
      3.77510679797217E+02,  5.26898325591498E+03,  1.95536055406304E+04,
      1.20431738098716E+04, -2.06482942053253E+04, -1.50863022876672E+04,
     -1.51383183411507E+03, -1.03770165173298E+04, -9.82710228142049E+05,
     -1.97183011586092E+07, -8.73167543823839E+07,  1.11938535429986E+08,
      4.81807710277363E+08, -2.44832176903288E+08, -2.40798698017337E+08,
      8.06588089900001E-04, -5.94997310888900E-04,  7.93650067542790E-04,
     -2.77777777688189E-03,  8.33333333333330E-02);

  Q: array [1..24] of float =
    ( 1.00000000000000E+00,  4.56467718758591E+01,  3.77837248482394E+02,
      9.51323597679706E+02,  8.46075536202078E+02,  2.62308347026946E+02,
      2.44351966250631E+01,  4.09779292109262E-01,  1.00000000000000E+00,
      1.28909318901296E+02,  3.03990304143943E+03,  2.20295621441566E+04,
      5.71202553960250E+04,  5.26228638384119E+04,  1.44020903717009E+04,
      6.98327414057351E+02,  1.00000000000000E+00, -2.01527519550048E+03,
     -3.11406284734067E+05, -1.04857758304994E+07, -1.11925411626332E+08,
     -4.04435928291436E+08, -4.35370714804374E+08, -7.90261111418763E+07);

begin
  Algval := MaxRealNumber;
  Scale := 1.0;
  Alinc := 0.0;
  Frac := 0.0;
  Rarg := Arg;
  Iof := 1;
  Qminus := False;
  Qdoit := True;
  if (Rarg < 0.0)
    then
      begin
        Qminus := True;
        Rarg := -Rarg;
        Top := Int(Rarg);
        Bot := 1.0;
        if ((INT(Top / 2.0) * 2.0) = 0.0)
          then Bot := -1.0;
        Top := Rarg - Top;
        if (Top = 0.0)
          then Qdoit := False
          else
            begin
              Frac := Bot * Const_pi / SIN(Top * Const_pi);
              Rarg := Rarg + 1.0;
              Frac := LN(ABS(Frac));
            end;
      end;
  if (Rarg = 0.0)
    then Qdoit := False
    else if (Rarg <= 0.5)
           then
             begin
               Alinc := -LN(Rarg);
               Scale := Rarg;
               Rarg := Rarg + 1.0;
               if (Scale < MachineEpsilon)
                 then
                   begin
                     Algval := Alinc;
                     Qdoit := False;
                   end;
             end
           else
             if (Rarg <= 1.5)
               then Scale := Rarg - 1.0
               else if (Rarg <= 4.0)
                      then
                        begin
                          Scale := Rarg - 2.0;
                          Iof := 9;
                        end
                      else if (Rarg <= 12.0)
                             then Iof := 17
                             else if (Rarg <= MaxRealNumber)
                                    then
                                      begin
                                        Alinc := (Rarg - 0.5) * LN(Rarg) - Rarg + Xln2sp;
                                        Scale := 1.0 / Rarg;
                                        Rarg := Scale * Scale;
                                        Top := P[25];
                                        for i := 26 to 29 do
                                          Top := Top * Rarg + P[I];
                                        Algval := Scale * Top + Alinc;
                                        Qdoit := False;
                                    end;
  if Qdoit
    then
      begin
        Ilo := Iof + 1;
        Ihi := Iof + 7;
        Top := P[Iof];
        Bot := Q[Iof];
        for I := Ilo to Ihi do
          begin
            Top := Top * Rarg + P[I];
            Bot := Bot * Rarg + Q[I];
          end;
        Algval := Scale * (Top / Bot) + Alinc;
      end;
  if (Qminus)
    then Algval := Frac - Algval;
  Result := Algval;
end;   // ALGama


function CDBeta(X, Alpha, Beta: float; var Cprec: float) : float;

var
  Epsz, A, B, C, F, Fx, Apb, Zm, Alo, Ahi, Aev,
  Aod, Blo, Bhi, Bod, Bev, Zm1, D1               : float;
  Ntries, iter                                   : integer;
  Qswap, Qdoit, Qconv                            : boolean;
  ch                                             : char;

begin
  Cprec := ValidFigures;
  Epsz := pot(10, -ValidFigures);
  A := Alpha;
  B := Beta;
  QSwap := False;
  CDBeta := -1.0;
  Qdoit := True;
  (* Check arguments. Error if: X <= 0  or x >= 1 or  A <= 0  or  B <= 0     *)
  if ((X <= 0.0) or (X >= 1.0) or (A <= 0.0) or (B <= 0.0))
    then
      begin
        ch := WriteErrorMessage('CD beta: illegal parameter');
        StatError := true;
        exit;
      end;
  CDBeta := 1.0;
  (* If X > A / ( A + B ) then swap A, B for more efficient eval.  *)
  if (X > (A / (A + B)))
    then
      begin
        X := 1.0 - X;
        A := Beta;
        B := Alpha;
        QSwap := True;
      end;
  (* Check for extreme values *)
  C := ALGama(A + B) + A*LN(X) + B*LN(1.0 - X) - ALGama(A) -  ALGama(B) - LN(A - X*(A + B));
  if (((C < -36.0) and QSwap) or (C < -180.0))
    then
      begin
        ch := WriteErrorMessage('CD beta: Outside range');
        StatError := true;
        exit;
      end;
  (*  Set up continued fraction expansion evaluation. *)
  CDBeta := 0.0;
  Apb := A + B;
  Zm := 0.0;
  Alo := 0.0;
  Bod := 1.0;
  Bev := 1.0;
  Bhi := 1.0;
  Blo := 1.0;
  Ahi := EXP(ALGama(Apb) + A * LN(X) + B * LN(1.0 - X) - ALGama(A + 1.0) - ALGama(B));
  F := Ahi;
  Iter := 0;
  Qconv := False;
  repeat   // Continued fraction loop
    Fx := F;
    Zm1 := Zm;
    Zm := Zm + 1.0;
    D1 := A + Zm + Zm1;
    Aev := -(A + Zm1) * (Apb + Zm1) * X / D1 / (D1 - 1.0);
    Aod := Zm * (B - Zm) * X / D1 / (D1 + 1.0);
    Alo := Bev * Ahi + Aev * Alo;
    Blo := Bev * Bhi + Aev * Blo;
    Ahi := Bod * Alo + Aod * Ahi;
    Bhi := Bod * Blo + Aod * Bhi;
    if ABS(Bhi) <= Zero then Bhi := 0.0;
    if (Bhi <> 0.0)
      then
        begin
          F := Ahi / Bhi;
          Qconv := (ABS((F - Fx) / F) < Epsz);
        end;
    inc(Iter);
  until ((Iter > MaxIter) or Qconv);
  if (iter > MaxIter)
    then
      begin
        ch := WriteErrorMessage('CD beta: Did not converge');
        StatError := true;
        exit;
      end;
   if (Qswap)
     then Result := 1.0 - F
     else Result := F;
end;   (* CDBeta *)


function BetaInv(P, Alpha, Beta: float): float;

var
  Eps, Xim1, Xi, Xip1, Fim1, Fi, W, Cmplbt, Adj,
  Sq, R, S, T, G, A, B, PP, H, A1, B1, Eprec      : float;
  Done                                            : boolean;
  Jter, iter                                      : integer;
  ch                                              : char;

label 10, 30;

begin
  Result := P;
  if ((Alpha <= 0.0) or (Beta <= 0.0) or (P >= 1.0) or (P <= 0.0))
  then
    begin
      ch := WriteErrorMessage('beta inv: illegal parameter');
      StatError := true;
      exit;
    end;
  (* Set precision *)
  Eps := pot(10, -2 * ValidFigures);
  (* Flip params if needed so that P for evaluation is <= .5     *)
  if (P > 0.5)
    then
      begin
        A := Beta;
        B := Alpha;
        PP := 1.0 - P;
      end
    else
      begin
        A := Alpha;
        B := Beta;
        PP := P;
      end;
  (* Generate initial approximation. Several different ones used, depending upon parameter values. *)
  Cmplbt := ALGama(A) + ALGama(B) - ALGama(A + B);
  Fi := Ninv(1.0 - PP);
  if ((A > 1.0) and (B > 1.0))
    then
      begin
        R := (Fi * Fi - 3.0) / 6.0;
        S := 1.0 / (A + A - 1.0);
        T := 1.0 / (B + B - 1.0);
        H := 2.0 / (S + T);
        W := Fi * SQRT(H + R) / H - (T - S) * (R + 5.0 / 6.0 -
          2.0 / (3.0 * H));
        Xi := A / (A + B * EXP(W + W));
      end
    else
      begin
        R := B + B;
        T := 1.0 / (9.0 * B);
        T := R * Pot((1.0 - T + Fi * SQRT(T)), 3);
        if (T <= 0.0)
          then Xi := 1.0 - EXP((LN((1.0 - PP) * B) + Cmplbt) / B)
          else
            begin
              T := (4.0 * A + R - 2.0) / T;
              if (T <= 1.0)
                then Xi := EXP((LN(PP * A) + Cmplbt) / PP)
                else Xi := 1.0 - 2.0 / (T + 1.0);
            end;
      end;
  (* Force initial estimate to reasonable range.         *)
  if (Xi < 0.0001) then Xi := 0.0001;
  if (Xi > 0.9999) then Xi := 0.9999;
  (* Set up Newton-Raphson loop *)
  A1 := 1.0 - A;
  B1 := 1.0 - B;
  Fim1 := 0.0;
  Sq := 1.0;
  Xim1 := 1.0;
  Iter := 0;
  Done := False;
  (* Begin Newton-Raphson loop  *)
  repeat
    inc(Iter);
    Done := Done or (Iter > MaxIter);
    Fi := CDBeta(Xi, A, B, ValidFigures + 1, MaxIter, Eprec, Jter, Ierr);
    if (Ierr <> 0)
      then
        begin
          Ierr := 2;
          Done := True;
        end
      else
        begin
          Fi := (Fi - PP) * EXP(Cmplbt + A1 * LN(Xi) + B1 * LN(1.0 - Xi));
          if ((Fi * Fim1) <= 0.0) then Xim1 := Sq;
          G := 1.0;
10:       repeat
              Adj := G * Fi;
              Sq := Adj * Adj;
              if (Sq >= Xim1) then G := G / 3.0;
          until (Sq < Xim1);
          Xip1 := Xi - Adj;
          if ((Xip1 < 0.0) or (Xip1 > 1.0))
            then
              begin
                G := G / 3.0;
                goto 10;
              end;
          if ((Xim1 <= Eps) or (Fi * Fi <= Eps)) then goto 30;
          if ((Xip1 = 0.0) or (Xip1 = 1.0))
            then
              begin
                G := G / 3.0;
                goto 10;
              end;
          if (Xip1 <> Xi)
            then
              begin
                Xi := Xip1;
                Fim1 := Fi;
              end
            else
              Done := True;
        end;
  until (Done);
30:
  if (P > 0.5)
    then Result := 1.0 - Xi
    else Result := Xi;
end   (* BetaInv *);


function GammaIn(Y, P: float): float;

const
  Oflo = 1.0E+37;
  MinExp = -87.0;

var
  F, C, A, B, Term, An, Gin, Rn, Dif, Eps : float;
  Pn                                      : array[1..6] of float;
  Done                                    : boolean;
  iter                                    : word;
  ch                                      : char;

begin
  (* Check arguments *)
  Ifault := 1;
  GammaIn := 1.0;
  if ((Y <= 0.0) or (P <= 0.0))
    then
      begin
        ch := WriteErrorMessage('gamma inv: illegal parameter');
        StatError := true;
        exit;
      end;
  (* Check value of F *)
  Ifault := 0;
  F := P * LN(Y) - ALGama(P + 1.0) - Y;
  if (F < MinExp)
    then
      begin
        ch := WriteErrorMessage('beta inv: illegal parameter');
        StatError := true;
        exit;
      end;
  F := EXP(F);
  if (F = 0.0)
    then
      begin
        ch := WriteErrorMessage('beta inv: illegal parameter');
        StatError := true;
        exit;
      end;
  (* Set precision *)
  Eps := Pot(10, -ValidFigures);
  (* Choose infinite series or continued fraction.       *)
  if ((Y > 1.0) and (Y >= P))
    then
      begin (* Continued Fraction *)
        A := 1.0 - P;
        B := A + Y + 1.0;
        Term := 0.0;
        Pn[1] := 1.0;
        Pn[2] := Y;
        Pn[3] := Y + 1.0;
        Pn[4] := Y * B;
        Gin := Pn[3] / Pn[4];
        Done := False;
        Iter := 0;
        repeat
          inc(Iter);
          A := A + 1.0;
          B := B + 2.0;
          Term := Term + 1.0;
          An := A * Term;
          Pn[5] := B * Pn[3] - An * Pn[1];
          Pn[6] := B * Pn[4] - An * Pn[2];
          if (Pn[6] <> 0.0)
            then
              begin
                Rn := Pn[5] / Pn[6];
                Dif := ABS(Gin - Rn);
                if (Dif <= Eps)
                  then if (Dif <= (Eps * Rn))
                         then Done := True;
                Gin := Rn;
              end;
          Pn[1] := Pn[3];
          Pn[2] := Pn[4];
          Pn[3] := Pn[5];
          Pn[4] := Pn[6];
          if (ABS(Pn[5]) >= Oflo)
            then
              begin
                Pn[1] := Pn[1] / Oflo;
                Pn[2] := Pn[2] / Oflo;
                Pn[3] := Pn[3] / Oflo;
                Pn[4] := Pn[4] / Oflo;
              end;
        until (Iter > MaxIter) or Done;
        Result := 1.0 - (F * Gin * P);
      end   (* Continued Fraction *)
    else
      begin (* Infinite series *)
        Iter := 0;
        Term := 1.0;
        C := 1.0;
        A := P;
        Done := False;
        repeat
          A := A + 1.0;
          Term := Term * Y / A;
          C := C + Term;
          inc(Iter);
        until (Iter > MaxIter) or ((Term / C) <= Eps);
        Result := C * F;
      end   (* Infinite series *);
end;    (* GammaIn *)


function Erf(Z: float): float;

const
  A: array[1..14] of float =
    (1.1283791670955,      0.34197505591854,     0.86290601455206E-1,
     0.12382023274723E-1,  0.11986242418302E-2,  0.76537302607825E-4,
     0.25365482058342E-5, -0.99999707603738,    -1.4731794832805,
    -1.0573449601594,     -0.44078839213875,    -0.10684197950781,
    -0.12636031836273E-1, -0.1149393366616E-8);

  B: array[1..12] of float =
    (-0.36359916427762,     0.52205830591727E-1, -0.30613035688519E-2,
     -0.46856639020338E-4,  0.15601995561434E-4, -0.62143556409287E-6,
      2.6015349994799,      2.9929556755308,      1.9684584582884,
      0.79250795276064,     0.18937020051337,     0.22396882835053E-1);

var
  U, X, S: float;

begin
  X := ABS(Z);
  if Z >= 0.0
    then S := 1.0
    else S := -1.0;
  if (abs(Z) < Zero)
    then Erf := 0.0
    else if (X >= 5.5)
           then Erf := S
           else
             begin
               U := X * X;
               if (X <= 1.5)
                 then Erf :=
                    (X*EXP(-U)*(A[1] + U*(A[2] + U*(A[3] + U*(A[4] +
                     U*(A[5] + U*(A[6] + U*A[7])))))) / (1.0 + U * (B[1] +
                     U*(B[2] + U*(B[3] + U*(B[4] + U*(B[5] + U*B[6]))))))) * S
                 else Erf :=
                    (EXP(-U)*(A[8] + X*(A[9] + X*(A[10] + X*(A[11] + X*(A[12] +
                     X*(A[13] + X*A[14])))))) / (1.0 + X*(B[7] + X*(B[8] +
                     X*(B[9] + X*(B[10] + X*(B[11] + X*B[12])))))) + 1.0) * S;
              end;
end;   (* Erf *)


function SigNorm(X: float): float;

begin
  if X >= 0.0
    then SigNorm := 1.0 - (1.0 + Erf(X / Const_sqrt2)) / 2.0
    else SigNorm := 1.0 - (1.0 - Erf(-X / Const_sqrt2)) / 2.0;
end   (* SigNorm *);


function Ninv(P: float): float;

const
  Lim = 1.0E-20;

var
  Y: float;
  Pr: float;
  Nv: float;

const
  PN: array [1..5] of float =
    (-0.322232431088,  -1.0,               -0.342242088547,
     -0.0204231210245, -0.453642210148E-4);

  QN: array [1..5] of float =
    ( 0.0993484626060,  0.588581570495,     0.531103462366,
      0.103537752850,   0.38560700634E-2);

begin (* Ninv *)
  Ninv := 0.0;
  if (P > 0.5)
    then Pr := 1.0 - P
    else Pr := P;
  if ((Pr >= Lim) and (Pr <> 0.5))
    then
      begin
        Y := SQRT(LN(1.0 / Pr / Pr));
        Nv := Y + ((((Y * PN[5] + PN[4]) * Y + PN[3]) * Y + PN[2]) * Y + PN[1]) /
          ((((Y * QN[5] + QN[4]) * Y + QN[3]) * Y + QN[2]) * Y + QN[1]);
        if (P < 0.5)
          then Result := -Nv
          else Result := Nv;
      end;
end   (* Ninv *);


function Ninv2(P: float): float;

var
  Xp, P1, Z, X3, X2, X1, Phi: float;

begin
  Xp := Ninv(P);
  P1 := SigNorm(Xp);
  Phi := SQRT(1.0 / (2.0 * Const_pi)) * EXP(-(Xp * Xp) / 2.0);
  Z := (P - P1) / Phi;
  X3 := (2.0 * (Xp * Xp) + 1.0) * Z / 3.0;
  X2 := (X3 + Xp) * Z / 2.0;
  X1 := ((X2 + 1.0) * Z);
  Result := Xp + X1;
end   (* Ninv2 *);


function CDNorm(X: float): float;

begin
  if X >= 0.0
    then CDNorm := (1.0 + Erf(X / Const_sqrt2)) / 2.0
    else CDNorm := (1.0 - Erf(-X / Const_sqrt2)) / 2.0;
end;


function SigF(F, Dfn, Dfd: float): float;

var
  Iter, Ifault: integer;
  Cprec, Pval: float;

begin
  Pval := -1.0;
  if (Dfn > 0.0) and (Dfd > 0.0)
    then
      begin
        Pval := CDBeta(Dfd / (Dfd + F * Dfn), Dfd / 2.0, Dfn /
          2.0, ValidFigures, MaxIter, Cprec, Iter, Ifault);
        if Ifault <> 0 then Pval := -1.0;
      end;
  Result := Pval;
end   (* SigF *);


function Finv(Alpha, Dfn, Dfe: float): float;

var
  Fin, Cprec: float;
  Iter, Ierr: integer;

begin
  Fin := -1.0;
  if ((Dfn > 0.0) and (Dfe > 0.0))
    then if ((Alpha >= 0.0) and (Alpha <= 1.0))
           then
             begin
               Fin := BetaInv(1.0 - Alpha, Dfn / 2.0, Dfe / 2.0, MaxIter,
                ValidFigures, Iter, Cprec, Ierr);
               if ((Fin >= 0.0) and (Fin < 1.0) and (Ierr = 0))
                 then Result := Fin * Dfe / (Dfn * (1.0 - Fin));
             end;
end   (* Finv *);


function SigChi(Chisq, Df: float): float;

var
  Ierr, Iter: integer;
  Cprec: float;

begin
  SigChi := 1.0 - GammaIn(Chisq / 2.0, Df / 2.0, ValidFigures, MaxIter,
    Cprec, Iter, Ierr);
  if (Ierr <> 0) then Result := -1.0;
end   (* SigChi *);


function Cinv(P, V: float): float;

const
  E = 1.0E-8;

var
  XX, C, Ch, Q, P1, P2, T, X, B, A, G, S1, S2, S3, S4, S5, S6, Cprec: float;
  Iter: integer;

begin
  (* Test arguments for validity *)
  Cinv := -1.0;
  if ((P < E) or (P > (1.0 - E)) or (V <= 0.0))
    then
      begin
        ch := WriteErrorMessage('C inv: illegal parameter');
        StatError := true;
        exit;
      end;
  Ifault := 2;
  (* Initialize *)
  P := 1.0 - P;
  XX := V / 2.0;
  G := ALGama(XX);
  C := XX - 1.0;
  (* Starting approx. for small chi-square *)
  if (V < (-1.24 * LN(P)))
    then
      begin
        Ch := Pot(P * XX * EXP(G + XX * Const_ln2), (1.0 / XX));
        if Ch < E
          then
            begin
              Result := Ch;
              exit;
            end;
      end
    else if (V <= 0.32)  (* Starting approx. for v <= .32 *)
           then
             begin
               Ch := 0.4;
               A := LN(1.0 - P);
               repeat
                 Q := Ch;
                 P1 := 1.0 + Ch * (4.67 + Ch);
                 P2 := Ch * (6.73 + Ch * (6.66 + Ch));
                 T := -0.5 + (4.67 + 2.0 * Ch) / P1 -
                  (6.73 + Ch * (13.32 + 3.0 * Ch)) / P2;
                 Ch := Ch - (1.0 - EXP(A + G + 0.5 * Ch + C * Const_ln2) *
                  P2 / P1) / T;
               until (ABS(Q / Ch - 1.0) <= 0.01);
             end
           else
             begin
               (* Starting approx. using Wilson and Hilferty estimate                 *)
               X := Ninv(P);
               P1 := 2.0 / (9.0 * V);
               Ch := V * pot((X * SQRT(P1) + 1.0 - P1), 3);
               (* Starting approx. for P ---> 1     *)
               if (Ch > (2.2 * V + 6.0))
                 then Ch := -2.0 * (LN(1.0 - P) - C * LN(0.5 * Ch) + G);
            end;
  (* We have starting approximation. Begin improvement loop. *)
  repeat
    (* Get probability of current approx. to percentage point. *)
    Q := Ch;
    P1 := 0.5 * Ch;
    P2 := P - GammaIn(P1, XX, ValidFigures, MaxIter, Cprec, Iter, Ifault);
    if (Ifault <> 0) or (Iter > MaxIter)
      then Ifault := 3
      else
        begin
          (* Calculate seven-term Taylor series *)
          T := P2 * EXP(XX * Const_ln2 + G + P1 - C * LN(Ch));
          B := T / Ch;
          A := 0.5 * T - B * C;
          S1 := (210.0 + A*(140.0 + A*(105.0 + A*(84.0 + A*(70.0 + 60.0*A))))) / 420.0;
          S2 := (420.0 + A*(735.0 + A*(966.0 + A*(1141.0 + 1278.0*A)))) / 2520.0;
          S3 := (210.0 + A*(462.0 + A*(707.0 + 932.0*A))) / 2520.0;
          S4 := (252.0 + A*(672.0 + 1182.0*A) + C*(294.0 + A*(889.0 + 1740.0*A))) / 5040.0;
          S5 := (84.0 + 264.0*A + C*(175.0 + 606.0*A)) / 2520.0;
          S6 := (120.0 + C*(346.0 + 127.0*C)) / 5040.0;
          Ch := Ch + T*(1.0 + 0.5*T*S1 - B*C*(S1 - B*(S2 - B*(S3 - B*(S4 - B*(S5 - B*S6))))));
        end;
  until (ABS((Q / Ch) - 1.0) <= E) or (Ifault <> 0);
  if Ifault = 0
    then Result := Ch
    else Result := -1.0;
end;   (* Cinv *)


function Sigt(t, Df: float): float;

var
  Iter, Ifault: integer;
  Cprec, Pval: float;

begin
  Pval := -1.0;
  if (Df > 0.0)
    then
      begin
        Pval := CDBeta(Df / (Df + t * t), Df / 2.0, 0.5, ValidFigures, MaxIter, Cprec, Iter, Ifault);
        if Ifault <> 0 then Pval := -1.0;
      end;
  Result := Pval;
end   (* Sigt *);


function tinv(Alpha, Df: float): float;

var
  tin, Cprec: float;
  Iter, Ierr: integer;

begin
  tin := -1.0;
  if (Df > 0.0)
    then if ((Alpha >= 0.0) and (Alpha <= 1.0))
           then
             begin
               tin := BetaInv(Alpha, 0.5, Df / 2.0, MaxIter,
                 ValidFigures, Iter, Cprec, Ierr);
               if ((tin >= 0.0) and (tin < 1.0) and (Ierr = 0))
                 then Result := SQRT(tin * Df / (1.0 - tin));
             end;
end   (* tinv *);


end. { StatFunc }
{ ************************************************************************* }
program Test;
procedure TesSigF;
var F, Dfh, Dfe, P: float;
Done: boolean;
begin
Done := False;
ClrScr;
repeat
write('Enter F-value, numerator DF, and denominator DF (-99 to stop): ');
READLN(F, Dfh, Dfe);
if (F <> -99.0)
then
begin
P := SigF(F, Dfh, Dfe);
WRITELN('Significance value = ', P: 8: 5);
end
else
Done := True;
until Done;
end; (* TestSigF *)procedure TestInvf;
var  F, Dfh, Dfe, P: float;
Done: boolean;
begin
Done := False;
ClrScr;
repeat
write('Enter tail prob., num. and den. DF (-99 to stop):  ');
READLN(P, Dfh, Dfe);
if (P > 0.0)
then
begin
F := Finv(P, Dfh, Dfe);
WRITELN('F percentage point = ', F: 12: 5);
end
else
Done := True;
until Done;
end;
procedure TestSigc;
var Chisq, Df, P: float;
Done: boolean;
begin
Done := False;
ClrScr;
repeat
write('Enter Chi-square value and degrees of freedom (-99 to stop): ');
READLN(Chisq, Df);
if (Chisq > 0.0)
then
begin
P := Sigchi(Chisq, Df);
WRITELN('Significance of chi-square = ', P: 8: 5);
end
else
Done := True;
until Done;
end;   (* TestSigc *)procedure TestInvC;
var Chisq, Df, P: float;
Done: boolean;
Ierr: integer;
begin
Done := False;
ClrScr;
repeat
write('Enter tail probability value and degrees of freedom: ');
READLN(P, Df);
if (P > 0.0)
then
begin
Chisq := Cinv(P, Df, Ierr);
WRITELN('Chi-square percentage point = ', Chisq: 12: 5);
end
else
Done := True;
until Done;
end;   (* TestInvC *)procedure TestSign;
var Z, P: float;
Done: boolean;
begin
Done := False;
repeat
write('Enter Z: ');
READLN(Z);
if (Z <> -99.0)
then
begin
P := SigNorm(Z);
WRITELN('P = ', P: 8: 6);
end
else
Done := True;
until Done;
end;   (* TestSign *)procedure TestInvN;
var Z, P: float;
Done: boolean;
begin
Done := False;
ClrScr;
repeat
write('Enter tail probability: ');
READLN(P);
if (P > 0.0)
then
begin
Z := Ninv2(P);
WRITELN('Normal distribution percentage point = ', Z: 12: 5);
end
else
Done := True;
until Done;
end;   (* TestInvN *)procedure TestSigt;
var t, Df, P: float;
Done: boolean;
begin
Done := False;
ClrScr;
repeat
write('Enter t-value and degrees of freedom (-99 to stop): ');
READLN(t, Df);
if (t > 0.0)
then
begin
P := Sigt(t, Df);
WRITELN('Two-tailed significance = ', P: 8: 5);
WRITELN('One-tailed significance - ', (P / 2.0): 8: 5);
end
else
Done := True;
until Done;
end;   (* TestSigt *)procedure TestInvt;
var t1, t2, Df, P: float;
Done: boolean;
begin
Done := False;
ClrScr;
repeat
write('Enter tail probability value and degrees of freedom: ');
READLN(P, Df);
if (P > 0.0)
then
begin
t1 := Tinv(P, Df);
t2 := Tinv(P / 2.0, Df);
WRITELN('One-tailed percentage point = ', t1: 12: 5);
WRITELN('Two-tailed percentage point = ', t2: 12: 5);
end
else
Done := True;
until Done;
end;   (* TestInvt *)begin { test main }  TestSigF;
TestInvF;
TestSigT;
TestInvT;
TestSigC;
TestInvC;
TestSigN;
TestInvN;
end.
















































