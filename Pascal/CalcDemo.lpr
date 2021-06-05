PROGRAM CalcDemo;

USES Calc,              // formula compiler
     Mathfunc,          // standard math functions
     Vector;            // vector arithmetic

VAR a, b, dx               : REAL;
    dummy                  : REAL;
    Id                     : Calc_IdStr;
    Formula                : Calc_String;
    FormProg, Deriv        : Calc_Prog;
    x                      : Calc_VarTab;
    xv, yv, dy             : VectorTyp;
    i, n                   : word;

BEGIN
  x := NewVarTab;
  dummy := AddToVarTab(x, 'X');
  CalcDecMod := TRUE;
  REPEAT
    WriteLn;
    WriteLn;
    WriteLn('Please enter function to be evaluated (H: help, <CR> finish): ');
    WriteLn;
    Write('f(x) = '); ReadLn(Formula); WriteLn;
    if upcase(Formula) = 'H'
      then
        HelpFormula
      else
        begin
          CompileExpression(Formula, x, FormProg);
          IF CalcResult
            THEN
              BEGIN
                WriteLn('Expression "', Formula, '" compiled correctly...');
                WriteLn;
                Write('Evaluate f(x) between  a = '); Read(a);
                Write('                   and b = '); Read(b);
                Write('      with step width dx = '); ReadLn(dx);
                WriteLn;
                deriv := CalcDerivation(FormProg, x, ID);
                write('y´(x) = ');
                CalcAOS(deriv, x);
                n := Round((abs(a-b)/dx));
                CreateVector(xv, succ(n), 0.0);
                CreateVector(yv, succ(n), 0.0);
                CreateVector(dy, succ(n), 0.0);
                FOR i := 1 TO n DO
                  BEGIN
                    SetVectorElement(xv, i, a);
                    AssignVar(x, 'X', a);
                    SetVectorElement(yv, i, CalcExpression(FormProg, x));
                    SetVectorElement(dy, i, CalcExpression(deriv, x));
                    a := a + dx;
                  END;
                for i := 1 to n do
                  writeln(FloatStr(GetVectorElement(xv, i), 10), ' ', FloatStr(GetVectorElement(yv, i), 10), ' ', FloatStr(GetVectorElement(dy, i), 10));
                KillExpression(FormProg);
                DestroyVector(xv);
                DestroyVector(yv);
            END;
        end;
  UNTIL Formula = '';
  KillVarTab(x);
  WriteLn('Demo finished...');
END.
