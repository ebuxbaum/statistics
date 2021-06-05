unit Calc;
{ calculates arithmetic expressions during run-time
  Gieselmann, K. and Ceol, M.: CALC - ein mathematischer Compiler,
  PASCAL Int. 1987 No. 8, 52-60                                     }

interface

uses Crt, math, MathFunc;

const
  Calc_IdLen = 10;
  Calc_MaxVar = 10;
  Calc_OpSize = 6;

type
  ErrString = string[80];
  Calc_IdStr = string[Calc_IdLen];
  Calc_String = string[255];
  Calc_Operand = double;

  Calc_VarType = record
                   VarId: Calc_IdStr;
                   Value: Calc_Operand;
                 end;

  Calc_VarTab = ^Calc_VarTable;
  Calc_VarTable = array[0..Calc_MaxVar] of Calc_VarType;
  Calc_Symbols = (Calc_Err, Calc_EOE, Calc_Const, Calc_Var, Calc_Pi,
                  Calc_E, Calc_lp, Calc_rp, Calc_Neg, CAlc_Add, Calc_Sub,
                  CAlc_Mul, Calc_Dvd, Calc_Div, Calc_Mod, Calc_ggT,
                  Calc_kgV, Calc_Pow, Calc_sqr, Calc_Sqrt, Calc_Exp,
                  Calc_Ln, CAlc_Lg, Calc_Ld, CAlc_Sin, Calc_Cos, Calc_Tan,
                  Calc_Cot, Calc_ArcSin, CAlc_ArcCos, Calc_ArcTan,
                  Calc_ArcCot, Calc_Sinh, Calc_Cosh, Calc_Tanh, CAlc_Coth,
                  Calc_ArcSinh, Calc_ArcCosh, Calc_ArcTanh, Calc_ArcCoth,
                  Calc_Abs, Calc_Deg, CAlc_Rad, Calc_Rez, CAlc_Fak,
                  Calc_Sign, Calc_Int, Calc_End);

  Calc_Prog = ^Calc_Instruct;

  Calc_Instruct = record
                    NextInst: Calc_Prog;
                    Instruct: Calc_Symbols;
                    case Calc_Symbols of
                      Calc_Var  : (VarIndex: integer);
                      Calc_Const: (Operand: Calc_Operand);
                  end;

var   CalcDecMod, CalcResult: boolean;


procedure CompileExpression(Expr: Calc_String; var VarTable: Calc_VarTab;
  var ExprPtr: Calc_Prog);
{ turn arithmetic expressions into UPN }

function CalcExpression(ExprPtr: Calc_Prog; VarTable: Calc_VarTab): Calc_Operand;
{ calculate the result of an expression }

function CalcDerivation(pptr: Calc_Prog; VarTab: Calc_VarTab;
          Nach: Calc_IDStr): Calc_Prog;
{ symbolic derivative of expression }

procedure CalcAOS(pptr: Calc_Prog; VarTable: Calc_VarTab);
{ turn UPN-notation into AOS-formula and display }

procedure CalcError(ErrNum: integer; Message: ErrString);
{catch errors }

procedure KillExpression(var ExprPtr: Calc_Prog);
{ delete Calc_prog if no longer needed }

function NewVarTab: Calc_VarTab;
{ create new, empty vatiable table}

procedure KillVarTab(var VarTab: Calc_VarTab);
{ delete variable table if no longer needed }

function SearchVarTab(VarTab: Calc_VarTab; ID: Calc_String): integer;
{ search a variable in the variable table and return its index }

function AddToVarTab(VarTab: Calc_VarTab; ID: Calc_String): integer;
{ enter a new variable into the variable table and return its index }

procedure AssignVarI(VarTab: Calc_VarTab; i: integer; x: Calc_Operand);
{ assign the i-th variable in the table the value x }

procedure AssignVar(VarTab: Calc_VarTab; ID: Calc_String; x: Calc_Operand);
{ assign the value x to the variable ID in the variable table}

procedure HelpFormula;
{ produces a short help text }

implementation

const
  Calc_Ids: array [Calc_Symbols] of Calc_IdStr =
    ('ERR', ';', 'CONST', 'VAR', 'PI', 'E', '(', ')', 'NEG',
    '+', '-', '*', '/', 'DIV', 'MOD', 'GGT', 'KGV', '^',
    'SQR', 'SQRT', 'EXP', 'LN', 'LG', 'LD', 'SIN', 'COS',
    'TAN', 'COT', 'ARCSIN', 'ARCCOS', 'ARCTAN', 'ARCCOT',
    'SINH', 'COSH', 'TANH', 'COTH', 'ARCSINH', 'ARCCOSH',
    'ARCTANH', 'ARCCOTH', 'ABS', 'DEG', 'RAD', 'REZ', 'FAK',
    'SGN', 'INT', 'END');


procedure CalcError(ErrNum: integer; Message: ErrString);

var
  Meldung: string;
  ch: char;

begin
  case ErrNum of
    1: Meldung := ' *** Run time error: Floating point overflow' + Message;
    2: Meldung := ' *** Run time error: Division durch Null' + Message;
    3: Meldung := ' *** Run time error: Argument error in ' + Message;
    else
      Meldung := ' *** Run time error: ' + Message;
  end;
  ch := MathFunc.WriteErrorMessage(Meldung);
  CalcResult := False;
end;


procedure KillExpression(var ExprPtr: Calc_Prog);

var
  NextPtr: Calc_Prog;

begin
  while ExprPtr <> nil do
  begin
    NextPtr := ExprPtr^.NextInst;
    Dispose(ExprPtr);
    ExprPtr := NextPtr;
  end;
  ExprPtr := nil;
end;


function NewVarTab: Calc_VarTab;

var
  VarTab: Calc_VarTab;

begin
  Result := nil;
  try
    begin
      New(VarTab);
      VarTab^[0].Value := 0.0;
      Result := VarTab;
    end
  except
    CalcError(0, ' not enough memory');
  end;
end;

procedure KillVarTab (var VarTab: Calc_VarTab);

begin
  if vartab <> nil then
    Dispose(VarTab);
  VarTab := nil;
end;


function SearchVarTab(VarTab: Calc_VarTab; Id: Calc_String): integer;

var
  i: integer;

begin
  if vartab <> nil
    then
      begin
        for i := 1 to Length(id) do
          id[i] := upcase(id[i]);
        i := Trunc(VarTab^[0].Value);
        VarTab^[0].VarId := Copy(Id, 1, Calc_IdLen);
        while VarTab^[i].VarId <> Id do
          i := Pred(i);
        Result := i;
      end
    else
      Result := 0;
end;

function AddToVarTab(VarTab: Calc_VarTab; Id: Calc_String): integer;

var
  i: integer;

begin
  if VarTab <> nil
    then
      begin
        for i := 1 to Length(id) do
          id[i] := upcase(id[i]);
        i := Trunc(VarTab^[0].Value);
        if i < Calc_MaxVar
          then
            begin
              i := Succ(i);
              VarTab^[0].Value := i;
              VarTab^[i].VarId := Id;
              VarTab^[i].Value := 0;
            end
          else
            i := -1;
        Result := i;
      end
    else
      Result := -1;
end;

procedure AssignVarI(VarTab: Calc_VarTab; i: integer; x: Calc_Operand);

begin
  if vartab <> nil
    then
      begin
        if (i > 0) and (i <= Trunc(VarTab^[0].Value))
          then VarTab^[i].Value := x
          else CalcError(0, 'value assigned to unknown variable');
      end
    else
      CalcError(0, 'value assigned to unknown variable');
end;



procedure AssignVar(VarTab: Calc_VarTab; Id: Calc_String; x: Calc_Operand);

begin
  AssignVarI(VarTab, SearchVarTab(VarTab, Id), x);
end;



procedure invert(pptrstart: calc_prog);

var
  pptr, pptr1, pptr2: calc_prog;
  max, i: integer;
  dummy: calc_instruct;

begin
  if pptrstart <> nil
    then
      begin
        pptr := pptrstart^.nextinst;
        max := 0;
        while pptr^.nextinst <> nil do
          begin
            pptr := pptr^.nextinst;
            max := Succ(max);
          end;
        pptr := pptrstart;
        repeat
          pptr := pptr^.nextinst;
          pptr1 := pptr;
          for i := 1 to max do
            pptr1 := pptr1^.nextinst;
          dummy := pptr^;
          pptr^ := pptr1^;
          pptr1^ := dummy;
          pptr2 := pptr^.nextinst;
          pptr^.nextinst := pptr1^.nextinst;
          pptr1^.nextinst := pptr2;
          max := max - 2;
        until max <= 0;
      end; { then }
end;


function endof(pptr: calc_prog): calc_prog;

var
  help: calc_prog;
  op: integer;

begin
  if pptr <> nil
    then
      begin
        op := 1;
        repeat
          help := pptr;
          if pptr^.instruct in [calc_var, calc_const]
            then op := Pred(op)
            else if not (pptr^.instruct in [calc_neg, calc_sqr..calc_fak])
               then
                 op := Succ(op);
          pptr := pptr^.nextinst
        until op = 0;
        Result := help;
      end
    else
      Result := nil;
end; // EndOf


procedure CalcSimplify(var pptr: calc_prog);

var
  pptr1, help1, help2     : calc_prog;
  op                      : integer;
  dummy                   : calc_operand;
  arg1, arg2, SimpleError : boolean;
  helpinstruct            : calc_symbols;


  function Equal(pptr1, pptr2: calc_prog): boolean;

  var
    help1, help2 : calc_prog;
    check        : boolean;

  begin
    help1 := endof(pptr1);
    help2 := endof(pptr2);
    if (pptr1 <> nil) and (pptr2 <> nil) and (help1 <> nil) and (help2 <> nil)
      then
        begin
          check := True;
          repeat
            check := check and (pptr1^.instruct = pptr2^.instruct);
            case pptr1^.instruct of
              calc_const : check := check and (pptr1^.operand = pptr2^.operand);
              calc_var   : check := check and (pptr1^.varindex = pptr2^.varindex)
            end;
            pptr1 := pptr1^.nextinst;
            pptr2 := pptr2^.nextinst
          until not check or (pptr1 = help1^.nextinst) and (pptr2 = help2^.nextinst);
          Result := check;
        end
      else
        Result := False;
  end; // Equal


  function compute(pptr, pptr1, pptr2: calc_prog): calc_operand;

  var
    exptr, a, b, c : calc_prog;
    vardummy       : calc_vartab;

  begin
    try
      begin
        vardummy := newvartab;
        New(a);
        New(b);
        New(c);
        a^ := pptr^;
        b^ := pptr1^;
        if pptr2 <> nil
          then c^ := pptr2^;
        a^.nextinst := nil;
        New(exptr);
        exptr^.nextinst := b;
        if pptr2 <> nil
          then
            begin
              b^.nextinst := c;
              c^.nextinst := a;
            end
          else
            b^.nextinst := a;
        Result := calcexpression(exptr, vardummy);
        SimpleError := SimpleError or not calcresult;
        Dispose(a);
        Dispose(b);
        Dispose(c);
        Dispose(exptr);
        killvartab(vardummy);
      end
    except
        begin
          Result := 0.0;
          SimpleError := True;
        end;
    end;
  end; // Compute


  procedure simple(pptr: calc_prog);

  var
    pptra, pptrb : calc_prog;


    procedure restoreptr;

    begin
      pptra := pptr^.nextinst;
      pptrb := endof(pptra);
      if pptrb <> nil
        then pptrb := pptrb^.nextinst;
    end;


    procedure erase_entry;

    begin
      while help1 <> pptrb do
        begin
          help2 := help1;
          help1 := help1^.nextinst;
          Dispose(help2);
        end;
    end;


    procedure pusha;

    begin
      pptr^ := pptra^;
      help1 := pptr;
      while help1^.nextinst <> pptrb do
        help1 := help1^.nextinst;
      help1^.nextinst := pptrb^.nextinst;
      Dispose(pptra);
      Dispose(pptrb);
      restoreptr;
    end;


    procedure skipa;

    begin
      help1 := pptr^.nextinst;
      erase_entry;
      pptr^ := pptrb^;
      Dispose(pptrb);
      restoreptr;
    end;


    procedure setconst(dummy: calc_operand);

    begin
      pptr^.instruct := calc_const;
      pptr^.operand := dummy;
      pptr^.nextinst := pptrb^.nextinst;
      help1 := pptra;
      erase_entry;
      Dispose(pptrb);
      restoreptr;
    end;

  begin // Simple
    if pptr <> nil
      then
        begin
          restoreptr;
          if pptr^.instruct in [calc_neg, calc_sqr..calc_fak]
            then
              begin
                simple(pptra);
                if pptra^.instruct = calc_const // ausrechnen !
                  then
                    begin
                      dummy := compute(pptr, pptra, nil);
                      if dummy = -0.0
                        then dummy := 0.0;
                      pptr^.instruct := calc_const;
                      pptr^.operand := dummy;
                      pptr^.nextinst := pptra^.nextinst;
                      Dispose(pptra);
                      restoreptr;
                    end;
                if pptr^.instruct = calc_neg
                  then
                    begin
                      if pptra^.instruct = calc_neg
                        then
                          begin
                            pptr^ := pptra^.nextinst^;
                            Dispose(pptra^.nextinst);
                            Dispose(pptra);
                            restoreptr;
                          end;
                    end;
                if (pptr^.instruct = calc_neg) and (pptra^.instruct in
                  [calc_mul.. calc_div])
                  then
                    begin
                      help1 := endof(pptra^.nextinst);
                      if pptra^.nextinst^.instruct = calc_const
                        then
                          begin
                            pptra^.nextinst^.operand := -pptra^.nextinst^.operand;
                            pptr^ := pptra^;
                            Dispose(pptra);
                            restoreptr;
                          end
                        else
                          if help1^.nextinst^.instruct = calc_const
                            then
                              begin
                                help1^.nextinst^.operand := -help1^.nextinst^.operand;
                                pptr^ := pptra^;
                                Dispose(pptra);
                                restoreptr;
                              end;
                    end;
              end
            else  // jetzt werden Operationen mit Konstanten vereinfacht
              if pptr^.instruct in [calc_add..calc_pow]
                then
                  begin
                    simple(pptra);
                    simple(pptrb);
                    arg1 := pptra^.instruct = calc_const;
                    arg2 := pptrb^.instruct = calc_const;
                    if arg1 and arg2
                      then
                        begin
                          dummy := compute(pptr, pptrb, pptra);
                          pptr^.instruct := calc_const;
                          pptr^.operand := dummy;
                          pptr^.nextinst := pptrb^.nextinst;
                          Dispose(pptra);
                          Dispose(pptrb);
                          restoreptr;
                        end
                      else
                        if arg2
                          then
                            begin
                              if pptrb^.operand = 0.0
                                then
                                  begin
                                    if pptr^.instruct in [calc_mul.. calc_div, calc_pow]
                                      then
                                        setconst(0.0)
                                      else
                                        if pptr^.instruct = calc_add
                                          then
                                            pusha
                                          else
                                            if pptr^.instruct = calc_sub
                                              then
                                                begin
                                                  pptr^.instruct := calc_neg;
                                                  help1 := endof(pptra);
                                                  help1^.nextinst := pptrb^.nextinst;
                                                  Dispose(pptrb);
                                                  restoreptr;
                                                end;
                                  end
                                else
                                  if (pptrb^.operand = 1.0) and (pptr^.instruct in
                                [calc_mul, calc_pow])
                                    then
                                      begin
                                        if pptr^.instruct = calc_mul
                                          then pusha
                                          else setconst(1.0);
                                      end;
                            end
                          else
                            if arg1
                              then
                                begin
                                  if pptra^.operand = 0.0
                                    then
                                      begin
                                        if pptr^.instruct in [calc_mul, calc_pow]
                                          then
                                            begin
                                              if pptr^.instruct = calc_mul
                                                then dummy := 0.0
                                                else dummy := 1.0;
                                              pptr^.instruct := calc_const;
                                              pptr^.operand := dummy;
                                              help1 := pptrb;
                                              op := 1;
                                              repeat
                                                help2 := help1;
                                                if help1^.instruct in [calc_add.. calc_pow]
                                                  then op := Succ(op)
                                                  else
                                                    if help1^.instruct in [calc_const, calc_var]
                                                      then op := Pred(op);
                                                help1 := help1^.nextinst;
                                                Dispose(help2);
                                              until op = 0;
                                              pptr^.nextinst := help1;
                                              Dispose(pptra);
                                              restoreptr;
                                            end
                                          else
                                            if pptr^.instruct in [calc_add, calc_sub]
                                              then skipa;
                                      end
                                    else
                                      if (pptra^.operand = 1.0) and (pptr^.instruct in
                                    [calc_mul..calc_div, calc_pow])
                                        then skipa;
                                end;
                                  if (pptr^.instruct = calc_mul) and (pptra^.instruct in [calc_div, calc_dvd])
                                    then
                                      begin
                                        help1 := endof(pptra^.nextinst);
                                        if (help1^.nextinst^.instruct = calc_const)
                                          then
                                            if help1^.nextinst^.operand = 1.0
                                              then
                                                begin
                                                  pptr^ := pptra^;
                                                  Dispose(pptra);
                                                  Dispose(help1^.nextinst);
                                                  help1^.nextinst := pptrb;
                                                  restoreptr;
                                                end
                                              else
                                                if help1^.nextinst^.operand = -1.0 then
                                                  begin
                                                    pptr^.instruct := calc_neg;
                                                    Dispose(help1^.nextinst);
                                                    help1^.nextinst := pptrb;
                                                    restoreptr;
                                                  end;
                                      end;
                                  if pptr^.instruct in [calc_mul..calc_div]
                                    then
                                      begin            // Negationen vereinfachen
                                        if pptra^.instruct = calc_neg
                                          then
                                            if pptrb^.instruct = calc_neg
                                              then
                                                begin
                                                  pptr^.nextinst := pptra^.nextinst;
                                                  Dispose(pptra);
                                                  pptra := pptr^.nextinst;
                                                  help1 := endof(pptra);
                                                  help1^.nextinst := pptrb^.nextinst;
                                                  Dispose(pptrb);
                                                  restoreptr;
                                                end
                                              else
                                                begin
                                                  if ((pptrb^.instruct = calc_const) and (pptrb^.operand < 0.0))
                                                    then
                                                      begin
                                                        pptr^.nextinst := pptra^.nextinst;
                                                        Dispose(pptra);
                                                        pptrb^.operand := Abs(pptrb^.operand);
                                                        restoreptr;
                                                      end;
                                                end
                                          else
                                            if ((pptra^.instruct = calc_const) and (pptra^.operand < 0.0))
                                              then
                                                if pptrb^.instruct = calc_neg
                                                  then
                                                    begin
                                                      pptra^.nextinst := pptrb^.nextinst;
                                                      Dispose(pptrb);
                                                      pptra^.operand := Abs(pptra^.operand);
                                                      restoreptr;
                                                    end;
                                        if (pptra^.instruct = calc_const) and (pptra^.operand = -1.0)
                                          then
                                            begin
                                              pptr^.instruct := calc_neg;
                                              pptr^.nextinst := pptrb;
                                              Dispose(pptra);
                                              restoreptr;
                                            end;
                                        if ((pptrb^.instruct = calc_const) and (pptrb^.operand = -1.0) and
                                          (pptr^.instruct = calc_mul))
                                          then
                                            begin
                                              help1 := endof(pptra);
                                              help1^.nextinst := pptrb^.nextinst;
                                              pptr^.instruct := calc_neg;
                                              Dispose(pptrb);
                                              restoreptr;
                                            end;
                                      end;
                                  if (pptr^.instruct = calc_add) and (pptra^.instruct = calc_neg)
                                    then
                                      begin
                                        pptr^.instruct := calc_sub;
                                        pptr^.nextinst := pptra^.nextinst;
                                        Dispose(pptra);
                                        restoreptr;
                                      end;
                                  if (pptr^.instruct = calc_sub) and (pptra^.instruct = calc_neg)
                                    then
                                      begin
                                        pptr^.instruct := calc_add;
                                        pptr^.nextinst := pptra^.nextinst;
                                        Dispose(pptra);
                                        restoreptr;
                                      end;
                                            // difficult : the kommutativ law
                                  if (((pptr^.instruct = calc_mul) and (pptra^.instruct in
                                    [calc_mul..calc_div])) or ((pptr^.instruct = calc_add) and
                                    (pptra^.instruct in [calc_add, calc_sub]))) and
                                    (pptrb^.instruct = calc_const)
                                    then
                                      begin
                                        help1 := endof(pptra^.nextinst);
                                        help2 := endof(help1^.nextinst);
                                        if pptra^.instruct in [calc_mul, calc_add]
                                          then
                                            begin
                                              if help1^.nextinst^.instruct = calc_const
                                                then
                                                  begin
                                                    help2^.nextinst := pptra^.nextinst;
                                                    pptra^.nextinst := pptrb;
                                                    help2 := pptrb^.nextinst;
                                                    pptrb^.nextinst := help1^.nextinst;
                                                    help1^.nextinst := help2;
                                                  end
                                                else
                                                  if pptra^.nextinst^.instruct = calc_const
                                                    then
                                                      begin
                                                        help2^.nextinst := pptrb^.nextinst;
                                                        pptrb^.nextinst := help1^.nextinst;
                                                        help1^.nextinst := pptrb;
                                                      end;
                                            end
                                          else
                                            begin
                                              if help1^.nextinst^.instruct = calc_const
                                                then
                                                  begin
                                                    helpinstruct := pptr^.instruct;
                                                    pptr^.instruct := pptra^.instruct;
                                                    pptra^.instruct := helpinstruct;
                                                    pptr^.nextinst := pptra^.nextinst;
                                                    pptra^.nextinst := help1^.nextinst;
                                                    help1^.nextinst := pptra;
                                                  end;
                                            end;
                                        restoreptr;
                                        simple(pptra);
                                        simple(pptrb);
                                      end
                                    else
                                      if (((pptr^.instruct = calc_mul) and (pptrb^.instruct in
                                         [calc_mul..calc_div])) or ((pptr^.instruct = calc_add) and
                                         (pptrb^.instruct in [calc_add, calc_sub]))) and
                                         (pptra^.instruct = calc_const)
                                        then
                                          begin
                                            help1 := endof(pptrb^.nextinst);
                                            help2 := endof(help1^.nextinst);
                                            if pptrb^.instruct in [calc_add, calc_mul]
                                              then
                                                begin
                                                  if pptrb^.nextinst^.instruct = calc_const
                                                    then
                                                      begin
                                                        pptr^.nextinst := help1^.nextinst;
                                                        help1^.nextinst := pptra;
                                                        pptra^.nextinst := help2^.nextinst;
                                                        help2^.nextinst := pptrb;
                                                      end
                                                    else
                                                      if help1^.nextinst^.instruct = calc_const
                                                        then
                                                          begin
                                                            pptr^.nextinst := pptrb^.nextinst;
                                                            pptra^.nextinst := help1^.nextinst;
                                                            help1^.nextinst := pptrb;
                                                            pptrb^.nextinst := pptra;
                                                          end;
                                              end
                                            else
                                              begin
                                                if help1^.nextinst^.instruct = calc_const
                                                  then
                                                    begin
                                                      helpinstruct := pptr^.instruct;
                                                      pptr^.instruct := pptrb^.instruct;
                                                      pptrb^.instruct := helpinstruct;
                                                      pptr^.nextinst := pptrb^.nextinst;
                                                      pptrb^.nextinst := pptra;
                                                      pptra^.nextinst := help1^.nextinst;
                                                      help1^.nextinst := pptrb;
                                                    end;
                                              end;
                                            restoreptr;
                                            simple(pptra);
                                            simple(pptrb);
                                          end;
                                            if pptr^.instruct = calc_mul then
                                            begin
                                              if pptra^.instruct = calc_pow then
                                              begin
                                                help2 := pptra^.nextinst;
                                                help1 := endof(help2);
                                                help1 := help1^.nextinst;
                                                if (help2^.instruct = calc_const) and equal(help1, pptrb) then
                                                begin
                                                  help2^.operand := help2^.operand + 1.0;
                                                  pptr^ := pptra^;
                                                  Dispose(pptra);
                                                  help2^.nextinst := pptrb;
                                                  erase_entry;
                                                  restoreptr;
                                                end;
                                              end;
                                            end;
                                  if pptr^.instruct in [calc_add, calc_sub, calc_dvd, calc_div, calc_mul]
                                    then
                                      if equal(pptra, pptrb)
                                        then           // sind die Operanden gleich ?
                                          begin
                                            case pptr^.instruct of
                                              calc_add: begin
                                                          pptr^.instruct := calc_mul;
                                                          help2 := endof(pptra);
                                                          help1 := pptra^.nextinst;
                                                          pptra^.instruct := calc_const;
                                                          pptra^.operand := 2.0;
                                                          pptra^.nextinst := help2^.nextinst;
                                                          erase_entry;
                                                          restoreptr;
                                                        end;
                                              calc_sub: setconst(0.0);
                                              calc_div,
                                              calc_dvd: setconst(1.0);
                                              calc_mul: begin
                                                          pptr^.instruct := calc_pow;
                                                          help2 := endof(pptra);
                                                          help1 := pptra^.nextinst;
                                                          pptra^.nextinst := help2^.nextinst;
                                                          pptra^.instruct := calc_const;
                                                          pptra^.operand := 2.0;
                                                          erase_entry;
                                                          restoreptr;
                                                        end;
                                            end; // case
                                          end;  // if equal
                  end; // pptr^.instruct in [calc_add..calc_pow]
        end; // pptr <> nil
  end; // Simple

begin // CalcSimplify
  if pptr <> nil
   then
     begin
       SimpleError := False;
       CalcResult := True;
       Invert(pptr);
       pptr1 := pptr^.nextinst;
       Simple(pptr1);
       if SimpleError
         then
           begin
             KillExpression(pptr);
             CalcResult := False;
           end
         else
           invert(pptr);
     end
   else
    calcresult := False;
end; // CalcSimlify


function CalcDerivation(pptr: calc_prog; vartab: calc_vartab;
  nach: calc_idstr): calc_prog;

const
  pi_durch_180 = 1.745329252e-2;

var
  pptrstart, pptr1 : calc_prog;
  ok               : boolean;


  procedure UpperCase(var varid: calc_idstr);

  var
    i: integer;

  begin
    for i := 1 to Length(varid) do
      varid[i] := Upcase(varid[i]);
  end;


  procedure NewConst(x: calc_operand);

  var
    pptr : calc_prog;

  begin
    try
      begin
        New(pptr);
        pptr^.instruct := calc_const;
        pptr^.operand := x;
        pptr^.nextinst := pptrstart^.nextinst;
        pptrstart^.nextinst := pptr;
      end
    except
        calcresult := False;
    end;
  end;


  procedure NewOp(id: calc_symbols);

  var
    pptr: calc_prog;

  begin
    try
      begin
        New(pptr);
        pptr^.instruct := id;
        pptr^.nextinst := pptrstart^.nextinst;
        pptrstart^.nextinst := pptr;
      end
    except
      calcresult := False;
    end;
  end;


  procedure push(pptr: calc_prog);

  var
    pptr1: calc_prog;
    op: integer;

  begin
    op := 1;
    repeat
      try
        begin
          if pptr^.instruct in [calc_add..calc_pow]
            then op := op + 1
            else
              if not (pptr^.instruct in [calc_neg, calc_sqr..calc_fak])
                then
                  op := op - 1;
          New(pptr1);
          pptr1^ := pptr^;
          pptr1^.nextinst := pptrstart^.nextinst;
          pptrstart^.nextinst := pptr1;
          pptr := pptr^.nextinst;
        end
      except
        calcresult := False
      end;
    until (op = 0) or not calcresult;
  end; // Push


  procedure derive(pptr : calc_prog);

  var
    pptra, pptrb : calc_prog;

  begin
    if calcresult
      then
        begin
          pptra := pptr^.nextinst;
          if (pptra <> nil)
            then
              begin
                pptrb := endof(pptra);
                pptrb := pptrb^.nextinst;
              end;
          case pptr^.instruct of
            calc_neg:  begin
                         newop(calc_neg);
                         derive(pptra);
                       end;
            calc_const,
            calc_div..calc_kgv,
            calc_fak : begin
                         newconst(0.0);
                       end;
            calc_var : begin
                         if nach = vartab^[pptr^.varindex].varid
                           then newconst(1.0)
                           else newconst(0.0);
                       end;
            calc_add : begin
                         if calc_const in [pptra^.instruct, pptrb^.instruct]
                           then
                             if pptra^.instruct = calc_const
                               then derive(pptrb)
                               else derive(pptra)
                           else
                             begin
                               newop(calc_add);
                               derive(pptra);
                               derive(pptrb);
                             end;
                       end;
            calc_sub : begin
                         if calc_const in [pptra^.instruct, pptrb^.instruct]
                           then
                             if pptra^.instruct = calc_const
                               then derive(pptrb)
                               else
                                 begin
                                   newop(calc_neg);
                                   derive(pptra);
                                 end
                           else
                             begin
                               newop(calc_sub);
                               derive(pptra);
                               derive(pptrb);
                             end;
                       end;
            calc_mul : begin
                         if calc_const in [pptra^.instruct, pptrb^.instruct]
                           then
                             if pptra^.instruct = calc_const
                               then
                                 begin
                                   newop(calc_mul);
                                   push(pptra);
                                   derive(pptrb);
                                 end
                               else
                                 begin
                                   newop(calc_mul);
                                   push(pptrb);
                                   derive(pptra);
                                 end
                           else
                             begin
                               newop(calc_add);
                               newop(calc_mul);
                               derive(pptra);
                               push(pptrb);
                               newop(calc_mul);
                               push(pptra);
                               derive(pptrb);
                             end;
                       end;
            calc_dvd: begin
                        if pptra^.instruct = calc_const
                          then
                            begin
                              newop(calc_dvd);
                              push(pptra);
                              derive(pptrb);
                            end
                          else
                            begin
                              newop(calc_dvd);
                              newop(calc_sqr);
                              push(pptra);
                              newop(calc_sub);
                              newop(calc_mul);
                              derive(pptra);
                              push(pptrb);
                              newop(calc_mul);
                              push(pptra);
                              derive(pptrb);
                            end;
                      end;
            calc_pow : begin
                         if (pptrb^.instruct = calc_const) and (pptrb^.operand < 0.0)
                           then calcresult := False
                           else
                             begin
                               ok := False;
                               case pptra^.instruct of
                                 calc_const : ok := True;
                                 calc_var   : ok := nach <> vartab^[pptra^.varindex].varid
                               end;
                               if ok then
                                 begin
                                   newop(calc_mul);
                                   newop(calc_mul);
                                   newop(calc_pow);
                                   newop(calc_sub);
                                   newconst(1.0);
                                   push(pptra);
                                   push(pptrb);
                                   push(pptra);
                                   derive(pptrb);
                                 end
                               else
                                 begin
                                   newop(calc_mul);
                                   newop(calc_pow);
                                   push(pptra);
                                   push(pptrb);
                                   newop(calc_add);
                                   newop(calc_dvd);
                                   push(pptrb);
                                   newop(calc_mul);
                                   push(pptra);
                                   derive(pptrb);
                                   newop(calc_mul);
                                   derive(pptra);
                                   newop(calc_ln);
                                   push(pptrb);
                                 end;
                            end;
                       end;
            calc_abs : begin
                         newop(calc_mul);
                         newop(calc_sign);       // calc_sig ???
                         push(pptra);
                         derive(pptra);
                       end;
            calc_int,
            calc_sign: newconst(0.0);
            calc_sqr : begin
                         newop(calc_mul);
                         newop(calc_mul);
                         push(pptra);
                         derive(pptra);
                         newconst(2.0);
                       end;
            calc_sqrt: begin
                         newop(calc_dvd);
                         newop(calc_mul);
                         newconst(2.0);
                         newop(calc_sqrt);
                         push(pptra);
                         derive(pptra);
                       end;
            calc_exp : begin
                         newop(calc_mul);
                         newop(calc_exp);
                         push(pptra);
                         derive(pptra);
                       end;
            calc_ln : begin
                        newop(calc_dvd);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_lg : begin
                        newop(calc_dvd);
                        newop(calc_mul);
                        newop(calc_ln);
                        newconst(10.0);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_ld : begin
                        newop(calc_dvd);
                        newop(calc_mul);
                        newop(calc_ln);
                        newconst(2.0);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_sin : begin
                        newop(calc_mul);
                        newop(calc_cos);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_cos : begin
                        newop(calc_mul);
                        newop(calc_neg);
                        newop(calc_sin);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_tan : begin
                        newop(calc_dvd);
                        newop(calc_sqr);
                        newop(calc_cos);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_cot: begin
                        newop(calc_neg);
                        newop(calc_dvd);
                        newop(calc_sqr);
                        newop(calc_sin);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_arcsin : begin
                        newop(calc_dvd);
                        newop(calc_sqrt);
                        newop(calc_sub);
                        newop(calc_sqr);
                        push(pptra);
                        newconst(1.0);
                        derive(pptra);
                      end;
            calc_arccos : begin
                        newop(calc_neg);
                        newop(calc_dvd);
                        newop(calc_sqrt);
                        newop(calc_sub);
                        newop(calc_sqr);
                        push(pptra);
                        newconst(1.0);
                        derive(pptra);
                      end;
            calc_arctan: begin
                        newop(calc_dvd);
                        newop(calc_add);
                        newop(calc_sqr);
                        push(pptra);
                        newconst(1.0);
                        derive(pptra);
                      end;
            calc_arccot : begin
                        newop(calc_neg);
                        newop(calc_dvd);
                        newop(calc_add);
                        newop(calc_sqr);
                        push(pptra);
                        newconst(1.0);
                        derive(pptra);
                      end;
            calc_sinh : begin
                        newop(calc_mul);
                        newop(calc_cosh);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_cosh : begin
                        newop(calc_mul);
                        newop(calc_sinh);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_tanh : begin
                        newop(calc_dvd);
                        newop(calc_sqr);
                        newop(calc_cosh);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_coth : begin
                        newop(calc_neg);
                        newop(calc_dvd);
                        newop(calc_sqr);
                        newop(calc_sinh);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_arcsinh : begin
                        newop(calc_dvd);
                        newop(calc_sqrt);
                        newop(calc_add);
                        newconst(1.0);
                        newop(calc_sqr);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_arccosh : begin
                        newop(calc_dvd);
                        newop(calc_sqrt);
                        newop(calc_sub);
                        newconst(1.0);
                        newop(calc_sqr);
                        push(pptra);
                        derive(pptra);
                      end;
            calc_arctanh,
            calc_arccoth : begin
                        newop(calc_dvd);
                        newop(calc_sub);
                        newop(calc_sqr);
                        push(pptra);
                        newconst(1.0);
                        derive(pptra);
                      end;
            calc_deg : begin
                        newop(calc_dvd);
                        newconst(pi_durch_180);
                        derive(pptra);
                      end;
            calc_rad: begin
                        newop(calc_mul);
                        newconst(pi_durch_180);
                        derive(pptra);
                      end;
            else      calcresult := False
          end; // case
      end; // if CalcResult
  end; // Derive

begin // CalcDerivation
  if pptr <> nil
    then
      begin
        uppercase(nach);
        invert(pptr);
        pptr1 := pptr;
        New(pptrstart);
        pptrstart^.nextinst := nil;
        pptr := pptr^.nextinst;
        calcresult := True;
        derive(pptr);
        if calcresult then
        begin
          calcsimplify(pptrstart);
          Result := pptrstart;
        end
        else
        begin
          killexpression(pptrstart);
          Result := nil;
          CalcError(4, 'Derivate of function not known');
        end;
        invert(pptr1);
      end
    else
      begin
        Result := nil;
        CalcResult := False;
      end;
end; // CalcDerivation


procedure CompileExpression(Expr: Calc_String; var VarTable: Calc_VarTab;
  var ExprPtr: Calc_Prog);

var
  VarTabFlag, ParsError, EndOfExpr : boolean;
  ch                               : string[1];    // akt. Zeichen aus String
  LastPos, StrPos                  : integer;      // zaehlt String-Position mit
  TempIdent, Ident                 : Calc_String;  // enth. aktuellen Bezeichner
  Symbol,                                          // akt. Symbol des Bezeichners
  LastSymbol                       : Calc_Symbols; // vorheriges Symbol
  Number                           : Calc_Operand; // akt. Zahl/Index aus String
  ProgPtr                          : Calc_Prog;


  procedure Error(ErrPos: integer; ErrMsg: Calc_String);

  begin
    if not ParsError
      then
        begin
          WriteLn;
          Write('*** ', ' Error compiling this expression: ');
          ClrEol;
          WriteLn;
          Write(' ', Expr);
          ClrEol;
          Writeln;
          Write(' ': ErrPos, '^');
          ClrEol;
          Writeln;
          Write(ErrMsg, '!');
          Clreol;
          Writeln;
          ClrEol;
          WriteLn;
          ClrEol;
        end;
    ParsError := True;
    Symbol := Calc_Err;
  end;


  procedure Add_To_Queue(op: Calc_Symbols; x: Calc_Operand);
  { add next operand to que }

  var
    UPN_Entry: Calc_Prog;

  begin
    try  // is enough memory available?
      begin
        New(UPN_Entry);
        with UPN_Entry^ do
          begin
            NextInst := nil;
            Instruct := op;
            case op of
              Calc_Var   : VarIndex := Trunc(x);
              Calc_Const : Operand  := x
            end;
          end;
        ProgPtr^.NextInst := UPN_Entry;
        ProgPtr := ProgPtr^.NextInst;
      end
    except
      Error(1, 'not enough free memory');
    end;
  end;


  procedure GetSymbol;
  { get next symbol from string }


    procedure GetChar;
    { get next character from string }

    begin
      ch := ' ';
      StrPos := Succ(StrPos);
      EndOfExpr := (StrPos > Length(Expr));
      if not EndOfExpr
        then ch := UpCase(Expr[StrPos]);
    end;


    procedure GetNumber;
    { Get the next number from string. The Turbo-Pascal val-procedure wants
      to see only valid characters of a floating point number, NumberEnd
      points to the first invalid character. Hence:                          }

    var
      NumberStr: Calc_String;
      NumberEnd, posi: integer;

    begin
      NumberStr := Copy(Expr, StrPos, 255); // everything from first figure
      NumberStr := NumberStr + '     ';
      posi := 1;
      while (not (numberstr[posi] in ['e', 'E'])) and (posi < length(numberstr)) do
        posi := succ(posi);
      if numberstr[posi] in ['e', 'E']
        then
          begin
            if numberstr[posi + 1] in ['+', '-']
              then posi := succ(posi);
            if not (numberstr[posi + 1] in ['0'..'9'])
              then error(StrPos + posi, 'incomplete expression')
              else
                if ((numberstr[posi + 1] = '3') and (numberstr[posi + 2] in
                    ['7'..'9'])) or ((numberstr[posi + 1] > '3') and
                    (numberstr[posi + 2] in ['0'..'9']))
                  then error(StrPos + posi, 'number out of range');
          end;
      Val(NumberStr, Number, NumberEnd);
      if NumberEnd > 0
        then                        // invalid character
          begin                     // number not at the end of expression
           StrPos := StrPos + NumberEnd - 2;
           NumberStr := Copy(NumberStr, 1, Pred(NumberEnd));
           Val(NumberStr, Number, NumberEnd);
         end
       else                         // worked, at the end of expression
        StrPos := Length(Expr);
    end;


    procedure searchsymtab;
    { look up symbol in symbol table }

    var
      symok: boolean;

    begin
      symok := False;
      symbol := calc_err;
      while (symbol < calc_end) and not symok do
        begin
          symbol := Succ(symbol);
          symok := (ident = calc_ids[symbol]);
        end;
      if not symok
        then symbol := calc_err;
    end;

  begin  // GetSymbol
    LastPos := StrPos;
    LastSymbol := Symbol;
    Ident := '';
    while (ch = ' ') and not EndOfExpr do
      GetChar;
    case ch[1] of
      'A'..'Z': repeat      // Name of operator, function or variable
                  Ident := Concat(Ident, ch);
                  GetChar
                until not (ch[1] in ['A'..'Z', '0'..'9']);
      '0'..'9', '.':
                begin
                  GetNumber;
                  Ident := 'CONST';
                  GetChar;
                end
      else      begin          // +, -, etc.
                  Ident := ch;
                  GetChar;
                end
    end;
    SearchSymtab;
    if Symbol = Calc_Err
      then  // operator not identified
        if Ident[1] in ['A'..'Z']
          then Symbol := Calc_Var  // use as variable
          else
            if Ident <> ' '
              then Error(LastPos, 'unknown symbol');
    if Symbol = Calc_Var
      then
        if LastSymbol <> Calc_Var
          then
            begin
             Number := SearchVarTab(VarTable, Ident);
             if CalcDecMod and (Number = 0.0)
               then  // refuse unidentified string
                 Error(LastPos, 'unknown symbol')
               else  // or accept it as new var
                 if Number = 0.0
                   then
                     begin
                       Number := AddToVarTab(VarTable, Ident);
                       if Number < 0.0
                         then Error(LastPos, 'too many variables');
                     end;
            end
      else
        if not EndOfExpr
          then Error(LastPos, 'operator expected');
  end; // GetSymbol


  procedure Expression;
  { expression contains several parts connected by operators }

  var
    ExprOp: Calc_Symbols;

    procedure Term;
    { expression contains several factors }

    var
      TermOp: Calc_Symbols;

      procedure Factor(fparen: boolean);
      { Factors can be variables, constants, functions or an expression in
        brackets, and all of them can be raised to a power.
        The parameter 'fparen' indicates whether brackets are for an expression
        or a function. In the latter case, the function needs to be evaluated
        before raising to power                                                }

       var
        FacOp: Calc_Symbols;

      begin
        if Symbol <> Calc_Err
          then
            begin
              case Symbol of
                Calc_Var   : begin
                               Add_To_Queue(Calc_Var, Number);
                               GetSymbol;
                             end;
                Calc_Const : begin
                               Add_To_Queue(Calc_Const, Number);
                               GetSymbol;
                             end;
                Calc_Pi    : begin
                               Add_To_Queue(Calc_Const, Const_pi);
                               GetSymbol;
                             end;
                Calc_E     : begin
                               Add_To_Queue(Calc_Const, Const_e);
                               GetSymbol;
                             end;
                Calc_lp    : begin                 // expression in brackets
                               GetSymbol;
                               Expression;
                               if (Symbol <> Calc_rp)
                                 then Error(StrPos - Ord(ch[1] > ' '),
                                            Concat(Calc_Ids[Calc_rp], ' expected'));
                               GetSymbol;
                             end;
                Calc_Pow   : ;      // power, dealt with later
                Calc_Sqr..
                Calc_Fak   : begin  // funktion
                               FacOp := Symbol;
                               GetSymbol;
                               if (Symbol = Calc_lp)
                                 then Factor(True)
                                 else Error(LastPos, Concat(Calc_Ids[Calc_lp], ' expected'));
                               Add_To_Queue(FacOp, 0);
                             end
                else         Error(LastPos, 'here unexpected')
              end; // CASE
              if not fparen
                then // brackets notfor function, therefore..
                  if Symbol = Calc_Pow
                    then  // power
                      if LastSymbol in [Calc_Const, Calc_Var, Calc_rp, Calc_PI, Calc_e]
                        then
                          begin  // evaluate
                            GetSymbol;
                            Factor(False);
                            Add_To_Queue(Calc_Pow, 0);
                          end
                        else      // power here not possible
                          Error(Pred(StrPos), 'here unexpected');
            end // IF Symbol <> Calc_Err
          else
            Error(StrPos - Ord(EndOfExpr), ' incomplete expression');
      end; // Faktor

    begin  // Term
      Factor(False);
      if Symbol in [Calc_Mul..Calc_kgV]
        then  // term contains several factors
          begin
            TermOp := Symbol;
            GetSymbol;
            Term;
            Add_To_Queue(TermOp, 0);
          end;
    end; // Term

  begin  // Expression
    if Symbol in [Calc_Add..Calc_Sub]
      then // expression starts with + or -
        begin
          ExprOp := Symbol;
          GetSymbol;
          Term;
          if ExprOp = Calc_Sub
            then Add_To_Queue(Calc_Neg, 0);      // negate everything
        end
      else
        Term;
    while (Symbol in [Calc_Add..Calc_Sub]) do
      begin                              // additional terms follow
        ExprOp := Symbol;
        GetSymbol;
        Term;
        Add_To_Queue(ExprOp, 0);
      end;
  end; // Expression

begin  // CompileExpression
  Symbol := Calc_Err;
  ParsError := False;
  EndOfExpr := False;
  StrPos := 0;
  ch := ' ';
  VarTabFlag := (VarTable = nil);
  if VarTabFlag
    then VarTable := NewVarTab;
  New(ExprPtr);
  ExprPtr^.NextInst := nil;
  ProgPtr := ExprPtr;
  GetSymbol;
  Expression;
  if Symbol <> Calc_EOE
    then Error(LastPos, '";" ' + ' expected');
  CalcResult := not ParsError;
  if ParsError
    then
      begin
        KillExpression(ExprPtr);
        if VarTabFlag
          then KillVarTab(VarTable);
      end;
end; // CompileExpression


function CalcExpression(ExprPtr: Calc_Prog; VarTable: Calc_VarTab): Calc_Operand;
{ evaluate an RPN-expression  }

const
  StackSize = 50;

var
  x        : Calc_Operand;
  StackPtr : integer;
  Stack    : array[1..StackSize] of Calc_Operand;


  procedure Push;                        // pushes number onto the stack

  begin
    Stack[StackPtr] := x;
    StackPtr := Succ(StackPtr);
  end;


  function Pop: Calc_Operand;            // gets a number from stack

  begin
    StackPtr := Pred(StackPtr);
    Result := Stack[StackPtr];
  end;

begin  // CalcExpression
  CalcResult := True;
  if (ExprPtr <> nil) and (VarTable <> nil)
    then
      begin
        ExprPtr := ExprPtr^.NextInst;
        StackPtr := 1;
        x := 0.0;
        while ExprPtr <> nil do
          begin
            with ExprPtr^ do
              case Instruct of
                Calc_Const : begin
                               Push;
                               x := Operand;
                             end;
                Calc_Var   : begin
                               Push;
                               x := VarTable^[VarIndex].Value;
                             end;
                else         begin
                              case Instruct of
                                Calc_Neg   : x := -x;
                                Calc_Add   : x := Pop + x;
                                Calc_Sub   : x := Pop - x;
                                Calc_Mul   : x := Pop * x;
                                Calc_Dvd   : if x <> 0.0
                                               then x := Pop / x
                                               else CalcError(2, '');
                                Calc_Div   : if Trunc(x) <> 0
                                             then x := Trunc(Pop) div Trunc(x)
                                             else CalcError(2, '');
                                Calc_Mod   : if Trunc(x) <> 0
                                               then x := Trunc(Pop) mod Trunc(x)
                                               else CalcError(2, '');
                                Calc_ggT   : x := GCD(Trunc(Pop), Trunc(x));
                                Calc_kgV   : x := SCM(Trunc(Pop), Trunc(x));
                                Calc_Pow   : x := pot(Pop, x);
                                Calc_Sqr   : x := Sqr(x);
                                Calc_Sqrt  : if x >= 0.0
                                               then x := Sqrt(x)
                                               else CalcError(3, 'Sqrt(x): x <= 0');
                                Calc_Exp   : x := Exp(x);
                                Calc_Ln    : if x > 0.0
                                               then x := Ln(x)
                                               else CalcError(3, 'ln(x): x <= 0');
                                Calc_Lg    : x := log(x, 10);
                                Calc_Ld    : x := log(x, 2);
                                Calc_Sin   : x := Sin(x);
                                Calc_Cos   : x := Cos(x);
                                Calc_Tan   : x := tan(x);
                                Calc_Cot   : x := cot(x);
                                Calc_ArcSin: x := arcsin(x);
                                Calc_ArcCos: x := arccos(x);
                                Calc_ArcTan: x := ArcTan(x);
                                Calc_ArcCot: x := arccot(x);
                                Calc_Sinh  : x := sinh(x);
                                Calc_Cosh  : x := cosh(x);
                                Calc_Tanh  : x := tanh(x);
                                Calc_Coth  : x := coth(x);
                                Calc_ArcSinh: x := arsinh(x);
                                Calc_ArcCosh: x := arcosh(x);
                                Calc_ArcTanh: x := artanh(x);
                                Calc_ArcCoth: x := arcoth(x);
                                Calc_Abs   : x := abs(x);
                                Calc_Deg   : x := grad(x);
                                Calc_Rad   : x := rad(x);
                                Calc_Rez   : if x <> 0
                                               then x := 1 / x
                                               else CalcError(3, '1 / 0');
                                Calc_Fak   : x := fak(Round(x));
                                Calc_Int   : x := Int(x);
                                Calc_Sign  : x := Signum(x);
                                else         CalcError(0, 'Function not known')
                               end; // CASE
                             end; // ELSE
              end; // CASE
            ExprPtr := ExprPtr^.NextInst;
            if StackPtr > StackSize then CalcError(0, 'Stack overflow');
            if not CalcResult then ExprPtr := nil;
          end; // WHILE
        Result := x;
      end // THEN
  else
    CalcError(0, 'Function not known');
end;


procedure CalcAOS(pptr: Calc_Prog; VarTable: Calc_VarTab);

var
  Value   : string[50];
  len     : byte ABSOLUTE Value;
  key     : char;
  dummy   : calc_operand;
  pptr1   : calc_prog;


  procedure writeaos(pptr: calc_prog);

  var
    pptra, pptrb : calc_prog;
    paren        : boolean;

  begin
    if (pptr <> nil)
      then
        begin
          if keypressed
            then
              begin
                key := ReadKey;
                if key = ^s then Key := readkey;
              end;
          pptra := pptr^.nextinst;
          if pptra <> nil
            then
              begin
                pptrb := endof(pptra);
                pptrb := pptrb^.nextinst;
              end;
          case pptr^.instruct of
            calc_const: begin
                          dummy := pptr^.operand;
                          Str(dummy: 0: 10, Value);
                          while Value[len] = '0' do
                            len := Pred(len);
                          if Value[len] = '.' then len := Pred(len);
                          paren := dummy < 0.0;
                          if paren then Write('(');
                          Write(Value);
                          if paren then Write(')');
                        end;
            calc_var  : Write(vartable^[pptr^.varindex].varid);
            calc_add..
            calc_pow  : begin
                          paren := (pptr^.instruct in [calc_mul..calc_pow]) and
                                   (pptrb^.instruct in [calc_add, calc_sub]);
                          paren := paren or (pptr^.instruct = calc_pow) and
                                   (pptrb^.instruct in [calc_add..calc_pow]);
                          if paren then Write('(');
                          writeaos(pptrb);
                          if paren then Write(')');
                          if pptr^.instruct in [calc_div..Calc_Kgv] then Write(' ');
                          Write(calc_ids[pptr^.instruct]);
                          if pptr^.instruct in [calc_div..Calc_Kgv] then Write(' ');
                          paren := (pptr^.instruct in [calc_mul..calc_pow]) and
                                   (pptra^.instruct in [calc_add, calc_sub]);
                          paren := paren or ((pptr^.instruct in [calc_dvd..calc_pow])
                                   and (pptra^.instruct in [calc_add..calc_pow]));
                          paren := paren or ((pptr^.instruct = calc_sub) and
                                   (pptra^.instruct in [calc_add, calc_sub]));
                          if paren then Write('(');
                          writeaos(pptra);
                          if paren then Write(')');
                        end;
            calc_neg  : begin
                          Write('(-');
                          paren := pptra^.instruct in [calc_add..calc_pow];
                          if paren then Write('(');
                          writeaos(pptra);
                          if paren then Write(')');
                          Write(')');
                        end;
            calc_sqr..
            calc_fak  : begin
                          Write(calc_ids[pptr^.instruct], '(');
                          writeaos(pptra);
                          Write(')');
                        end
            end; // case
        end;     // then
  end;           // WriteAOS

begin            // CalcAOS
  if pptr <> nil
    then
      begin
        pptr1 := pptr;
        invert(pptr);
        pptr := pptr^.nextinst;
        writeaos(pptr);
        invert(pptr1);
        Writeln;
      end
    else
      Writeln('Function not defined');
end; // CalcAOS


procedure HelpFormula;

begin
    writeln('The formula is entered in Pascal-syntax and must end with a semicolon.');
    writeln;
    writeln('The compiler ''knows'' the following constants and functions, which ');
    writeln('can not be redefined:');
    writeln('Constants: e, pi                        Basic operators: +, -, *, /, ^');
    writeln('Integer: div, mod, ggt, kgv             Logarithms: ln, lg, ld, exp');
    writeln('sin, cos, tan, cot and the equivalent hyperbolic and arcus functions');
    writeln('Various Functions: abs, deg, rad, fak, sgn');
    writeln;
end;



end. // Unit Calc
