PROGRAM ReadTest;

USES ReadExt;

VAR s  : STRING;
    c  : CHAR;
    i  : INTEGER;
    r  : real;
    w  : WORD;
    b  : BYTE;
    SI : SHORTINT;
    li : LONGINT;

BEGIN
    Write('Geben sie einn Buchstaben ein: ');
    c := ReadChar;
    Writeln;
    Writeln('Das Zeichen war ', c);
    Writeln;
    Write('Geben Sie eine Zeichenkette ein: ');
    s := ReadStr(10);
    Writeln;
    Writeln('Die Zeichenkette war ', s);
    Writeln;
    Write('Geben Sie eine Integer-Zahl ein: ');
    i := ReadInt(5);
    Writeln;
    Writeln('Die Zahl war ', i);
    Writeln('Geben Sie eine Real-Zahl ein: ');
    r := ReadReal(12, 6);
    Writeln;
    Writeln('Die Zahl war ', r:12:6);
    Writeln;
    Write('Geben Sie eine Word-Zahl ein: ');
    w := ReadWord(5);
    Writeln;
    Writeln('Die Zahl war ', w);
    Writeln;
    Writeln('Geben Sie eine Byte-Zahl ein: ');
    b := ReadByte(3);
    Writeln;
    Writeln('Die Zahl war ', b);
    Writeln;
    Write('Geben Sie eine ShortInt-Zahl ein: ');
    SI := ReadShortInt(4);
    Writeln;
    Writeln('Die Zahl war ', SI);
    Writeln;
    Write('Geben Sie eine LongInt-Zahl ein: ');
    li := ReadLongInt(10);
    Writeln;
    Writeln('Die Zahl war ', li);
    Writeln;
END.

