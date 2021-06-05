unit Haeufig;

interface

uses crt, graph, Graphik, Print_Gr, Printer, Mathfunc, Windows,
     ReadExt, Vektor, Matrix;

const MaxKlassen   = 20;
      HaeufigError : boolean = false;

type Verteilung = record
                       AnzKlassen    : byte;
                       AnzDaten      : word;
                       UntereGrenze,
                       ObereGrenze,
                       KlassenBreite : float;
                       Anzahl        : array [1..MaxKlassen] of word;
                       Darueber,
                       Darunter      : word;
                  end;

{ *************** Grundfunktionen auf den Datentyp Verteilung ************* }

procedure Klassifiziere (Daten : VektorTyp; var Frequenzen : Verteilung);

procedure WriteVerteilung (MedStr : string; Frequenzen : Verteilung);

procedure ReadVerteilung (MedStr : string; var Frequenzen : Verteilung);

function AnzahlKlassen (Frequenzen : Verteilung) : byte;

function Haeufigkeit (Frequenzen : Verteilung; KlNr : byte) : word;

function GetLower (Frequenzen : Verteilung; KlNr : byte) : float;

function GetUpper (Frequenzen : Verteilung; KlNr : byte) : float;

function GetAverage (Frequenzen : Verteilung; KlNr : byte) : float;

function AnzDarueber (Frequenzen : Verteilung) : word;

function AnzDarunter (Frequenzen : Verteilung) : word;

function Breite (Frequenzen : Verteilung) : float;

function AnzahlDaten (Frequenzen : Verteilung) : word;

{ *************** Statistische Funktionen mit Verteilungen **************** }

function ArithmetischesMittel (Frequenzen : Verteilung) : float;

function Varianz (Frequenzen : Verteilung) : float;

function GeometrischesMittel (Frequenzen : Verteilung) : float;

function HarmonischesMittel (Frequenzen : Verteilung) : float;

function ShepphardKorrigierteVarianz (Frequenzen : Verteilung) : float;

function Quantile (Frequenzen : Verteilung; QuantilWert : float) : float;

function Entropie (Frequenzen : Verteilung) : float;

procedure Beschreibe (Frequenzen : Verteilung);

procedure KumulPlot (Frequenzen : Verteilung; var xVektor, yVektor : VektorTyp);

function LorenzMuenzner (Frequenzen : Verteilung; var xKoord, yKoord : VektorTyp) : float;

procedure SaeulenDiagramm (Frequenzen : Verteilung; NameStr, LizStr : string);

procedure ZeichneTorte (Frequenzen : Verteilung; NameStr, LizStr : string);

{ ************************************************************************* }

implementation

var ch : char;

procedure Klassifiziere (Daten : VektorTyp; var Frequenzen : Verteilung);

{$IFDEF Deutsch}
const Meldung1 = ' Anzahl der Klassen (max. 20, empf. ';
      Meldung2 = ' Minimum der Daten ist: ';
      Meldung3 = ' Untere Grenze der Klassen-Bildung: ';
      Meldung4 = ' Maximum der Daten ist: ';
      Meldung5 = ' Obere Grenze der Klassen-Bildung: ';
      Meldung6 = ' Anzahl der Daten: ';
{$ENDIF}
{$IFDEF Englisch}
const Meldung1 = ' Number of classes (max. 20, recomended ';
      Meldung2 = ' Minimum of data: ';
      Meldung3 = ' Lower border of class formation: ';
      Meldung4 = ' Maximum of data: ';
      Meldung5 = ' Higher border of class formation: ';
      Meldung6 = ' Number of data: ';
{$ENDIF}

var i     : byte;
    Datum : float;

begin
    with Frequenzen do
       begin
           AnzDaten := VektorLaenge(Daten);
           OpenWindowType(5, 7, 75, 17, SubMenueWindow);
           writeln(Meldung6, AnzDaten:4);
           repeat
                write(Meldung1, round(sqrt(AnzDaten)):3, '): ');
                AnzKlassen := ReadByte(3);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         CloseWindow;
                         exit
                     end;
                writeln;
           until (AnzKlassen <= MaxKlassen) and (AnzKlassen > 0);
           writeln(Meldung2, HoleVektorElement(Daten, 1):Stellen);
           write(Meldung3);
           UntereGrenze := ReadReal(11, 5);
           if ReadError
             then
                begin
                    ReadError := false;
                    CloseWindow;
                    exit
                end;
           writeln;
           writeln(Meldung4, HoleVektorElement(Daten, AnzDaten):Stellen);
           write(Meldung5);
           ObereGrenze := ReadReal(11, 5);
           if ReadError
             then
                begin
                    ReadError := false;
                    CloseWindow;
                    exit
                end;
           writeln;
           CloseWindow;
           KlassenBreite := (ObereGrenze - UntereGrenze) / AnzKlassen;
           for i := 1 to AnzKlassen do
              Anzahl[i] := 0;
           Darueber := 0;
           Darunter := 0;
           for i := 1 to AnzDaten do
              begin
                  Datum := HoleVektorElement(Daten, i);
                  if (Datum < ObereGrenze) and (Datum >= UntereGrenze)
                    then
                       inc(Anzahl[succ(trunc((Datum - UntereGrenze) / KlassenBreite))])
                    else
                       if Datum < UntereGrenze
                         then
                            inc(Frequenzen.Darunter)
                         else
                            inc(Frequenzen.Darueber);
              end; { for }
       end; { with }
end; { Klassifiziere }


procedure WriteVerteilung (MedStr : string; Frequenzen : Verteilung);

{$IFDEF Deutsch}
const Meldung1 = ' H„ufigkeits-Verteilung ';
      Meldung2 = ' konnte nicht ge”ffnet werden';
      Meldung3 = ' Anzahl der Klassen: ';
      Meldung4 = ' Anzahl der Daten:   ';
      Meldung5 = ' Untere Grenze:      ';
      Meldung6 = ' Obere Grenze:       ';
      Meldung7 = ' Breite der Klassen: ';
      Meldung8 = ' Darunter:           ';
      Meldung9 = ' Darber:            ';
{$ENDIF}
{$IFDEF Englisch}
const Meldung1 = ' Frequency-Distribution ';
      Meldung2 = ' could not be opend';
      Meldung3 = ' Number of classes:  ';
      Meldung4 = ' Number of data:     ';
      Meldung5 = ' Lower border:       ';
      Meldung6 = ' Higher border:      ';
      Meldung7 = ' Class width:        ';
      Meldung8 = ' Below:              ';
      Meldung9 = ' Above:              ';
{$ENDIF}

var j, k    : word;
    Medium  : text;
    Fenster : WindowType;

begin
    for k := 1 to length(MedStr) do
       MedStr[k] := UpCase(MedStr[k]);
    if (MedStr = 'CON')
      then
         begin
             Fenster.TextAtt := 23;
             Fenster.RahmenAtt := 23;
             Fenster.Rahmen := einfach;
             Fenster.Titel := Meldung1;
             OpenWindowType(5, 1, 75, 25, Fenster);
         end;
    if (MedStr <> 'CON') and (MedStr <> 'LST')
      then
         begin
             assign(Medium, MedStr);
             rewrite(Medium);
             if IOResult <> 0
               then
                  begin
                      ch := Fehler(MedStr + Meldung2);
                      HaeufigError := true;
                      exit;
                  end;
         end;
    with Frequenzen do
       begin
           if (MedStr = 'CON')
             then
                begin
                    writeln(Meldung3, AnzKlassen:3);
                    writeln(Meldung4, AnzDaten:4);
                    writeln(Meldung5, FloatStr(UntereGrenze, Stellen));
                    writeln(Meldung6, FloatStr(ObereGrenze, Stellen));
                    writeln(Meldung7, FloatStr(KlassenBreite, Stellen));
                    for j := 1 to AnzKlassen do
                       writeln(FloatStr(UntereGrenze+pred(j)*KlassenBreite, Stellen), ' - ',
                               FloatStr(UntereGrenze+j*Klassenbreite, Stellen), ' : ', Anzahl[j]:3);
                    writeln(Meldung8, Darunter:3);
                    writeln(Meldung9, Darueber:3);
                    repeat until keypressed;
                    CloseWindow;
                    exit;
                end;
           if MedStr = 'LST'
             then
                begin
                    writeln(lst, Meldung3, AnzKlassen:3);
                    writeln(lst, Meldung4, AnzDaten:4);
                    writeln(lst, Meldung5, FloatStr(UntereGrenze, Stellen));
                    writeln(lst, Meldung6, FloatStr(ObereGrenze, Stellen));
                    writeln(lst, Meldung7, FloatStr(KlassenBreite, Stellen));
                    for j := 1 to AnzKlassen do
                       writeln(lst, FloatStr(UntereGrenze+pred(j)*KlassenBreite, Stellen), ' - ',
                               FloatStr(UntereGrenze+j*Klassenbreite, Stellen), ' : ', Anzahl[j]:3);
                    writeln(lst, Meldung8, Darunter:3);
                    writeln(lst, Meldung9, Darueber:3);
                    exit;
                end;
           writeln(Medium, AnzKlassen:3);
           writeln(Medium, AnzDaten:4);
           writeln(Medium, FloatStr(UntereGrenze, Stellen));
           writeln(Medium, FloatStr(ObereGrenze, Stellen));
           writeln(Medium, FloatStr(KlassenBreite, Stellen));
           for j := 1 to AnzKlassen do
              writeln(Medium, Anzahl[j]:3);
           writeln(Medium, Darunter:3);
           writeln(Medium, Darueber:3);
           close(Medium);
       end;
end;


procedure ReadVerteilung (MedStr : string; var Frequenzen : Verteilung);

{$IFDEF Deutsch}
const Meldung1 = ' Anzahl der Klassen (max. 20, empf. ';
      Meldung2 = ' Darunter: ';
      Meldung3 = ' Untere Grenze der Klassen-Bildung: ';
      Meldung4 = ' Darber: ';
      Meldung5 = ' Obere Grenze der Klassen-Bildung: ';
      Meldung6 = ' Anzahl der Daten: ';
      Meldung7 = ' konnte nicht ge”ffnet werden';
{$ENDIF}
{$IFDEF Englisch}
const Meldung1 = ' Number of classes (max. 20, recomended ';
      Meldung2 = ' Below: ';
      Meldung3 = ' Lower border of class formation: ';
      Meldung4 = ' Above: ';
      Meldung5 = ' Higher border of class formation: ';
      Meldung6 = ' Number of data: ';
      Meldung7 = ' could not be opend';
{$ENDIF}

var j, k    : word;
    Medium  : text;
    Fenster : WindowType;
    c       : char;

begin
    for k := 1 to length(MedStr) do
       MedStr[k] := UpCase(MedStr[k]);
    with Frequenzen do
       if MedStr = 'CON'
         then
            begin
                OpenWindowType(5, 7, 75, 23, SubMenueWindow);
                write(Meldung1);
                AnzKlassen := ReadByte(3);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         HaeufigError := true;
                         CloseWindow;
                         exit
                     end;
                Writeln;
                write(Meldung6);
                AnzDaten := ReadWord(4);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         HaeufigError := true;
                         CloseWindow;
                         exit
                     end;
                writeln;
                write(Meldung3);
                UntereGrenze := ReadReal(11, 5);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         HaeufigError := true;
                         CloseWindow;
                         exit;
                     end;
                writeln;
                write(Meldung5);
                ObereGrenze := ReadReal(11, 5);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         Haeufigerror := true;
                         CloseWindow;
                         exit
                     end;
                writeln;
                KlassenBreite := (ObereGrenze - UntereGrenze) / AnzKlassen;
                for j := 1 to AnzKlassen do
                   begin
                       write((UntereGrenze+pred(j)*KlassenBreite):Stellen, ' - ',
                             (UntereGrenze+j*Klassenbreite):Stellen, ' : ');
                       Anzahl[j] := ReadByte(3);
                       if ReadError
                         then
                            begin
                                ReadError := false;
                                HaeufigError := true;
                                CloseWindow;
                                exit
                            end;
                       writeln;
                   end;
                write(Meldung2);
                Darunter := ReadByte(3);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         HaeufigError := true;
                         CloseWindow;
                         exit
                     end;
                writeln;
                write(Meldung4);
                Darueber := ReadByte(3);
                if ReadError
                  then
                     begin
                         ReadError := false;
                         HaeufigError := true;
                         CloseWindow;
                         exit
                     end;
                repeat until keypressed;
                CloseWindow;
            end
         else
            begin
                assign(Medium, MedStr);
                reset(Medium);
                if IOResult <> 0
                  then
                     begin
                         ch := Fehler(MedStr + Meldung7);
                         HaeufigError := true;
                         exit;
                     end;
                readln(Medium, AnzKlassen);
                readln(Medium, AnzDaten);
                readln(Medium, UntereGrenze);
                readln(Medium, ObereGrenze);
                readln(Medium, KlassenBreite);
                for j := 1 to AnzKlassen do
                   readln(Medium, Anzahl[j]);
                readln(Medium, Darunter);
                readln(Medium, Darueber);
                close(Medium);
            end;
end;


function Haeufigkeit (Frequenzen : Verteilung; KlNr : byte) : word;

begin
    Haeufigkeit := Frequenzen.Anzahl[KlNr];
end;


function GetLower (Frequenzen : Verteilung; KlNr : byte) : float;

begin
    GetLower := Frequenzen.UntereGrenze + pred(KlNr) * Frequenzen.KlassenBreite;
end;


function GetUpper (Frequenzen : Verteilung; KlNr : byte) : float;

begin
    GetUpper := Frequenzen.UntereGrenze + KlNr * Frequenzen.KlassenBreite;
end;


function GetAverage (Frequenzen : Verteilung; KlNr : byte) : float;

begin
    GetAverage := (GetUpper(Frequenzen, KlNr) + GetLower(Frequenzen, KlNr)) / 2;
end;


function AnzDarueber (Frequenzen : Verteilung) : word;

begin
    AnzDarueber := Frequenzen.Darueber;
end;


function AnzDarunter (Frequenzen : Verteilung) : word;

begin
    AnzDarunter := Frequenzen.Darunter;
end;


function Breite (Frequenzen : Verteilung) : float;

begin
    Breite := Frequenzen.Klassenbreite;
end;


function AnzahlKlassen (Frequenzen : Verteilung) : byte;

begin
    AnzahlKlassen := Frequenzen.AnzKlassen;
end;


function AnzahlDaten (Frequenzen : Verteilung) : word;

begin
    AnzahlDaten := Frequenzen.AnzDaten;
end;

{ ************************************************************************** }

function ArithmetischesMittel (Frequenzen : Verteilung) : float;

var i            : word;
    Summe        : float;

begin
    if AnzahlDaten(Frequenzen) = 0
      then
         begin
             HaeufigError := true;
             exit;
         end;
    Summe := 0;
    for i := 1 to AnzahlKlassen(Frequenzen) do
       begin
           Summe := Summe + GetAverage(Frequenzen, i) * Haeufigkeit(Frequenzen, i);
       end;
    ArithmetischesMittel := Summe / AnzahlDaten(Frequenzen);
end;


function GeometrischesMittel (Frequenzen : Verteilung) : float;

var Produkt : float;
    i       : word;

begin
    Produkt := 1;
    for i := 1 to AnzahlKlassen(Frequenzen) do
       begin
           if GetAverage(Frequenzen, i) > 0
             then
                Produkt := Produkt + ln(GetAverage(Frequenzen, i)) * Haeufigkeit(Frequenzen, i)
             else
                begin
                    HaeufigError := true;
                    exit;
                end;
       end;
    GeometrischesMittel := exp(Produkt / AnzahlDaten(Frequenzen));
end;


function HarmonischesMittel (Frequenzen : Verteilung) : float;

var Summe : float;
    i     : word;

begin
    Summe := 0;
    for i := 1 to AnzahlKlassen(Frequenzen) do
       begin
           if GetAverage(Frequenzen, i) <> 0
             then
                Summe := Summe + Haeufigkeit(Frequenzen,i) / GetAverage(Frequenzen, i)
             else
                begin
                    HaeufigError := true;
                    exit;
                end;
       end;
    HarmonischesMittel := AnzahlDaten(Frequenzen) / Summe;
end;


function Varianz (Frequenzen : Verteilung) : float;

var Summe1, Summe2 : float;
    i              : word;

begin
    Summe1 := 0;
    Summe2 := 0;
    for i := 1 to AnzahlKlassen(Frequenzen) do
       begin
           Summe1 := Summe1 + GetAverage(Frequenzen, i) * Haeufigkeit(Frequenzen, i);
           Summe2 := Summe2 + sqr(GetAverage(Frequenzen, i)) * Haeufigkeit(Frequenzen, i);
       end;
    Varianz := (Summe2 - sqr(Summe1)/AnzahlDaten(Frequenzen)) / AnzahlDaten(Frequenzen);
end;


function ShepphardKorrigierteVarianz (Frequenzen : Verteilung) : float;

begin
    ShepphardKorrigierteVarianz := Varianz(Frequenzen) - sqr(Breite(Frequenzen)) / 12;
end;


function Quantile (Frequenzen : Verteilung; QuantilWert : float) : float;

var Summe, Nummer,
    Davor         : word;
    Grenze        : float;

begin
    Grenze := AnzahlDaten(Frequenzen) * QuantilWert;
    Summe := 0;
    Nummer := 0;
    while (Summe < Grenze) do
        begin
            inc(Nummer);
            Summe := Summe + Haeufigkeit(Frequenzen, Nummer);
        end;
    Davor := Summe - Haeufigkeit(Frequenzen, Nummer);
    Quantile := (Grenze - Davor) / Haeufigkeit(Frequenzen, Nummer) * Breite(Frequenzen)
                + GetLower(Frequenzen, Nummer);
end;


function Entropie (Frequenzen : Verteilung) : float;

var Summe : float;
    i     : word;

begin
    Summe := 0;
    for i := 1 to AnzahlKlassen(Frequenzen) do
       if Haeufigkeit(Frequenzen, i) > 0
         then
            Summe := Summe + Haeufigkeit(Frequenzen, i) / AnzahlDaten(Frequenzen)
                     * log(2, AnzahlDaten(Frequenzen) / Haeufigkeit(Frequenzen, i));
    Entropie := Summe;
end;


procedure Beschreibe (Frequenzen : Verteilung);

{$IFDEF Deutsch}
const Meldung1  = ' Anzahl der Klassen:   ';
      Meldung2  = ' Anzahl der Daten:     ';
      Meldung3  = ' Arthmetisches Mittel: ';
      Meldung4  = ' Varianz:              ';
      Meldung5  = ' kor. nach Shepphard:  ';
      Meldung6  = ' Geometrisches Mittel: ';
      Meldung7  = ' Harmonisches Mittel:  ';
      Meldung8  = ' Entropie:             ';
      Meldung9  = ' Q1:                   ';
      Meldung10 = ' Q2:                   ';
      Meldung11 = ' Q3:                   ';
{$ENDIF}
{$IFDEF Englisch}
const Meldung1  = ' Number of classes:   ';
      Meldung2  = ' Number of data:      ';
      Meldung3  = ' Arithmetic mean:     ';
      Meldung4  = ' Variance:            ';
      Meldung5  = ' Shepphard-corrected: ';
      Meldung6  = ' Geometric mean:      ';
      Meldung7  = ' Harmonic mean:       ';
      Meldung8  = ' Entropy:             ';
      Meldung9  = ' Q1:                  ';
      Meldung10 = ' Q2:                  ';
      Meldung11 = ' Q3:                  ';
{$ENDIF}

begin
    OpenWindowType(5, 7, 75, 23, SubMenueWindow);
    writeln(Meldung2, Frequenzen.AnzDaten:3);
    writeln(Meldung1, Frequenzen.AnzKlassen:2);
    writeln(Meldung3, FloatStr(ArithmetischesMittel(Frequenzen), Stellen));
    writeln(Meldung4, FloatStr(Varianz(Frequenzen), Stellen));
    writeln(Meldung5, FloatStr(ShepphardKorrigierteVarianz(Frequenzen), Stellen));
    writeln(Meldung6, FloatStr(GeometrischesMittel(Frequenzen), Stellen));
    writeln(Meldung7, FloatStr(HarmonischesMittel(Frequenzen), Stellen));
    writeln(Meldung8, FloatStr(Entropie(Frequenzen), Stellen));
    writeln(Meldung9, FloatStr(Quantile(Frequenzen, 0.25), Stellen));
    writeln(Meldung10, FloatStr(Quantile(Frequenzen, 0.50), Stellen));
    writeln(Meldung11, FloatStr(Quantile(Frequenzen, 0.75), Stellen));
    repeat until keypressed;
    CloseWindow;
end;


procedure KumulPlot (Frequenzen : Verteilung; var xVektor, yVektor : VektorTyp);

var Summe, i : word;
    xKoord   : float;

begin
    with Frequenzen do
       begin
           InitVektor(xVektor, AnzKlassen);
           InitVektor(yVektor, AnzKlassen);
           Summe := 0;
           for i := 1 to AnzKlassen do
              begin
                  xKoord := UntereGrenze + i * Klassenbreite;
                  Summe := Summe + Anzahl[i];
                  SetzeVektorElement(xVektor, i, xKoord);
                  SetzeVektorElement(yVektor, i, Summe);
              end;
       end;
end;


function LorenzMuenzner (Frequenzen : Verteilung; var xKoord, yKoord : VektorTyp) : float;

var Summe,
    gSumme,
    xSumme,
    ySumme,
    v, vSum,
    x, y, Alt: float;
    i        : word;

begin
    with Frequenzen do
       begin
           gSumme := 0;
           InitVektor(xKoord, AnzKlassen);
           InitVektor(yKoord, AnzKlassen);
           for i := 1 to AnzahlKlassen(Frequenzen) do
              GSumme := gSumme + Haeufigkeit(Frequenzen, i) * GetAverage(Frequenzen, i);
           ySumme := 0;
           xSumme := 0;
           vSum   := 0;
           Alt := 0;
           for i := 1 to AnzahlKlassen(Frequenzen) do
              begin
                  x := Haeufigkeit(Frequenzen, i) / AnzahlDaten(Frequenzen);
                  xSumme := xSumme + x;
                  SetzeVektorElement(xKoord, i, xSumme);
                  y := Haeufigkeit(Frequenzen, i) * GetAverage(Frequenzen, i) / gSumme;
                  ySumme := ySumme + y;
                  SetzeVektorElement(yKoord, i, ySumme);
                  vSum := vSum + (Alt + ySumme) * x;
                  Alt := ySumme;
(*                write(i:2, '  ', xSumme:5:3, '  ', ySumme:5:3, '  ', vSum:6:3);
                  readln;
*)            end;
           LorenzMuenzner := 1 - vSum;
       end;
end;


procedure SaeulenDiagramm (Frequenzen : Verteilung; NameStr, LizStr : string);

{$IFDEF Deutsch}  const Meldung = ' Bezeichnung fr x-Achse: '; {$ENDIF}
{$IFDEF Englisch} const Meldung = ' Label for x-axis: ';        {$ENDIF}

var xMaximum, yMaximum,
    xMax, yMax         : word;
    xFaktor            : float;
    xname              : string[8];
    c                  : char;


     function GetYMaximum : word;

     var i, temp : word;

     begin
         temp := 0;
         for i := 1 to AnzahlKlassen(Frequenzen) do
           if Haeufigkeit(Frequenzen, i) > temp then temp := Haeufigkeit(Frequenzen, i);
         if AnzDarunter(Frequenzen) > temp then temp := AnzDarunter(Frequenzen);
         if AnzDarueber(Frequenzen) > temp then temp := AnzDarueber(Frequenzen);
         if odd(temp) then inc(temp);
         GetYMaximum := temp;
     end;

     procedure Koordinaten;

     var xPos, yPos, i : word;
         ZahlStr       : string[11];

     begin
         SetColor(White);
         Line(25, yMax-20, xMax-60, yMax-20);      { x-Achse zeichnen }
         OutTextXY(xMax-40, yMax-20, xName);
         Line(30, yMax-15, 30, 20);                { y-Achse zeichnen }
         OutTextXY(10, 10, 'n');
         with Frequenzen do                        { x-Achse beschriften }
            begin
                xPos := xMax - 90;
                yPos := yMax-20;
                Line(xPos, yPos, xPos, yPos+5);
                OutTextXY(xPos-8, yPos+10, FloatStr(ObereGrenze, Stellen));
                xPos := 60;
                Line(xPos, yPos, xPos, yPos+5);
                OutTextXY(xPos-8, yPos+10, FloatStr(UntereGrenze, Stellen));
                xFaktor := (xMax - 150) / (ObereGrenze - UntereGrenze);
                for i := 1 to pred(AnzKlassen) do
                   begin
                       xPos := 60 + round((GetUpper(Frequenzen, i) - UntereGrenze) * xFaktor);
                       Line(xPos, yPos, xPos, yPos+5);
                   end;
                yPos := yMax-20;                   { y-Achse beschriften }
                xPos := 30;
                OutTextXY(0, yPos-5, '0');
                yPos := 30;
                str(yMaximum:3, ZahlStr);
                Line(xPos, yPos, xPos-5, yPos);
                OutTextXY(0, yPos-5, ZahlStr);
                str((yMaximum div 2):3, ZahlStr);
                yPos := (yMax-50) div 2 + 20;
                Line(xPos, yPos, xPos-5, yPos);
                OutTextXY(0, yPos-5, ZahlStr);
            end;
     end;


     procedure DatenEinzeichnen;

     var i,
         Oben, Unten,
         Rechts, Links,
         Tiefe, xPos    : word;
         yFaktor        : float;
         ZahlStr        : string[3];

     begin
         Unten := yMax - 21;
         yFaktor := (yMax - 40) / yMaximum;
         SetLineStyle(SolidLn, LightGreen, NormWidth);
         SetFillStyle(LtSlashFill, Lightgreen);
         Tiefe := 4;
         for i := 1 to AnzahlKlassen(Frequenzen) do
            begin
                Rechts := round((GetLower(Frequenzen, i) - Frequenzen.UntereGrenze)
                                 * xFaktor) + 63;
                Links := round((GetUpper(Frequenzen, i) - Frequenzen.UntereGrenze )
                                 * xFaktor) + 57;
                Oben := yMax - (round(Haeufigkeit(Frequenzen, i) * yFaktor) + 20);
                Bar3D(Links, Oben, Rechts, Unten, Tiefe, true);
                xPos := (Rechts + Links) div 2 - 20;
                str(Haeufigkeit(Frequenzen, i):3, ZahlStr);
                OutTextXY(xPos, Oben-20, ZahlStr);
            end;
         SetLineStyle(DottedLn, LightGreen, NormWidth);
         SetFillStyle(CloseDotFill, LightGreen);
         if AnzDarunter(Frequenzen) > 0
           then
              begin
                  Rechts := 26;
                  Links := 52;
                  Oben := yMax - (round(AnzDarunter(Frequenzen) * yFaktor) + 20);
                  Bar3D(Links, Oben, Rechts, Unten, Tiefe, true);
                  xPos := (Rechts + Links) div 2 - 10;
                  str(AnzDarunter(Frequenzen):3, ZahlStr);
                  OutTextXY(xPos, Oben-20, ZahlStr);
              end;
         if AnzDarueber(Frequenzen) > 0
           then
              begin
                  Rechts := xMax - 60;
                  Links := xMax - 90;
                  Oben := yMax - (round(AnzDarueber(Frequenzen) * yFaktor) + 20);
                  Bar3D(Links, Oben, Rechts, Unten, Tiefe, true);
                  xPos := (Rechts + Links) div 2 - 10;
                  str(AnzDarueber(Frequenzen):3, ZahlStr);
                  OutTextXY(xPos, Oben-20, ZahlStr);
              end;
     end;


begin
    OpenWindowType(5, 7, 75, 20, SubMenueWindow);
    write(Meldung);
    xName := ReadStr(8);
    CloseWindow;
    if ReadError
      then
         begin
             ReadError := false;
             exit;
         end;
    OpenWindow(1, 1, 80, 25);
    yMaximum := GetYMaximum;
    OpenGraphik;
    if GraphikError
      then
         begin
             GraphikError := false;
	     HaeufigError := true;
             CloseWindow;
             exit;
         end;
    xMax := GetMaxX;
    yMax := GetMaxY;
    Koordinaten;
    DatenEinzeichnen;
    repeat until KeyPressed;
    if (Upcase(ReadKey) = 'P') then GraphikDrucken(NameStr+'  fr '+LizStr);
    CloseGraphik;
    CloseWindow;
end;


procedure ZeichneTorte (Frequenzen : Verteilung; NameStr, LizStr : string);

var xMax, yMax,
    xMitte, yMitte, n,
    Radius, i, Kumulativ,
    AlterWinkel, NeuerWinkel : word;
    Multiplikator            : float;


begin
    OpenWindow(1, 1, 80, 25);
    OpenGraphik;
    if GraphikError
      then
         begin
             GraphikError := false;
	     HaeufigError := true;
             CloseWindow;
             exit;
         end;
    xMax := GetMaxX;
    yMax := GetMaxY;
    xMitte := xMax div 2;
    yMitte := yMax div 2;
    Radius := round(0.5 * yMax * 0.9);
    Multiplikator := 360 / AnzahlDaten(Frequenzen);
    AlterWinkel := 0;
    Kumulativ := 0;
    n := AnzahlKlassen(Frequenzen);
    for i := 1 to n do
       begin
           SetFillStyle(i mod 11, succ(i mod 15));
           Kumulativ := Kumulativ + Haeufigkeit(Frequenzen, i);
           NeuerWinkel := round(Kumulativ * Multiplikator);
           if NeuerWinkel < 360
             then
                PieSlice(xMitte, yMitte, AlterWinkel, NeuerWinkel, Radius)
             else
                Arc(xMitte, yMitte, AlterWinkel, 360, Radius);
           AlterWinkel := NeuerWinkel;
       end;
    repeat until KeyPressed;
     if (Upcase(ReadKey) = 'P') then GraphikDrucken(NameStr+'  fr '+LizStr);
    CloseGraphik;
    CloseWindow;
end;

end. { haeufig }

