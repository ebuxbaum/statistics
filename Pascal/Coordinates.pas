program Coordinates;

uses math, MathFunc;

const Codes : array[1..20] of char = ('2', '3', '4', '5', '6', '7', '8', '9',
                                      'C', 'F', 'G', 'H', 'J', 'M', 'P', 'Q',
                                      'R', 'V', 'W', 'X');
      Resolution : array[1..5] of double = (20, 1, 1/20, 1/400, 1/8000);
      GridCols         = 4;
      GridRows         = 5;
      GridSizeDeg      = 1/400;
      GridRowSize      = GridSizeDeg / GridRows;
      GridColSize      = GridSizeDeg / GridCols;
      CodeLengthNormal = 10;
      CodeLengthExtra  = 11;
      Separator        = '+';
      SeparatorPos     = 8;

var Latitude, Longitude : double;

function Encode (Latitude, Longitude : double; CodeLength : byte) : string;

var Code : string;
    len, idx, PairCount, enc  : byte;

begin
  CodeLength := Min(CodeLengthExtra, Max(CodeLength, 2));
  Latitude = Min(90, Max(-90, Latitude));
  while (Longitude < -180) do Longitude := Longitude + 360;
  while (Longitude > 180) do Longitude := Longitude - 360;
  if (Latitude = 90)
    then
      if (CodeLength <= CodeLengthNormal)
        then Latitude := Latitude - pot(Divisor, floor(2 - CodeLength/2))
        else Latitude := Latitude - pot(Divisor, -3);
  Latitude := Latitude + 90;
  Longitude := Longitude + 180;
  Code := '';
  len := min(CodeLengthNormal, CodeLength);
  for idx = 0 to len do
    begin
      PairCount = floor(idx/2);
      enc :=
    end;
end;

begin
end.

