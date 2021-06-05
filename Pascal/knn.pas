UNIT kNN;

INTERFACE

USES Math,
     MathFunc,
     Vector,
     Matrix,
     Deskript,
     Zufall,
     CrossValidation;

CONST
  kNNError: BOOLEAN = FALSE;
  KNNmax = 15;
  pMax = 100;

TYPE
  GroupTyp = ARRAY[1..MaxVector] OF WORD;

FUNCTION kMeansCluster(CONST Data: MatrixTyp;       // data matrix
  k: WORD;                                          // number OF clusters
  VAR Centroids: MatrixTyp;                         // centroids FOR all k clusters AND p variables
  VAR Group: GroupTyp): float;                      // group no FOR each item
{ splits a data matrix into k groups by the k-nearest neighbour method and
  returns the sum of Euklidian distances of the data from their centroids. }

PROCEDURE AssignTestData(CONST TestData, Centroids: MatrixTyp;
  kValidate : WORD;                                         // number OF cross-validation groups
  VAR Count : WORD;                                         // no OF results calculated already
  VAR ValidationResult: MatrixTyp);                         // Result
{ assigns all rows of Data to the closest Centroid and
  returns the study number, assigned group and squared distance to the centroid
  in ValidationResult }

FUNCTION WithinSumOfSquares(CONST ValidationResult: MatrixTyp): float;
{ calculates the sum of the squared distance of all data from their centroid,
  which is stored in the third column of ValidationResult }

IMPLEMENTATION

PROCEDURE StartingCentoids(CONST Data: MatrixTyp; k: WORD; VAR Centroids: MatrixTyp);
// select the starting centroids using knn++ METHOD (distance-controlled Random)

VAR
  v, vi, row, Distance : VectorTyp;
  i, j, l, n, p        : WORD;
  r                    : LONGINT;
  Used                 : ARRAY[1..KNNmax] OF WORD;
  u                    : BOOLEAN;
  c                    : CHAR;
  d, Sum               : float;

BEGIN
  n := MatrixRows(Data);
  p := MatrixColumns(Data);
  CreateMatrix(Centroids, k, p, 0.0);
  i := RandomLaplace(1, n);          // randomly select first centroid
  Writeln('first centroid = ', i: 3);
  GetRow(Data, i, v);
  SetRow(Centroids, v, 1);
  Used[1] := i;
  FOR j := 2 TO k DO                        // select the remaining centroids
    BEGIN
      CreateVector(Distance, n, 0.0);
      IF VectorError
        THEN
          BEGIN
            c := WriteErrorMessage('program terminated');
            HALT;
          END;
      Sum := 0;
      FOR i := 1 TO n DO // calculate distance between all items AND last centroid
        BEGIN
          GetRow(Data, i, vi);
          d := SquaredEuklidianDistance(v, vi, TRUE);
          DestroyVector(vi);
          SetVectorElement(Distance, i, d);
          FOR l := 1 TO Pred(j) DO          // ignore IF datum i IS already a centroid
            IF Used[l] = i THEN u := TRUE;
          IF NOT (u) THEN Sum := Sum + d;
        END;
      l := ceil(Sum);
      r := RandomLaplace(1, l); // distance-controlled Random selection OF NEW centroid
      Sum := 0;
      i := 0;
      REPEAT
        INC(i);
        u := FALSE;
        FOR l := 1 TO Pred(j) DO
          IF Used[l] = i THEN u := TRUE;
        IF NOT (u) THEN Sum := sum + Round(GetVectorElement(Distance, i));
      UNTIL (Sum >= r);
      GetRow(Data, i, vi);
      SetRow(Centroids, vi, j);
      FOR i := 1 TO p DO
        SetVectorElement(v, i, GetVectorElement(vi, i));
      DestroyVector(vi);
      DestroyVector(Distance);
    END; // FOR j
  DestroyVector(v);
END;  // StartingCentroids

FUNCTION AssignItems(CONST Data, Centroids: MatrixTyp;
  k: WORD;
  VAR Group: GroupTyp): float;
  // Assign each datum TO the cluster from the centre OF which it has minimal
  // distance. Returns sum OF the minimal distances over all data

VAR
  i, j, l, n, p : WORD;
  d, min, Sum   : float;
  vj, vi        : VectorTyp;

BEGIN
  n := MatrixRows(Data);
  p := MatrixColumns(Data);
  Sum := 0;
  FOR i := 1 TO n DO             // FOR all data rows
    BEGIN
      min := MaxRealNumber;
      GetRow(Data, i, vi);
      FOR j := 1 TO k DO         // find centroid WITH minimal distance TO the datum
      BEGIN
        GetRow(Centroids, j, vj);
        d := SquaredEuklidianDistance(vj, vi, TRUE);
        IF (d < min)
          THEN
            BEGIN
              min := d;
              l := j;
            END; // THEN
        DestroyVector(vj);
      END; // FOR j
      Group[i] := l;
      Sum := Sum + min;
      DestroyVector(vi);
    END; // FOR i
  Result := Sum;
END; // AssignItems


PROCEDURE CalculateNewCentroids(CONST Data: MatrixTyp;
  k: WORD;
  VAR Group: GroupTyp;
  VAR Centroids: MatrixTyp
  );

VAR
  i, j, l, n, p, nk : WORD;
  Sums              : ARRAY [1..pMax] OF float;
  x                 : float;

BEGIN
  n := MatrixRows(Data);
  p := MatrixColumns(Data);
  FOR l := 1 TO k DO
    BEGIN
      FOR j := 1 TO p DO
        Sums[j] := 0.0;
      nk := 0;
      FOR i := 1 TO n DO
        BEGIN
          IF Group[i] = l
            THEN
              BEGIN
                FOR j := 1 TO p DO
                  BEGIN
                    x := GetMatrixElement(Data, i, j);
                    Sums[j] := Sums[j] + x;
                  END;
              END;
          INC(nk);
        END;
      FOR j := 1 TO p DO
        IF (nk = 0)
          THEN
            SetMatrixElement(Centroids, k, j, 0) // no data IN group, shouldn't happen
          ELSE
            SetMatrixElement(Centroids, k, j, Sums[j] / nk);
    END; // for l
END; // CalculateNewCentroids


FUNCTION kMeansCluster(CONST Data: MatrixTyp;
  k: WORD;
  VAR Centroids: MatrixTyp;
  VAR Group: GroupTyp): float;

VAR
  p, n, Iter         : WORD;
  TotalDist, OldDist : float;

BEGIN  // kMeansCluster
  n := MatrixRows(Data);
  p := MatrixColumns(Data);
  StartingCentoids(Data, k, Centroids);
  TotalDist := MaxRealNumber;
  Iter := 0;
  REPEAT
    OldDist := TotalDist;
    Inc(Iter);
    TotalDist := AssignItems(Data, Centroids, k, Group);
    CalculateNewCentroids(Data, k, Group, Centroids);
  UNTIL (abs(OldDist - TotalDist) < MaxError) OR (Iter > MaxIter);
  Result := TotalDist;
END;   // kMeansCluster

PROCEDURE AssignTestData(CONST TestData, Centroids: MatrixTyp;
  kValidate : WORD;
  VAR Count : WORD;
  VAR ValidationResult: MatrixTyp);

VAR
  i, j, k, n, p, Opt : WORD;
  vi, vj             : VectorTyp;
  d, MinD            : float;

BEGIN
  n := MatrixRows(TestData);
  p := MatrixColumns(TestData);
  k := MatrixRows(Centroids);
  CreateVector(vi, p, 0.0);
  FOR i := 1 TO n DO       // for all data in test data matrix
    BEGIN
      GetRow(TestData, i, vi);
      MinD := MaxRealNumber;
      Opt := 0;
      FOR j := 1 TO k DO  // find centroid of minimal distance
        BEGIN
          GetRow(Centroids, j, vj);
          d := SquaredEuklidianDistance(vi, vj, True);
          IF d < MinD
            THEN
              BEGIN
                MinD := d;
                Opt := j;
              END; // then
          DestroyVector(vj);
        END; // for j
      SetMatrixElement(ValidationResult, Count, 1, GetMatrixElement(TestData, i, 1));
      // study number
      SetMatrixElement(ValidationResult, Count, 2, Opt);  // optimal centroid
      SetMatrixElement(ValidationResult, Count, 3, MinD); // distance from centroid
      DestroyVector(vi);
      INC(Count);
    END; // for i
END; // AssignTestData


FUNCTION WithinSumOfSquares(CONST ValidationResult: MatrixTyp): float;

VAR
  i, n: WORD;
  WSS: float;

BEGIN
  n := MatrixRows(ValidationResult);
  WSS := 0;
  FOR i := 1 TO n DO
    WSS := WSS + GetMatrixElement(ValidationResult, i, 3);
  Result := WSS;
END;  // WithinSumOfSquares

end. //kNN
