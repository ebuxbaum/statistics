PROGRAM Principal;

USES
  Math,              // Lazarus math UNIT
  mathfunc,          // basic mathematical routines FOR real AND INTEGER
  Complex,           // routines FOR complex numbers
  Vector,            // vector algebra
  Matrix,            // basic matrix algebra
  stat,              // statistical TEST distributions
  Correlations,      // various types OF correlation coefficients
  PCA,               // principal component analysis
  PerformShrinkage,  // Ledoit-Wolf shrinkage OF correlation matrix
  Rotation,          // rotation OF loading matrix
  SystemSolve,       // linear equations
  EigenValues,       // calculation OF eigenvalues
  Zufall;            // generate Random numbers

VAR
  Types                                                : TypeVector;
  EV, Variance, CumVariance, ImputVector, Acc,
  Communalities, Uniqueness, VarianceAccounted         : VectorTyp;
  Data, Cor, EigenVectoren, Scores, Identity, Loadings : MatrixTyp;
  Freqs                                                : FreqsType;
  x                                                    : float;
  Cases, i, iter                                       : WORD;


BEGIN
  ProblemName := 'Anonym';
  SigVecs := 40;
  ReadCSV(Types, Data);
  CalculateCorrelationMatrix(Data, Types, Cor);
  WriteCorrelations(Cor, FALSE);
  Cases := MatrixRows(Data);
  Average(Cor, 0.5);               // perform shrinkage by average correlation
  WriteCorrelations(Cor, TRUE);
  AnalyseFrequencies(Cor, Types, Freqs);
  WriteCorrelations(Cor, TRUE);
  x := Bartlett(Cor, Cases);
  i := Jacobi(Cor, EV, EigenVectoren, Iter); // eigenanalysis
  Writeln('Jacobi: ', Iter: 3, ' iterations, Result = ', i:2, ': ', EigenError[i]);
  DestroyMatrix(Cor);
  SortEigenValues(EV, EigenVectoren);
  ExplainedVariance(EV, Acc, Variance, CumVariance);
  WriteEigenVector(EigenVectoren);
  Writeln('Eigenvectors and -values written');
  MaximumLikelihood(Data, Types, ImputVector);       // Scores
  RobustProduct(Data, EigenVectoren, ImputVector, Scores);
  WriteComponentScores(Scores);
  Writeln('scores written');
  CalculateLoadings(Data, Scores, Types, Loadings);  // unrotated Loadings
  CalculateCommunalities(Loadings, Communalities, Uniqueness, VarianceAccounted);
  WriteLoadings(Loadings, Communalities, Uniqueness, VarianceAccounted, FALSE);
  Writeln('Loadings and communalities calculated and written');
  DestroyVector(Variance);
  DestroyVector(CumVariance);
  DestroyVector(Communalities);
  CreateIdentityMatrix(Identity, SigVecs);
  GradProjAlgOrth(Loadings, Identity, Varimax);      // rotation
  CalculateCommunalities(Loadings, Communalities, Uniqueness, VarianceAccounted);
  WriteLoadings(Loadings, Communalities, Uniqueness, VarianceAccounted, TRUE);
  DestroyMatrix(Identity);
  DestroyMatrix(Data);
  DestroyMatrix(EigenVectoren);
  DestroyMatrix(Scores);
  DestroyMatrix(Loadings);
  DestroyVector(EV);
  DestroyVector(ImputVector);
  DestroyVector(Communalities);
  Writeln('All done, press <CR> to finish:');
  ReadLn;
END.
