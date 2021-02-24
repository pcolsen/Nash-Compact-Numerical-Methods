{$I turbo.cnm}
Program dr13(input,output);
{dr13.pas == run Nash svd for eigenvalue computations (Alg13)

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I calceps}
{$I tdstamp.pas}  {time and date stamp}
{$I resids}       {compute residuals -- used in alg2.pas}
{$I Matrixin}     {input or generate a matrix of reals}
{$I alg01.pas}    {Nashsvd -- singular value decomposition}
{$I alg13.pas}    {evsvd -- symmetric matrix eigenproblem}
{$I startup.pas}

var
  i,j,k, nRow, nCol : integer;
  A, V, ACOPY : rmatrix;
  Bvec, Y, Z : rvector;
  W : wmatrix; {to store the working array}
  ch : char;
  t1: real;
  tvec : smatvec; {needed only for Matrixin}
  initev, sym : boolean;

begin
  banner:='dr13.pas -- driver for svd eigensolutions of a symmetric matrix';
  startup;
  write('Order of problem (n): '); readln(infile,nRow);
  writeln(confile,'Order of problem (n): ',nRow);
  if infname<>'con' then writeln(nRow);
  nCol := nRow;
  Matrixin(nRow, nCol, A, tvec, sym);
  if not sym then halt;
  for j := 1 to nRow do
  begin
    for i := 1 to nRow do
    begin
      write(A[i,j]:10:5,' ');
      write(confile,A[i,j]:10:5,' ');
      ACOPY[i,j] := A[i,j];
      if (7 * (i div 7) = i) and (i<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  initev := true; {Here we want to get the eigenvectors of A, not some
            generalized problem.}
  evsvd( nRow, A, V, initev, W, Z);
  for j := 1 to nRow do
  begin
    t1 := Z[j];
    writeln;
    writeln(confile);
    writeln('Eigenvalue ',j,' = ',t1);
    writeln(confile,'Eigenvalue ',j,' = ',t1);
    for i := 1 to nRow do
    begin
      write(V[i,j]:10:7,' ');
      write(confile,V[i,j]:10:7,' ');
      if (i = 7 * (i div 7)) and (i<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
      Bvec[i] := V[i,j]; {to initialize residual test}
      Y[i] := t1*Bvec[i];
    end;
    writeln;
    writeln(confile);
    t1 := resids(nRow, nCol, ACOPY, Y, Bvec,true);
  end; {loop on solutions j}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr13.pas}
