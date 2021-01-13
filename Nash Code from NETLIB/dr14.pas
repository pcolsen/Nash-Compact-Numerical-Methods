{$I turbo.cnm}
Program dr14(input, output);
{dr14.pas == driver for Jacobi method (Alg14) for eigensolutions of a real
          symmetric matrix

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I resids}       {compute residuals -- used in alg2.pas}
{$I Matrixin}     {input or generate a matrix of reals}
{$I alg14.pas}    {evJacobi -- symmetric matrix eigenproblem}
{$I startup.pas}

{main program}
var
  i, j, k, nRow, nCol : integer;
  A, V, ACOPY : rmatrix;
  Bvec, Y : rvector; {to test residuals}
  t1: real;
  tvec : smatvec; {needed only for Matrixin}
  initev, sym : boolean;
  ch : char;

begin
  banner:='dr14.pas -- driver for Jacobi eigensolution method';
  startup;
  write('Order of problem (n) = ');
  readln(infile, nRow); if infname<>'con' then writeln(nRow);
  writeln(confile, 'Order of problem (n) = ', nRow);
  nCol:=nRow;
  Matrixin(nRow, nCol, A, tvec, sym);
  if not sym then halt; {program only designed for symmetric matrices}
  for j:=1 to nRow do
  begin
    for i:=1 to nRow do
    begin
      write(A[i, j]:10:5, ' ');
      write(confile, A[i, j]:10:5, ' ');
      ACOPY[i, j]:=A[i, j];
      if (7 * (i div 7) = i) and (i<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  writeln;
  writeln(confile);
  initev:=true; {Here we want to get the eigenvectors of A, not some
            generalized problem.}
  evJacobi( nRow, A, V, initev);
  for j:=1 to nRow do
  begin
    t1:=A[j, j];
    writeln('Eigenvalue ', j, ' = ', t1);
    writeln(confile, 'Eigenvalue ', j, ' = ', t1);
    for i:=1 to nRow do
    begin
      write(V[i, j]:10:7, ' ');
      write(confile, V[i, j]:10:7, ' ');
      if (i = 7 * (i div 7)) and (i<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
      Bvec[i]:=V[i, j]; {to initialize residual test}
      Y[i]:=t1*Bvec[i];
    end;
    writeln;
    writeln(confile);
    t1 := resids(nRow, nCol, ACOPY, Y, Bvec, true);
    writeln;
    writeln(confile);
  end; {loop on solutions j}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr14.pas}
