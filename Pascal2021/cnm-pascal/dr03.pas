{$I turbo.cnm}
program givrun(input, output);
{dr03.PAS ==  driver for Givens' reduction of a matrix

          Copyright 1988 J.C.Nash
}
{$I constype.def}   {definitions of constants and types}
{$I tdstamp.pas}    {time and date stamp}
{$I calceps.pas}    {compute machine precision}
{$I matrixin.pas}   {input or generate a matrix of reals}
{$I alg03.pas }     {Givens' reduction of a matrix}
{$I startup.pas}
{Replace the last include file with alg03a.pas to get a
  row-wise version of the Givens' reduction.}

var
  A, Q: rmatrix;
  i, j, k, nRow, nCol : integer;
  Acopy : rmatrix;
  avec: smatvec; {for matrixin only}
  sym : boolean; {for matrixin only}
  s : real;

begin
  banner:='dr03.pas -- driver for Givens'+chr(39)+' reduction';
  startup;
  write('Size of problem (rows, columns) ');
  readln(infile,nRow,nCol); if infname<>'con' then writeln(nRow,' ',nCol);
  writeln(confile,'Size of problem (rows, columns) ',nRow,' ',nCol);
  matrixin(nRow,nCol,A,avec,sym);
  writeln('Matrix A');
  writeln(confile,'Matrix A');
  for i:=1 to nRow do
  begin
    for j:=1 to nCol do
    begin
      Acopy[i,j]:=A[i,j];
      write(A[i,j]:10:5,' ');
      write(confile,A[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nCol) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  givens(nRow,nCol,A,Q);
  writeln('Decomposition');
  writeln(confile,'Decomposition');
  writeln('Q');
  writeln(confile,'Q');
  for i:=1 to nRow do
  begin
    for j:=1 to nRow do
    begin
      write(Q[i,j]:10:5,' ');
      write(confile,Q[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  writeln('R');
  writeln(confile,'R');
  for i:=1 to nRow do
  begin
    for j:=1 to nCol do
    begin
      write(A[i,j]:10:5,' ');
      write(confile,A[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nCol) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  writeln('Q*R - Acopy');
  writeln(confile,'Q*R - Acopy');
  for i:=1 to nRow do
  begin
    for j:=1 to nCol do
    begin
      s:=-Acopy[i,j];
      for k:=1 to nRow do s:=s+Q[i,k]*A[k,j];
      write(s:10,' ');
      write(confile,s:10,' ');
      if (7 * (j div 7) = j) and (j<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr03.pas == Givens' reduction driver}
