{$I turbo.cnm}
program dr09(input,output);
{dr09.pas == driver program to test procedure for the Bauer-Reinsch
          inversion of a symmetric positive definite real matrix stored
          in row-wise vector form

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I matrixin.pas}
{$I vectorin.pas}
{$I alg09.pas}
{$I startup.pas}

var
  A, Ainverse : rmatrix;
  avector : smatvec;
  i, imax, j, jmax, k, nCol, nRow : integer;
  errmax, s : real;
  singmat, sym : boolean;

begin
  banner:='dr09.pas -- test Bauer Reinsch inversion of a';
  startup;
  write('Order of matrix = '); readln(infile,nCol);
  if infname<>'con' then writeln(nCol);
  writeln(confile,'Order of matrix = ',nCol);
  nRow := nCol; {use symmetric matrix in Choleski}
  Matrixin(nRow,nCol,A,avector,sym);
  writeln;
  writeln(confile);
  writeln('returned matrix of order ',nRow);
  writeln(confile,'returned matrix of order ',nRow);
  if not sym then halt; {must have symmetric matrix}
  begin
    writeln('Symmetric matrix -- Vector form');
    writeln(confile,'Symmetric matrix -- Vector form');
    k := 0;
    for i := 1 to nRow do
    begin
      for j := 1 to i do
      begin
        k := k+1;
        write(avector[k]:10:5,' ');
        write(confile,avector[k]:10:5,' ');
        {Note: A[i,j] from matrixin is unaltered.}
        if (7 * (j div 7) = j) and (j<i) then
        begin
          writeln;
          writeln(confile);
        end;
      end;
      writeln;
      writeln(confile);
    end;
  end;
  brspdmi(nRow, avector,singmat);
  if singmat then halt; {safety check}
  writeln('Computed inverse');
  writeln(confile,'Computed inverse');
  k := 0; {initialize index to smatvec elements}
  for i := 1 to nRow do
  begin
    for j := 1 to i do
    begin
      k := k+1;
      write(avector[k]:10:5,' ');
      write(confile,avector[k]:10:5,' ');
      Ainverse[i,j] := avector[k]; {save square form of inverse}
      Ainverse[j,i] := avector[k];
      if (7 * (j div 7) = j) and (j<i) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  {Compute maximum error in A * Ainverse and note where it occurs.}
  errmax := 0.0; imax := 0; jmax := 0;
  for i := 1 to nRow do
  begin
    for j := 1 to nRow do
    begin
      s := 0.0; if i=j then s := -1.0;
      for k := 1 to nRow do s := s + Ainverse[i,k]*A[k,j];
      {Note: A has not been altered, since avector was used.}
      if abs(s)>abs(errmax) then
      begin
        errmax := s; imax := i; jmax := j; {save maximum error, indices}
      end;
    end; {loop on j}
  end; {loop on i}
  writeln('Maximum element in Ainverse * A - 1(n) = ',errmax,
          '  position ',imax,',',jmax);
  writeln(confile,'Maximum element in Ainverse * A - 1(n) = ',errmax,
          '  position ',imax,',',jmax);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr09.pas == Bauer Reinsch inversion}
