{$I turbo.cnm}
Program dr15(input,output);
{dr15.pas == driver for generalized symmetric eigenvalue problem (alg15)

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I matrixin.pas}
{$I genevres.pas}
{$I rayquo.pas}
{$I alg14.pas}
{$I alg15.pas}
{$I startup.pas}

{main program}
var
  i, j, k, nRow, nCol : integer;
  s, s2, t1 : real;
  A, V, Acopy, B, Bcopy : rmatrix;
  Z, Y : rvector; {to test residuals}
  sym : boolean;
  ch : char;
  tvec : smatvec; {needed only for Matrixin}

begin
  banner:='dr15.pas -- generalized symmetric matrix eigenvalue problem';
  startup;
  write('Order of matrices = '); readln(infile,nCol);
  if infname<>'con' then writeln(nCol);
  writeln(confile,'Order of matrices = ',nCol);
  nRow := nCol;
  writeln('Provide matrix A');
  writeln(confile,'Provide matrix A');
  Matrixin(nRow, nCol, A, tvec, sym);
  if not sym then halt;
  for i := 1 to nRow do
  begin
    for j := 1 to nRow do
    begin
    write(A[i,j]:10:5,' ');
    write(confile,A[i,j]:10:5,' ');
    Acopy[i,j] := A[i,j];
    if (7 * (j div 7) = j) and (j<nRow) then
    begin
      writeln;
      writeln(confile);
    end;
    end;
    writeln;
    writeln(confile);
  end;
  writeln('Provide matrix B');
  writeln(confile,'Provide matrix B');
  Matrixin(nRow, nCol, B, tvec, sym);
  if not sym then halt;
  for i := 1 to nRow do
  begin
    for j := 1 to nRow do
    begin
    write(B[i,j]:10:5,' ');
    write(confile,B[i,j]:10:5,' ');
    Bcopy[i,j] := B[i,j];
    if (7 * (j div 7) = j) and (j<nRow) then
    begin
      writeln;
      writeln(confile);
    end;
    end;
    writeln;
    writeln(confile);
  end;
  genevJac( nRow, A, B, V);
  writeln;
  writeln(confile);
  for j := 1 to nRow do
  begin
    t1 := B[j,j];
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
    Z[i] := V[i,j]; {to initialize residual test}
    end;
    writeln;
    writeln(confile);
    s2 := genevres(nRow,Acopy,Bcopy,t1,Z,true);
    t1 := rayquo(nRow,Acopy,Bcopy,Z);
    writeln('Rayleigh quotient = ',t1);
    writeln(confile,'Rayleigh quotient = ',t1);
    s2 := genevres(nRow,Acopy,Bcopy,t1,Z,true);
    writeln;
    writeln(confile);
  end; {loop on solutions j}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr15.pas}
