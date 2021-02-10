{$I turbo.cnm}
program dr25(input, output);
{dr25.pas == eigensolutions by minimisation of the Rayleigh quotient

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}
{$I matrixin.pas}
{$I vectorin.pas}
{$I matmul.pas}
{$I alg25.pas}
{$I resids.pas}
{$I startup.pas}

var
  A, B : rmatrix;
  X : rvector; {eigenvector}
  Y : rvector; {for residuals}
  avec : smatvec; {for matrixin only}
  sym : boolean; {to tell if matrix symmetric}
  ch : char;
  i, j, n, itcount : integer;
  ev, t, s : real;

begin
  banner:='dr25.pas -- minimise Rayleigh quotient';
  startup;
  write('Order of problem =');
  readln(infile,n); if infname<>'con' then writeln(n);
  writeln(confile,'Order of problem = ',n);
  writeln('Matrix A');
  writeln(confile,'Matrix A');
  matrixin(n, n, A, avec, sym);
  if not sym then
  begin
    writeln('Matrix not symmetric -- halting');
    writeln(confile,'Matrix not symmetric -- halting');
    halt;
  end;
  writeln('Metric matrix B');
  writeln(confile,'Metric matrix B');
  matrixin(n, n, B, avec, sym);
  if not sym then
  begin
    writeln('Matrix not symmetric -- halting');
    writeln(confile,'Matrix not symmetric -- halting');
    halt;
  end;
  writeln('Initial eigenvector approximation');
  writeln(confile,'Initial eigenvector approximation');
  vectorin(n, X);
  itcount:=100*n; {safety setting}
  rqmcg( n, A, B, X, itcount, ev);
  writeln('Solution after ',itcount,' products. Est. eigenvalue =',ev);
  writeln(confile,
        'Solution after ',itcount,' products. Est. eigenvalue =',ev);
  for i:=1 to n do
  begin
    write(X[i]:10:7,' ');
    write(confile,X[i]:10:7,' ');
    if (7 * (i div 7) = i) and (i<n) then
    begin
    writeln;
    writeln(confile);
    end;
    t:=0.0;
    for j:=1 to n do t:=t+B[i,j]*X[j];
    Y[i]:=ev*t; {to save eigenvalue * matrix-vector product for residuals}
  end;
  writeln;
  writeln(confile);
  s := resids(n, n, A, Y, X, true);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr25.pas}

