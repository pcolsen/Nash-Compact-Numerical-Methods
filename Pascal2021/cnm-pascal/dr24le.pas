{$I turbo.cnm}
program dr24le(input, output);
{dr24le.pas == linear equations by conjugate gradients

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}
{$I matrixin.pas}
{$I vectorin.pas}
{$I matmul.pas}
{$I alg24.pas}
{$I resids.pas}
{$I startup.pas}

var
  A : rmatrix;
  Y : rvector; {RHS}
  Bvec : rvector; {solution}
  avec : smatvec; {for matrixin only}
  sym : boolean; {to tell if matrix symmetric}
  ch : char;
  i, j, n, itcount : integer;
  ssmin, t, s : real;

begin
  banner:='dr24le.pas -- linear equations by conjugate gradients';
  startup;
  write('Order of problem = ');
  readln(infile,n);
  if infname<>'con' then writeln(n);
  writeln(confile,'order of matrix = ',n);
  writeln('Coefficient matrix');
  matrixin(n, n, A, avec, sym);
  if not sym then halt;
  writeln('RHS vector');
  writeln(confile,'RHS vector');
  vectorin(n, Y);
  writeln('Initial guess for solution');
  writeln(confile,'Initial guess for solution');
  vectorin(n, Bvec);
  itcount:=10*n; {safety setting}
  lecg( n, A, Y, Bvec, itcount, ssmin);
  writeln('Solution after ',itcount,' iterations. Est. sumsquares ',ssmin);
  writeln(confile,
        'Solution after ',itcount,' iterations. Est. sumsquares ',ssmin);
  for i:=1 to n do
  begin
    write(Bvec[i]:10:5,' ');
    write(confile,Bvec[i]:10:5,' ');
    if (7 * (i div 7) = i) and (i<n) then
    begin
    writeln;
    writeln(confile);
    end;
  end;
  writeln;
  writeln(confile);
  s:=resids(n, n, A, Y, Bvec, true);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr24le.pas}
