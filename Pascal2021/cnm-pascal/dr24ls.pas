{$I turbo.cnm}
program dr24ls(input, output);
{dr24ls.pas == linear least squares by conjugate gradients

Note: this implementation uses the normal equations, which
are not recommended as a general approach.

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
  A, AtransA : rmatrix;
  Y, AtransY : rvector; {RHS}
  Bvec : rvector; {solution}
  avec : smatvec; {for matrixin only}
  sym : boolean; {to tell if matrix symmetric}
  ch : char;
  i, j, k, nRow, nCol, itcount : integer;
  ssmin, t, s : real;

begin
  banner:='dr24ls.pas -- linear least squares by conjugate gradients';
  startup;
  write('Number of rows in coefficient matrix = ');
  readln(infile,nRow);
  if infname<>'con' then writeln(nRow);
  writeln(confile,'Number of rows in coefficient matrix = ',nRow);
  write('Number of columns in coefficient matrix = ');
  readln(infile,nCol); {This is the order of the conjugate gradients problem.}
  if infname<>'con' then writeln(nCol);
  writeln(confile,'Number of columns in coefficient matrix = ',nCol);
  writeln('Coefficient matrix');
  matrixin(nRow, nCol, A, avec, sym);
  writeln('RHS vector');
  writeln(confile,'RHS vector');
  vectorin(nRow, Y);
  writeln('Initial guess for solution');
  writeln(confile,'Initial guess for solution');
  vectorin(nCol, Bvec);
  {Now form the normal equations.}
  writeln('Normal equations -- coefficient matrix');
  writeln(confile,'Normal equations -- coefficient matrix');
  for i:=1 to nCol do
  begin
    t:=0.0;
    for k:=1 to nRow do t:=t+A[k,i]*Y[k];
    AtransY[i]:=t;
    for j:=1 to nCol do
    begin
    s:=0.0; for k:=1 to nRow do s:=s+A[k,i]*A[k,j];
    AtransA[i,j]:=s;
    write(s:10:5,' ');
    write(confile,s:10:5,' ');
    if (7 * (j div 7) = j) and (j<nCol) then
    begin
      writeln;
      writeln(confile);
    end;
    end;
    writeln;
    writeln(confile);
  end; {loop on i and normal equations build}
  writeln('Normal equations - RHS');
  writeln(confile,'Normal equations - RHS');
  for j:=1 to nCol do
  begin
    write(AtransY[j]:10:5,' ');
    write(confile,AtransY[j]:10:5,' ');
    if (7 * (j div 7) = j) and (j<nCol) then
    begin
    writeln;
    writeln(confile);
    end;
  end;
  writeln;
  writeln(confile);
  {***WARNING*** this is NOT a good way to solve this problem generally}
  itcount:=10*nRow; {safety setting}
  lecg( nCol, AtransA, AtransY, Bvec, itcount, ssmin);
  writeln('Solution after ',itcount,
      ' iterations. Est. normal eqn. sumsquares ',ssmin);
  writeln(confile,'Solution after ',itcount,
      ' iterations. Est. normal eqn. sumsquares ',ssmin);
  for i:=1 to nCol do
  begin
    write(Bvec[i]:10:5,' ');
    write(confile,Bvec[i]:10:5,' ');
    if (7 * (i div 7) = i) and (i<nCol) then
    begin
    writeln;
    writeln(confile);
    end;
  end;
  writeln;
  writeln(confile);
  write('For original least squares problem -- ');
  write(confile,'For original least squares problem -- ');
  s:=resids(nRow, nCol, A, Y, Bvec, true);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr24ls.pas}
