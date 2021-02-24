function resids(nRow, nCol: integer; A : rmatrix;
          Y: rvector; Bvec : rvector; print : boolean):real;
{resids.pas
  == Computes residuals and , if print is TRUE, displays them 7
    per line for the linear least squares problem. The sum of
    squared residuals is returned.

    residual vector = A * Bvec - Y
}
var
i, j: integer;
t1, ss : real;

begin
  if print then
  begin
    writeln('Residuals');
    writeln(confile,'Residuals');
  end;
  ss:=0.0;
  for i:=1 to nRow do
  begin
    t1:=-Y[i]; {note form of residual is residual = A * B - Y }
    for j:=1 to nCol do
      t1:=t1+A[i,j]*Bvec[j];
    ss:=ss+t1*t1;
    if print then
    begin
      write(t1:10,' ');
      write(confile,t1:10,' ');
      if (i = 7 * (i div 7)) and (i<nRow) then writeln;
      if (i = 7 * (i div 7)) and (i<nRow) then
      writeln(confile);
    end;
  end; {loop on i}
  if print then
  begin
    writeln;
    writeln(confile);
    writeln('Sum of squared residuals =',ss);
    writeln(confile,'Sum of squared residuals =',ss);
  end;
  resids:=ss
end; {resids.pas == residual calculation for linear least squares}
