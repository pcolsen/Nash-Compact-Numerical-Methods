{$I turbo.cnm}
program dr26(input,output);

{dr26.pas == eigensolutions of a complex matrix by Eberlein's
          complex Jacobi procedure

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I calceps.pas}
{$I matrixin.pas}
{$I alg26.pas}  {comeig.pas}
{$I alg11.pas}  {standardize complex eigensolutions}
{$I alg12.pas}  {complex residual calculation}
{$I startup.pas}

{Main program}
var
  A, Z, Acopy, Zcopy, T, U : rmatrix;
  i, it, j, k, n : integer;
  sym : boolean;
  avec : smatvec; {for compatibility of Matrixin only}

begin
  banner:='dr26.pas -- Eigensolutions of a general complex matrix';
  startup;
  write(' Order of matrix = '); readln(infile,n);

  writeln(confile,' Order of matrix = ',n);
  if infname<>'con' then writeln(n);
  writeln('Provide real part of matrix (A)');
  writeln(confile,'Provide real part of matrix (A)');
  matrixin(n,n,A,avec,sym);
  writeln('Provide imaginary part of matrix (Z)');
  writeln(confile,'Provide imaginary part of matrix (Z)');
  matrixin(n,n,Z,avec,sym);
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
    Acopy[i,j]:=A[i,j]; Zcopy[i,j]:=Z[i,j];
    write('(',A[i,j]:10:5,',',Z[i,j]:10:5,') ');
    write(confile,'(',A[i,j]:10:5,',',Z[i,j]:10:5,') ');
    if (3 * (j div 3) = j) and (j<n) then
    begin
      writeln;
      writeln(confile);
    end;
    end; {copy loop j}
    writeln;
    writeln(confile);
  end; {copy loop i}
  it:=50; {allow a maximum of 50 iterations}
  comeig( n, it, A, Z, T, U);
  if it>0 then writeln('Converged in ',it,' iterations')
    else writeln('Not converged after ',it,' iterations');
  if it>0 then
  writeln(confile,'Converged in ',it,' iterations')
    else
    writeln(confile,'Not converged after ',it,' iterations');
  stdceigv(n, T, U); {standardize the eigensolutions -- alg11.pas}
  for i:=1 to n do
  begin
    writeln;
    writeln(confile);
    writeln('EIGENVALUE ',i,'=(',A[i,i],',',Z[i,i],')');
    writeln(confile,'EIGENVALUE ',i,'=(',A[i,i],',',Z[i,i],')');
    writeln('VECTOR');
    writeln(confile,'VECTOR');
    for k:=1 to n do
    begin
    writeln('(',T[k,i],',',U[k,i],')');
    writeln(confile,'(',T[k,i],',',U[k,i],')');
    end; { loop on  k}
    comres( i, n, A, Z, T, U, Acopy, Zcopy); {residuals -- alg12.pas}
  end; {loop on i}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr26.pas == eigensolutions of a complex matrix}
