function genevres(n: integer; {order of matrices}
              A, B: rmatrix;
              evalue : real; {eigenvalue}
              X : rvector; {presumed eigenvector}
              print: boolean) {true for printing}
              : real; {returns sum of squared residuals}

{genevres.pas

  -- to compute residuals of generalized (symmetric) matrix
  eigenvalue problem

        A x = evalue B x
}

var
  t, ss : real;
  i,j : integer;

begin
  if print then
  begin
    writeln('Generalized matrix eigensolution residuals');
    writeln(confile,'Generalized matrix eigensolution residuals');
  end;
  ss:=0.0; {to accumulate the sum of squared residuals}
  for i:=1 to n do
  begin
    t:=0.0;
    for j:=1 to n do
    t:=t+(A[i,j]-evalue*B[i,j])*X[j];
    if print then
    begin
      write(t:10,' ');
      write(confile,t:10,' ');
    if (7 * (i div 7) = i) and (i<n) then
    begin
      writeln;
      writeln(confile);
    end;
    end; {if print}
    ss:=ss+t*t;
  end;
  if print then
  begin
    writeln;
    writeln(confile);
    writeln('Sum of squared residuals =',ss);
    writeln(confile,'Sum of squared residuals =',ss);
  end; {if print}
  genevres:=ss; {return sum of squared residuals}
end; {genevres.pas == residuals for generalized eigenproblem}
