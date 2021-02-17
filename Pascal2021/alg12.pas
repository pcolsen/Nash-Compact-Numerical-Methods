procedure comres( i, n: integer;
                  A, Z, T, U, Acopy, Zcopy : rmatrix);

var
  j, k: integer;
  g, s, ss : real;

begin
  writeln('alg12.pas -- complex eigensolution residuals');
  ss := 0.0;
  for j := 1 to n do
  begin
    s := -A[i,i]*T[j,i]+Z[i,i]*U[j,i]; g := -Z[i,i]*T[j,i]-A[i,i]*U[j,i];

    for k := 1 to n do
    begin
      s := s+Acopy[j,k]*T[k,i]-Zcopy[j,k]*U[k,i];
      g := g+Acopy[j,k]*U[k,i]+Zcopy[j,k]*T[k,i];
    end;
    writeln('(',s,',',g,')');
    ss := ss+s*s+g*g;
  end;
  writeln('Sum of squares = ',ss);
  writeln;
end;
