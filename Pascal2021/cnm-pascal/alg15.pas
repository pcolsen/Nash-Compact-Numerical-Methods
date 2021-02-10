procedure genevJac( n : integer;
                var A, B, V : rmatrix);

var
 i,j,k,m : integer;
 s : real;
 initev : boolean;

begin
  writeln('alg15.pas -- generalized symmetric matrix eigenproblem');
  writeln(confile,
         'alg15.pas -- generalized symmetric matrix eigenproblem');
  initev := true;
  writeln('Eigensolutions of metric B');
  writeln(confile,'Eigensolutions of metric B');
  evJacobi(n, B, V, initev);



  for i := 1 to n do
  begin
    if B[i,i]<=0.0 then halt;
    s := 1.0/sqrt(B[i,i]);
    for j := 1 to n do V[j,i] := s * V[j,i];
  end;


  for i := 1 to n do
  begin
    for j := i to n do
    begin
      s := 0.0;
      for k := 1 to n do
        for m := 1 to n do
          s := s+V[k,i]*A[k,m]*V[m,j];
      B[i,j] := s; B[j,i] := s;
    end;
  end;

  initev := false;
  writeln('Eigensolutions of general problem');
  writeln(confile,'Eigensolutions of general problem');
  evJacobi( n, B, V, initev);
end;
