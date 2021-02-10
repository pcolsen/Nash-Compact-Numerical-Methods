procedure gebacksub(n, p:integer;
                     var A : rmatrix);


var
  s : real;
  i, j, k: integer;

begin
  writeln('alg06.pas -- Gauss elimination back-substitution');
  writeln(confile,'alg06.pas -- Gauss elimination back-substitution');
  for i:=(n+1) to (n+p) do
  begin
    A[n,i]:=A[n,i]/A[n,n];
    for j:=(n-1) downto 1 do
    begin
      s:=A[j,i];
      for k:=(j+1) to n do
      begin
        s:=s-A[j,k]*A[k,i];
      end;
      A[j,i]:=s/A[j,j];
    end;
  end;
end;
