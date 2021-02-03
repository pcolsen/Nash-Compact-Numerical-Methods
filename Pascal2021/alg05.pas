Procedure gelim( n : integer;
                 p : integer;
                var A : rmatrix;
                 tol : real);
var
  det, s : real;
  h,i,j,k: integer;

begin
  det := 1.0;
  writeln('alg05.pas -- Gauss elimination with partial pivoting');
  for j := 1 to (n-1) do
  begin
    s := abs(A[j,j]); k := j;
    for h := (j+1) to n do
    begin
      if abs(A[h,j])>s then
      begin
        s := abs(A[h,j]); k := h;
      end;
    end;
    if k<>j then
    begin
      writeln('Interchanging rows ',k,' and ',j);
      for i := j to (n+p) do
      begin
        s := A[k,i]; A[k,i] := A[j,i]; A[j,i] := s;
      end;
      det := -det;
    end;
    det := det*A[j,j];
    if abs(A[j,j])<tol then
    begin
      writeln('Matrix computationally singular -- pivot < ',tol);
      halt;
    end;
    for k := (j+1) to n do
    begin
      A[k,j] := A[k,j]/A[j,j];
      for i := (j+1) to (n+p) do
          A[k,i] := A[k,i]-A[k,j]*A[j,i];
    end;
    det := det*A[n,n];
    if abs(A[n,n])<tol then
    begin
      writeln('Matrix computationally singular -- pivot < ',tol);
      halt;
    end;
  end;
  writeln('Gauss elimination complete -- determinant = ',det);
end;

