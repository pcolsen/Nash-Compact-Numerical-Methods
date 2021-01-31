Procedure evsvd(n: integer; var A,V : rmatrix; initev: boolean;
             W : wmatrix; var Z: rvector);

var
  i, j: integer;
  shift, t : real ;
  
begin
  writeln('alg13.pas -- symmetric matrix eigensolutions via svd');
  shift:=0.0;
  for i:=1 to n do
  begin
    t:=A[i,i];
    for j:=1 to n do
      if i<>j then t:=t-abs(A[i,j]);
    if t<shift then shift:=t;
  end;
  shift:=-shift;
  if shift<0.0 then shift:=0.0;
  writeln('Adding a shift of ',shift,' to diagonal of matrix.');
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
      W[i,j]:=A[i,j];
      if i=j then W[i,i]:=A[i,i]+shift;
      if initev then
      begin
        if i=j then W[i+n,i]:=0.0
        else
        begin
          W[i+n,j]:=0.0;
        end;
      end;
    end;
  end;
  if (n > 1) then 
     NashSVD(n, n, W, Z)
  else
  begin { order 1 matrix }
     Z[1] := A[1,1]*A[1,1];
     W[2,1]:= 1.0; {Eigenvector!}
  end;
  for i:=1 to n do
  begin
    Z[i]:=sqrt(Z[i])-shift;
    for j:=1 to n do
      V[i,j]:=W[n+i,j];
  end;
end;

