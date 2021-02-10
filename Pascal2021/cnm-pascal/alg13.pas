Procedure evsvd(n: integer;
             var A,V : rmatrix;
             initev: boolean;
             W : wmatrix;
             var Z: rvector);

var
  count, i, j, k, limit, skipped : integer;
  c, p, q, s, shift, t : real ;
  oki, okj, rotn : boolean;
  ch : char;

begin
  writeln('alg13.pas -- symmetric matrix eigensolutions via svd');
  writeln(confile,'alg13.pas -- symmetric matrix eigensolutions via svd');

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
  writeln(confile,'Adding a shift of ',shift,' to diagonal of matrix.');
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
  NashSVD(n, n, W, Z);
  for i:=1 to n do
  begin
    Z[i]:=sqrt(Z[i])-shift;
    for j:=1 to n do
      V[i,j]:=W[n+i,j];
  end;
end;
