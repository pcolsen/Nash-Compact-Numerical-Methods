procedure gii(nRow : integer;
             var A : rmatrix;
             var Y : rvector;
             var shift : real;
             var itcount: integer);

var
  i, itlimit, j, m, msame, nRHS :integer;
  ev, s, t, tol : real;
  X : rvector;

begin
  itlimit:=itcount;
  nRHS:=nRow;
  tol:=Calceps;
  s:=0.0;
  for i:=1 to nRow do
  begin
    X[i]:=Y[i];
    Y[i]:=0.0;
    for j:=1 to nRow do
    begin
      A[i,j]:=A[i,j]-shift*A[i,j+nRow];
      s:=s+abs(A[i,j]);
    end;
  end;
  tol:=tol*s;
  gelim(nRow, nRHS, A, tol);
  itcount:=0;
  msame :=0;
  while (msame<nRow) and (itcount<itlimit) do
  begin
    itcount:=itcount+1;
    m:=nRow; s:=X[nRow];
    X[nRow]:=Y[nRow];
    if abs(A[nRow,nRow])<tol then Y[nRow]:=s/tol
                             else Y[nRow]:=s/A[nRow,nRow];
    t:=abs(Y[nRow]);
    for i:=(nRow-1) downto 1 do
    begin
      s:=X[i]; X[i]:=Y[i];
      for j:=(i+1) to nRow do
      begin
        s:=s-A[i,j]*Y[j];
      end;
      if abs(A[i,i])<tol then Y[i]:=s/tol else  Y[i]:=s/A[i,i];
      if abs(Y[i])>t then
      begin
        m:=i; t:=abs(Y[i]);
      end;
    end;
    ev:=shift+X[m]/Y[m];
(*    writeln('Iteration ',itcount,'  approx. ev=',ev);*)

    t:=Y[m]; msame:=0;
    for i:=1 to nRow do
    begin
      Y[i]:=Y[i]/t;
      if reltest+Y[i] = reltest+X[i] then msame:=msame+1;

    end;

    if msame<nRow then
    begin
      for i:=1 to nRow do
      begin
        s:=0.0;
        for j:=1 to nRow do s:=s+A[i,j+nRow]*Y[j];
        X[i]:=s;
      end;
    end;
  end;
  if itcount>=itlimit then itcount:=-itcount;
  shift:=ev;
end;

