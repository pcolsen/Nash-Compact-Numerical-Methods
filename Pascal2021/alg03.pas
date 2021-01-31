procedure givens( nRow,nCol : integer;
                  var A, Q: rmatrix);
var
 i, j, k, mn: integer;
 b, c, eps, p, s : real;

begin
  writeln('alg03.pas -- Givens',chr(39),' reduction -- column-wise');
  mn  :=  nRow; if nRow>nCol then mn  :=  nCol;
  for i  :=  1 to nRow do
  begin
    for j := 1 to nRow do Q[i,j] := 0.0;
    Q[i,i] := 1.0;
  end;
  eps := calceps;
  for j := 1 to (mn-1) do
  begin
    for k := (j+1) to nRow do
    begin
      c := A[j,j]; s := A[k,j];
      b := abs(c); if abs(s)>b then b := abs(s);
      if b>0 then
      begin
        c := c/b; s := s/b;
        p := sqrt(c*c+s*s);
        s := s/p;
        if abs(s)>=eps then
        begin
          c := c/p;
          for i := 1 to nCol do
          begin
            p := A[j,i]; A[j,i] := c*p+s*A[k,i]; A[k,i] := -s*p+c*A[k,i];
          end;
          for i := 1 to nRow do
          begin
            p := Q[i,j]; Q[i,j] := c*p+s*Q[i,k]; Q[i,k] := -s*p+c*Q[i,k];
          end;
        end;
      end;
    end;
  end;
end;


