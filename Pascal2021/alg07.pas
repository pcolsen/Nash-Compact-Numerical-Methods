procedure choldcmp(n: integer;
                   var a: smatvec;
                   var singmat: boolean);
var
  i,j,k,m,q: integer;
  s : real;

begin
  singmat := false;
  for j := 1 to n do
  begin
    q := j*(j+1) div 2;
    if j>1 then
    begin
      for i := j to n do
      begin
        m := (i*(i-1) div 2)+j; s := a[m];
        for k := 1 to (j-1) do s := s-a[m-k]*a[q-k];
        a[m] := s;
      end;
    end;
    if a[q]<=0.0 then
    begin
      singmat := true;
      a[q] := 0.0;
    end;
    s := sqrt(a[q]);
    for i := j to n do
    begin
      m := (i*(i-1) div 2)+j;
      if s=0.0 then a[m] := 0
          else a[m] := a[m]/s;
    end;
  end;
end;

