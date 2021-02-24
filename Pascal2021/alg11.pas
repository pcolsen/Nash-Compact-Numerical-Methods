procedure stdceigv(n: integer;
                var T, U: rmatrix);

var
  i, k, m : integer;
  b, e, g, s : real;

begin
  writeln('alg11.pas -- standardized eigensolutions');
  for i := 1 to n do
  begin
    g := T[1,i]*T[1,i]+U[1,i]*U[1,i];

    k := 1;
    if n>1 then
    begin
      for m := 2 to n do
      begin
        b := T[m,i]*T[m,i]+U[m,i]*U[m,i];
        if b>g then
        begin
          k := m;
          g := b;
        end;
      end;
    end;
    e := T[k,i]/g;
    s := -U[k,i]/g;
    for k := 1 to n do
    begin
      g := T[k,i]*e-U[k,i]*s; U[k,i] := U[k,i]*e+T[k,i]*s; T[k,i] := g;
    end;
  end;
end;
