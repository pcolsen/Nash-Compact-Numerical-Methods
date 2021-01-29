procedure brspdmi(n : integer;
                var avector : smatvec;
                var singmat : boolean);

var
  i,j,k,m,q : integer;
  s,t : real;
  X : rvector;

begin
  writeln('alg09.pas -- Bauer Reinsch inversion');
  writeln(confile,'alg09.pas -- Bauer Reinsch inversion');
  singmat  :=  false;
  for k  :=  n downto 1 do
  begin
    if (not singmat) then
    begin
      s  :=  avector[1];
      if s>0.0 then
      begin
        m  :=  1;
        for i := 2 to n do
        begin
          q := m; m := m+i; t := avector[q+1]; X[i] := -t/s;

          if i>k then X[i] := -X[i];
          for j := (q+2) to m do
          begin
            avector[j-i] := avector[j]+t*X[j-q];
          end;
        end;
        q := q-1; avector[m] := 1.0/s;
        for i := 2 to n do avector[q+i] := X[i];
      end
      else
        singmat := true;
    end;
  end;
end;
