procedure brspdmi(n : integer;
                var avector : smatvec;
                var singmat : boolean);

var
  i,j,k,m,q : integer;
  s,t : real;
  X : rvector;

begin
  writeln('alg09.pas -- Bauer Reinsch inversion');
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

var
  A, Ainverse : rmatrix;
  avector : smatvec;
  i, imax, j, jmax, k, n : integer;
  errmax, s : real;
  singmat: boolean;


BEGIN { main program }
  banner:='dr09.pas -- test Bauer Reinsch sym, posdef matrix inversion';
  writeln(banner);
  n:=4; {Fixed example size 20210113}
  Frank(n,A,avector);
  writeln;
  writeln('returned matrix of order ',n);
  begin
    for i:=1 to n do
    begin
        for j:=1 to n do
        begin
            write(A[i,j],' ');
        end;
        writeln;
    end;
  end;
  mat2vec(n, A, avector);
  begin
    writeln('Symmetric matrix -- Vector form');
    k := 0;
    for i := 1 to n do
    begin
      for j := 1 to i do
      begin
        k := k+1;
        write(avector[k]:10:5,' ');
      end;
      writeln;
    end;
  end;
  brspdmi(n, avector,singmat);
  if singmat then halt; {safety check}
  writeln('Computed inverse');
  k := 0; {initialize index to smatvec elements}
  for i := 1 to n do
  begin
    for j := 1 to i do
    begin
      k := k+1;
      write(avector[k]:10:5,' ');
      Ainverse[i,j] := avector[k]; {save square form of inverse}
      Ainverse[j,i] := avector[k];
      if (7 * (j div 7) = j) and (j<i) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
  {Compute maximum error in A * Ainverse and note where it occurs.}
  errmax := 0.0; imax := 0; jmax := 0;
  for i := 1 to n do
  begin
    for j := 1 to n do
    begin
      s := 0.0; if i=j then s := -1.0;
      for k := 1 to n do s := s + Ainverse[i,k]*A[k,j];
      {Note: A has not been altered, since avector was used.}
      if abs(s)>abs(errmax) then
      begin
        errmax := s; imax := i; jmax := j; {save maximum error, indices}
      end;
    end; {loop on j}
  end; {loop on i}
  writeln('Maximum element in Ainverse * A - 1(n) = ',errmax,
          '  position ',imax,',',jmax);
end. {dr09.pas == Bauer Reinsch inversion}
