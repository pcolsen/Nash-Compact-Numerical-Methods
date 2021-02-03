Procedure evJacobi(n: integer;
               var A,V : rmatrix;
                   initev: boolean);
var
  count, i, j, k, limit, skipped : integer;
  c, p, q, s, t : real;
  oki, okj, rotn : boolean;

begin
  writeln('alg14.pas -- eigensolutions of a real symmetric');
  writeln('             matrix via a Jacobi method');
  if initev then
  begin
    for i := 1 to n do
    begin
      for j := 1 to n do V[i,j] := 0.0;
      V[i,i] := 1.0;
    end;
  end;
  count := 0;
  limit := 30;
  skipped := 0;

  while (count<=limit) and (skipped<((n*(n-1)) div 2) ) do

  begin
    count := count+1;
    write('sweep ',count,' '); 
    skipped := 0;
    for i := 1 to (n-1) do
    begin
      for j := (i+1) to n do
      begin
        rotn := true;
        p := 0.5*(A[i,j]+A[j,i]);
        q := A[i,i]-A[j,j];
        t := sqrt(4.0*p*p+q*q);
        if t=0.0 then
        begin
          rotn := false;
        end
        else
        begin
          if q>=0.0 then
          begin

            oki := (abs(A[i,i])=abs(A[i,i])+100.0*abs(p));
            okj := (abs(A[j,j])=abs(A[j,j])+100.0*abs(p));
            if oki and okj then rotn := false
                           else rotn := true;

            if rotn then
            begin
               c := sqrt((t+q)/(2.0*t)); s := p/(t*c);
            end;
          end
          else
          begin
            rotn := true;
            s := sqrt((t-q)/(2.0*t));
            if p<0.0 then s := -s;
            c := p/(t*s);
          end;
          if 1.0=(1.0+abs(s)) then rotn := false;
        end;
        if rotn then
        begin
          for k := 1 to n do
          begin
            q := A[i,k]; A[i,k] := c*q+s*A[j,k]; A[j,k] := -s*q+c*A[j,k];
          end;

          for k := 1 to n do
          begin
            q := A[k,i]; A[k,i] := c*q+s*A[k,j]; A[k,j] := -s*q+c*A[k,j];

            q := V[k,i]; V[k,i] := c*q+s*V[k,j]; V[k,j] := -s*q+c*V[k,j];
          end;
        end
        else

           skipped := skipped+1;
      end;
    end;
    writeln('  ',skipped,' / ',n*(n-1) div 2,'  rotations skipped');
  end;
end;

