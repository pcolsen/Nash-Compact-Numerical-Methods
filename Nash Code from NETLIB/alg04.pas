procedure GivSVD( n : integer;
                nRHS: integer;
                var B: rmatrix;
                var rss: rvector;
                var svs: rvector;
                var W: rmatrix;
                var nobs : integer);


var
  count, EstRowRank, i, j, k, m, slimit, sweep, tcol : integer;
  bb, c, e2, eps, p, q, r, s, tol, trss, vt : real;
  enddata : boolean;
  endflag : real;

procedure rotnsub;
var
  i: integer;
  r: real;

begin
  for i := m to tcol do
  begin
    r := W[j,i];
    W[j,i] := r*c+s*W[k,i];
    W[k,i] := -r*s+c*W[k,i];
  end;
end;

begin
  writeln('alg04.pas -- Givens',chr(39),
         ' reduction, svd, and least squares solution');
  writeln(confile,'alg04.pas -- Givens',chr(39),
         ' reduction, svd, and least squares solution');
  write('Order of ls problem and no. of right hand sides: ');
  readln(infile,n,nRHS);
  if infname<>'con' then writeln(n,' ',nRHS);
  writeln(confile,'Order of ls problem and no. of right hand sides = ',
                      n,' ',nRHS);
  write('Enter a number to indicate end of data ');
  readln(infile,endflag);
  if infname<>'con' then writeln(endflag);
  writeln(confile,'Enter a number to indicate end of data ',endflag);
  tcol := n+nRHS;
  k := n+1;
  for i := 1 to n do
    for j := 1 to tcol do
      W[i,j] := 0.0;
  for i := 1 to nRHS do rss[i] := 0.0;

  eps := calceps;
  tol := n*n*eps*eps;
  nobs := 0;

  enddata := false;
  while (not enddata) do
  begin
    getobsn( n, nRHS, W, k, endflag, enddata);
    if (not enddata) then
    begin
      nobs := nobs+1;
      write('Obsn ',nobs,' ');
      for j := 1 to (n+nRHS) do
      begin
        write(W[k,j]:10:5,' ');
        if (7 * (j div 7) = j) and (j<n+nRHS) then writeln;
      end;
      writeln;
      write(confile,'Obsn ',nobs,' ');
      for j := 1 to (n+nRHS) do
      begin
        write(confile,W[k,j]:10:5,' ');
        if (7 * (j div 7) = j) and (j<n+nRHS) then writeln(confile);
      end;
      writeln(confile);
      for j := 1 to n do
      begin
        m := j; s := W[k,j]; c := W[j,j];
        bb := abs(c); if abs(s)>bb then bb := abs(s);
        if bb>0.0 then
        begin
          c := c/bb; s := s/bb; p := sqrt(c*c+s*s);
          s := s/p;
          if abs(s)>=tol then
          begin
            c := c/p;
            rotnsub;
          end;
        end;
      end;

      write('    Uncorrelated residual(s):');
      write(confile,'    Uncorrelated residual(s):');
      for j := 1 to nRHS do
      begin
        rss[j] := rss[j]+sqr(W[k,n+j]);
        write(W[k,n+j]:10,' ');
        write(confile,W[k,n+j]:10,' ');
        if (7 * (j div 7) = j) and (j < nRHS) then
        begin
          writeln; writeln(confile);
        end;
      end;
      writeln; writeln(confile);

    end;
  end;


  m := 1;
  slimit := n div 4;  if slimit<6 then slimit := 6;

  sweep := 0;
  e2 := 10.0*n*eps*eps;
  tol := eps*0.1;
  EstRowRank := n; ;
  repeat
    count := 0;
    for j := 1 to (EstRowRank-1) do
    begin
      for k := (j+1) to EstRowRank do
      begin
        p := 0.0; q := 0.0; r := 0.0;
        for i := 1 to n do
        begin
          p := p+W[j,i]*W[k,i]; q := q+sqr(W[j,i]); r := r+sqr(W[k,i]);
        end;
        svs[j] := q; svs[k] := r;



        IF q >= r then
        begin

          if not ((q<=e2*svs[1]) or (abs(p)<=tol*q)) then

          begin

            p := p/q; r := 1-r/q; vt := sqrt(4*p*p + r*r);
            c := sqrt(0.5*(1+r/vt)); s := p/(vt*c);
            rotnsub;
            count := count+1;
          end;
        end
        else

        begin
          p := p/r; q := q/r-1; vt := sqrt(4*p*p + q*q);
          s := sqrt(0.5*(1-q/vt));
          if p<0 then s := -s;
          c := p/(vt*s);
          rotnsub;
          count := count+1;
        end;


      end;
    end;
    sweep := sweep +1;
    writeln('Sweep ',sweep,' ',count,' rotations performed');
    writeln(confile,'Sweep ',sweep,' ',count,' rotations performed');

      {JN941126}
      {Remove the next "while" statement to force higher precision
       computation in cases of severe multicollinearity. Otherwise do
       NOT use results in cases where svs[n]/svs[1] is very small.}

     while (EstRowRank >= 3) and (svs[EstRowRank] <= svs[1]*tol+tol*tol)
          do EstRowRank := EstRowRank-1;
  until (sweep>slimit) or (count=0);

  writeln('Singular values and principal components');
  writeln(confile,'Singular values and principal components');
  for j := 1 to n do
  begin
    s := svs[j];
    s := sqrt(s); svs[j] := s;
    writeln('Singular value [',j,']= ',s);
    writeln(confile,'Singular value [',j,']= ',s);
    if s>=tol then
    begin
      for i := 1 to n do W[j,i] := W[j,i]/s;
      for i := 1 to n do
      begin
        write(W[j,i]:8:5,' '); write(confile,W[j,i]:8:5,' ');
        if (8 * (i div 8) = i) and (i<n) then
        begin
          writeln; writeln(confile);
        end;
      end;

      writeln; writeln(confile);
    end;

  end;

  q := 0.0;
  while q>=0.0 do
  begin
    write('Enter a tolerance for zero (<0 to exit) ');
    readln(infile,q);
    writeln(confile,'Enter a tolerance for zero (<0 to exit) ',q);
    if infname<>'con' then writeln(q);
    if q>=0.0 then
    begin

      for i := 1 to nRHS do
      begin
        trss := rss[i];
        for j := 1 to n do
        begin
          p := 0.0;
          for k := 1 to n do
          begin
            if svs[k]>q then p := p+W[k,j]*W[k,n+i]/svs[k];
          end;
          B[j,i] := p;
          writeln('Solution component [',j,']= ',p);
          writeln(confile,'Solution component [',j,']= ',p);
          if svs[j]<=q then trss := trss+sqr(W[j,n+i]);
        end;
        writeln('Residual sum of squares =',trss);
        writeln(confile,'Residual sum of squares =',trss);
      end;
    end;
  end;
end;
