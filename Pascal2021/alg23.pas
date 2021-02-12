procedure modmrt( n : integer;
          var Bvec : rvector;
          var X : rvector;
          var Fmin  : real;
            Workdata : probdata);

{modified 1991 - 01 - 13}
var
  a, c: smatvec;
  delta, v : rvector;
  dec, eps, inc, lambda, p, phi, res : real;
  count, i, ifn, igrad, j, k, nn2, q : integer;
  notcomp, singmat, calcmat: boolean;

begin
  writeln('alg23.pas -- Nash Marquardt nonlinear least squares');
  with Workdata do
  begin
    if nlls = false then halt;
    Fmin:=big;
    inc:=10.0;
    dec:=0.4;
    eps:=calceps;
    lambda:=0.0001;
    phi:=1.0;
    ifn:=0;  igrad:=0;
    calcmat:=true;
    nn2:=(n*(n+1)) div 2;
    p:=0.0;
    for i:=1 to m do
    begin
      res:=nlres(i, n, Bvec, notcomp, Workdata);

      if notcomp then halt;
      p:=p+res*res;
    end;
    ifn:=ifn+1;
    Fmin:=p;
    count:=0;

    while count<n do
    begin

      if calcmat then
      begin
        writeln(igrad,' ',ifn,'  sum of squares=',Fmin);
        for i:=1 to n do
        begin
          write(Bvec[i]:10:5,' ');
          if (7 * (i div 7) = i) and (i<n) then writeln;
        end;
        writeln;
        igrad:=igrad+1;
        for j:=1 to nn2 do a[j]:=0.0;
        for j:=1 to n do v[j]:=0.0;
        for i:=1 to m do
        begin
          nljac(i, n, Bvec, X, workdata);
          res:=nlres(i, n, Bvec, notcomp, Workdata);
          for j:=1 to n do
          begin
            v[j]:=v[j]+X[j]*res;
            q:=(j*(j-1)) div 2;
            for k:=1 to j do a[q+k]:=a[q+k]+X[j]*X[k];
          end;
        end;
        for j:=1 to nn2 do c[j]:=a[j];
        for j:=1 to n do X[j]:=Bvec[j];
      end;
      writeln('LAMDA =',lambda:8);
      for j:=1 to n do
      begin
        q:=(j*(j+1)) div 2;
        a[q]:=c[q]*(1.0+lambda)+phi*lambda;
        delta[j]:=-v[j];
        if j>1 then
          for i:=1 to (j-1) do a[q-i]:=c[q-i];
      end;
      notcomp:=false;
      Choldcmp(n, a, singmat);
      if (not singmat) then
      begin
        Cholback(n, a, delta);
        count:=0;
        for i:=1 to n do
        begin
          Bvec[i]:=X[i]+delta[i];
          if (reltest + Bvec[i])=(reltest+X[i]) then count:=count+1;
        end;
        if count<n then
        begin
          p:=0.0; i:=0;
          repeat
            i:=i+1; res:=nlres(i,n,Bvec,notcomp, Workdata);
            if (not notcomp) then  p:=p+res*res;
          until notcomp or (i>=m);  {MODIFICATION m replaces n 1991-01-13}
          ifn:=ifn+1;
        end;
      end;
      if count<n then
        if (not singmat) and (not notcomp) and (p<Fmin) then
        begin
          lambda:=lambda*dec;
          Fmin:=p;
          calcmat:=true;
        end
      else
      begin
        lambda:=lambda*inc;
        if lambda<eps*eps then lambda:=eps;
        calcmat:=false;
      end;
    end;
  end;
end;
