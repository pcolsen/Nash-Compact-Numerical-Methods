procedure rqmcg( n : integer;
            A, B : rmatrix;
          var X : rvector;
          var ipr : integer;
          var rq : real);

var
  count, i, itn, itlimit : integer;
  avec, bvec, yvec, zvec, g, t : rvector;
  beta, d, eps, gg, pa, pn, step : real;
  ta, tabt, tat, tbt, tol, u, v, w, xat, xax, xbt, xbx : real;
  conv: boolean;

begin
  writeln('alg25.pas -- Rayleigh quotient minimisation');
  itlimit := ipr;
  conv := false;
  ipr := 0;
  eps := calceps;
  tol := n*n*eps*eps;

  pa := big;
  while (ipr<=itlimit) and (not conv) do
  begin
    matmul(n, A, X, avec);
    matmul(n, B, X, bvec);
    ipr := ipr+1;

    xax := 0.0; xbx := 0.0;
    for i := 1 to n do
    begin
      xax := xax+X[i]*avec[i]; xbx := xbx+X[i]*bvec[i];
    end;
    if xbx<=tol then halt;
    rq := xax/xbx;
    write(ipr,' products -- ev approx. =',rq:18);
    if rq<pa then
    begin
      pa := rq;
      gg := 0.0;
      for i := 1 to n do
      begin
        g[i] := 2.0*(avec[i]-rq*bvec[i])/xbx; gg := gg+g[i]*g[i];
      end;
      writeln(' squared gradient norm =',gg:8);
      if gg>tol then

      begin

        for i := 1 to n do t[i] := -g[i];
        itn := 0;
        repeat
          itn := itn+1;
          matmul(n, A, t, yvec);
          matmul(n, B, t, zvec); ipr := ipr+1;
          tat := 0.0; tbt := 0.0; xat := 0.0; xbt := 0.0;
          for i := 1 to n do
          begin
            xat := xat+X[i]*yvec[i]; tat := tat+t[i]*yvec[i];
            xbt := xbt+X[i]*zvec[i]; tbt := tbt+t[i]*zvec[i];
          end;

          u := tat*xbt-xat*tbt; v := tat*xbx-xax*tbt;
          w := xat*xbx-xax*xbt; d := v*v-4.0*u*w;
          if d<0.0 then halt;

          d := sqrt(d);
          if v>0.0 then step := -2.0*w/(v+d) else step := 0.5*(d-v)/u;

          count := 0;
          xax := 0.0; xbx := 0.0;
          for i := 1 to n do
          begin
            avec[i] := avec[i]+step*yvec[i];
            bvec[i] := bvec[i]+step*zvec[i];
            w := X[i]; X[i] := w+step*t[i];
            if (reltest+w)=(reltest+X[i]) then count := count+1;
            xax := xax+X[i]*avec[i]; xbx := xbx+X[i]*bvec[i];
          end;
          if xbx<=tol then halt
                  else pn := xax/xbx;
          if (count<n) and (pn<rq) then
          begin
            rq := pn; gg := 0.0;
            for i := 1 to n do
            begin
              g[i] := 2.0*(avec[i]-pn*bvec[i])/xbx; gg := gg+g[i]*g[i];
            end;
            if gg>tol then
            begin
              xbt := 0.0; for i := 1 to n do xbt := xbt+X[i]*zvec[i];

              tabt := 0.0; beta := 0.0;
              for i := 1 to n do
              begin
                w := yvec[i]-pn*zvec[i]; tabt := tabt+t[i]*w;
                beta := beta+g[i]*(w-g[i]*xbt);
              end;
              beta := beta/tabt;

              for i := 1 to n do t[i] := beta*t[i]-g[i];
            end;
          end

          else
          begin
            if itn=1 then conv := true;
            itn := n+1;
          end;
        until (itn>=n) or (count=n) or (gg<=tol) or conv;
      end
      else conv := true;
    end
    else
    begin
      conv := true;
    end;
    ta := 0.0;
    for i := 1 to n do ta := ta+sqr(X[i]); ta := 1.0/sqrt(ta);
    for i := 1 to n do X[i] := ta*X[i];
  end;
  if ipr>itlimit then ipr := -ipr;
  writeln;
end;
