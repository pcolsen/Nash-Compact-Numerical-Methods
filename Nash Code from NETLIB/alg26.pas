procedure comeig( n : integer;
          var itcount: integer;
          var A, Z, T, U : rmatrix);

var
  Rvec : rvector;

  i, itlimit, j, k, k1, m, m9, n1 : integer;
  aki, ami, bv, br, bi : real;
  c, c1i, c1r, c2i, c2r, ca, cb, ch, cos2a, cot2x, cotx, cx : real;
  d, de, di, diag, dr, e, ei, er, eps, eta, g, hi, hj, hr : real;
  isw, max, nc, nd, root1, root2, root : real;
  s, s1i, s1r, s2i, s2r, sa, sb, sh, sig, sin2a, sx : real;
  tanh, tau, te, tee, tem, tep ,tse, zki, zmi : real;
  mark : boolean;

begin

  writeln('alg26.pas -- comeig');
  writeln(confile,'alg26.pas -- comeig');
  eps := Calceps;

  mark := false; n1 := n-1;
  for i := 1 to n do
  begin
    for j := 1 to n do
    begin
      T[i,j] := 0.0; U[i,j] := 0.0; if i=j then T[i,i] := 1.0;
    end;
  end;
  itlimit := itcount;
  itcount := 0;
  while (itcount<=itlimit) and (not mark) do
  begin
    itcount := itcount+1;
    tau := 0.0;
    diag := 0.0;
    for k := 1 to n do
    begin
      tem := 0;
      for i := 1 to n do if i<>k then tem := tem+ABS(A[i,k])+ABS(Z[i,k]);
      tau := tau+tem; tep := abs(A[k,k])+abs(Z[k,k]);
      diag := diag+tep;
      Rvec[k] := tem+tep;
    end;
    writeln('TAU=',tau,'  AT ITN ',itcount);
    writeln(confile,'TAU=',tau,'  AT ITN ',itcount);
    for k := 1 to n1 do
    begin
      max := Rvec[k]; i := k; k1 := k+1;
      for j := k1 to n do
      begin
        if max<Rvec[j] then
        begin
          max := Rvec[j]; i := j;
        end;
      end;
      if i<>k then
      begin
        Rvec[i] := Rvec[k];
        for j := 1 to n do
        begin
          tep := A[k,j]; A[k,j] := A[i,j]; A[i,j] := tep; tep := Z[k,j];
          Z[k,j] := Z[i,j]; Z[i,j] := tep;
        end;
        for j := 1 to n do
        begin
          tep := A[j,k]; A[j,k] := A[j,i]; A[j,i] := tep; tep := Z[j,k];
          Z[j,k] := Z[j,i]; Z[j,i] := tep; tep := T[j,k]; T[j,k] := T[j,i];
          T[j,i] := tep; tep := U[j,k]; U[j,k] := U[j,i]; U[j,i] := tep;
        end;
      end;
    end;
    if tau>=100.0*eps then
    begin
      mark := true;
      for k := 1 to n1 do
      begin
        k1 := k+1;
        for m := k1 to n do
        begin
          hj := 0.0; hr := 0.0; hi := 0.0; g := 0.0;
          for i := 1 to n do
          begin
            if (i<>k) and (i<>m) then
            begin
              hr := hr+A[k,i]*A[m,i]+Z[k,i]*Z[m,i];
              hr := hr-A[i,k]*A[i,m]-Z[i,k]*Z[i,m];
              hi := hi+Z[k,i]*A[m,i]-A[k,i]*Z[m,i];
              hi := hi-A[i,k]*Z[i,m]+Z[i,k]*A[i,m];
              te := A[i,k]*A[i,k]+Z[i,k]*Z[i,k]+A[m,i]*A[m,i]+Z[m,i]*Z[m,i];
              tee := A[i,m]*A[i,m]+Z[i,m]*Z[i,m]+A[k,i]*A[k,i]+Z[k,i]*Z[k,i];
              g := g+te+tee; hj := hj-te+tee;
            end;
          end;
          br := A[k,m]+A[m,k]; bi := Z[k,m]+Z[m,k]; er := A[k,m]-A[m,k];
          ei := Z[k,m]-Z[m,k]; dr := A[k,k]-A[m,m]; di := Z[k,k]-Z[m,m];
          te := br*br+ei*ei+dr*dr; tee := bi*bi+er*er+di*di;
          if te>=tee then
          begin
            isw := 1.0; c := br; s := ei; d := dr; de := di;
            root2 := sqrt(te);
          end
          else
          begin
            isw := -1.0; c := bi; s := -er; d := di; de := dr;
            root2 := sqrt(tee);
          end;
          root1 := sqrt(s*s+c*c); sig := -1.0; if d>=0.0 then sig := 1.0;
          sa := 0.0; ca := -1.0; if c>=0.0 then ca := 1.0;
          if root1<=eps then
          begin
            sx := 0.0; sa := 0.0; cx := 1.0; ca := 1.0;
            if isw<=0.0 then
            begin
              e := ei; bv := -br;
            end
            else
            begin
              e := er; bv := bi;
            end;
            nd := d*d+de*de;
          end
          else
          begin
            if abs(s)>eps then
            begin
              ca := c/root1; sa := s/root1;
            end;
            cot2x := d/root1; cotx := cot2x+(sig*sqrt(1.0+cot2x*cot2x));
            sx := sig/sqrt(1.0+cotx*cotx); cx := sx*cotx;

            eta := (er*br+ei*bi)/root1; tse := (br*bi-er*ei)/root1;
            te := sig*(tse*d-de*root1)/root2; tee := (d*de+root1*tse)/root2;
            nd := root2*root2+tee*tee; tee := hj*cx*sx; cos2a := ca*ca-sa*sa;
            sin2a := 2.0*ca*sa; tem := hr*cos2a+hi*sin2a;
            tep := hi*cos2a-hr*sin2a; hr := hr*cx*cx-tem*sx*sx-ca*tee;
            hi := hi*cx*cx+tep*sx*sx-sa*tee;
            bv := isw*te*ca+eta*sa; e := ca*eta-isw*te*sa;
          end;

          s := hr-sig*root2*e; c := hi-sig*root2*bv; root := sqrt(c*c+s*s);
          if root<eps then
          begin
            cb := 1.0; ch := 1.0; sb := 0.0; sh := 0.0;
          end
          else
          begin
            cb := -c/root; sb := s/root; tee := cb*bv-e*sb; nc := tee*tee;
            tanh := root/(g+2.0*(nc+nd)); ch := 1.0/sqrt(1.0-tanh*tanh);
            sh := ch*tanh;
          end;
          tem := sx*sh*(sa*cb-sb*ca); c1r := cx*ch-tem; c2r := cx*ch+tem;
          c1i := -sx*sh*(ca*cb+sa*sb); c2i := c1i; tep := sx*ch*ca;
          tem := cx*sh*sb; s1r := tep-tem; s2r := -tep-tem; tep := sx*ch*sa;
          tem := cx*sh*cb; s1i := tep+tem; s2i := tep-tem;
          tem := sqrt(s1r*s1r+s1i*s1i); tep := sqrt(s2r*s2r+s2i*s2i);
          if tep>eps then mark := false;
          if (tep>eps) and (tem>eps) then
          begin
            for i := 1 to n do
            begin
              aki := A[k,i]; ami := A[m,i]; zki := Z[k,i]; zmi := Z[m,i];
              A[k,i] := c1r*aki-c1i*zki+s1r*ami-s1i*zmi;
              Z[k,i] := c1r*zki+c1i*aki+s1r*zmi+s1i*ami;
              A[m,i] := s2r*aki-s2i*zki+c2r*ami-c2i*zmi;
              Z[m,i] := s2r*zki+s2i*aki+c2r*zmi+c2i*ami;
            end;
            for i := 1 to n do
            begin
              aki := A[i,k]; ami := A[i,m]; zki := Z[i,k]; zmi := Z[i,m];
              A[i,k] := c2r*aki-c2i*zki-s2r*ami+s2i*zmi;
              Z[i,k] := c2r*zki+c2i*aki-s2r*zmi-s2i*ami;
              A[i,m] := -s1r*aki+s1i*zki+c1r*ami-c1i*zmi;
              Z[i,m] := -s1r*zki-s1i*aki+c1r*zmi+c1i*ami;
              aki := T[i,k]; ami := T[i,m]; zki := U[i,k]; zmi := U[i,m];
              T[i,k] := c2r*aki-c2i*zki-s2r*ami+s2i*zmi;
              U[i,k] := c2r*zki+c2i*aki-s2r*zmi-s2i*ami;
              T[i,m] := -s1r*aki+s1i*zki+c1r*ami-c1i*zmi;
              U[i,m] := -s1r*zki-s1i*aki+c1r*zmi+c1i*ami;
            end;
          end;
        end;
      end;
    end
    else mark := true;
  end;
  if itcount>itlimit then itcount := -itcount;
end;
