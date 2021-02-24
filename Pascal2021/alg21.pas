procedure vmmin(n: integer;
            var Bvec, X: rvector;
            var Fmin: real;
                Workdata: probdata;
            var fail: boolean;
            var intol: real);

const
  Maxparm = 25;
  stepredn = 0.2;
  acctol = 0.0001;
  reltest = 10.0;

var
  accpoint  : boolean;
  B         : array[1..Maxparm, 1..Maxparm] of real;

  c         : rvector;
  count     : integer;
  D1, D2     : real;
  f         : real;
  funcount  : integer;
  g         : rvector;
  gradcount : integer;
  gradproj  : real;
  i, j       : integer;
  ilast     : integer;
  notcomp   : boolean;
  s         : real;
  steplength: real;
  t         : rvector;

begin
  writeln('alg21.pas -- version 2 1988-03-24');
  writeln('  Variable metric function minimiser');
  fail:=false;
  f:=fminfn(n, Bvec, Workdata, notcomp);
  if notcomp then
  begin
    writeln('**** Function cannot be evaluated at initial parameters ****');
    fail := true;
  end
  else
  begin
    Fmin:=f;
    funcount:=1;
    gradcount:=1;
    fmingr(n, Bvec, Workdata, g);
    ilast:=gradcount;

    repeat
      if ilast=gradcount then
      begin
        for i:=1 to n do
        begin
          for j:=1 to n do B[i, j]:=0.0; B[i, i]:=1.0;
        end;
      end;
      writeln(gradcount,' ', funcount,' ', Fmin);
      write('parameters ');
      for i:=1 to n do write(Bvec[i]:10:5,' ');
      writeln;
      for i:=1 to n do
      begin
        X[i]:=Bvec[i];
        c[i]:=g[i];
      end;

      gradproj:=0.0;
      for i:=1 to n do
      begin
        s:=0.0;
        for j:=1 to n do s:=s-B[i, j]*g[j];
        t[i]:=s; gradproj:=gradproj+s*g[i];
      end;

      if gradproj<0.0 then {!! note change to floating point}
      begin
        steplength:=1.0;

        accpoint:=false;
        repeat
          count:=0;
          for i:=1 to n do
          begin
            Bvec[i]:=X[i]+steplength*t[i];
            if (reltest+X[i])=(reltest+Bvec[i]) then count:=count+1;
          end;
          if count<n then
          begin
            f:=fminfn(n, Bvec, Workdata, notcomp);
            funcount:=funcount+1;
            accpoint:=(not notcomp) and (f<=Fmin+gradproj*steplength*acctol);

            if not accpoint then
            begin
              steplength:=steplength*stepredn; write('*'); 
            end;
          end;
        until (count=n) or accpoint;
        if count<n then
        begin
          Fmin:=f;
          fmingr(n, Bvec, Workdata, g);
          gradcount:=gradcount+1;
          D1:=0.0;
          for i:=1 to n do
          begin
            t[i]:=steplength*t[i]; c[i]:=g[i]-c[i];
            D1:=D1+t[i]*c[i];
          end;
          if D1>0 then
          begin
            D2:=0.0;
            for i:=1 to n do
            begin
              s:=0.0;
              for j:=1 to n do s:=s+B[i, j]*c[j];
              X[i]:=s; D2:=D2+s*c[i];
            end;
            D2:=1.0+D2/D1;
            for i:=1 to n do
            begin
              for j:=1 to n do
              begin
                B[i, j]:=B[i, j]-(t[i]*X[j]+X[i]*t[j]-D2*t[i]*t[j])/D1;
              end;
            end;
          end
          else
          begin
            writeln(' UPDATE NOT POSSIBLE');
            ilast:=gradcount;
          end;
        end
        else
        begin
          if ilast<gradcount then
          begin
            count:=0;
            ilast:=gradcount;
          end;
        end;
      end
      else
      begin
          writeln('UPHILL SEARCH DIRECTION');
          count:=0; {!! order of statements}
          if ilast=gradcount then count:=n else ilast:=gradcount;
          {!! Resets Hessian inverse if it has not just been set,
              otherwise forces a convergence.}
      end;
    until (count=n) and (ilast=gradcount);
  end;
