procedure cgmin(n: integer;
          var Bvec, X: rvector;
          var Fmin: real;
            Workdata: probdata;
          var fail: boolean;
          var intol: real);

type
  methodtype= (Fletcher_Reeves, Polak_Ribiere, Beale_Sorenson);

const
  Maxparm = 25;
  stepredn = 0.2;
  acctol = 0.0001;
  reltest = 10.0;

var
  accpoint  : boolean;
  c         : rvector;
  count     : integer;
  cycle     : integer;
  cyclimit  : integer;
  f         : real;
  funcount  : integer;
  g         : rvector;
  G1, G2     : real;
  G3, gradproj     : real;
  gradcount : integer;
  i, j       : integer;
  method    : methodtype;
  newstep   : real;
  notcomp   : boolean;
  oldstep   : real;
  s         : real;
  setstep   : real;
  steplength: real;
  t         : rvector;
  tol       : real;

begin
  writeln('alg22.pas -- Nash Algorithm 22 version 2 1988-03-24');
  writeln('  Conjugate gradients function minimiser');
  writeln('Steplength saving factor multiplies best steplength found at the');
  writeln('  end of each iteration as a starting value for next search');
  write('Enter a steplength saving factor (sugg. 1.7) -- setstep ');
  readln(setstep);
  writeln(setstep);
  write('Choose method (1=FR, 2=PR, 3=BS) ');
  readln(i); writeln(i);
  case i of
    1: method:=Fletcher_Reeves;
    2: method:=Polak_Ribiere;
    3: method:=Beale_Sorenson;
    else halt;
  end;
  case method of
    Fletcher_Reeves: writeln('Method: Fletcher Reeves');
    Polak_Ribiere: writeln('Method: Polak Ribiere');
    Beale_Sorenson: writeln('Method: Beale Sorenson');
  end;
  fail:=false;
  cyclimit:=n;
  if intol<0.0 then intol:=Calceps;
  tol:=intol*n*sqrt(intol);

  writeln('tolerance used in gradient test=', tol);
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
    gradcount:=0;
    repeat
      for i:=1 to n do
      begin
        t[i]:=0.0;
        c[i]:=0.0;
      end;
      cycle:=0;
      oldstep:=1.0;
      count:=0;
      repeat
        cycle:=cycle+1;
        writeln(gradcount, ' ', funcount, ' ', Fmin);
        write('parameters ');
        for i:=1 to n do
        begin
          write(Bvec[i]:10:5, ' ');
          if (7 * (i div 7) = i) and (i<n) then writeln;
        end;
        writeln;
        gradcount:=gradcount+1;
        fmingr(n, Bvec, Workdata, g);
        G1:=0.0; G2:=0.0;
        for i:=1 to n do
        begin
          X[i]:=Bvec[i];
          case method of
            Fletcher_Reeves: begin
              G1:=G1+sqr(g[i]); G2:=G2+sqr(c[i]);
            end;
            Polak_Ribiere  : begin
              G1:=G1+g[i]*(g[i]-c[i]); G2:=G2+sqr(c[i]);
            end;
            Beale_Sorenson : begin
              G1:=G1+g[i]*(g[i]-c[i]); G2:=G2+t[i]*(g[i]-c[i]);
            end;
          end;
          c[i]:=g[i];
        end;
        if G1>tol then
        begin
          if G2>0.0 then G3:=G1/G2 else G3:=1.0;
          gradproj:=0.0;
          for i:=1 to n do
          begin
            t[i]:=t[i]*G3-g[i]; gradproj:=gradproj+t[i]*g[i];
          end;
          steplength:=oldstep;

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
                steplength:=steplength*stepredn;
                write('*');
              end;
            end;
          until (count=n) or accpoint;
          if count<n then
          begin
            newstep:=2*((f-Fmin)-gradproj*steplength);
            if newstep>0 then
            begin
              newstep:=-gradproj*sqr(steplength)/newstep;
              for i:=1 to n do
              begin
                Bvec[i]:=X[i]+newstep*t[i];
              end;
              Fmin:=f;
              f:=fminfn(n, Bvec, Workdata, notcomp);
              funcount:=funcount+1;
              if f<Fmin then
              begin
                Fmin:=f; write(' i< ');
              end
              else
              begin
                write(' i> ');
                for i:=1 to n do Bvec[i]:=X[i]+steplength*t[i];
              end;
            end;
          end;
        end;
        oldstep:=setstep*steplength;
        if oldstep>1.0 then oldstep:=1.0;
      until (count=n) or (G1<=tol) or (cycle=cyclimit);

    until (cycle=1) and ((count=n) or (G1<=tol));

  end;
  writeln('Exiting from Alg22.pas conjugate gradients minimiser');
  writeln('    ', funcount, ' function evaluations used');
  writeln('    ', gradcount, ' gradient evaluations used');
end;
