procedure nmmin(n: integer;
            var Bvec,X: rvector;
            var Fmin: real;
                Workdata: probdata;
            var fail: boolean;
            var intol: real);


const
  Pcol = 27;
  Prow = 26;
  alpha = 1.0;
  beta = 0.5;
  gamma = 2.0;

var
  action : string[15];
  C         : integer;
  calcvert  : boolean;
  convtol   : real;
  f         : real;
  funcount  : integer;
  H         : integer;
  i,j       : integer;
  L         : integer;
  notcomp   : boolean;
  n1        : integer;
  oldsize   : real;
  P         : array[1..Prow,1..Pcol] of real;
  shrinkfail: boolean;
  size      : real;
  step      : real;
  temp      : real;
  trystep   : real;
  tstr : string[5];
  VH,VL,VN  : real;
  VR        : real;

begin
  writeln('Nash Algorithm 19 version 2 1988-03-17');
  writeln('  Nelder Mead polytope direct search function minimiser');
  writeln(confile,'Nash Algorithm 19 version 2 1988-03-17');
  writeln(confile,'  Nelder Mead polytope direct search function minimiser');
  fail := false;
  f := fminfn(n,Bvec,Workdata,notcomp);
  if notcomp then
  begin
    writeln('**** Function cannot be evaluated at initial parameters ****');
    writeln(confile,
            '**** Function cannot be evaluated at initial parameters ****');
    fail := true;
  end
  else
  begin
    writeln('Function value for initial parameters = ',f);
    writeln(confile,'Function value for initial parameters = ',f);
    if intol<0.0 then intol := Calceps;
    funcount := 1;
    convtol := intol*(abs(f)+intol);
    writeln('  Scaled convergence tolerance is ',convtol);
    writeln(confile,'  Scaled convergence tolerance is ',convtol);
    n1 := n+1;  C := n+2;  P[n1,1] := f;
    for i := 1 to n do P[i,1] := Bvec[i];

    L := 1;
    size := 0.0;

    step := 0.0;
    for i := 1 to n do if 0.1*abs(Bvec[i])>step then step := 0.1*abs(Bvec[i]);
    writeln('Stepsize computed as ',step);
    writeln(confile,'Stepsize computed as ',step);
    for j := 2 to n1 do
    begin
      action := 'BUILD          ';
      for i := 1 to n do P[i,j] := Bvec[i];



      trystep := step;
      while P[j-1,j]=Bvec[j-1] do
      begin
        P[j-1,j] := Bvec[j-1]+trystep; trystep := trystep*10;
      end;
      size := size+trystep;
    end;
    oldsize := size;
    calcvert := true;
    shrinkfail := false;
    repeat
      if calcvert then
      begin
        for j := 1 to n1 do
        begin
          if j<>L then
          begin
            for i := 1 to n do Bvec[i] := P[i,j];
            f := fminfn(n,Bvec,Workdata,notcomp);
            if notcomp then f := big; funcount := funcount+1;  P[n1,j] := f;
          end;
        end;
        calcvert := false;
      end;

      VL := P[n1,L];
      VH := VL;
      H := L;

      for j := 1 to n1 do
      begin
        if j<>L then
        begin
          f := P[n1,j];
          if f<VL then
          begin
            L := j; VL := f;
          end;
          if f>VH then
          begin
            H := j; VH := f;
          end;
        end;
      end;

      if VH>VL+convtol then
      begin
        str(funcount:5,tstr);
        writeln(action,tstr,' ',VH,' ',VL);
        writeln(confile,action,tstr,' ',VH,' ',VL);
        VN := beta*VL+(1.0-beta)*VH;

        for i := 1 to n do
        begin
          temp := -P[i,H];
          for j := 1 to n1 do temp := temp+P[i,j];
          P[i,C] := temp/n;
        end;
        for i := 1 to n do
            Bvec[i] := (1.0+alpha)*P[i,C]-alpha*P[i,H];
        f := fminfn(n,Bvec,Workdata,notcomp);
        if notcomp then f := big;
        funcount := funcount+1;
        action := 'REFLECTION     ';
        VR := f;
        if VR<VL then
        begin
          P[n1,C] := f;
          for i := 1 to n do
          begin
            f := gamma*Bvec[i]+(1-gamma)*P[i,C];
            P[i,C] := Bvec[i];
            Bvec[i] := f;
          end;
          f := fminfn(n,Bvec,Workdata,notcomp);
          if notcomp then f := big; funcount := funcount+1;
          if f<VR then
          begin
            for i := 1 to n do P[i,H] := Bvec[i];
            P[n1,H] := f;
            action := 'EXTENSION      ';
          end
          else
          begin
            for i := 1 to n do P[i,H] := P[i,C];
            P[n1,H] := VR;
          end
        end
        else
        begin
          action := 'HI-REDUCTION    ';
          if VR<VH then
          begin
            for i := 1 to n do P[i,H] := Bvec[i];
            P[n1,H] := VR;
            action := 'LO-REDUCTION    ';
          end;

          for i := 1 to n do Bvec[i] := (1-beta)*P[i,H]+beta*P[i,C];
          f := fminfn(n,Bvec,Workdata,notcomp);
          if notcomp then f := big; funcount := funcount+1;

          if f<P[n1,H] then
          begin
            for i := 1 to n do P[i,H] := Bvec[i];
            P[n1,H] := f;
          end
          else

          if VR>=VH then
          begin
            action := 'SHRINK         ';
            calcvert := true;
            size := 0.0;
            for j := 1 to n1 do
            begin
              if j<>L then
              for i := 1 to n do
              begin
                P[i,j] := beta*(P[i,j]-P[i,L])+P[i,L];
                size := size+abs(P[i,j]-P[i,L]);
              end;
            end;
            if size<oldsize then
            begin
              shrinkfail := false;
              oldsize := size;
            end
            else
            begin
              writeln('Polytope size measure not decreased in shrink');
              writeln(confile,
                      'Polytope size measure not decreased in shrink');
              shrinkfail := true;
            end;
          end;
        end;
      end;

    until ((VH<=VL+convtol) or shrinkfail );

  end;

  writeln('Exiting from Alg19.pas Nelder Mead polytope minimiser');
  writeln('    ',funcount,' function evaluations used');
  writeln(confile,'Exiting from Alg19.pas Nelder Mead polytope minimiser');
  writeln(confile,'    ',funcount,' function evaluations used');
  Fmin := P[n1,L];
  for i := 1 to n do X[i] := P[i,L];
  if shrinkfail then fail := true;

end;
