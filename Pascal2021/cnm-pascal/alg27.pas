procedure hjmin(n: integer;
        var B,X: rvector;
        var Fmin: real;
            Workdata: probdata;
        var fail: boolean;
            intol: real);

var
  i, j: integer;
  stepsize: real;
  fold: real;
  fval: real;
  notcomp: boolean;
  temp: real;
  samepoint: boolean;
  ifn: integer;

begin
  if intol<0.0 then intol := calceps;
  ifn := 1;
  fail := false;

  stepsize := 0.0;
  for i := 1 to n do
    if stepsize < stepredn*abs(B[i]) then stepsize := stepredn*abs(B[i]);
  if stepsize=0.0 then stepsize := stepredn;

  for i := 1 to n do X[i] := B[i];

  fval := fminfn(n, B,Workdata,notcomp);
  if notcomp then
  begin
    writeln('*** FAILURE *** Function not computable at initial point');
    writeln(confile,
        '*** FAILURE *** Function not computable at initial point');
    fail := true;
  end
  else
  begin
    writeln('Initial function value =',fval);
    writeln(confile,'Initial function value =',fval);
    for i := 1 to n do
    begin
      write(B[i]:10:5,' ');
      write(confile,B[i]:10:5,' ');
      if (7 * (i div 7) = i) and (i<n) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
    fold := fval; Fmin := fval;
    while stepsize>intol do
    begin

      for i := 1 to n do
      begin
        temp := B[i]; B[i] := temp+stepsize;
        fval := fminfn(n, B,Workdata,notcomp); ifn := ifn+1;
        if notcomp then fval := big;
        if fval<Fmin then
          Fmin := fval
        else
        begin
          B[i] := temp-stepsize;
          fval := fminfn(n, B,Workdata,notcomp); ifn := ifn+1;
          if notcomp then fval := big;
          if fval<Fmin then
            Fmin := fval
          else
            B[i] := temp;
        end;
      end;
      if Fmin<fold then
      begin

        for i := 1 to n do
        begin
          temp := 2.0*B[i]-X[i];
          X[i] := B[i]; B[i] := temp;
        end;
        fold := Fmin;
      end
      else
      begin
        samepoint := true;
        i := 1;
        repeat
          if B[i]<>X[i] then samepoint := false;
          i := i+1;
        until (not samepoint) or (i>n);
        if samepoint then
        begin
          stepsize := stepsize*stepredn;

          write('stepsize now ',stepsize:10,'  Best fn value=',Fmin);
          write(confile,'stepsize now ',stepsize:10,'  Best fn value=',Fmin);
          writeln(' after ',ifn);
          writeln(confile,' after ',ifn);
          for i := 1 to n do
          begin
            write(B[i]:10:5,' ');
            write(confile,B[i]:10:5,' ');
            if (7 * (i div 7) = i) and (i<n) then
            begin
              writeln;
              writeln(confile);
            end;
          end;
          writeln;
          writeln(confile);
        end
        else
        begin
          for i := 1 to n do B[i] := X[i];

          writeln('Return to old base point');
          writeln(confile,'Return to old base point');
        end;
      end;
    end;
    writeln('Converged to Fmin=',Fmin,' after ',ifn,' evaluations');
    writeln(confile,'Converged to Fmin=',Fmin,' after ',ifn,' evaluations');

  end;
end;
