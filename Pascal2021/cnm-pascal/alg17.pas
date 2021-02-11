procedure min1d(var bb : real;
                var st: real;
                var ifn : integer;
                var fnminval : real );


var
  a1, a2, fii, s0, s1, s2, tt0, tt1, tt2, x0, x1, x2, xii : real;
  notcomp, tripleok: boolean;

begin
  writeln('alg17.pas -- One dimensional function minimisation');
  writeln(confile,'alg17.pas -- One dimensional function minimisation');

  ifn := 0;
  a1 := 1.5;
  a2 := -0.25;
  x1 := bb;
  notcomp := false;
  s0 := fn1d(x1,notcomp); ifn := ifn+1;
  if notcomp then
  begin
    writeln('*** FAILURE *** Function cannot be computed at initial point');
    writeln(confile,
         '*** FAILURE *** Function cannot be computed at initial point');
    halt;
  end;
  repeat
    x0 := x1;
    bb := x0;
    x1 := x0+st;
    s1 := fn1d(x1,notcomp); if notcomp then s1 := big; ifn := ifn+1;

    tripleok := false;
    if s1<s0 then
    begin
      repeat
        st := st*a1;
        x2 := x1+st;
        s2 := fn1d(x2,notcomp); if notcomp then s2 := big; ifn := ifn+1;
        if s2<s1 then
        begin
          s0 := s1; s1 := s2;
          x0 := x1; x1 := x2;
          write('Success1 '); write(confile,'Success1 ');
        end
        else
        begin
          tripleok := true;
          write('Failure1');  write(confile,'Failure1');
        end;
      until tripleok;
    end
    else
    begin
      st := a2*st;
      tt2 := s0; s0 := s1; s1 := tt2;
      tt2 := x0; x0 := x1; x1 := tt2;
      repeat
        x2 := x1+st;
        s2 := fn1d(x2,notcomp); if notcomp then s2 := big; ifn := ifn+1;
        if s2<s1 then
        begin
          s0 := s1; s1 := s2; x0 := x1; x1 := x2;
          st := st*a1;
          write('Success2 '); write(confile,'Success2 ');
        end
        else
        begin
          tripleok := true; write('Failure2'); write(confile,'Failure2');
        end;
      until tripleok;
    end;

    writeln; writeln(confile);  writeln('Triple (',x0,',',s0,')');
    writeln('       (',x1,',',s1,')');  writeln('       (',x2,',',s2,')');
    writeln(confile,'Triple (',x0,',',s0,')');
    writeln(confile,'       (',x1,',',s1,')');
    writeln(confile,'       (',x2,',',s2,')');
    tt0 := x0-x1;
    tt1 := (s0-s1)*st; tt2 := (s2-s1)*tt0;
    if tt1<>tt2 then
    begin
      st := 0.5*(tt2*tt0-tt1*st)/(tt2-tt1);
      xii := x1+st;
      writeln('Paramin step and argument :',st,' ',xii);
      writeln(confile,'Paramin step and argument :',st,' ',xii);
      if (reltest+xii)<>(reltest+x1) then
      begin
        fii := fn1d(xii,notcomp); ifn := ifn+1;
        if notcomp then fii := big;
        if fii<s1 then
        begin
          s1 := fii; x1 := xii;
          writeln('New min f(',x1,')=',s1);
          writeln(confile,'New min f(',x1,')=',s1);
        end;
      end;
    end;
    writeln(ifn,' evalns    f(',x1,')=',s1);
    writeln(confile,ifn,' evalns    f(',x1,')=',s1);
    s0 := s1;
  until (bb=x1);
  writeln('Apparent minimum is f(',bb,')=',s1);
  writeln('     after ',ifn,' function evaluations');
  writeln(confile,'Apparent minimum is f(',bb,')=',s1);
  writeln(confile,'     after ',ifn,' function evaluations');
  fnminval := s1;
end;
