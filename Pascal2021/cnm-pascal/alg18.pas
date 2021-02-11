procedure root1d(var lbound, ubound: real;
                 var ifn: integer;
                     tol : real;
                 var noroot: boolean );

var
 nbis: integer;
 b, fb, flow, fup : real;
 notcomp: boolean;

begin
  writeln('alg18.pas -- root of a function of one variable');
  writeln(confile,'alg18.pas -- root of a function of one variable');

  notcomp := false;
  ifn := 2;
  nbis := 5;
  fup := fn1d(ubound,notcomp);
  if notcomp then halt;
  flow := fn1d(lbound,notcomp);
  if notcomp then halt;
  writeln('f(',lbound:8:5,')=',flow,'  f(',ubound:8:5,')=',fup);
  writeln(confile,'f(',lbound:8:5,')=',flow,'  f(',ubound:8:5,')=',fup);
  if fup*flow>0 then noroot := true else noroot := false;
  while (not noroot) and ((ubound-lbound)>tol) do
  begin

    if (nbis * ((ifn - 2) div nbis) = (ifn - 2)) then
    begin
      write('Bisect  '); write(confile,'Bisect  ');
      b := lbound + 0.5*(ubound - lbound)
    end
    else
    begin
      write('False P '); write(confile,'False P ');
      b := (lbound*fup-ubound*flow)/(fup-flow);
    end;

    if b<=lbound then
    begin
      b := lbound;
      ubound := lbound;
    end;
    if b>=ubound then
    begin
      b := ubound; lbound := ubound;
    end;
    ifn := ifn+1;
    fb := fn1d(b, notcomp);
    if notcomp then halt;
    write(ifn,' evalns: f(',b:16,')=',fb:10);
    write(confile,ifn,' evalns: f(',b:16,')=',fb:10);
    writeln('  width interval= ',(ubound-lbound):10);
    writeln(confile,'  width interval= ',(ubound-lbound):10);
    if (ubound-lbound)>tol then
    begin
      if fb*flow<0.0 then
      begin
        fup := fb; ubound := b;
      end
      else
      begin
        flow := fb; lbound := b;
      end;
    end;
  end;
  writeln('Converged to f(',b,')=',fb);
  writeln('  Final interval width =',ubound-lbound);
  writeln(confile,'Converged to f(',b,')=',fb);
  writeln(confile,'  Final interval width =',ubound-lbound);
end;
