{$I turbo.cnm}
program dr1618(input, output);
{dr1618.PAS == this program is designed to allow one-dimensional root
          finding using

  alg16.pas -- grid search -- gridsrch
  alg18.pas -- one-dimensional root-finding by a bisection and
          regula falsi method -- root1d

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I calceps.pas}
{$I htanfn.pas}  {called fn1d}
{$I alg16.pas}
{$I alg18.pas}
{$I startup.pas}

{Main program}
var
  gfmin, lbound, tfmin, ubound, widthtol: real;
  changarg, ifn, minarg, nint: integer;
  noroot: boolean;

begin
  banner:='dr1618.pas -- One dimensional root-finding';
  startup;
  write('Enter lower bound for search ');readln(infile,lbound);
  write('Enter upper bound for search ');readln(infile,ubound);
  write('Enter a tolerance for root search interval width ');
  readln(infile,widthtol); if infname<>'con' then writeln(widthtol);
  writeln(confile,'Enter lower bound for search ',lbound);
  writeln(confile,'Enter upper bound for search ',ubound);
  writeln(confile,'Enter a tolerance for root search interval width ',
            widthtol);
  write('Enter the number of intervals for grid search (0 for none) ');
  readln(infile,nint);
  writeln(confile,
    'Enter the number of intervals for grid search (0 for none) ',nint);
  if infname<>'con' then writeln(nint);
  if nint>0 then gridsrch(lbound, ubound, nint, tfmin, minarg, changarg);
  writeln;
  writeln(confile);
  writeln('Now try rootfinder');
  writeln(confile,'Now try rootfinder');
  ubound := (ubound-lbound)/nint; {to temporarily save stepsize}
  lbound := lbound+(changarg-1)*ubound; {new lower bound}
  ubound := lbound+ubound; {new upper bound}
  root1d(lbound, ubound, ifn, widthtol, noroot);
  if noroot then writeln('Possibly no root in interval');
  writeln;
  writeln(confile);
  noroot := true; {set TRUE to display root of function}
  lbound := 0.0;
  tfmin := fn1d(lbound, noroot);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr1618.pas}
