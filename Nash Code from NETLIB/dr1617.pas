{$I turbo.cnm}
program dr1617(input, output);
{dr1617.PAS == this program is designed to allow one-dimensional minimum
          finding using

  ALG16.PAS -- grid search
  ALG17.PAS -- one-dimensional minimisation using success-failure
            search and parabolic inverse interpolation

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I calceps.pas}
{$I cubefn.pas}
{$I alg16.pas}
{$I alg17.pas}
{$I startup.pas}

{Main program}
var
  b, gfmin, lbound, st, tfmin, ubound, widthtol : real;
  changarg, ifn, minarg, nint: integer;

begin
  banner:='dr1617.pas -- One dimensional minimisation';
  startup;
  write('Enter lower bound for search ');readln(infile,lbound);
  if infname<>'con' then writeln(lbound);
  write('Enter upper bound for search ');readln(infile,ubound);
  if infname<>'con' then writeln(ubound);
  write('Enter a tolerance for search interval width ');
  readln(infile,widthtol); if infname<>'con' then writeln(widthtol);
  writeln(confile,'Enter lower bound for search ',lbound);
  writeln(confile,'Enter upper bound for search ',ubound);
  writeln(confile,'Enter a tolerance for search interval width ',
          widthtol);
  write('Enter the number of intervals per search (0 for no grid search) ');
  readln(infile,nint); if infname<>'con' then writeln(nint);
  writeln(confile,
    'Enter the number of intervals per search (0 for no grid search) ',nint);
  if nint>0 then
  begin
    gfmin := big;
    repeat
    gridsrch(lbound, ubound, nint, tfmin, minarg, changarg);
    if tfmin<gfmin then
    begin
      st := (ubound-lbound)/nint;
      ubound := lbound+(minarg+1)*st; lbound := ubound-2.0*st;
      {Note that we use one step either side of minimal value to
            define the new interval.}
      writeln('New lowest function value =',tfmin,' in [',lbound,',',
                  ubound,']');
      writeln(confile,'New lowest function value =',tfmin,' in [',lbound,',',
                  ubound,']');
      gfmin := tfmin;
    end
    else
    begin
      writeln('Unable to reduce function');
      writeln('lowest function value still in [',lbound,',',ubound,']');
      writeln(confile,'Unable to reduce function');
      writeln(confile,
          'lowest function value still in [',lbound,',',ubound,']');
    end;
    until ((ubound-lbound)<=widthtol) or (tfmin>=gfmin);{end of repeat loop}
  end; {grid search section -- nint>0 }
  writeln('Now call the minimiser');
  writeln(confile,'Now call the minimiser');
  {write('Enter a starting guess for function argument ');
  readln(b);
  write('Enter a starting stepsize '); readln(st);}
  {Above lines could be used for user entry of initial data.}
  {Following two lines give an alternative choice of starting
  value for function argument and step size:}
  b := (ubound+lbound)/2.0; {try middle of interval}
  st := (ubound-lbound)/10.0; {and 10% of interval as stepsize}
  min1d(b, st, ifn, gfmin);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr1617.pas}
