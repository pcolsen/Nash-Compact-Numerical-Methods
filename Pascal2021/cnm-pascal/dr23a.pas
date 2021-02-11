{$I turbo.cnm}
program dr23(input,output);
{dr23a.PAS == driver for Nash Marquardt nonlinear least squares

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file logistic.pas, which must be
  replaced with similar code for the user's problem.


          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}
{$I logistic.pas}
{$I alg07.pas}
{$I alg08.pas}
{$I alg23.pas}
{$I startup.pas}

{main program}
var
  n          : integer; {the order of the problem}
  Bvec       : rvector; {current set of parameters}
  X          : rvector; {"best" set of parameters}
  Workdata   : probdata;
  i          : integer;
  Fmin       : real;
  fail       : boolean;
  mytol      : real;

begin
  banner:='dr23.pas -- Marquardt Nash nonlinear least squares';
  startup;
  fminset(n,Bvec,Workdata);
            {Sets up problem and defines starting values of Bvec}
  mytol:=-1.0; {Note: set the tolerance negative to indicate that
            procedure must obtain an appropriate value.}
  Workdata.nlls:=true;
  modmrt( n, Bvec, X, Fmin, Workdata);
  writeln;
  writeln(confile);
  writeln(' Minimum function value found =',Fmin);
  writeln(confile,' Minimum function value found =',Fmin);
  writeln(' At parameters');
  writeln(confile,' At parameters');
  for i:=1 to n do
  begin
    writeln(' Bvec[',i,']=',X[i]);
    writeln(confile,' Bvec[',i,']=',X[i]);
  end;
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr23a.pas == Nash Marquardt nonlinear least squares}
