{$I turbo.cnm}
program dr21(input,output);
{dr21.PAS == driver for Variable Metric method

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file ROSEN.PAS, which must be
  replaced with similar code for the user's problem.

          Copyright 1988 J.C.Nash
}

{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I rosen.pas}
{$I alg21.pas}
{$I startup.pas}

var
  n          : integer; {the order of the problem}
  Bvec       : rvector; {current set of parameters}
  X          : rvector; {"best" set of parameters}
  Workdata   : probdata; { the problem data type from CONSTYPE.DEF}
  i          : integer;
  Fmin       : real;   {for the minimal function value found}
  fail       : boolean; {set TRUE if the method fails in some way}
  mytol      : real; {to store a convergence tolerance}

begin
  banner:='dr21.pas -- driver for variable metric minimisation';
  startup;
  fminset(n,Bvec,Workdata); {sets up problem and defines starting
                  values of Bvec}
  mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
            must obtain an appropriate value.}
  vmmin(n,Bvec,X,Fmin,Workdata,fail,mytol); {minimise the function}
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
  end; {loop to write out parameters}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr21.pas -- variable metric minimisation}

