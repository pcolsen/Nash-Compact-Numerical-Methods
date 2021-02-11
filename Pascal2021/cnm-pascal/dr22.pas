{$I turbo.cnm}
program dr22(input,output);
{dr22.PAS == driver for conjugate gradients method

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file ROSEN.PAS, which must be
  replaced with similar code for the user's problem.

          Copyright 1988 J.C.Nash
}

{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I Calceps.pas}
{$I rosen.pas}
{$I alg22.pas}
{$I startup.pas}
{
fnmin.pas
  -- a main program to run function minimisation procedures
}

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
  banner:='dr21.pas -- driver for conjugate gradients minimisation';
  startup;
  fminset(n,Bvec,Workdata); {sets up problem and defines starting
                  values of Bvec}
  mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
            must obtain an appropriate value.}
  cgmin(n,Bvec,X,Fmin,Workdata,fail,mytol); {minimise the function}
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
end. {dr22.pas -- conjugate gradients minimisation}

