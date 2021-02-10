{$I turbo.cnm}
program dr27(input,output);
{dr27.pas == driver for Hooke and Jeeves method

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file ROSEN.PAS, which must be
  replaced with similar code for the user's problem.

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I Calceps.pas}
{$I rosen.pas}
(* remove the comments and delete the inclusion of ROSEN.PAS
   to use the JJACF.PAS test with EX27R.CNM
   {$I JJACF.PAS}
   Note that we move the inclusion to the right just in case.
*)
{$I alg27.pas}
{$I startup.pas}

       {main program}
var
  n          : integer; {the order of the problem}
  B          : rvector; {current set of parameters}
  X          : rvector; {"best" set of parameters}
  Workdata   : probdata; { the problem data type from CONSTYPE.DEF}
  i          : integer;
  Fmin       : real;   {for the minimal function value found}
  fail       : boolean; {set TRUE if the method fails in some way}
  mytol      : real; {to store a convergence tolerance}

begin
  banner:='dr27.pas -- driver for Hooke & Jeeves minimisation';
  startup;
  fminset(n,B,Workdata); {sets up problem and defines starting
                  values of B}
  mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
            must obtain an appropriate value.}
  hjmin(n,B,X,Fmin,Workdata,fail,mytol); {minimise the function}
  writeln;
  writeln(confile);
  writeln(' Minimum function value found =',Fmin);
  writeln(confile,' Minimum function value found =',Fmin);
  writeln(' At parameters');
  writeln(confile,' At parameters');
  for i:=1 to n do
  begin
    writeln(' B[',i,']=',X[i]);
    writeln(confile,' B[',i,']=',X[i]);
  end; {loop to write out parameters}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr27.pas -- Hooke & Jeeves driver}
