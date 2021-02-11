{$I turbo.cnm}
program dr1920(input,output);
{dr1920.PAS == driver for Nelder-Mead method with axial search

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file ROSEN.PAS, which must be
  replaced with similar code for the user's problem.

  Modified here to allow for axial search (alg20.pas).

          Copyright 1988 J.C.Nash
}

{$I constype.def}
{$I tdstamp.pas} {time and date stamp}
{$I calceps.pas}
{$I rosen.pas}
(* remove the comments and delete the inclusion of ROSEN.PAS
   to use the JJACF.PAS test with EX1920J.CNM
   {$I jjacf.pas}
*)
{$I alg19.pas}
{$I alg20.pas}
{$I startup.pas}

var
  n          : integer; {the order of the problem}
  B          : rvector; {current set of parameters}
  X          : rvector; {"best" set of parameters}
  Workdata   : probdata; { the problem data type from CONSTYPE.DEF}
  i          : integer;
  Fmin       : real;   {for the minimal function value found}
  fail       : boolean; {set TRUE if the method fails in some way}
  mytol      : real; {to store a convergence tolerance}
  lowerfn    : boolean; {set TRUE if a lower function value is found
                  during the axial search}

begin
  banner:='dr1920.pas -- driver for Nelder-Mead minimisation';
  startup;
  fminset(n,B,Workdata); {sets up problem and defines starting
                  values of B}
  lowerfn:=false; {safety setting}
  repeat
    mytol:=-1.0; {Note: set the tolerance negative to indicate that
              procedure must obtain an appropriate value.}
    nmmin(n,B,X,Fmin,Workdata,fail,mytol); {minimise the function}
    writeln;
    writeln(confile);
    writeln(' Minimum function value found =',Fmin);
    writeln(' At parameters');
    writeln(confile,' Minimum function value found =',Fmin);
    writeln(confile,' At parameters');
    for i:=1 to n do
    begin
    writeln(' B[',i,']=',X[i]);
    writeln(confile,' B[',i,']=',X[i]);
    end; {loop to write out parameters}
    axissrch(n, B, Fmin, lowerfn, Workdata);  {alg20.pas}
    if lowerfn then
    begin
    writeln('Lower function value found');
    writeln(confile,'Lower function value found');
    end;
  until (not lowerfn);
  flush(confile); close(confile);
  if infname<>'con' then close(infile);
end. {dr1920.pas -- Nelder Mead minimisation with axial search}
