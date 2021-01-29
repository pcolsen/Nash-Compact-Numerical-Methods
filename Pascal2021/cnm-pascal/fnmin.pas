{
fnmin.pas
 -- a main program to run function minimisation procedures
}

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
  fminset(n,B,Workdata); {sets up problem and defines starting
                                values of B}
  mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
                must obtain an appropriate value.}
  minmeth(n,B,X,Fmin,Workdata,fail,mytol); {minimise the function}
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
end. {fnmin.pas}
