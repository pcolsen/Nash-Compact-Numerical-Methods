program dr27(input,output);
{dr27.pas == driver for Hooke and Jeeves method

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file ROSEN.PAS, which must be
  replaced with similar code for the user's problem.

          Copyright 1988 J.C.Nash
}
{constype.def ==
  This file contains various definitions and type statements which are
  used throughout the collection of "Compact Numerical Methods".  In many
  cases not all definitions are needed, and users with very tight memory
  constraints may wish to remove some of the lines of this file when
  compiling certain programs.

  Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}

const
  big = 1.0E+35;    {a very large number}
  Maxconst = 25;    {Maximum number of constants in data record}
  Maxobs = 100;     {Maximum number of observations in data record}
  Maxparm  = 25;    {Maximum number of parameters to adjust}
  Maxvars = 10;     {Maximum number of variables in data record}
  acctol = 0.0001;  {acceptable point tolerance for minimisation codes}
  maxm = 20;        {Maximum number or rows in a matrix}
  maxn = 20;        {Maximum number of columns in a matrix}
  maxmn = 40;       {maxn+maxm, the number of rows in a working array}
  maxsym = 210;     {maximum number of elements of a symmetric matrix
              which need to be stored = maxm * (maxm + 1)/2 }
  reltest = 10.0;   {a relative size used to check equality of numbers.
              Numbers x and y are considered equal if the
              floating-point representation of reltest+x equals
              that of reltest+y.}
  stepredn = 0.2;   {factor to reduce stepsize in line search}
  yearwrit = 1990;  {year in which file was written}

type
  str2 = string[2];
  rmatrix = array[1..maxm, 1..maxn] of real; {a real matrix}
  wmatrix = array[1..maxmn, 1..maxn] of real; {a working array, formed
                  as one real matrix stacked on another}
  smatvec = array[1..maxsym] of real; {a vector to store a symmetric matrix
              as the row-wise expansion of its lower triangle}
  rvector = array[1..maxm] of real;  {a real vector. We will use vectors
              of m elements always. While this is NOT space efficient,
              it simplifies program codes.}
  cgmethodtype= (Fletcher_Reeves,Polak_Ribiere,Beale_Sorenson);
    {three possible forms of the conjugate gradients updating formulae}
  probdata = record
          m     : integer; {number of observations}
          nvar  : integer; {number of variables}
          nconst: integer; {number of constants}
          vconst: array[1..Maxconst] of real;
          Ydata : array[1..Maxobs, 1..Maxvars] of real;
          nlls  : boolean; {true if problem is nonlinear least squares}
        end;
{
  NOTE: Pascal does not let us define the work-space for the function
  within the user-defined code.  This is a weakness of Pascal for this
  type of work.

}
var {global definitions}
  banner     : string[80]; {program name and description}

function calceps:real;
{calceps.pas ==
  This function returns the machine EPSILON or floating point tolerance,
  the smallest positive real number such that 1.0 + EPSILON > 1.0.
  EPSILON is needed to set various tolerances for different algorithms.
  While it could be entered as a constant, I prefer to calculate it, since
  users tend to move software between machines without paying attention to
  the computing environment. Note that more complete routines exist.
}
var
  e,e0: real;
  i: integer;
begin {calculate machine epsilon}
  e0 := 1; i:=0;
  repeat
    e0 := e0/2; e := 1+e0;  i := i+1;
  until (e=1.0) or (i=50); {note safety check}
  e0 := e0*2;
{ Writeln('Machine EPSILON =',e0);}
  calceps:=e0;
end; {calceps}

(* remove the comments and delete the inclusion of ROSEN.PAS
   to use the JJACF.PAS test with EX27R.CNM
   {$I JJACF.PAS}
   Note that we move the inclusion to the right just in case.
*)

{rosen.pas
  == suite of procedures and functions defining the Rosenbrock
    banana shaped valley problem.
}
procedure fminset(var n:integer;var Bvec: rvector; var Workdata: probdata);
{sets up problem and defines starting values of Bvec}
{setup for Rosenbrock problem from rosen.pas}
begin
  writeln('Function: Rosenbrock Banana Valley');
  n:=2;
  Workdata.m:=2; {for nonlinear least squares problems}
  Workdata.nvar:=0;
  Bvec[1]:=-1.2;
  Bvec[2]:=1.0;
  writeln('Classical starting point (-1.2,1)');
end; {fminset from rosen.pas}

function fminfn(n: integer; var Bvec: rvector; var Workdata:probdata;
            var nocomp:boolean):real;
{this is the Rosenbrock banana valley function from rosen.pas}
begin
  nocomp:=false; {never undefined here}
  fminfn:=sqr(Bvec[2]-sqr(Bvec[1]))*100.0+sqr(1.0-Bvec[1]);
end; {fminfn from rosen.pas}
procedure fmingr(n:integer;Bvec:rvector; var Workdata:probdata; 
                                                   var g:rvector);
{computes the gradient of the Rosenbrock banana valley at point Bvec
  from rosen.pas}
begin
  g[1]:=-400.0*Bvec[1]*(Bvec[2]-sqr(Bvec[1]))-2.0*(1.0-Bvec[1]);
  g[2]:=200.0*(Bvec[2]-sqr(Bvec[1]));
end; {fmingrad from rosen.pas}

function nlres(i, n : integer; Bvec: rvector; var nocomp: boolean;
                                           var Workdata: probdata): real;
{computes residuals for the nonlinear least squares form of the
  Rosenbrock function from rosen.pas}
var
  temp: real;
begin
  nocomp:=false; {never set here}
  case i of
    1: begin
      temp:=10.0*(Bvec[2]-sqr(Bvec[1]));
    end;
    2: begin
      temp:=1.0-Bvec[1];
    end;
    else halt; {safety stop}
  end; {case}
  nlres := temp; {assign residual}
end; {nlres from rosen.pas}
procedure nljac(i, n: integer; Bvec: rvector; var jacrow: rvector;
                                              var Workdata: probdata);
{computes derivatives of residuals for the nonlinear least squares
  form of the Rosenbrock function from rosen.pas}
begin
  case i of
    1: begin
      jacrow[1]:=-20.0*Bvec[1];
      jacrow[2]:=10.0;
    end;
    2: begin
      jacrow[1]:=-1.0;
      jacrow[2]:=0.0;
    end;
    else halt; {safety stop}
  end; {case}
end; {nljac from rosen.pas}

{end of rosen.pas test function code suite}
procedure hjmin(n: integer;
        var B,X: rvector;
        var Fmin: real;
            Workdata: probdata;
        var fail: boolean;
            intol: real);

var
  i: integer;
  stepsize: real;
  fold: real;
  fval: real;
  notcomp: boolean;
  temp: real;
  samepoint: boolean;
  ifn: integer;

begin
  if intol<0.0 then intol := calceps;
  ifn := 1;
  fail := false;

  stepsize := 0.0;
  for i := 1 to n do
    if stepsize < stepredn*abs(B[i]) then stepsize := stepredn*abs(B[i]);
  if stepsize=0.0 then stepsize := stepredn;

  for i := 1 to n do X[i] := B[i];

  fval := fminfn(n, B,Workdata,notcomp);
  if notcomp then
  begin
    writeln('*** FAILURE *** Function not computable at initial point');
    fail := true;
  end
  else
  begin
    writeln('Initial function value =',fval);
    for i := 1 to n do
    begin
      write(B[i]:10:5,' ');
      if (7 * (i div 7) = i) and (i<n) then writeln;
    end;
    writeln;
    fold := fval; Fmin := fval;
    while stepsize>intol do
    begin

      for i := 1 to n do
      begin
        temp := B[i]; B[i] := temp+stepsize;
        fval := fminfn(n, B,Workdata,notcomp); ifn := ifn+1;
        if notcomp then fval := big;
        if fval<Fmin then
          Fmin := fval
        else
        begin
          B[i] := temp-stepsize;
          fval := fminfn(n, B,Workdata,notcomp); ifn := ifn+1;
          if notcomp then fval := big;
          if fval<Fmin then
            Fmin := fval
          else
            B[i] := temp;
        end;
      end;
      if Fmin<fold then
      begin

        for i := 1 to n do
        begin
          temp := 2.0*B[i]-X[i];
          X[i] := B[i]; B[i] := temp;
        end;
        fold := Fmin;
      end
      else
      begin
        samepoint := true;
        i := 1;
        repeat
          if B[i]<>X[i] then samepoint := false;
          i := i+1;
        until (not samepoint) or (i>n);
        if samepoint then
        begin
          stepsize := stepsize*stepredn;

          write('stepsize now ',stepsize:10,'  Best fn value=',Fmin);
          writeln(' after ',ifn);
          for i := 1 to n do
          begin
            write(B[i]:10:5,' ');
            if (7 * (i div 7) = i) and (i<n) then  writeln;
          end;
          writeln;
        end
        else
        begin
          for i := 1 to n do B[i] := X[i];
          writeln('Return to old base point');
        end;
      end;
    end;
    writeln('Converged to Fmin=',Fmin,' after ',ifn,' evaluations');
  end;
end;

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
  fminset(n,B,Workdata); {sets up problem and defines starting
                  values of B}
  mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
            must obtain an appropriate value.}
  hjmin(n,B,X,Fmin,Workdata,fail,mytol); {minimise the function}
  writeln;
  writeln(' Minimum function value found =',Fmin);
  writeln(' At parameters');
  for i:=1 to n do
  begin
    writeln(' B[',i,']=',X[i]);
  end; {loop to write out parameters}
end. {dr27.pas -- Hooke & Jeeves driver}
