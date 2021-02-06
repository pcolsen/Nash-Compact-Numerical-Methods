program dr1617(input, output);
{dr1617.PAS == this program is designed to allow one-dimensional minimum
finding using
ALG16.PAS -- grid search
ALG17.PAS -- one-dimensional minimisation using success-failure
search and parabolic inverse interpolation
Copyright 1988 J.C.Nash
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
function fn1d(x: real;var nocomp: boolean): real;
  {real valued test function of x for [1D] minimisation}
{
This is a cubic function x*(x*x - 2.) - 5.
From Forsythe, Malcolm & Moler (1977) page 184.
Minimum is at 0.81650 .
}
begin
  nocomp:=false; {always set for this function}
  fn1d:=(x*x-2.0)*x-5.0;
end; {fn1d.pas == cubefn.pas}

procedure gridsrch( var lbound, ubound : real;
                    nint : integer;
                    var fmin: real;
                    var minarg: integer;
                    var changarg: integer  );

var
  j : integer;
  h, p, t : real;
  notcomp : boolean;

begin
  writeln('alg16.pas -- one-dimensional grid search');
  writeln('In gridsrch lbound=',lbound,'  ubound=',ubound);
  notcomp:=false;
  t:=fn1d(lbound, notcomp);
  writeln('  lb  f(',lbound,')=',t);
  if notcomp then halt;
  fmin:=t;
  minarg:=0;
  changarg:=0;
  h:=(ubound-lbound)/nint;
  for j:=1 to nint do

  begin
    p:=fn1d(lbound+j*h, notcomp);
    write('      f(',lbound+j*h,')=',p);
    if notcomp then halt;
    if p<fmin then
    begin
      fmin:=p; minarg:=j;
    end;
    if p*t<=0 then
    begin
      writeln(' *** sign change ***');
      changarg:=j;
    end
    else
    begin
      writeln; 
    end;
    t:=p;
  end;
  writeln('Minimum so far is f(',lbound+minarg*h,')=',fmin);
  if changarg>0 then
  begin
    writeln('Sign change observed last in interval ');
    writeln(' [',lbound+(changarg-1)*h,',',lbound+changarg*h,']');
  end
  else
  begin
    writeln('Apparently no sign change in [',lbound,',',ubound,']');
  end;
end;
procedure min1d(var bb : real;
                var st: real;
                var ifn : integer;
                var fnminval : real );


var
  a1, a2, fii, s0, s1, s2, tt0, tt1, tt2, x0, x1, x2, xii : real;
  notcomp, tripleok: boolean;

begin
  writeln('alg17.pas -- One dimensional function minimisation');

  ifn := 0;
  a1 := 1.5;
  a2 := -0.25;
  x1 := bb;
  notcomp := false;
  s0 := fn1d(x1,notcomp); ifn := ifn+1;
  if notcomp then
  begin
    writeln('*** FAILURE *** Function cannot be computed at initial point');
    halt;
  end;
  repeat
    x0 := x1;
    bb := x0;
    x1 := x0+st;
    s1 := fn1d(x1,notcomp); if notcomp then s1 := big; ifn := ifn+1;

    tripleok := false;
    if s1<s0 then
    begin
      repeat
        st := st*a1;
        x2 := x1+st;
        s2 := fn1d(x2,notcomp); if notcomp then s2 := big; ifn := ifn+1;
        if s2<s1 then
        begin
          s0 := s1; s1 := s2;
          x0 := x1; x1 := x2;
          write('Success1 ');
        end
        else
        begin
          tripleok := true;
          write('Failure1');
        end;
      until tripleok;
    end
    else
    begin
      st := a2*st;
      tt2 := s0; s0 := s1; s1 := tt2;
      tt2 := x0; x0 := x1; x1 := tt2;
      repeat
        x2 := x1+st;
        s2 := fn1d(x2,notcomp); if notcomp then s2 := big; ifn := ifn+1;
        if s2<s1 then
        begin
          s0 := s1; s1 := s2; x0 := x1; x1 := x2;
          st := st*a1;
          write('Success2 ');
        end
        else
        begin
          tripleok := true; write('Failure2');
        end;
      until tripleok;
    end;

    writeln; writeln('Triple (',x0,',',s0,')');
    writeln('       (',x1,',',s1,')');  writeln('       (',x2,',',s2,')');
    tt0 := x0-x1;
    tt1 := (s0-s1)*st; tt2 := (s2-s1)*tt0;
    if tt1<>tt2 then
    begin
      st := 0.5*(tt2*tt0-tt1*st)/(tt2-tt1);
      xii := x1+st;
      writeln('Paramin step and argument :',st,' ',xii);
      if (reltest+xii)<>(reltest+x1) then
      begin
        fii := fn1d(xii,notcomp); ifn := ifn+1;
        if notcomp then fii := big;
        if fii<s1 then
        begin
          s1 := fii; x1 := xii;
          writeln('New min f(',x1,')=',s1);
        end;
      end;
    end;
    writeln(ifn,' evalns    f(',x1,')=',s1);
    s0 := s1;
  until (bb=x1);
  writeln('Apparent minimum is f(',bb,')=',s1);
  writeln('     after ',ifn,' function evaluations');
  fnminval := s1;
end;
{Main program}
var
b, gfmin, lbound, st, tfmin, ubound, widthtol : real;
changarg, ifn, minarg, nint: integer;
begin
banner:='dr1617.pas -- One dimensional minimisation';
write('Enter lower bound for search ');readln(lbound);
writeln(lbound);
write('Enter upper bound for search ');readln(ubound);
writeln(ubound);
if (lbound >= ubound) then halt;
write('Enter a tolerance for search interval width ');
readln(widthtol); writeln(widthtol);
write('Enter the number of intervals per search (0 for no grid search) ');
readln(nint); writeln(nint);
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
gfmin := tfmin;
end
else
begin
writeln('Unable to reduce function');
writeln('lowest function value still in [',lbound,',',ubound,']');
end;
until ((ubound-lbound)<=widthtol) or (tfmin>=gfmin);{end of repeat loop}
end; {grid search section -- nint>0 }
writeln('Now call the minimiser');
{write('Enter a starting guess for function argument ');
readln(b);
write('Enter a starting stepsize '); readln(st);}
{Above lines could be used for user entry of initial data.}
{Following two lines give an alternative choice of starting
value for function argument and step size:}
b := (ubound+lbound)/2.0; {try middle of interval}
st := (ubound-lbound)/10.0; {and 10% of interval as stepsize}
min1d(b, st, ifn, gfmin);
end. {dr1617.pas}
