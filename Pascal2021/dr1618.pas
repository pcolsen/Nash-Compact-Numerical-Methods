program dr1618(input, output);
{dr1618.PAS == this program is designed to allow one-dimensional root
finding using
alg16.pas -- grid search -- gridsrch
alg18.pas -- one-dimensional root-finding by a bisection and
regula falsi method -- root1d
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
function fn1d(x: real;var nocomp: boolean): real;
  {real valued test function of x for [1D] minimisation}
{cubfn.pas}
begin
    fn1d:=x*(x*x-2.0)-5.0; {assign function value for return}
end; {fn1d.pas == cubfn.pas}
procedure root1d(var lbound, ubound: real;
                 var ifn: integer;
                     tol : real;
                 var noroot: boolean );

var
 nbis: integer;
 b, fb, flow, fup : real;
 notcomp: boolean;

begin
  writeln('alg18.pas -- root of a function of one variable');
  notcomp := false;
  ifn := 2;
  nbis := 5;
  fup := fn1d(ubound,notcomp);
  if notcomp then halt;
  flow := fn1d(lbound,notcomp);
  if notcomp then halt;
  writeln('f(',lbound:8:5,')=',flow,'  f(',ubound:8:5,')=',fup);
  if fup*flow>0 then noroot := true else noroot := false;
  while (not noroot) and ((ubound-lbound)>tol) do
  begin

    if (nbis * ((ifn - 2) div nbis) = (ifn - 2)) then
    begin
      write('Bisect  '); 
      b := lbound + 0.5*(ubound - lbound)
    end
    else
    begin
      write('False P ');
      b := (lbound*fup-ubound*flow)/(fup-flow);
    end;

    if b<=lbound then
    begin
      b := lbound;
      ubound := lbound;
    end;
    if b>=ubound then
    begin
      b := ubound; lbound := ubound;
    end;
    ifn := ifn+1;
    fb := fn1d(b, notcomp);
    if notcomp then halt;
    write(ifn,' evalns: f(',b:16,')=',fb:10);
    writeln('  width interval= ',(ubound-lbound):10);
    if (ubound-lbound)>tol then
    begin
      if fb*flow<0.0 then
      begin
        fup := fb; ubound := b;
      end
      else
      begin
        flow := fb; lbound := b;
      end;
    end;
  end;
  writeln('Converged to f(',b,')=',fb);
  writeln('  Final interval width =',ubound-lbound);
end;
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
{Main program}
var
  lbound, tfmin, ubound, widthtol: real;
  changarg, ifn, minarg, nint: integer;
  noroot: boolean;
begin
  banner:='dr1618.pas -- One dimensional root-finding';

readln(lbound);
readln(ubound);
writeln('Enter lower bound for search ',lbound);
writeln('Enter upper bound for search ',ubound);
readln(nint);
writeln('Enter the number of intervals for grid search (0 for none) ',nint);
readln(widthtol);
writeln('Enter a tolerance for root search interval width ',widthtol);
if nint>0 then 
begin 
   gridsrch(lbound, ubound, nint, tfmin, minarg, changarg);
   ubound := (ubound-lbound)/nint; {to temporarily save stepsize}
   lbound := lbound+(changarg-1)*ubound; {new lower bound}
   ubound := lbound+ubound; {new upper bound}
end; {adjustment of bounds from grid search}
writeln('Now try rootfinder');
root1d(lbound, ubound, ifn, widthtol, noroot);
if noroot then writeln('Possibly no root in interval');
writeln;
noroot := true; {set TRUE to display root of function}
lbound := 0.0;
tfmin := fn1d(lbound, noroot);
end. {dr1618.pas}
