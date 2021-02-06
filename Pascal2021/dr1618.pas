{$M 20000,0,512000}
{TURBO.CNM for Turbo Pascal 5.0}
program dr1618(input, output);
{dr1618.PAS == this program is designed to allow one-dimensional root
finding using
alg16.pas -- grid search -- gridsrch
alg18.pas -- one-dimensional root-finding by a bisection and
regula falsi method -- root1d
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
uses Dos, Crt; {Turbo Pascal 5.0 Modules}
{ 1. Interrupt, Unit, Interface, Implementation, Uses are reserved words now.}
{ 2. System,Dos,Crt are standard unit names in Turbo 5.0.}

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

  The following variables allow us to keep a copy of all screen
  information in a file for some of the codes.  Pascal requires a
  variable (confile in this case) for the file itself.  The string
  variable confname is used for the name of the file.  Similar variables
  allow problem data to be read from the file dfile named dfname.
}
var {global definitions}
  banner     : string[80]; {program name and description}
  confile    : text;       {file for output of console image}
  confname   : string[64]; {a name for confile}
  dfile      : text;       {file for output of console image}
  dfname     : string[64]; {a name for confile}
  infile     : text;       {an input file (keyboard image)}
  infname    : string[64]; {a name for the input file}
  con        : text;       {the console file}
{********TDSTAMP.PAS***********************************}
procedure leadzero(value:word; var outstr: str2);
{
          Copyright 1990 J.C.Nash
}
begin   {code needed by tdstamp.pas}
  if value<10 then
  begin
    str(value:1,outstr);
    outstr:='0'+outstr;
  end else str(value:2,outstr);
end; {leadzero}
{******************************************************}
procedure tdstamp(var outfile : text);
{tdstamp.pas == writes a time/date stamp to outfile

          Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}
 var
  tdstr: string[20];
  year: word;
  month: word;
  day: word;
  dayofweek: word;
  yrstr: string[4];
  mnstr: str2;
  daystr: str2;
  hrstr : str2;
  minstr: str2;
  secstr: str2;
  hour: word;
  minute: word;
  second: word;
  sec100: word;

begin
  getdate(year, month, day, dayofweek);
  gettime(hour, minute, second, sec100);
  writeln(year,month,day,hour,minute,second);{??}
  if year<yearwrit then
  begin
    writeln('Program senses that clock is not up to date.');
    readln;
    halt;
  end;
  str(year:4,yrstr);
  leadzero(month,mnstr);
  leadzero(day,daystr);
  leadzero(hour,hrstr);
  leadzero(minute,minstr);
  leadzero(second,secstr);
  tdstr:=yrstr+'/'+mnstr+'/'+daystr+' '+hrstr+':'+minstr+':'+secstr;
  writeln(outfile,tdstr);
end; {tdstamp.pas}
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
{htanfn.pas
  == This is the hyperbolic tangent Example 13.2 from
    Compact Numerical Methods.

If this function is called with nocomp TRUE, the root will
be displayed. Otherwise, the function at x will be computed.
}
var
  val,y: real;
  s,t,w,z: real; {constants which define the function}

begin
  s:=1; t:=0.5; w:=0.2; z:=100.0; {fix a particular function}
  {This function is well-behaved.}
  {The following altenative set of parameters gives a function which is
  nearly a step-function.}
{  s:=100; t:=0.5; w:=0.99; z:=100.0;}
  if nocomp then
  begin
    writeln('root is at ',t + ln((1-w)/(1+w))/(2*s)  );
    writeln(confile,'root is at ',t + ln((1-w)/(1+w))/(2*s)  );
  end
  else
  begin
    y:=s*(x - t);
    val:= exp(y);{ partial result for building tanh}
    val:=(val - 1.0/val)/(val + 1.0/val); {tanh}
    val:=z*(val+w);
    nocomp:=false; {never not computable}
    fn1d:=val; {assign function value for return}
  end;
end; {fn1d.pas == htanfn.pas}
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
  writeln(confile,'alg16.pas -- one-dimensional grid search');

  writeln('In gridsrch lbound=',lbound,'  ubound=',ubound);
  writeln(confile,'In gridsrch lbound=',lbound,'  ubound=',ubound);
  notcomp:=false;
  t:=fn1d(lbound, notcomp);
  writeln('  lb  f(',lbound,')=',t);
  writeln(confile,'  lb  f(',lbound,')=',t);
  if notcomp then halt;
  fmin:=t;
  minarg:=0;
  changarg:=0;
  h:=(ubound-lbound)/nint;
  for j:=1 to nint do

  begin
    p:=fn1d(lbound+j*h, notcomp);
    write('      f(',lbound+j*h,')=',p);
    write(confile,'      f(',lbound+j*h,')=',p);
    if notcomp then halt;
    if p<fmin then
    begin
      fmin:=p; minarg:=j;
    end;
    if p*t<=0 then
    begin
      writeln(' *** sign change ***');
      writeln(confile,' *** sign change ***');
      changarg:=j;
    end
    else
    begin
      writeln; writeln(confile);
    end;
    t:=p;
  end;
  writeln('Minimum so far is f(',lbound+minarg*h,')=',fmin);
  writeln(confile,'Minimum so far is f(',lbound+minarg*h,')=',fmin);
  if changarg>0 then
  begin
    writeln('Sign change observed last in interval ');
    writeln(' [',lbound+(changarg-1)*h,',',lbound+changarg*h,']');
    writeln(confile,'Sign change observed last in interval ');
    writeln(confile,' [',lbound+(changarg-1)*h,',',lbound+changarg*h,']');
  end
  else
  begin
    writeln('Apparently no sign change in [',lbound,',',ubound,']');
    writeln(confile,'Apparently no sign change in [',lbound,',',ubound,']');
  end;
end;
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
  writeln(confile,'alg18.pas -- root of a function of one variable');

  notcomp := false;
  ifn := 2;
  nbis := 5;
  fup := fn1d(ubound,notcomp);
  if notcomp then halt;
  flow := fn1d(lbound,notcomp);
  if notcomp then halt;
  writeln('f(',lbound:8:5,')=',flow,'  f(',ubound:8:5,')=',fup);
  writeln(confile,'f(',lbound:8:5,')=',flow,'  f(',ubound:8:5,')=',fup);
  if fup*flow>0 then noroot := true else noroot := false;
  while (not noroot) and ((ubound-lbound)>tol) do
  begin

    if (nbis * ((ifn - 2) div nbis) = (ifn - 2)) then
    begin
      write('Bisect  '); write(confile,'Bisect  ');
      b := lbound + 0.5*(ubound - lbound)
    end
    else
    begin
      write('False P '); write(confile,'False P ');
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
    write(confile,ifn,' evalns: f(',b:16,')=',fb:10);
    writeln('  width interval= ',(ubound-lbound):10);
    writeln(confile,'  width interval= ',(ubound-lbound):10);
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
  writeln(confile,'Converged to f(',b,')=',fb);
  writeln(confile,'  Final interval width =',ubound-lbound);
end;
procedure startup;
{startup -- startup code modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}
begin
  assign(con,'CON');
  rewrite(con);
  assign(input,''); reset(input);
  assign(output,''); rewrite(output); {needed for I/O redirection}
  writeln(banner);
  {display the program banner}
  tdstamp(CON); {displays a time and date stamp}
  write('File for input of control data ([cr] for keyboard) ');
  readln(infname); if length(infname)=0 then infname:='con';
  {Gets a filename for a file with data for running a problem.}
  assign(infile, infname); reset(infile);
  {Opens this file}
  write('File for console image ([cr] = nul) ');
  readln(infile,confname);
  {Gets a name for a file to which console output is repeated.}
  if length(confname)=0 then confname:='nul';
  {Changes a carriage return to a 'nul' so console copy is omitted.}
  if infname<>'con' then writeln(confname);
  {Writes out the console file name on the screen if it comes from a file.
  Also moves to new line when reading from a file rather than keyboard.}
  assign(confile,confname); rewrite(confile);
  {Opens to console image file.}
  writeln(confile,banner);
  tdstamp(confile); {Repeat banner and time stamp -- it is a little bit
    different from the console timestamp.}
  writeln(confile,'File for input of control data ([cr] for keyboard) ',
                infname);
  writeln(confile,'File for console image ([cr] = nul) ',confname);
  {Repeat the setup information to the console image file.}
  {Now start the program proper.}
end; {startup.pas}
{Main program}
var
gfmin, lbound, tfmin, ubound, widthtol: real;
changarg, ifn, minarg, nint: integer;
noroot: boolean;
begin
banner:='dr1618.pas -- One dimensional root-finding';
startup;
write('Enter lower bound for search ');readln(infile,lbound);
write('Enter upper bound for search ');readln(infile,ubound);
write('Enter a tolerance for root search interval width ');
readln(infile,widthtol); if infname<>'con' then writeln(widthtol);
writeln(confile,'Enter lower bound for search ',lbound);
writeln(confile,'Enter upper bound for search ',ubound);
writeln(confile,'Enter a tolerance for root search interval width ',
widthtol);
write('Enter the number of intervals for grid search (0 for none) ');
readln(infile,nint);
writeln(confile,
'Enter the number of intervals for grid search (0 for none) ',nint);
if infname<>'con' then writeln(nint);
if nint>0 then gridsrch(lbound, ubound, nint, tfmin, minarg, changarg);
writeln;
writeln(confile);
writeln('Now try rootfinder');
writeln(confile,'Now try rootfinder');
ubound := (ubound-lbound)/nint; {to temporarily save stepsize}
lbound := lbound+(changarg-1)*ubound; {new lower bound}
ubound := lbound+ubound; {new upper bound}
root1d(lbound, ubound, ifn, widthtol, noroot);
if noroot then writeln('Possibly no root in interval');
writeln;
writeln(confile);
noroot := true; {set TRUE to display root of function}
lbound := 0.0;
tfmin := fn1d(lbound, noroot);
flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr1618.pas}
