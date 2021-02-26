{$M 20000,0,512000}
{TURBO.CNM for Turbo Pascal 5.0}

program b21(input,output);
{b21.PAS == driver for Variable Metric method

  This program is designed to minimise functions of n parameters.

  Present example uses the problem file ROSEN.PAS, which must be
  replaced with similar code for the user's problem.

          Copyright 1988 J.C.Nash

MODIFICATIONS: 900624 for bounds constraints on parameters
}
{ bounds constraint mods marked with comment} {BC!!}


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
  bcind = (onupper, onlower, masked, free); {BC!!} {indicator type
                       to tell program the bounds constraint status}
  probdtax = record           {BC!! -- mod from probdtax}
          m     : integer; {number of observations}
          nvar  : integer; {number of variables}
          nconst: integer; {number of constants}
          vconst: array[1..Maxconst] of real;
          Ydata : array[1..Maxobs, 1..Maxvars] of real;
          nlls  : boolean; {true if problem is nonlinear least squares}
          bconst: array[1..Maxparm] of bcind; {BC!!} {bounds constraint
                  indicators}
          lower: array[1..Maxparm] of real; {BC!!} {lower bounds}
          upper: array[1..Maxparm] of real; {BC!!} {upper bounds}
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
 {time and date stamp}
{rosen.pas
  == suite of procedures and functions defining the Rosenbrock
    banana shaped valley problem.
}
procedure fminset(var n:integer;var Bvec: rvector; var Workdata: probdtax);
{sets up problem and defines starting values of Bvec}
{setup for Rosenbrock problem from rosen.pas}
begin
  writeln('Function: Rosenbrock Banana Valley');
  writeln(confile,'Function: Rosenbrock Banana Valley');
  n:=2;
  Workdata.m:=2; {for nonlinear least squares problems}
  Workdata.nvar:=0;
  Bvec[1]:=-1.2;
  Bvec[2]:=1.0;
  Workdata.lower[1]:=-1.5; {BC!!}
  Workdata.lower[2]:=-1.0; {BC!!}
  Workdata.upper[1]:=6.0;  {BC!!}
  Workdata.upper[2]:=2.0;  {BC!!}
  Workdata.bconst[1]:=free;  {BC!!}
  Workdata.bconst[2]:=free;  {BC!!}
  writeln('Classical starting point (-1.2,1)');
  writeln(confile,'Classical starting point (-1.2,1)');
end; {fminset from rosen.pas}

function fminfn(n: integer; var Bvec: rvector; Workdata:probdtax;
            var nocomp:boolean):real;
{this is the Rosenbrock banana valley function from rosen.pas}
begin
  nocomp:=false; {never undefined here}
  fminfn:=sqr(Bvec[2]-sqr(Bvec[1]))*100.0+sqr(1.0-Bvec[1]);
end; {fminfn from rosen.pas}
procedure fmingr(n:integer;Bvec:rvector; Workdata:probdtax; var g:rvector);
{computes the gradient of the Rosenbrock banana valley at point Bvec
  from rosen.pas}
begin
  g[1]:=-400.0*Bvec[1]*(Bvec[2]-sqr(Bvec[1]))-2.0*(1.0-Bvec[1]);
  g[2]:=200.0*(Bvec[2]-sqr(Bvec[1]));
end; {fmingrad from rosen.pas}

function nlres(i, n : integer; Bvec: rvector; var nocomp: boolean): real;
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
procedure nljac(i, n: integer; Bvec: rvector; var jacrow: rvector);
{computes derivatives of residuals for the nonlinear least squares
  form of the Rosenbrock function from rosen.pas}
var
  t1, t2: real;
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

procedure vmmin(n: integer;
            var Bvec, X: rvector;
            var Fmin: real;
                Workdata: probdtax;
            var fail: boolean;
            var intol: real);

const
  Maxparm = 25;
  stepredn = 0.2;
  acctol = 0.0001;
  reltest = 10.0;

var
  accpoint  : boolean;
  B         : array[1..Maxparm, 1..Maxparm] of real;

  c         : rvector;
  count     : integer;
  D1, D2     : real;
  f         : real;
  funcount  : integer;
  g         : rvector;
  gradcount : integer;
  gradproj  : real;
  i, j       : integer;
  ilast     : integer;
  nbound    : integer; {number of bounds constraits active} {BC!!}
  notcomp   : boolean;
  s         : real;
  steplength: real;
  t         : rvector;

begin
  writeln('alg21.pas -- version 2 1988-03-24');
  writeln('  Variable metric function minimiser');
  writeln(confile,'alg21.pas --  version 2 1988-03-24');
  writeln(confile,'  Variable metric function minimiser');
  fail:=false;
  f:=fminfn(n, Bvec, Workdata, notcomp);
  if notcomp then
  begin
    writeln('**** Function cannot be evaluated at initial parameters ****');
    writeln(confile,
            '**** Function cannot be evaluated at initial parameters ****');
    fail := true;
  end
  else
  begin
    Fmin:=f;
    funcount:=1;
    gradcount:=1;
    fmingr(n, Bvec, Workdata, g);
    ilast:=gradcount;
    repeat
      if ilast=gradcount then
      begin
        for i:=1 to n do
        begin
          for j:=1 to n do B[i, j]:=0.0; B[i, i]:=1.0;
        end;
      end;
      writeln(gradcount,' ', funcount,' ', Fmin);
      writeln(confile, gradcount,' ', funcount,' ', Fmin);
      write('parameters ');
      for i:=1 to n do write(Bvec[i]:10:5,' ');
      writeln; write(confile,'parameters ');
      for i:=1 to n do write(confile, Bvec[i]:10:5,' ');
      writeln(confile);
      for i:=1 to n do
      begin
        X[i]:=Bvec[i];
        c[i]:=g[i];
      end;
      gradproj:=0.0;
      for i:=1 to n do
      begin
        s:=0.0;
        for j:=1 to n do s:=s-B[i, j]*g[j];
        t[i]:=s; gradproj:=gradproj+s*g[i];
      end;
      if gradproj<0.0 then {!! note change to floating point}
      begin
        {BC!! start}
        steplength:=1.0; {start with a maximal steplength}
        nbound:=0; {to count the number of constraints active}
        for i:=1 to n do {loop to find correct steplength maximum}
        begin
          if Workdata.bconst[i]=free then
          begin
            if t[i]<>0.0 then {no limit on steplength for zero step element}
            begin
              if t[i]>0.0 then
                s:=(Workdata.upper[i]-X[i])/t[i]
                else s:=(Workdata.lower[i]-X[i])/t[i];
                {Note: since t[i]<0, this s is positive}
(*              write('parm ',i,' t = ',t[i],'  step =',s); readln; {??}  *)
              if s<steplength then steplength:=s; {reduce steplength}
            end; {if t[i]<>0.0}
          end {free parameter}
          else nbound:=nbound+1;
        end; {for loop to check steplength bounds}
        writeln('Number of bounds constraints active =',nbound);
        writeln(confile,'Number of bounds constraints active =',nbound);
        s:=steplength;
        if s<1.0 then
        begin
          writeln('Step reduced to ',s,' by bounds');
          writeln(confile,'Step reduced to ',s,' by bounds');
        end;
        {BC!! end}
        accpoint:=false;
        repeat
          count:=0;
          for i:=1 to n do
          begin
            if Workdata.bconst[i]=free then {BC!!}
            begin  {BC!!}
              Bvec[i]:=X[i]+steplength*t[i];
              if (reltest+X[i])=(reltest+Bvec[i]) then count:=count+1;
            end {free parameter} {BC!!}
            else Bvec[i]:=X[i]; {BC!! do we need this??}
          end;
          if count<n-nbound then {BC!!}
          begin
            f:=fminfn(n, Bvec, Workdata, notcomp);
            funcount:=funcount+1;
            accpoint:=(not notcomp) and (f<=Fmin+gradproj*steplength*acctol);

            if not accpoint then
            begin
              steplength:=steplength*stepredn; write('*'); write(confile,'*');
            end;
          end;
        until (count=n) or accpoint;
        if count<n then
        begin
          Fmin:=f;
          fmingr(n, Bvec, Workdata, g);
          gradcount:=gradcount+1;
          D1:=0.0;
          for i:=1 to n do
          begin
            if Workdata.bconst[i]=free then t[i]:=steplength*t[i]
               else t[i]:=0.0;
            c[i]:=g[i]-c[i];
            D1:=D1+t[i]*c[i];
          end;
          if D1>0 then
          begin
            D2:=0.0;
            for i:=1 to n do
            begin
              s:=0.0;
              for j:=1 to n do s:=s+B[i, j]*c[j];
              X[i]:=s; D2:=D2+s*c[i];
            end;
            D2:=1.0+D2/D1;
            for i:=1 to n do
            begin
              for j:=1 to n do
              begin
                B[i, j]:=B[i, j]-(t[i]*X[j]+X[i]*t[j]-D2*t[i]*t[j])/D1;
              end;
            end;
          end
          else
          begin
            writeln(' UPDATE NOT POSSIBLE');
            writeln(confile,' UPDATE NOT POSSIBLE');
            ilast:=gradcount;
          end;
        end
        else
        begin
          if ilast<gradcount then
          begin
            count:=0;
            ilast:=gradcount;
          end;
        end;
      end
      else
      begin
          writeln('UPHILL SEARCH DIRECTION');
          writeln(confile,'UPHILL SEARCH DIRECTION');
          count:=0; {!! order of statements}
          if ilast=gradcount then count:=n else ilast:=gradcount;
          {!! Resets Hessian inverse if it has not just been set,
              otherwise forces a convergence.}
      end;
    until (count=n) and (ilast=gradcount);
  end;

  writeln('Exiting from alg21.pas variable metric minimiser');
  writeln('    ', funcount,' function evaluations used');
  writeln('    ', gradcount,' gradient evaluations used');
  writeln(confile,'Exiting from alg21.pas variable metric minimiser');
  writeln(confile,'    ', funcount,' function evaluations used');
  writeln(confile,'    ', gradcount,' gradient evaluations used');
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


var
  n          : integer; {the order of the problem}
  Bvec       : rvector; {current set of parameters}

  X          : rvector; {"best" set of parameters}
  Workdata   : probdtax; { the problem data type from CONSTYPE.DEF}
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

