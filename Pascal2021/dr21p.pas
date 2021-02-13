program dr21(input,output);
{dr21.PAS == driver for Variable Metric method
This program is designed to minimise functions of n parameters.
Present example uses the problem file ROSEN.PAS, which must be
replaced with similar code for the user's problem.
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
procedure vmmin(n: integer;
            var Bvec, X: rvector;
            var Fmin: real;
                Workdata: probdata;
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
  notcomp   : boolean;
  s         : real;
  steplength: real;
  t         : rvector;

begin
  writeln('alg21.pas -- version 2 1988-03-24');
  writeln('  Variable metric function minimiser');
  fail:=false;
  f:=fminfn(n, Bvec, Workdata, notcomp);
  if notcomp then
  begin
    writeln('**** Function cannot be evaluated at initial parameters ****');
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
      write('parameters ');
      for i:=1 to n do write(Bvec[i]:10:5,' ');
      writeln;
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
        steplength:=1.0;

        accpoint:=false;
        repeat
          count:=0;
          for i:=1 to n do
          begin
            Bvec[i]:=X[i]+steplength*t[i];
            if (reltest+X[i])=(reltest+Bvec[i]) then count:=count+1;
          end;
          if count<n then
          begin
            f:=fminfn(n, Bvec, Workdata, notcomp);
            funcount:=funcount+1;
            accpoint:=(not notcomp) and (f<=Fmin+gradproj*steplength*acctol);

            if not accpoint then
            begin
              steplength:=steplength*stepredn; write('*'); 
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
            t[i]:=steplength*t[i]; c[i]:=g[i]-c[i];
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
end;

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
fminset(n,Bvec,Workdata); {sets up problem and defines starting
values of Bvec}
mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
must obtain an appropriate value.}
vmmin(n,Bvec,X,Fmin,Workdata,fail,mytol); {minimise the function}
writeln;
writeln(' Minimum function value found =',Fmin);
writeln(' At parameters');
for i:=1 to n do writeln(' Bvec[',i,']=',X[i]);
end. {dr21.pas -- variable metric minimisation}
