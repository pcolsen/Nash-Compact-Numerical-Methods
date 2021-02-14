program dr22(input,output);
{dr22.PAS == driver for conjugate gradients method
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
  cgmethodtype= (Fletcher_Reeves,Polak_Ribiere,Beale_Sorenson);
    {three possible forms of the conjugate gradients updating formulae}
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

procedure cgmin(n: integer;
          var Bvec, X: rvector;
          var Fmin: real;
            Workdata: probdata;
          var fail: boolean;
          var intol: real);

type
  methodtype= (Fletcher_Reeves, Polak_Ribiere, Beale_Sorenson);

const
  Maxparm = 25;
  stepredn = 0.2;
  acctol = 0.0001;
  reltest = 10.0;

var
  accpoint  : boolean;
  c         : rvector;
  count     : integer;
  cycle     : integer;
  cyclimit  : integer;
  f         : real;
  funcount  : integer;
  g         : rvector;
  G1, G2     : real;
  G3, gradproj     : real;
  gradcount : integer;
  i, j       : integer;
  method    : methodtype;
  newstep   : real;
  notcomp   : boolean;
  oldstep   : real;
  s         : real;
  setstep   : real;
  steplength: real;
  t         : rvector;
  tol       : real;

begin
  writeln('alg22.pas -- Nash Algorithm 22 version 2 1988-03-24');
  writeln('  Conjugate gradients function minimiser');
  writeln('Steplength saving factor multiplies best steplength found at the');
  writeln('  end of each iteration as a starting value for next search');
  write('Enter a steplength saving factor (sugg. 1.7) -- setstep ');
  readln(setstep);
  writeln(setstep);
  write('Choose method (1=FR, 2=PR, 3=BS) ');
  readln(i);  writeln(i);
  case i of
    1: method:=Fletcher_Reeves;
    2: method:=Polak_Ribiere;
    3: method:=Beale_Sorenson;
    else halt;
  end;
  case method of
    Fletcher_Reeves: writeln('Method: Fletcher Reeves');
    Polak_Ribiere: writeln('Method: Polak Ribiere');
    Beale_Sorenson: writeln('Method: Beale Sorenson');
  end;
  fail:=false;
  cyclimit:=n;
  if intol<0.0 then intol:=1e-14; {artificial eps}
  tol:=intol*n*sqrt(intol);
  writeln('tolerance used in gradient test=', tol);
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
    gradcount:=0;
    repeat
      for i:=1 to n do
      begin
        t[i]:=0.0;
        c[i]:=0.0;
      end;
      cycle:=0;
      oldstep:=1.0;
      count:=0;
      repeat
        cycle:=cycle+1;
        writeln(gradcount, ' ', funcount, ' ', Fmin);
        write('parameters ');
        for i:=1 to n do
        begin
          write(Bvec[i]:10:5, ' ');
          if (7 * (i div 7) = i) and (i<n) then writeln;
        end;
        writeln;
        gradcount:=gradcount+1;
        fmingr(n, Bvec, Workdata, g);
        G1:=0.0; G2:=0.0;
        for i:=1 to n do
        begin
          X[i]:=Bvec[i];
          case method of
            Fletcher_Reeves: begin
              G1:=G1+sqr(g[i]); G2:=G2+sqr(c[i]);
            end;
            Polak_Ribiere  : begin
              G1:=G1+g[i]*(g[i]-c[i]); G2:=G2+sqr(c[i]);
            end;
            Beale_Sorenson : begin
              G1:=G1+g[i]*(g[i]-c[i]); G2:=G2+t[i]*(g[i]-c[i]);
            end;
          end;
          c[i]:=g[i];
        end;
        if G1>tol then
        begin
          if G2>0.0 then G3:=G1/G2 else G3:=1.0;
          gradproj:=0.0;
          for i:=1 to n do
          begin
            t[i]:=t[i]*G3-g[i]; gradproj:=gradproj+t[i]*g[i];
          end;
          steplength:=oldstep;
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
                steplength:=steplength*stepredn;
                write('*');
              end;
            end;
          until (count=n) or accpoint;
          if count<n then
          begin
            newstep:=2*((f-Fmin)-gradproj*steplength);
            if newstep>0 then
            begin
              newstep:=-gradproj*sqr(steplength)/newstep;
              for i:=1 to n do
              begin
                Bvec[i]:=X[i]+newstep*t[i];
              end;
              Fmin:=f;
              f:=fminfn(n, Bvec, Workdata, notcomp);
              funcount:=funcount+1;
              if f<Fmin then
              begin
                Fmin:=f; write(' i< ');
              end
              else
              begin
                write(' i> ');
                for i:=1 to n do Bvec[i]:=X[i]+steplength*t[i];
              end;
            end;
          end;
        end;
        oldstep:=setstep*steplength;
        if oldstep>1.0 then oldstep:=1.0;
      until (count=n) or (G1<=tol) or (cycle=cyclimit);

    until (cycle=1) and ((count=n) or (G1<=tol));

  end;
  writeln('Exiting from Alg22.pas conjugate gradients minimiser');
  writeln('    ', funcount, ' function evaluations used');
  writeln('    ', gradcount, ' gradient evaluations used');
end;
{fnmin.pas -- a main program to run function minimisation procedures}
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
fminset(n,Bvec,Workdata); {sets up problem and defines starting values of Bvec}
mytol:=-1.0; {Note: set the tolerance negative to indicate that procedure
must obtain an appropriate value.}
cgmin(n,Bvec,X,Fmin,Workdata,fail,mytol); {minimise the function}
writeln;
writeln(' Minimum function value found =',Fmin);
writeln(' At parameters');
for i:=1 to n do
begin
writeln(' Bvec[',i,']=',X[i]);
end; {loop to write out parameters}
end. {dr22.pas -- conjugate gradients minimisation}
