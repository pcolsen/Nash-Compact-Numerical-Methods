program dr23(input,output);
{dr23.PAS == driver for Nash Marquardt nonlinear least squares
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
{*** replace 'rosen' with 'logistic' in above line to run 3 parameter
logistic function fit ***}
procedure choldcmp(n: integer;
                   var a: smatvec;
                   var singmat: boolean);
var
  i,j,k,m,q: integer;
  s : real;

begin
  singmat := false;
  for j := 1 to n do
  begin
    q := j*(j+1) div 2;
    if j>1 then
    begin
      for i := j to n do
      begin
        m := (i*(i-1) div 2)+j; s := a[m];
        for k := 1 to (j-1) do s := s-a[m-k]*a[q-k];
        a[m] := s;
      end;
    end;
    if a[q]<=0.0 then
    begin
      singmat := true;
      a[q] := 0.0;
    end;
    s := sqrt(a[q]);
    for i := j to n do
    begin
      m := (i*(i-1) div 2)+j;
      if s=0.0 then a[m] := 0
          else a[m] := a[m]/s;
    end;
  end;
end;
procedure cholback(n: integer;
                   a: smatvec;
                   var x: rvector);


var
 i,j,q  : integer;

begin
  if a[1]=0.0 then x[1]:=0.0
              else x[1]:=x[1]/a[1];
  if n>1 then
  begin
    q:=1;
    for i:=2 to n do
    begin
      for j:=1 to (i-1) do
      begin
        q:=q+1; x[i]:=x[i]-a[q]*x[j];
      end;
      q:=q+1;
      if a[q]=0.0 then x[i]:=0.0
                  else x[i]:=x[i]/a[q];
    end;
  end;

  if a[n*(n+1) div 2]=0.0 then x[n]:=0.0
                           else x[n]:=x[n]/a[n*(n+1) div 2];
  if n>1 then
  begin
    for i:=n downto 2 do
    begin
      q:=i*(i-1) div 2;
      for j:=1 to (i-1) do x[j]:=x[j]-x[i]*a[q+j];
      if a[q]=0.0 then x[i-1]:=0.0
                  else x[i-1]:=x[i-1]/a[q];
    end;
  end;
end;
procedure modmrt( n : integer;
          var Bvec : rvector;
          var X : rvector;
          var Fmin  : real;
            Workdata : probdata);

{modified 1991 - 01 - 13}
var
  a, c: smatvec;
  delta, v : rvector;
  dec, eps, inc, lambda, p, phi, res : real;
  count, i, ifn, igrad, j, k, nn2, q : integer;
  notcomp, singmat, calcmat: boolean;

begin
  writeln('alg23.pas -- Nash Marquardt nonlinear least squares');
  with Workdata do
  begin
    if nlls = false then halt;
    Fmin:=big;
    inc:=10.0;
    dec:=0.4;
    eps:=1.0e-14; {Use a fixed value of double float machine precision}
    lambda:=0.0001;
    phi:=1.0;
    ifn:=0;  igrad:=0;
    calcmat:=true;
    nn2:=(n*(n+1)) div 2;
    p:=0.0;
    for i:=1 to m do
    begin
      res:=nlres(i, n, Bvec, notcomp, Workdata);

      if notcomp then halt;
      p:=p+res*res;
    end;
    ifn:=ifn+1;
    Fmin:=p;
    count:=0;

    while count<n do
    begin

      if calcmat then
      begin
        writeln(igrad,' ',ifn,'  sum of squares=',Fmin);
        for i:=1 to n do
        begin
          write(Bvec[i]:10:5,' ');
          if (7 * (i div 7) = i) and (i<n) then
          begin
            writeln;
          end;
        end;
        writeln;
        igrad:=igrad+1;
        for j:=1 to nn2 do a[j]:=0.0;
        for j:=1 to n do v[j]:=0.0;
        for i:=1 to m do
        begin
          nljac(i, n, Bvec, X, workdata);
          res:=nlres(i, n, Bvec, notcomp, Workdata);
          for j:=1 to n do
          begin
            v[j]:=v[j]+X[j]*res;
            q:=(j*(j-1)) div 2;
            for k:=1 to j do a[q+k]:=a[q+k]+X[j]*X[k];
          end;
        end;
        for j:=1 to nn2 do c[j]:=a[j];
        for j:=1 to n do X[j]:=Bvec[j];
      end;
      writeln('LAMDA =',lambda:8);
      for j:=1 to n do
      begin
        q:=(j*(j+1)) div 2;
        a[q]:=c[q]*(1.0+lambda)+phi*lambda;
        delta[j]:=-v[j];
        if j>1 then
          for i:=1 to (j-1) do a[q-i]:=c[q-i];
      end;
      notcomp:=false;
      Choldcmp(n, a, singmat);
      if (not singmat) then
      begin
        Cholback(n, a, delta);
        count:=0;
        for i:=1 to n do
        begin
          Bvec[i]:=X[i]+delta[i];
          if (reltest + Bvec[i])=(reltest+X[i]) then count:=count+1;
        end;
        if count<n then
        begin
          p:=0.0; i:=0;
          repeat
            i:=i+1; res:=nlres(i,n,Bvec,notcomp, Workdata);
            if (not notcomp) then  p:=p+res*res;
          until notcomp or (i>=m);  {MODIFICATION m replaces n 1991-01-13}
          ifn:=ifn+1;
        end;
      end;
      if count<n then
        if (not singmat) and (not notcomp) and (p<Fmin) then
        begin
          lambda:=lambda*dec;
          Fmin:=p;
          calcmat:=true;
        end
      else
      begin
        lambda:=lambda*inc;
        if lambda<eps*eps then lambda:=eps;
        calcmat:=false;
      end;

    end;
  end;
end;
{main program}
var
n          : integer; {the order of the problem}
Bvec       : rvector; {current set of parameters}
X          : rvector; {"best" set of parameters}
Workdata   : probdata;
i          : integer;
Fmin       : real;
fail       : boolean;
mytol      : real;
begin
banner:='dr23.pas -- Marquardt Nash nonlinear least squares';
fminset(n,Bvec,Workdata);
{Sets up problem and defines starting values of Bvec}
mytol:=-1.0; {Note: set the tolerance negative to indicate that
procedure must obtain an appropriate value.}
Workdata.nlls:=true;
modmrt( n, Bvec, X, Fmin, Workdata);
writeln;
writeln(' Minimum function value found =',Fmin);
writeln(' At parameters');
for i:=1 to n do
begin
writeln(' Bvec[',i,']=',X[i]);
end;
{END dr23.pas == Nash Marquardt nonlinear least squares}
end.

