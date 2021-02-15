program dr24le(input, output);
{dr24le.pas == linear equations by conjugate gradients
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

Procedure Frank(var n: integer; var A: rmatrix);
var
  i,j: integer;
begin {Frank symmetric}
  for i:=1 to n do
  begin
    for j:=1 to i do
    begin
      A[i,j]:=j;
      A[j,i]:=j;
    end;
  end;
end; {Frank}
procedure matmul(nn : integer; {order of matrix}
            matrix: rmatrix;{the matrix or order nn}
            vectorin: rvector;{vector which is multiplied}
          var vectorout: rvector); {product vector}
{matmul.pas == Here we use an explicit matrix multiplication. This may
  be replaced by implicit forms as appropriate.}
var
  ii, jj : integer;
  tt : real;
begin
  for ii := 1 to nn do
  begin
    tt := 0.0;
    for jj := 1 to nn do tt := tt+matrix[ii,jj]*vectorin[jj];
    vectorout[ii] := tt;
  end; {loop on ii}
end; {matmul.pas}

procedure lecg( n : integer;
            H : rmatrix;
            C : rvector;
        var Bvec : rvector;
        var itcount : integer;
        var ssmin : real);
var
  count, i, itlimit : integer;
  eps, g2, oldg2, step, steplim, t2, tol : real;
  g, t, v : rvector;

begin
  itlimit := itcount;
  itcount := 0;
  eps := calceps;
  steplim := 1.0/sqrt(eps);
  g2 := 0.0;
  for i := 1 to n do g2 := g2+abs(C[i]);  tol := g2*eps*eps*n;
  matmul(n, H, Bvec, g);
  for i := 1 to n do g[i] := g[i]-C[i];
  g2 := 0.0;
  for i := 1 to n do
  begin
    g2 := g2+g[i]*g[i]; t[i] := -g[i];
  end;
  ssmin := big;
  while (g2>tol) and (itcount<itlimit) and (ssmin>0.0) do
  begin
    itcount := itcount+1;
    matmul( n, H, t, v);
    t2 := 0.0;
    for i := 1 to n do t2 := t2+t[i]*v[i];
    step := g2/t2; oldg2 := g2;
    if abs(step)>steplim then
    begin
      writeln('Step too large -- coefficient matrix indefinite?');
      ssmin := -big;
    end
    else
    begin
      g2 := 0.0; count := 0;
      for i := 1 to n do
      begin
        g[i] := g[i]+step*v[i];
        t2 := Bvec[i]; Bvec[i] := t2+step*t[i];
        if Bvec[i]=t2 then count := count+1;
        g2 := g2+g[i]*g[i];
      end;
      if count<n then
      begin
        if g2>tol then
        begin
          t2 := g2/oldg2;
          for i := 1 to n do t[i] := t2*t[i]-g[i];
        end;
      end;
      ssmin := g2;
    end;
  end;
  if itcount>=itlimit then itcount := -itcount;
end;
function resids(nRow, nCol: integer; A : rmatrix;
          Y: rvector; Bvec : rvector; print : boolean):real;
{resids.pas
  == Computes residuals and , if print is TRUE, displays them 7
    per line for the linear least squares problem. The sum of
    squared residuals is returned.
    residual vector = A * Bvec - Y
}
var
i, j: integer;
t1, ss : real;

begin
  if print then
  begin
    writeln('Residuals');
  end;
  ss:=0.0;
  for i:=1 to nRow do
  begin
    t1:=-Y[i]; {note form of residual is residual = A * B - Y }
    for j:=1 to nCol do
      t1:=t1+A[i,j]*Bvec[j];
    ss:=ss+t1*t1;
    if print then
    begin
      write(t1:10,' ');
      if (i = 7 * (i div 7)) and (i<nRow) then writeln;
    end;
  end; {loop on i}
  if print then
  begin
    writeln;
    writeln('Sum of squared residuals =',ss);
  end;
  resids:=ss
end; {resids.pas == residual calculation for linear least squares}

var
A : rmatrix;
Y : rvector; {RHS}
Bvec : rvector; {solution}
i, n, itcount : integer;
ssmin, s : real;
begin
banner:='dr24le.pas -- linear equations by conjugate gradients';
write('Order of problem = '); { 10 }
readln(n);
writeln(n);
writeln('Coefficient matrix');
Frank(n, A);  {4 == Frank}
(* RHS vector and initial guess *)
for i:=1 to n do
begin
  Y[i]:=1.0; Bvec[i]:=1.0;
end;
itcount:=10*n; {safety setting}
lecg( n, A, Y, Bvec, itcount, ssmin);
writeln('Solution after ',itcount,' iterations. Est. sumsquares ',ssmin);
for i:=1 to n do
begin
  write(Bvec[i]:10:5,' ');
  if (7 * (i div 7) = i) and (i<n) then writeln;
end;
writeln;
s:=resids(n, n, A, Y, Bvec, true);
end. {dr24le.pas}
