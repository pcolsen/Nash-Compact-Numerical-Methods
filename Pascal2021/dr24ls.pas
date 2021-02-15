program dr24ls(input, output);
{dr24ls.pas == linear least squares by conjugate gradients
Note: this implementation uses the normal equations, which
are not recommended as a general approach.
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

Procedure Frankm(var m, n: integer; var A: rmatrix);
var
  i,j: integer;
begin {Frank symmetric}
  for i:=1 to m do
  begin
    for j:=1 to n do
    begin
      A[i,j]:=j;
      if (i <= j) then A[j,i]:=j;
    end;
  end;
end; {Frank}

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
  eps := calceps; {artificial machine epsilon}
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
var
A, AtransA : rmatrix;
Y, AtransY : rvector; {RHS}
Bvec : rvector; {solution}
i, j, k, nRow, nCol, itcount : integer;
ssmin, t, s : real;
begin
banner:='dr24ls.pas -- linear least squares by conjugate gradients';
write('Number of rows in coefficient matrix = ');
readln(nRow);
writeln(nRow);
write('Number of columns in coefficient matrix = ');
readln(nCol); {This is the order of the conjugate gradients problem.}
writeln(nCol);
writeln('Coefficient matrix');
Frankm(nRow, nCol, A);
writeln('RHS vector and initial guess all 1s');
for i:=1 to nRow do 
begin
  Y[i]:=1.0; Bvec[i]:=1.0;
end;
{Now form the normal equations.}
writeln('Normal equations -- coefficient matrix');
for i:=1 to nCol do
begin
  t:=0.0;
  for k:=1 to nRow do t:=t+A[k,i]*Y[k];
  AtransY[i]:=t;
  for j:=1 to nCol do
  begin
    s:=0.0; for k:=1 to nRow do s:=s+A[k,i]*A[k,j];
    AtransA[i,j]:=s;
    write(s:10:5,' ');
    if (7 * (j div 7) = j) and (j<nCol) then writeln;
  end;
  writeln;
end; {loop on i and normal equations build}
writeln('Normal equations - RHS');
for j:=1 to nCol do
begin
write(AtransY[j]:10:5,' ');
if (7 * (j div 7) = j) and (j<nCol) then writeln;
end;
writeln;
{***WARNING*** this is NOT a good way to solve this problem generally}
itcount:=10*nRow; {safety setting}
lecg( nCol, AtransA, AtransY, Bvec, itcount, ssmin);
writeln('Solution after ',itcount,
' iterations. Est. normal eqn. sumsquares ',ssmin);
for i:=1 to nCol do
begin
write(Bvec[i]:10:5,' ');
if (7 * (i div 7) = i) and (i<nCol) then writeln;
end;

writeln;
write('For original least squares problem -- ');
s:=resids(nRow, nCol, A, Y, Bvec, true);
end. {dr24ls.pas}
