program dr24ii(input, output);
{dr24ii.pas == modification of conjugate gradients method to solve
generalized symmetric matrix eigenvalue problems by inverse
iteration
Note, this approach ignores the definiteness of the matrices
which is supposedly needed for proper operation of the conjugate
gradients method.
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

function rayquo(n :integer;
            A,B : rmatrix;
            Y : rvector): real;
{rayquo.pas
  == compute Rayleigh quotient. If the denominator of
      the quotient is zero, set the quotient to a large
      negative number (-big)
}
var
  s,t : real;
  i,j : integer;

begin
  s:=0.0; t:=0.0;
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
      s:=s+Y[i]*A[i,j]*Y[j];
      t:=t+Y[i]*B[i,j]*Y[j];
    end; {loop on j}
  end; {loop on i}
  if t>0.0 then rayquo:=s/t else rayquo:=-big;
  {note the safety value for the result}
end; {rayquo.pas == compute Rayleigh quotient}


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

function genevres(n: integer; {order of matrices}
              A, B: rmatrix;
              evalue : real; {eigenvalue}
              X : rvector; {presumed eigenvector}
              print: boolean) {true for printing}
              : real; {returns sum of squared residuals}

{genevres.pas

  -- to compute residuals of generalized (symmetric) matrix
  eigenvalue problem

        A x = evalue B x
}
var
  t, ss : real;
  i,j : integer;

begin
  if print then
  begin
    writeln('Generalized matrix eigensolution residuals');
  end;
  ss:=0.0; {to accumulate the sum of squared residuals}
  for i:=1 to n do
  begin
    t:=0.0;
    for j:=1 to n do
    t:=t+(A[i,j]-evalue*B[i,j])*X[j];
    if print then
    begin
      write(t:10,' ');
    if (7 * (i div 7) = i) and (i<n) then
    begin
      writeln;
    end;
    end; {if print}
    ss:=ss+t*t;
  end;
  if print then
  begin
    writeln;
    writeln('Sum of squared residuals =',ss);
  end; {if print}
  genevres:=ss; {return sum of squared residuals}
end; {genevres.pas == residuals for generalized eigenproblem}
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
  eps, g2, oldg2, s2, step, steplim, t2, tol : real;
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

var
A, B, H : rmatrix;
Y : rvector; {RHS in solution of linear equations}
X : rvector; {solution}
RHS : rvector; {used for solution of linear equations sub-problem}
avec : smatvec; {for matrixin only}
shift : real; {for the eigenvalue shift}
evalue : real; {for the eigenvalue}
sym : boolean; {to tell if matrix symmetric}
ch : char;
i, j, m, n, itcount : integer;
iicount: integer; {to count iterations}
ssmin, t, s, ymax : real;
tol : real;
mchoice : integer; {to control flow in program}
snum : real; {to enter numbers as strings}
evcon : boolean;
oldev : real;
iilimit : integer;
rq : real; {Rayleigh quotient value}
begin
banner:='dr24ii.pas -- cg inverse iteration for matrix eigenproblem';
write('Order of problem = '); readln(n);
writeln(n);
for i:=1 to n do 
begin
  for j:=1 to i do
  begin
    A[i,j]:=j-2.0;
    A[j,i]:=j-2.0;
    B[i,j]:=j;
    B[j,i]:=j;
  end;
  A[i,i]:=i;
end;
writeln('Matrix A = Moler; Matrix B = Frank');
shift:=0.0; {safety setting}
write('Shift value for eigenvalues ([cr] =',shift,') = ');
readln(snum);
writeln(snum);
if snum <= -9999.0 then halt; {useful only for this program}
evalue:=shift; {initial guess}
writeln('Initial guess for eigenvector = 1s');
for i:=1 to n do Y[i]:=1.0;
{This is where inverse iteration loop begins.}
ssmin:=genevres(n,A,B,evalue,Y,false);
writeln('Residual sumsquares from trial solution = ',ssmin);
iicount:=0; iilimit:=4*n; {loose limit here}
tol:=calceps; tol:=tol*(abs(shift)+tol); {to be used as a small number}
oldev:=big; {to ensure no convergence in one iteration}
writeln('   ev itns':10,'   cg itns':10,'   cgss':14,
'   ev approx':16,'   RQ':16);
repeat {inverse iteration loop}
iicount:=iicount+1;
for i:=1 to n do {store shifted matrix in H}
for j:=1 to n do H[i,j]:=A[i,j]-shift*B[i,j];
for i:=1 to n do
begin {form RHS:=B * X as rhs of local equations}
t:=0.0;
for j:=1 to n do t:=t+B[i,j]*Y[j];
RHS[i]:=t; X[i]:=Y[i]; {save eigenvector elements for comparison}
end;
itcount:=n+1; {safety setting -- short iteration}
lecg( n, H, RHS, Y, itcount, ssmin);
{call to alg24.pas -- linear equations by conjugate gradients}
write(iicount:8,'  ',itcount:8,'    ',ssmin:12,'  ');
s:=abs(Y[1]); m:=1; {index of current maximum}
for i:=2 to n do {find largest eigenvector component}
begin
if abs(Y[i]) > s then
begin
m:=i; s:=abs(Y[i]);
end;
end; {loop on i}
evalue:=shift+X[m]/Y[m];
write(evalue:14,'  ');
{Normalize eigenvector.}
ymax:=Y[m];
for i:=1 to n do Y[i]:=Y[i]/ymax; {to normalize the eigenvector}
rq:=rayquo(n,A,B,Y);    {Compute Rayleigh quotient.}
writeln(rq:14); 
evcon:= ( (reltest+evalue) = (reltest+oldev) ); {check if evalue converged}
oldev:=evalue; {to save eigenvalue}
until ((ssmin<=tol) and (itcount>0) and evcon)
or (iicount>iilimit) or (ssmin=-big);
{loop until convergence or failure}
{Note this is a strict criterion. It could be relaxed somewhat and
still be satisfactory for many purposes.}
if iicount>iilimit then
begin
writeln('Iteration limit exceeded');
end;
ssmin:=genevres(n,A,B,evalue,Y,false);
writeln('Residual sumsquares of normalized vector = ',ssmin);
write('Enter a new shift, or -10000 to end ');
readln(snum); writeln(snum);
if (snum<-999.0) then halt;
end; {while mchoice<>=0}
writeln('Eigensolution for eigenvalue =',evalue);
writeln('           Rayleigh quotient =',rq);
for i:=1 to n do
begin
write(Y[i]:10:7,' '); 
if (7 * (i div 7) = i) and (i<n) then writeln; 
end;
writeln;
ssmin:=genevres(n,A,B,evalue,Y,true);
end. {dr24ii.pas == inverse iteration by conjugate gradients}
