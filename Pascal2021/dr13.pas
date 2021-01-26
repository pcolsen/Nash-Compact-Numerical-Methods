Program dr13(input,output);
{dr13.pas == run Nash svd for eigenvalue computations (Alg13)

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

}
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

function resids(nRow, nCol: integer; A : rmatrix;
          Y: rvector; Bvec : rvector):real;
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
  writeln('Residuals');
  ss:=0.0;
  for i:=1 to nRow do
  begin
    t1:=-Y[i]; {note form of residual is residual = A * B - Y }
    for j:=1 to nCol do
      t1:=t1+A[i,j]*Bvec[j];
    ss:=ss+t1*t1;
    write(t1:10,' ');
    if (i = 7 * (i div 7)) and (i<nRow) then writeln;
  end; {loop on i}
  writeln;
  writeln('Sum of squared residuals =',ss);
  resids:=ss
end; {resids.pas == residual calculation for linear least squares}

procedure NashSVD(nRow, nCol: integer;
               var W: wmatrix;
               var Z: rvector);

var
  i, j, k, EstColRank, RotCount, SweepCount, slimit : integer;
  eps, e2, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2 : real;

procedure rotate;
var
  ii : integer;

begin
  for ii := 1 to nRow+nCol do
  begin
    D1 := W[ii,j]; D2 := W[ii,k];
    W[ii,j] := D1*c0+D2*s0; W[ii,k] := -D1*s0+D2*c0
  end;
end;

begin

  writeln('alg01.pas -- NashSVD');
  eps := Calceps;
  slimit := nCol div 4;  if slimit<6 then slimit := 6;

  SweepCount := 0;
  e2 := 10.0*nRow*eps*eps;
  tol := eps*0.1;

  EstColRank := nCol; ;

  for i := 1 to nCol do
    begin
    for j := 1 to nCol do
      W[nRow+i,j] := 0.0;
    W[nRow+i,i] := 1.0;
  end;

  repeat
    RotCount := EstColRank*(EstColRank-1) div 2;
    SweepCount := SweepCount+1;

    for j := 1 to EstColRank-1 do
    begin
      for k := j+1 to EstColRank do
      begin
        p := 0.0;  q := 0.0; r := 0.0;
        for i := 1 to nRow do
        begin
          x0 := W[i,j]; y0 := W[i,k];
          p := p+x0*y0; q := q+x0*x0;  r := r+y0*y0;
        end;
        Z[j] := q; Z[k] := r;

        if q >= r then
        begin
          if (q<=e2*Z[1]) or (abs(p)<= tol*q) then RotCount := RotCount-1

          else
          begin
            p := p/q; r := 1-r/q; vt := sqrt(4*p*p + r*r);
            c0 := sqrt(0.5*(1+r/vt)); s0 := p/(vt*c0);
            rotate;
          end
        end
        else
        begin

          p := p/r; q := q/r-1; vt := sqrt(4*p*p + q*q);
          s0 := sqrt(0.5*(1-q/vt));
          if p<0 then s0 := -s0;
          c0 := p/(vt*s0);
          rotate;
        end;

      end;
    end;
    writeln('End of Sweep #', SweepCount,
            '-  no. of rotations performed =', RotCount);
    while (EstColRank >= 3) and (Z[EstColRank] <= Z[1]*tol + tol*tol)
          do EstColRank := EstColRank-1;
  until (RotCount=0) or (SweepCount>slimit);
  if (SweepCount > slimit) then writeln('**** SWEEP LIMIT EXCEEDED');
end;

Procedure evsvd(n: integer; var A,V : rmatrix; initev: boolean;
             W : wmatrix; var Z: rvector);

var
  i, j: integer;
  shift, t : real ;
  
begin
  writeln('alg13.pas -- symmetric matrix eigensolutions via svd');
  shift:=0.0;
  for i:=1 to n do
  begin
    t:=A[i,i];
    for j:=1 to n do
      if i<>j then t:=t-abs(A[i,j]);
    if t<shift then shift:=t;
  end;
  shift:=-shift;
  if shift<0.0 then shift:=0.0;
  writeln('Adding a shift of ',shift,' to diagonal of matrix.');
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
      W[i,j]:=A[i,j];
      if i=j then W[i,i]:=A[i,i]+shift;
      if initev then
      begin
        if i=j then W[i+n,i]:=0.0
        else
        begin
          W[i+n,j]:=0.0;
        end;
      end;
    end;
  end;
  if (n > 1) then 
     NashSVD(n, n, W, Z)
  else
  begin { order 1 matrix }
     Z[1] := A[1,1]*A[1,1];
     W[2,1]:= 1.0; {Eigenvector!}
  end;
  for i:=1 to n do
  begin
    Z[i]:=sqrt(Z[i])-shift;
    for j:=1 to n do
      V[i,j]:=W[n+i,j];
  end;
end;

Procedure Frank2(var m, n: integer; var A: rmatrix);
var
  i,j: integer;
begin
    for i:=1 to m do
    begin
        for j:=1 to n do
        begin
          if (i <= j) then
             A[i,j]:=i
          else
             A[i,j]:=j;
        end;
    end;
end;

var
  i, j, nRow, nCol : integer;
  A, V, ACOPY : rmatrix;
  Bvec, Y, Z : rvector;
  W : wmatrix; {to store the working array}
  t1: real;
  initev : boolean;

begin
  banner:='dr13.pas -- driver for svd eigensolutions of a symmetric matrix';
  nRow := 1; {To get loop going}
  while (nRow > 0) do
  begin
  write('Order of problem (n): '); readln(nRow);
  if (nRow <= 0) then halt;
  nCol := nRow;
  Frank2(nRow, nCol, A);
  writeln('Initial matrix of order ', nRow);
  for j := 1 to nRow do
  begin
    for i := 1 to nRow do
    begin
      write(A[i,j]:10:5,' ');
      ACOPY[i,j] := A[i,j];
      W[i,j]:=0.0; {to avoid warning 'uninitialized' from fpc}
      if (7 * (i div 7) = i) and (i<nRow) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
  initev := true; {Here we want to get the eigenvectors of A, not some
            generalized problem.}
  writeln('Calling evsvd');
  evsvd( nRow, A, V, initev, W, Z);
  for j := 1 to nRow do
  begin
    t1 := Z[j];
    writeln;
    writeln('Eigenvalue ',j,' = ',t1);
    for i := 1 to nRow do
    begin
      write(V[i,j]:10:7,' ');
      if (i = 7 * (i div 7)) and (i<nRow) then
      begin
        writeln;
      end;
      Bvec[i] := V[i,j]; {to initialize residual test}
      Y[i] := t1*Bvec[i];
    end;
    writeln;
    t1 := resids(nRow, nCol, ACOPY, Y, Bvec);
  end; {loop on solutions j}
  end; {main while loop}
end. {dr13.pas}
