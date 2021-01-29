Program runsvd(input,output);
{dr0102.pas == Calculation of Singular values and vectors of an arbitrary
          real matrix, solution of linear least squares approximation
          problem.

  Modifies a method due to Kaiser.  See Nash and Shlien (1987): Simple
  algorithms for the partial singular value decomposition.  Computer
  Journal, vol.30, pp.268-275.

          Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
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


Procedure matcopy(nRow ,nCol: integer; A: rmatrix; var B:wmatrix);
{matcopy.pas
  -- copies matrix A, nRow by nCol, into matrix B }
var i,j: integer;
begin
  for i:=1 to nRow do
    for j:=1 to nCol do
      B[i,j]:=A[i,j];
end;{matcopy.pas}


Procedure PrtSVDResults( nRow, nCol:integer;
                U, V: rmatrix; Z: rvector);
{psvdres.pas
  == routine to display svd results and print them to confile
}
var
  i, j : integer;

begin
  writeln(' Singular values and vectors:');
  for j := 1 TO nCol do
  begin
    writeln('Singular value (',j,') =', Z[j]);
    writeln('Principal coordinate (U):');
    for i := 1 to nRow do
    begin
      write(U[i,j]:10:7);
      if (7 * (i div 7) = i) and (i<nRow) then writeln;
    end;
    writeln;
    writeln('Principal component (V):');
    for i:=1 to nCol do
    begin
      write(V[i,j]:10:7);
      if (7 * (i div 7) = i) and (i<nCol) then writeln;
    end;
    writeln;
  end;
end; {psvdres == print svd results via procedure PrtSVDResults }


Procedure svdtst( A, U, V: rmatrix; Z: rvector;
              nRow, UCol, VCol: integer);
{svdtst.pas
  == This routine tests the results of a singular value
  decomposion calculation.  The matrix A is presumed to contain
  the matrix of which the purported decomposition is

        U  Z  V-transpose

  This routine tests column orthogonality of U and V,
  row orthogonality of V, and the reconstruction suggested
  by the decomposition. It does not carry out the tests of
  the Moore-Penrose inverse A+, which can be computed as

      A+ :=  V  Z  U-transpose.

  FORTRAN codes for the conditions

          A+ A A+  = ? =  A+
          A  A+ A  = ? =  A
          (A+  A)-transpose  = ? =  A+ A
          (A  A+)-transpose  = ? =  A  A+

  are given in Nash, J.C. and Wang, R.L.C. (1986)
}

var
  i,j,k:integer;
  t1: real;
  imax, jmax: integer;
  valmax: real;

begin
  writeln('Column orthogonality of U');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to UCol do
  begin
    for j:=i to UCol do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1;
      for k:=1 to nRow do t1:=t1+U[k,i]*U[k,j];
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,'=',valmax);
  writeln('Row orthogonality of U (NOT guaranteed in svd)');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to nRow do
  begin
    for j:=i to nRow do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1;
      for k:=1 to UCol do t1:=t1+U[i,k]*U[j,k];
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,'=',valmax);
  writeln('Column orthogonality of V');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to VCol do
  begin
    for j:=i to VCol do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1.0;
      for k:=1 to VCol do t1:=t1+V[k,i]*V[k,j];
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,'=',valmax);
  writeln('Row orthogonality of V');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to VCol do
  begin
    for j:=i to VCol do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1;
      for k:=1 to VCol do t1:=t1+V[i,k]*V[j,k];
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,'=',valmax);
  writeln('Reconstruction of initial matrix');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to nRow do
  begin
    for j:=1 to VCol do
    begin
      t1:=0;
      for k:=1 to VCol do
        t1:=t1+U[i,k]*Z[k]*V[j,k]; { U * S * V-transpose}
      {writeln('A[',i,',',j,']=',A[i,j],' Recon. =',t1,'  error=',A[i,j]-t1);}
      if abs(A[i,j]-t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=A[i,j]-t1;
      end;
    end;
  end;
  writeln('Largest error is ',imax,',',jmax,'=',valmax);
end; {svdtst.pas}


{I matrixin.pas} {input or generate a matrix of reals}
{I vectorin.pas} {input or generate a vector of reals}


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


procedure svdlss(nRow, nCol: integer;
                 W : wmatrix;
                 Y: rvector;
                 Z : rvector;
                 A : rmatrix;
                 var Bvec: rvector;
                 q : real);

var
 i, j, k : integer;
 s : real;

begin
  writeln('alg02.pas == svdlss');
{  write('Y:');
  for i := 1 to nRow do 
  begin
    write(Y[i],' ');
  end;
  writeln;

  for i := 1 to (nRow+nCol) do
  begin
     write('W row ',i,':');
     for j:= 1 to nCol do
     begin
       write(W[i,j],' ');
     end;
     writeln;
   end;
}
{    writeln('Singular values');
    for j := 1 to nCol do
    begin
      write(Z[j]:18,' ');
      if j=4 * (j div 4) then writeln;
    end;
    writeln;
}
    if q>=0.0 then
    begin
    q := q*q; 
      for i := 1 to nCol do
      begin
        s := 0.0;
        for j := 1 to nCol do
        begin
          for k := 1 to nRow do
          begin
            if Z[j]>q then
              s := s + W[i+nRow,j]*W[k,j]*Y[k]/Z[j];
                       { V   S+   U'  y }

          end;
        end;
        Bvec[i] := s;
      end;
      writeln('Least squares solution');
      for j := 1 to nCol do
      begin
        write(Bvec[j]:12,' ');
        if j=5 * (j div 5) then writeln;
      end;
      writeln;
      s := resids(nRow, nCol, A, Y, Bvec, true);
    end;
end;

{main program}
var
  nRow, nCol : integer;
  A, V, U : rmatrix;
  W : wmatrix; {a working matrix which will contain U Zd in the
    upper nRow rows, and V in the bottom nCol rows, where Zd
    is the diagonal matrix of singular values.  That is, W
    becomes

                  (  U  Zd )
                  (        )
                  (    V   )

  }
  Z, Zsq : rvector; {Z will contain either the squares of singular
          values or the singular values themselves}
  Y : rvector; {Y will contain the 'right hand side' of the
          least squares problem, i.e.  the vector to be
          approximated }
  Bvec : rvector; {the least squares solution }
  inchar : char;
  i,j,k, imax, jmax : integer;
  t1, t2: real;

begin
  banner:='dr0102.pas -- driver for svd and least squares solution';
  {Test matrix from CNM pg 34}
  nRow:=4;
  nCol:=3;
  {Read in matrix the hard way!}
  A[1,1]:=5; A[1,2]:=1.0E-6; A[1,3]:=1; Y[1]:=1;
  A[2,1]:=6; A[2,2]:=0.999999; A[2,3]:=1; Y[2]:=2;
  A[3,1]:=7; A[3,2]:=2.00001; A[3,3]:=1; Y[3]:=3;
  A[4,1]:=8; A[4,2]:=2.9999; A[4,3]:=1; Y[4]:=4;
  
  Matcopy(nRow,nCol, A, W); {The matrix A is copied into working array W.}
  NashSVD( nRow, nCol, W, Z); {The singular value decomposition is
        computed for matrix A by columnwise orthogonalization of the
        working array W, to which a unit matrix of order nCol is added
        in order to form the matrix V in the bottom nCol rows of W.}
  begin
    for j:=1 to nCol do
    begin
      Zsq[j] := Z[j];
      Z[j]:= sqrt(Z[j]);
      for i:=1 to nRow do U[i,j]:=W[i,j]/Z[j];
      for i:=1 to nCol do V[i,j]:=W[i+nRow,j];
    end;
    PrtSVDResults( nRow, nCol, U, V,Z);
    begin
      svdtst(A,U,V,Z,nRow,nCol,nCol);
      writeln('Reconstruction of initial matrix from Nash working form');
      t2:=0.0; {to store largest error in reconstruction}
      for i:=1 to nRow do
      begin
        for j:=1 to nCol do
        begin
          t1:=0.0;
          for k:=1 to nCol do
            t1:=t1+W[i,k]*W[j+nRow,k]; { U * S * V-transpose}
          t1:=A[i,j]-t1; {to compute the residual}
          if abs(t1)>t2 then
          begin
            t2:=abs(t1); imax:=i; jmax:=j; {to save biggest element}
          end;
        end; {loop over columns}
      end; {loop over rows}
      writeln('Largest error is ',imax,',',jmax,'=',t2);
    end; {test svd results}
  end; {print results}
  svdlss(nRow, nCol, W, Y, Zsq, A, Bvec, 1.0e-16); 
end. {dr0102.pas == svd and least squares solution}
