program dr09(input,output);
{dr09.pas == driver program to test procedure for the Bauer-Reinsch
          inversion of a symmetric positive definite real matrix stored
          in row-wise vector form

          Copyright 1988 J.C.Nash
}
{I constype.def}
{constype.def ==
  This file contains various definitions and type statements which are
  used throughout the collection of "Compact Numerical Methods".  In many
  cases not all definitions are needed, and users with very tight memory
  constraints may wish to remove some of the lines of this file when
  compiling certain programs.

  Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}
{uses Dos, Crt;} {Turbo Pascal 5.0 Modules}
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

Procedure Frank(var n: integer; var A: rmatrix; var avector: smatvec);
var
  i,j: integer;
begin
  writeln('Frank symmetric');
    for i:=1 to n do
    begin
        for j:=1 to i do
        begin
          A[i,j]:=j;
          A[j,i]:=j;
        end;
    end;
end;

Procedure mat2vec(var n: integer; var A: rmatrix; var avector: smatvec);
var
  i,j,k: integer;

begin {convert to vector form}
    k:=0; {index for vector element}
    for i:=1 to n do
    begin
      for j:=1 to i do
      begin
        k:=k+1;
        avector[k]:=A[i,j];
      end;
    end;
end; {matrixin}

Procedure vec2mat(var n: integer; var A: rmatrix; var avector: smatvec);
var
  i,j,k: integer;

  begin {convert to matrix form}
    k:=0; {index for vector element}
    for i:=1 to n do
    begin
      for j:=1 to i do
      begin
        k:=k+1;
        A[i,j]:=avector[k];
      end;
    end;
end; {matrixin}

{ I alg09.pas}
procedure brspdmi(n : integer;
                var avector : smatvec;
                var singmat : boolean);

var
  i,j,k,m,q : integer;
  s,t : real;
  X : rvector;

begin
  writeln('alg09.pas -- Bauer Reinsch inversion');
  singmat  :=  false;
  for k  :=  n downto 1 do
  begin
    if (not singmat) then
    begin
      s  :=  avector[1];
      if s>0.0 then
      begin
        m  :=  1;
        for i := 2 to n do
        begin
          q := m; m := m+i; t := avector[q+1]; X[i] := -t/s;

          if i>k then X[i] := -X[i];
          for j := (q+2) to m do
          begin
            avector[j-i] := avector[j]+t*X[j-q];
          end;
        end;
        q := q-1; avector[m] := 1.0/s;
        for i := 2 to n do avector[q+i] := X[i];
      end
      else
        singmat := true;
    end;
  end;
end;

var
  A, Ainverse : rmatrix;
  avector : smatvec;
  i, imax, j, jmax, k, n : integer;
  errmax, s : real;
  singmat: boolean;


BEGIN { main program }
  banner:='dr09.pas -- test Bauer Reinsch sym, posdef matrix inversion';
  writeln(banner);
  n:=4; {Fixed example size 20210113}
  Frank(n,A,avector);
  writeln;
  writeln('returned matrix of order ',n);
  begin
    for i:=1 to n do
    begin
        for j:=1 to n do
        begin
            write(A[i,j],' ');
        end;
        writeln;
    end;
  end;
  mat2vec(n, A
  begin
    writeln('Symmetric matrix -- Vector form');
    k := 0;
    for i := 1 to n do
    begin
      for j := 1 to i do
      begin
        k := k+1;
        write(avector[k]:10:5,' ');
      end;
      writeln;
    end;
  end;
  brspdmi(n, avector,singmat);
  if singmat then halt; {safety check}
  writeln('Computed inverse');
  k := 0; {initialize index to smatvec elements}
  for i := 1 to n do
  begin
    for j := 1 to i do
    begin
      k := k+1;
      write(avector[k]:10:5,' ');
      Ainverse[i,j] := avector[k]; {save square form of inverse}
      Ainverse[j,i] := avector[k];
      if (7 * (j div 7) = j) and (j<i) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
  {Compute maximum error in A * Ainverse and note where it occurs.}
  errmax := 0.0; imax := 0; jmax := 0;
  for i := 1 to n do
  begin
    for j := 1 to n do
    begin
      s := 0.0; if i=j then s := -1.0;
      for k := 1 to n do s := s + Ainverse[i,k]*A[k,j];
      {Note: A has not been altered, since avector was used.}
      if abs(s)>abs(errmax) then
      begin
        errmax := s; imax := i; jmax := j; {save maximum error, indices}
      end;
    end; {loop on j}
  end; {loop on i}
  writeln('Maximum element in Ainverse * A - 1(n) = ',errmax,
          '  position ',imax,',',jmax);
end. {dr09.pas == Bauer Reinsch inversion}
