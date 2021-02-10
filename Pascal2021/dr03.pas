program givrun(input, output);
{dr03.PAS ==  driver for Givens' reduction of a matrix

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



procedure givens( nRow,nCol : integer;
                  var A, Q: rmatrix);
var
 i, j, k, mn: integer;
 b, c, eps, p, s : real;

begin
  writeln('alg03.pas -- Givens',chr(39),' reduction -- column-wise');
  mn  :=  nRow; if nRow>nCol then mn  :=  nCol;
  for i  :=  1 to nRow do
  begin
    for j := 1 to nRow do Q[i,j] := 0.0;
    Q[i,i] := 1.0;
  end;
  eps := calceps;
  for j := 1 to (mn-1) do
  begin
    for k := (j+1) to nRow do
    begin
      c := A[j,j]; s := A[k,j];
      b := abs(c); if abs(s)>b then b := abs(s);
      if b>0 then
      begin
        c := c/b; s := s/b;
        p := sqrt(c*c+s*s);
        s := s/p;
        if abs(s)>=eps then
        begin
          c := c/p;
          for i := 1 to nCol do
          begin
            p := A[j,i]; A[j,i] := c*p+s*A[k,i]; A[k,i] := -s*p+c*A[k,i];
          end;
          for i := 1 to nRow do
          begin
            p := Q[i,j]; Q[i,j] := c*p+s*Q[i,k]; Q[i,k] := -s*p+c*Q[i,k];
          end;
        end;
      end;
    end;
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
          write(i,' ',j,';');
          if (i <= j) then
             A[i,j]:=i
          else
             A[i,j]:=j;
          writeln(A[i,j]);
        end;
    end;
end;



var
  A, Q: rmatrix;
  i, j, k, nRow, nCol : integer;
  Acopy : rmatrix;
  s : real;

begin
  banner:='dr03.pas -- driver for Givens'+chr(39)+' reduction';
  nRow := 5;
  nCol := 3; {Specific to this example.}
  writeln('Size of problem (rows, columns)  (',nRow,', ',nCol,')');
  writeln('Frank matrix example');
  Frank2(nRow, nCol, A);
  writeln('Matrix A');
  for i:=1 to nRow do
  begin
    for j:=1 to nCol do
    begin
      Acopy[i,j]:=A[i,j];
      write(A[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nCol) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
  givens(nRow,nCol,A,Q);
  writeln('Decomposition');
  writeln('Q');
  for i:=1 to nRow do
  begin
    for j:=1 to nRow do
    begin
      write(Q[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nRow) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
  writeln('R');
  for i:=1 to nRow do
  begin
    for j:=1 to nCol do
    begin
      write(A[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nCol) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
  writeln('Q*R - Acopy');
  for i:=1 to nRow do
  begin
    for j:=1 to nCol do
    begin
      s:=-Acopy[i,j];
      for k:=1 to nRow do s:=s+Q[i,k]*A[k,j];
      write(s:10,' ');
      if (7 * (j div 7) = j) and (j<nRow) then
      begin
        writeln;
      end;
    end;
    writeln;
  end;
end. {dr03.pas == Givens' reduction driver}
