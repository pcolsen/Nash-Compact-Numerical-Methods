program dr0506(input,output);
{dr05606.pas == driver for gelim (alg05) and gebacksub (alg06),
    Gauss elimination and linear equations solution

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

Procedure gelim( n : integer;
                 p : integer;
                var A : rmatrix;
                 tol : real);


var
  det, s : real;
  h,i,j,k: integer;

begin
  det := 1.0;
  writeln('alg05.pas -- Gauss elimination with partial pivoting');
  for j := 1 to (n-1) do
  begin
    s := abs(A[j,j]); k := j;
    for h := (j+1) to n do
    begin
      if abs(A[h,j])>s then
      begin
        s := abs(A[h,j]); k := h;
      end;
    end;
    if k<>j then
    begin
      writeln('Interchanging rows ',k,' and ',j);
      for i := j to (n+p) do
      begin
        s := A[k,i]; A[k,i] := A[j,i]; A[j,i] := s;
      end;
      det := -det;
    end;
    det := det*A[j,j];
    if abs(A[j,j])<tol then
    begin
      writeln('Matrix computationally singular -- pivot < ',tol);
      halt;
    end;
    for k := (j+1) to n do
    begin
      A[k,j] := A[k,j]/A[j,j];
      for i := (j+1) to (n+p) do
          A[k,i] := A[k,i]-A[k,j]*A[j,i];
    end;
    det := det*A[n,n];
    if abs(A[n,n])<tol then
    begin
      writeln('Matrix computationally singular -- pivot < ',tol);
      halt;
    end;
  end;
  writeln('Gauss elimination complete -- determinant = ',det);
end;

procedure gebacksub(n, p:integer;
                     var A : rmatrix);


var
  s : real;
  i, j, k: integer;

begin
  writeln('alg06.pas -- Gauss elimination back-substitution');
  for i:=(n+1) to (n+p) do
  begin
    A[n,i]:=A[n,i]/A[n,n];
    for j:=(n-1) downto 1 do
    begin
      s:=A[j,i];
      for k:=(j+1) to n do
      begin
        s:=s-A[j,k]*A[k,i];
      end;
      A[j,i]:=s/A[j,j];
    end;
  end;
end;

var
  banner     : string[80]; {program name and description}
  i,j,k,n, nRHS,nRow,ncase : integer;
  rss : real; {to accumulate residual sum of squares}
  s : real;
  tol : real; {tolerance for zero pivot}
  A, Acopy : rmatrix;

begin
  banner:='dr0506.pas -- Gauss elimination and linear equations';

  for ncase:=1 to 2 do
  begin  {create A}
     n:=4*ncase;
     nRow:=n;
     nRHS:=1;
     for i:=1 to n do
     begin
        for j:=1 to i do
        begin
            A[i, j] := i;
            A[j, i] := i;
        end;
     end;
     for i:=1 to n do
     begin
        s:=0;
        for j:=1 to n do s:=s+A[i,j];
        A[i, n+1] := s;
     end;
     writeln('Data matrix :',nRow,' by ',nRow+nRHS);
     for i:=1 to nRow do
     begin
       writeln('Row ',i);
       for j:=1 to (nRow+nRHS) do
       begin
          write(A[i,j]:10:5,' ');
          if (7 * (j div 7) = j) and (j<nRow+nRHS) then writeln;
          Acopy[i,j]:=A[i,j];
       end;
       writeln;
     end;
     writeln;
     begin
       tol:=0.0;
       for i:=1 to nRow do
       begin
         s:=0.0;
         for j:=1 to nRow do s:=s+abs(A[i,j]);
         if s>tol then tol:=s;
       end;
       tol:=tol*calceps; {tol has maximum row norm * EPSILON}
       writeln('tol for pivod = ',tol);
     end;

    Gelim(nRow,nRHS,A,tol);
     {writeln('after Gauss elimination');}
     writeln('returned matrix ',nRow,' by ',nRow+nRHS);
     for i:=1 to nRow do
     begin
       writeln('Row ',i);
       for j:=1 to (nRow+nRHS) do
       begin
          write(A[i,j]:10:5,' ');
          if (7 * (j div 7) = j) and (j<nRow+nRHS) then writeln;
       end;
     end;

     writeln;
     GEBacksub(nRow,nRHS,A);
     for i:=1 to nRHS do
     begin
        writeln('Solution ',i);
        for j:=1 to nRow do
        begin
          write(A[j,nRow+i]:10:5,' ');
          if (7 * (j div 7) = j) and (j<nRow) then writeln;
        end;
        writeln;

        writeln('Residuals');
        rss:=0.0; {initialize sum of squared residuals}
        {Rather than use resids.pas, this simple code has been placed in line.}

        for j:=1 to nRow do
        begin
          s:=Acopy[j,nRow+i]; {right hand side}
          for k:=1 to nRow do s:=s-Acopy[j,k]*A[k,nRow+i];
          write(s:10,' ');
          rss:=rss+s*s;
          if (7 * (j div 7) = j) and (j<nRow) then writeln;
        end;
        writeln;
        writeln('Sum of squared residuals = ',rss);
     end; {loop over solutions}

  end; { loop over n }
end. {dr0506.pas}
