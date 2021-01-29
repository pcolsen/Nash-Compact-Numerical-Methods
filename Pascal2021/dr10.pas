program dr10(input,output);
{dr10.pas == driver to use Gauss elimination for inverse iteration
calculation of matrix eigensolutions

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
    if (7 * (i div 7) = i) and (i<n) then writeln;
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
procedure gii(nRow : integer;
             var A : rmatrix;
             var Y : rvector;
             var shift : real;
             var itcount: integer);

var
  i, itlimit, j, m, msame, nRHS :integer;
  ev, s, t, tol : real;
  X : rvector;

begin
  itlimit:=itcount;
  nRHS:=nRow;
  tol:=Calceps;
  s:=0.0;
  for i:=1 to nRow do
  begin
    X[i]:=Y[i];
    Y[i]:=0.0;
    for j:=1 to nRow do
    begin
      A[i,j]:=A[i,j]-shift*A[i,j+nRow];
      s:=s+abs(A[i,j]);
    end;
  end;
  tol:=tol*s;
  gelim(nRow, nRHS, A, tol);
  itcount:=0;
  msame :=0;
  while (msame<nRow) and (itcount<itlimit) do
  begin
    itcount:=itcount+1;
    m:=nRow; s:=X[nRow];
    X[nRow]:=Y[nRow];
    if abs(A[nRow,nRow])<tol then Y[nRow]:=s/tol
                             else Y[nRow]:=s/A[nRow,nRow];
    t:=abs(Y[nRow]);
    for i:=(nRow-1) downto 1 do
    begin
      s:=X[i]; X[i]:=Y[i];
      for j:=(i+1) to nRow do
      begin
        s:=s-A[i,j]*Y[j];
      end;
      if abs(A[i,i])<tol then Y[i]:=s/tol else  Y[i]:=s/A[i,i];
      if abs(Y[i])>t then
      begin
        m:=i; t:=abs(Y[i]);
      end;
    end;
    ev:=shift+X[m]/Y[m];
(*    writeln('Iteration ',itcount,'  approx. ev=',ev);*)

    t:=Y[m]; msame:=0;
    for i:=1 to nRow do
    begin
      Y[i]:=Y[i]/t;
      if reltest+Y[i] = reltest+X[i] then msame:=msame+1;

    end;

    if msame<nRow then
    begin
      for i:=1 to nRow do
      begin
        s:=0.0;
        for j:=1 to nRow do s:=s+A[i,j+nRow]*Y[j];
        X[i]:=s;
      end;
    end;
  end;
  if itcount>=itlimit then itcount:=-itcount;
  shift:=ev;
end;

var
  A, Acopy, B, Bcopy : rmatrix;
  Y : rvector;
  i, itcount, j, n, nRow, nCol : integer;
  rq, ss, Shift : real;
  vectorOK : boolean;

begin
  banner:='dr10.pas -- inverse iteration via Gauss elimination';
  repeat
    write('order of problem (n) = '); readln(n);
    writeln(n);
    if (n > 0) then 
    begin {main loop over examples}
      nRow:=n; nCol:=2*n; {store matrices in an array n by 2n}
      writeln('Provide the A matrix');
      Frank2(nRow, nCol, Acopy);
      writeln('A matrix');
      for i:=1 to n do
      begin
        for j:=1 to n do
          begin
          write(Acopy[i,j]:10:5,' ');
          if (7 * (j div 7) = j) and (j<n) then writeln;
        end;
        writeln;
      end;
      writeln('B matrix set to unit matrix');
      for i:=1 to n do 
      begin 
        for j:=1 to n do B[i,j]:=0.0;
        Bcopy[i,i]:=1.0;
      end;
      writeln('B matrix');
      for i:=1 to n do 
      begin
        for j:=1 to n do
        begin
          write(Bcopy[i,j]:10:5,' ');
          if (7 * (j div 7) = j) and (j<n) then writeln;
        end;
        writeln;
      end;
      shift:=0.0; {rem initial and safety value for the eigenvalue shift}
      vectorOK:=true; {approximate eigenvector will be all 1s for example}
      for i:=1 to n do Y[i]:=1.0; {Set starting vector}
      repeat
      if (not vectorOK) then
      begin
        writeln('Need a starting vector for inverse iteration');
        halt;
      end;
      vectorOK:=true; {set flag to indicate a vector is now in place}
(*
    writeln('Enter a shift for eigenvalues ([cr] = ',shift,') ');
    write(' A value > 1E+30 halts execution. Entry = ');
*)
      write('shift=?');
      readln(shift);
      writeln(shift);
      if (not (shift > 1e30)) then 
      begin
        for i:=1 to n do  {copy matrices into working matrices}
        begin
          for j:=1 to n do
          begin
            A[i,j]:=Acopy[i,j]; B[i,j]:=Bcopy[i,j];
            A[i,j+n]:=B[i,j]; {to provide work matrix for ALG10}
          end; {loop on j}
        end; {loop on i}
        itcount:=100; {rem fairly liberal bound}
        gii(n,  A , Y, shift, itcount);
        writeln;
        if itcount > 0 then
        begin
          writeln(
          'Converged to eigenvalue =',shift,'  in  ',itcount,' iterations');
        end
        else
        begin
          writeln('Not converged. Approximate eigenvalue=',shift,
                  ' after ',-itcount,' iterations');
        end; {else not converged}
        writeln('Eigenvector');
        for i:=1 to n do
        begin
          write(Y[i]:10:5,' ');
          if (7* (i div 7) = i) and (i<n) then  writeln;
        end; {loop on i}
        writeln;
        ss:=genevres(n, Acopy, Bcopy, shift, Y, false);
        rq:=rayquo(n, Acopy, Bcopy, Y);
        writeln('Rayleigh quotient = ',rq,'   sumsquared err=',ss:12:9);
      end; {if shift <=1e30}
    until (shift>1e30); {end of loop over possible shifts}
    end; {main loop over n block}
  until (n <= 0);
end. {dr10.pas}
