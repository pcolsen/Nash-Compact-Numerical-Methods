program dr0708(input,output);
{dr0708.pas == driver program to test procedures for Choleski (Alg07)
      and Choleski back-substitution (Alg08)

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
    
Procedure matrixin(var m, n: integer; var A: rmatrix;
              var avector: smatvec; var sym :boolean);

{matrixin.pas --

  This procedure generates an m by n real matrix in both
  A or avector.

  A is of type rmatrix, an array[1..nmax, 1..nmax] of real
  where nmax >= n for all possible n to be provided.

  avector is of type rvector, an array[1..nmax*(nmax+1)/2]
  of real, with nmax as above.

  sym is set true if the resulting matrix is symmetric.
}

var
  temp : real;
  i,j,k: integer;
  inchar: char;
  mtype: integer;
  mn : integer;

begin
  if (m=0) or (n=0) then
  begin
    writeln;
    writeln('******* Matrix dimensions zero *********');
    halt;
  end;
  writeln('Matrixin.pas -- generate or input a real matrix ',m,' by ',n);
  writeln('Possible matrices to generate:');
  writeln('0) Keyboard or console file input');
  writeln('1) Hilbert segment');
  writeln('2) Ding Dong');
  writeln('3) Moler');
  writeln('4) Frank symmetric');
  writeln('5) Bordered symmetric');
  writeln('6) Diagonal');
  writeln('7) Wilkinson W+');
  writeln('8) Wilkinson W-');
  writeln('9) Constant');
  writeln('10) Unit');
{ Note: others could be added.}
  mn:=n;
  if m>mn then mn:=m; {mn is maximum of m and n}
  write('Enter type to generate ');
  readln(mtype);
  writeln(mtype);
  case mtype of
    0: begin
      sym:=false;
      if m=n then
      begin
        write('Is matrix symmetric? ');  
        readln(inchar);
        writeln(inchar);
        if (inchar='y') or (inchar='Y') then sym:=true else sym:=false;
      end; {ask if symmetric}
      if sym then
      begin
        for i:=1 to n do
        begin
          writeln('Row ',i,' lower triangle elements');
          for j:=1 to i do
          begin
            read(A[i,j]);
            write(A[i,j]:10:5,' ');
            A[j,i]:=A[i,j];
            if (7*(j div 7) = j) and (j<i) then writeln;
          end;
          writeln;
        end;
      end {symmetric matrix}
      else
      begin {not symmetric}
        for i:=1 to m do
        begin
          writeln('Row ',i);
          for j:=1 to n do
          begin
            read(A[i,j]);
            write(A[i,j]:10:5,' ');
          end; {loop on j}
          writeln;
        end; {loop on i}
      end; {else not symmetric}
    end; {case 0 -- input of matrix}
    1: begin {Hilbert}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=1.0/(i+j-1.0);
      if m=n then sym:=true;
    end;
    2: begin {Ding Dong}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=0.5/(1.5+n-i-j);
      if m=n then sym:=true;
    end;
    3: begin {Moler}
      for i:=1 to mn do
      begin
        for j:=1 to i do
        begin
          A[i,j]:=j-2.0;
          A[j,i]:=j-2.0;
        end;
        A[i,i]:=i;
        if m=n then sym:=true;
      end;
    end;
    4: begin {Frank symmetric}
      for i:=1 to mn do
        for j:=1 to i do
        begin
          A[i,j]:=j;
          A[j,i]:=j;
        end;
        if m=n then sym:=true;
    end;
    5: begin {Bordered}
      temp:=2.0;
      for i:=1 to (mn-1) do
      begin
        temp:=temp/2.0; {2^(1-i)}
        for j:=1 to mn do
          A[i,j]:=0.0;
        A[i,mn]:=temp;
        A[mn,i]:=temp;
        A[i,i]:=1.0;
      end;
      A[mn,mn]:=1.0;
      if m=n then sym:=true;
    end;
    6: begin {Diagonal}
      for i:=1 to mn do
      begin
        for j:=1 to mn do
          A[i,j]:=0.0;
        A[i,i]:=i;
      end;
      if m=n then sym:=true;
    end;
    7: begin {W+}
      k:=mn div 2; {[n/2]}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=0.0;
      if m=n then sym:=true;
      for i:=1 to k do
      begin
        A[i,i]:=k+1-i;
        A[mn-i+1,mn-i+1]:=k+1-i;
      end;
      for i:=1 to mn-1 do
      begin
        A[i,i+1]:=1.0;
        A[i+1,i]:=1.0;
      end;
    end;
    8: begin {W-}
      k:=mn div 2; {[n/2]}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=0.0;
      if m=n then sym:=true;
      for i:=1 to k do
      begin
        A[i,i]:=k+1-i;
        A[mn-i+1,mn-i+1]:=i-1-k;
      end;
      for i:=1 to mn-1 do
      begin
        A[i,i+1]:=1.0;
        A[i+1,i]:=1.0;
      end;
      if m=n then sym:=true;
    end;
    9: begin {Constant}
      write('Set all elements to a constant value = ');
      readln(temp);
      writeln(temp);
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=temp;
      if m=n then sym:=true;
    end;
    10: begin {Unit}
      for i:=1 to mn do
      begin
        for j:=1 to mn do A[i,j]:=0.0;
        A[i,i]:=1.0;
      end;
      if m=n then sym:=true;
    end;
    else {case statement else} {!!!! Note missing close bracket here}
    begin
      writeln;
      writeln('*** ERROR *** unrecognized option');
      halt;
    end; {else of case statement}
  end; {case statement}
  if sym then
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
  end;
end; {matrixin}

procedure vectorin(n: integer; var Y: rvector);
{vectorin.pas
  == enter or generate a vector of n real elements
}
var
  i, j, k, m, nt : integer;
  x : real;

begin
  write('vectorin.pas');
  writeln('   -- enter or generate a real vector of ',n,' elements');
  writeln('Options:');
  writeln('   1) constant');
  writeln('   2) uniform random in [0,user_value) ');
  writeln('   3) user entered from console ');
  writeln('   4) entered from RHS columns in matrix file ');
  write('   Choose option :');
  readln(i);
  writeln(i);
  Case i of
    1 : begin
        write('Enter constant value ='); 
        readln(x);
        writeln(x);
        for j:=1 to n do Y[j]:=x;
    end;
    2 : begin
        write('Enter the upper bound to the generator =');
        readln(x);
        writeln(x);
        for j:=1 to n do Y[j]:=Random;
        {According to the Turbo Pascal manual, version 3.0, Random
          returns a number in [0,1). My experience is that most such
          pseudo-random number generators leave a lot to be desired
          in terms of statistical properties. I do NOT recommend it
          for serious use in Monte Carlo calculations or other situations
          where a quality generator is required. For a better generator
          in Pascal, see Wichman B. and Hill, D., (1987)}
    end;
    3 : begin
        writeln('Enter elements of vector one by one');
        for j:=1 to n do
        begin
          write('Y[',j,']=');
          readln(Y[j]);
          writeln(Y[j]);
        end;
    end;
  end {case};
end {vectorin.pas};

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

var
  A : rmatrix;
  avector : smatvec;
  i, j, k, nCol, nRow : integer;
  sym : boolean;
  Y, Ycopy : rvector; {to store the right hand side of the equations}
  singmat : boolean; {set true if matrix discovered to be computationally
              singular during alg07.pas}
  s : real; {an accumulator}

begin
  banner:='dr0708 -- Choleski decomposition and back-substitution';
  write('order of problem = ');  
  readln(nRow);
  writeln(nRow);
  nCol:=nRow; {use symmetric matrix in Choleski}
  Matrixin(nRow,nCol,A,avector,sym);
  writeln;
  writeln('returned matrix of order ',nRow);
  if not sym then halt; {must have symmetric matrix}
  begin
    writeln('Symmetric matrix -- Vector form');
    k:=0;
    for i:=1 to nRow do
    begin
      for j:=1 to i do
      begin
         k:=k+1;
        write(avector[k]:10:5,' ');
        if (7 * (j div 7) = j) and (j<i) then writeln;
      end;
      writeln;
    end;
  end;
  writeln('Enter right hand side of equations');
  vectorin(nRow, Y);
  for i:=1 to nRow do Ycopy[i]:=Y[i];
  writeln;
  choldcmp(nRow,avector, singmat); {decompose matrix}
  begin
    writeln('Decomposed matrix -- Vector form');
    k:=0;
    for i:=1 to nRow do
    begin
      for j:=1 to i do
      begin
        k:=k+1;
        write(avector[k]:10:5,' ');
        if (7 * (j div 7) = j) and (j<i) then writeln;
      end;
      writeln;
    end;
  end;
  if not singmat then
  begin
    Cholback(nRow,avector,Y);
    writeln('Solution');
    for i:=1 to nRow do
    begin
      write(Y[i]:10:5,' ');
      if (7 * (i div 7) = i) and (i<nRow) then writeln;
      writeln;
    end;
    s:=resids(nRow,nCol,A,Ycopy,Y,true);
  end {non-singular case}
  else
  begin
    writeln('Matrix computationally singular -- solution not possible');
  end;
end. {dr0708.pas}
