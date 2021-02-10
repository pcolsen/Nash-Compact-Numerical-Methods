program dr24le(input, output);
{dr24le.pas == linear equations by conjugate gradients
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
    writeln(confile);
    writeln(confile,'******* Matrix dimensions zero *********');
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
  readln(infile,mtype);
  if infname<>'con' then writeln(mtype);
  writeln(confile,'Enter type to generate ',mtype);
  case mtype of
    0: begin
      sym:=false;
      if m=n then
      begin
        write('Is matrix symmetric? ');  readln(infile,inchar);
        if infname<>'con' then writeln(inchar);
        writeln(confile,'Is matrix symmetric? ',inchar);
        if (inchar='y') or (inchar='Y') then sym:=true else sym:=false;
      end; {ask if symmetric}
      if sym then
      begin
        for i:=1 to n do
        begin
          writeln('Row ',i,' lower triangle elements');
          writeln(confile,'Row ',i,' lower triangle elements');
          for j:=1 to i do
          begin
            read(infile,A[i,j]);
            if infname<>'con' then write(A[i,j]:10:5,' ') else write(' ');
            write(confile,A[i,j]:10:5,' ');
            A[j,i]:=A[i,j];
            if (7*(j div 7) = j) and (j<i) then
            begin
              writeln;
              writeln(confile);
            end;
          end;
          writeln;
          writeln(confile);
        end;
      end {symmetric matrix}
      else
      begin {not symmetric}
        for i:=1 to m do
        begin
          writeln('Row ',i);
          writeln(confile,'Row ',i);
          for j:=1 to n do
          begin
            read(infile,A[i,j]);
            if infname<>'con' then write(A[i,j]:10:5,' ') else write(' ');
            write(confile,A[i,j]:10:5,' ');
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
      readln(infile,temp);
      if infname<>'con' then writeln(temp);
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
      writeln(confile);
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
  readln(infile,i);
  if infname<>'con' then writeln(i);
  write(confile,'vectorin.pas');
  writeln(confile,'   -- enter or generate a real vector of ',n,' elements');
  writeln(confile,'Options:');
  writeln(confile,'   1) constant');
  writeln(confile,'   2) uniform random in [0,user_value) ');
  writeln(confile,'   3) user entered ');
  writeln(confile,'   4) entered from RHS columns in matrix file ');
  writeln(confile,'   Choose option :',i);
  Case i of
    1 : begin
        write('Enter constant value ='); readln(infile,x);
        if infname<>'con' then writeln(x);
        writeln(confile,'Enter constant value =',x);
        for j:=1 to n do Y[j]:=x;
    end;
    2 : begin
        write('Enter the upper bound to the generator =');
        readln(infile,x);
        if infname<>'con' then writeln(x);
        writeln(confile,'Enter the upper bound to the generator =',x);
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
        writeln(confile,'Enter elements of vector one by one');
        for j:=1 to n do
        begin
          write('Y[',j,']=');
          readln(infile,Y[j]);
          if infname<>'con' then writeln(Y[j]);
          writeln(confile,'Y[',j,']=',Y[j]);
        end;
    end;
    4 : begin  {Get elements of RHS from a matrix + vectors file}
        write('Datafile ');
        readln(infile,dfname);
        if infname<>'con' then writeln(dfname);
        writeln(confile,'Datafile ',dfname);
        write('Which RHS vector should be retrieved? ');
        readln(j);
        if infname<>'con' then writeln(j);
        writeln(confile,'Which RHS vector should be retrieved? ',j);
        {We will use the least squares data file a a form of input}
        if length(dfname)>0 then
        begin
          assign(dfile, dfname);
          reset(dfile);
          read(dfile, nt, i);{reading i rhs elements}
          writeln('Number of columns =',nt);
          writeln(confile,'Number of columns =',nt);
          m:=0; {to count the number of rows}
          while (not eof(dfile)) do
          begin
            m:=m+1;
            for k:=1 to nt do read(dfile,x); {ignore coefficient matrix}
            for k:=1 to j do
            begin
              read(dfile,x); if k=j then Y[m]:=x;
            end;
          end; {while}
          close(dfile);
          writeln('Found ',m,' elements in vector');
          writeln(confile,'Found ',m,' elements in vector');
          if m<>n then
          begin
            writeln('*** ERROR *** not in agreement with procedure call');
            writeln(confile,
                '*** ERROR *** not in agreement with procedure call');
            halt;
          end;
        end {if length(dfname)}
      end; {case 4}
  end {case};
end {vectorin.pas};


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
  count, i, itn, itlimit : integer;
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
      writeln(confile,'Step too large -- coefficient matrix indefinite?');
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
    writeln(confile,'Residuals');
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
      write(confile,t1:10,' ');
      if (i = 7 * (i div 7)) and (i<nRow) then writeln;
      if (i = 7 * (i div 7)) and (i<nRow) then
      writeln(confile);
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
avec : smatvec; {for matrixin only}
sym : boolean; {to tell if matrix symmetric}
ch : char;
i, j, n, itcount : integer;
ssmin, t, s : real;
begin
banner:='dr24le.pas -- linear equations by conjugate gradients';
write('Order of problem = '); { 10 }
readln(n);
writeln(n);
writeln('Coefficient matrix');
matrixin(n, n, A, avec, sym);  {4 == Frank}
if not sym then halt;
writeln('RHS vector');
vectorin(n, Y);
writeln('Initial guess for solution');
vectorin(n, Bvec);
itcount:=10*n; {safety setting}
lecg( n, A, Y, Bvec, itcount, ssmin);
writeln('Solution after ',itcount,' iterations. Est. sumsquares ',ssmin);
for i:=1 to n do
begin
write(Bvec[i]:10:5,' ');
if (7 * (i div 7) = i) and (i<n) then
begin
writeln;
end;
end;
writeln;
s:=resids(n, n, A, Y, Bvec, true);
end. {dr24le.pas}
