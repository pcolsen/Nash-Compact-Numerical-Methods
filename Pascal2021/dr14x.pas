Program dr14(input, output);
{dr14.pas == driver for Jacobi method (Alg14) for eigensolutions of a real
symmetric matrix
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
        write('Is matrix symmetric? ');  readln(inchar);
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






Procedure evJacobi(n: integer;
               var A,V : rmatrix;
                   initev: boolean);


var
  count, i, j, k, limit, skipped : integer;
  c, p, q, s, t : real;
  oki, okj, rotn : boolean;

begin
  writeln('alg14.pas -- eigensolutions of a real symmetric');
  writeln('             matrix via a Jacobi method');
  if initev then
  begin
    for i := 1 to n do
    begin
      for j := 1 to n do V[i,j] := 0.0;
      V[i,i] := 1.0;
    end;
  end;
  count := 0;
  limit := 30;
  skipped := 0;

  while (count<=limit) and (skipped<((n*(n-1)) div 2) ) do

  begin
    count := count+1;
    write('sweep ',count,' '); 
    skipped := 0;
    for i := 1 to (n-1) do
    begin
      for j := (i+1) to n do
      begin
        rotn := true;
        p := 0.5*(A[i,j]+A[j,i]);
        q := A[i,i]-A[j,j];
        t := sqrt(4.0*p*p+q*q);
        if t=0.0 then
        begin
          rotn := false;
        end
        else
        begin
          if q>=0.0 then
          begin

            oki := (abs(A[i,i])=abs(A[i,i])+100.0*abs(p));
            okj := (abs(A[j,j])=abs(A[j,j])+100.0*abs(p));
            if oki and okj then rotn := false
                           else rotn := true;

            if rotn then
            begin
               c := sqrt((t+q)/(2.0*t)); s := p/(t*c);
            end;
          end
          else
          begin
            rotn := true;
            s := sqrt((t-q)/(2.0*t));
            if p<0.0 then s := -s;
            c := p/(t*s);
          end;
          if 1.0=(1.0+abs(s)) then rotn := false;
        end;
        if rotn then
        begin
          for k := 1 to n do
          begin
            q := A[i,k]; A[i,k] := c*q+s*A[j,k]; A[j,k] := -s*q+c*A[j,k];
          end;

          for k := 1 to n do
          begin
            q := A[k,i]; A[k,i] := c*q+s*A[k,j]; A[k,j] := -s*q+c*A[k,j];

            q := V[k,i]; V[k,i] := c*q+s*V[k,j]; V[k,j] := -s*q+c*V[k,j];
          end;
        end
        else

           skipped := skipped+1;
      end;
    end;
    writeln('  ',skipped,' / ',n*(n-1) div 2,'  rotations skipped');
  end;
end;
{main program}
var
i, j, nRow, nCol : integer;
A, V, ACOPY : rmatrix;
Bvec, Y : rvector; {to test residuals}
t1: real;
tvec : smatvec; {needed only for Matrixin}
initev, sym : boolean;
begin
banner:='dr14.pas -- driver for Jacobi eigensolution method';

nRow := 1; {To get loop started}

while (nRow > 0) do
begin
  write('Order of problem (n) = ');
  readln(nRow); 
  writeln(nRow);
  if (nRow <= 0) then halt;
  nCol:=nRow;
  Matrixin(nRow, nCol, A, tvec, sym);
  if not sym then halt; {program only designed for symmetric matrices}
  for j:=1 to nRow do
  begin
    for i:=1 to nRow do
    begin
      write(A[i, j]:10:5, ' ');
      ACOPY[i, j]:=A[i, j];
      if (7 * (i div 7) = i) and (i<nRow) then writeln;
    end;
    writeln;
  end;
  initev:=true; {Here we want to get the eigenvectors of A, not some
  generalized problem.}
  evJacobi( nRow, A, V, initev);
  for j:=1 to nRow do
  begin
    t1:=A[j, j];
    writeln('Eigenvalue ', j, ' = ', t1);
    for i:=1 to nRow do
    begin
      write(V[i, j]:10:7, ' ');
      if (i = 7 * (i div 7)) and (i<nRow) then writeln;
      Bvec[i]:=V[i, j]; {to initialize residual test}
      Y[i]:=t1*Bvec[i];
    end;
    writeln;
    t1 := resids(nRow, nCol, ACOPY, Y, Bvec, true);
    writeln;
  end; {loop on solutions j}
end; {while}
end. {dr14.pas}
