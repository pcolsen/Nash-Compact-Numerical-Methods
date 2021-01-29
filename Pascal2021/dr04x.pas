Program runsvd(input,output);
{dr04.pas == driver for singular value decomposition and least squares
solution using row-wise entry of data
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

  The following variables allow us to keep a copy of all screen
  information in a file for some of the codes.  Pascal requires a
  variable (confile in this case) for the file itself.  The string
  variable confname is used for the name of the file.  Similar variables
  allow problem data to be read from the file dfile named dfname.
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
procedure Getobsn( n: integer;
            nRHS: integer;
            var W : rmatrix;
            k: integer; {row into which data is loaded}
            endflag: real;
            var enddata: boolean); {set true when end of data reached}

{Getobsn.pas
  -- reads one row of data for a least squares problem from console
    into row k of a working array W. The boolean variable enddata is
    set true when an end of date flag is encountered.
    This procedure assumes that the order of the problem, number of
    right hand sides, etc. are already defined, and that the datafile
    is already open.
}
var
  i : integer;

begin
  begin
    for i:=1 to (n+nRHS) do read(W[k,i]);
    if W[k,1]=endflag then enddata:=true else enddata:=false;
  end;
end; {Getobsn.pas}
procedure GivSVD( n : integer;
                nRHS: integer;
                var B: rmatrix;
                var rss: rvector;
                var svs: rvector;
                var W: rmatrix;
                var nobs : integer);


var
  count, EstRowRank, i, j, k, m, slimit, sweep, tcol : integer;
  bb, c, e2, eps, p, q, r, s, tol, trss, vt : real;
  enddata : boolean;
  endflag : real;

procedure rotnsub;
var
  i: integer;
  r: real;

begin
  for i := m to tcol do
  begin
    r := W[j,i];
    W[j,i] := r*c+s*W[k,i];
    W[k,i] := -r*s+c*W[k,i];
  end;
end;

begin
  writeln('alg04.pas -- Givens',chr(39),
         ' reduction, svd, and least squares solution');
  write('Order of ls problem and no. of right hand sides: ');
  readln(n,nRHS);
  writeln(n,' ',nRHS);
  write('Enter a number to indicate end of data ');
  readln(endflag);
  writeln(endflag);
  tcol := n+nRHS;
  k := n+1;
  for i := 1 to n do
    for j := 1 to tcol do
      W[i,j] := 0.0;
  for i := 1 to nRHS do rss[i] := 0.0;

  eps := calceps;
  tol := n*n*eps*eps;
  nobs := 0;

  enddata := false;
  while (not enddata) do
  begin
    getobsn( n, nRHS, W, k, endflag, enddata);
    if (not enddata) then
    begin
      nobs := nobs+1;
      write('Obsn ',nobs,' ');
      for j := 1 to (n+nRHS) do
      begin
        write(W[k,j]:10:5,' ');
        if (7 * (j div 7) = j) and (j<n+nRHS) then writeln;
      end;
      writeln;
      for j := 1 to n do
      begin
        m := j; s := W[k,j]; c := W[j,j];
        bb := abs(c); if abs(s)>bb then bb := abs(s);
        if bb>0.0 then
        begin
          c := c/bb; s := s/bb; p := sqrt(c*c+s*s);
          s := s/p;
          if abs(s)>=tol then
          begin
            c := c/p;
            rotnsub;
          end;
        end;
      end;

      write('    Uncorrelated residual(s):');
      for j := 1 to nRHS do
      begin
        rss[j] := rss[j]+sqr(W[k,n+j]);
        write(W[k,n+j]:10,' ');
        if (7 * (j div 7) = j) and (j < nRHS) then writeln;
      end;
      writeln; 

    end;
  end;


  m := 1;
  slimit := n div 4;  if slimit<6 then slimit := 6;

  sweep := 0;
  e2 := 10.0*n*eps*eps;
  tol := eps*0.1;
  EstRowRank := n; ;
  repeat
    count := 0;
    for j := 1 to (EstRowRank-1) do
    begin
      for k := (j+1) to EstRowRank do
      begin
        p := 0.0; q := 0.0; r := 0.0;
        for i := 1 to n do
        begin
          p := p+W[j,i]*W[k,i]; q := q+sqr(W[j,i]); r := r+sqr(W[k,i]);
        end;
        svs[j] := q; svs[k] := r;
        IF q >= r then
        begin
          if not ((q<=e2*svs[1]) or (abs(p)<=tol*q)) then
          begin
            p := p/q; r := 1-r/q; vt := sqrt(4*p*p + r*r);
            c := sqrt(0.5*(1+r/vt)); s := p/(vt*c);
            rotnsub;
            count := count+1;
          end;
        end
        else

        begin
          p := p/r; q := q/r-1; vt := sqrt(4*p*p + q*q);
          s := sqrt(0.5*(1-q/vt));
          if p<0 then s := -s;
          c := p/(vt*s);
          rotnsub;
          count := count+1;
        end;


      end;
    end;
    sweep := sweep +1;
    writeln('Sweep ',sweep,' ',count,' rotations performed');
     while (EstRowRank >= 3) and (svs[EstRowRank] <= svs[1]*tol+tol*tol)
          do EstRowRank := EstRowRank-1;
  until (sweep>slimit) or (count=0);

  writeln('Singular values and principal components');
  for j := 1 to n do
  begin
    s := svs[j];
    s := sqrt(s); svs[j] := s;
    writeln('Singular value [',j,']= ',s);
    if s>=tol then
    begin
      for i := 1 to n do W[j,i] := W[j,i]/s;
      for i := 1 to n do
      begin
        write(W[j,i]:8:5,' '); 
        if (8 * (i div 8) = i) and (i<n) then writeln;
      end;

      writeln;
    end;

  end;

  q := 0.0;
  while q>=0.0 do
  begin
    write('Enter a tolerance for zero (<0 to exit) ');
    readln(q);
    writeln(q);
    if q>=0.0 then
    begin

      for i := 1 to nRHS do
      begin
        trss := rss[i];
        for j := 1 to n do
        begin
          p := 0.0;
          for k := 1 to n do
          begin
            if svs[k]>q then p := p+W[k,j]*W[k,n+i]/svs[k];
          end;
          B[j,i] := p;
          writeln('Solution component [',j,']= ',p);
          if svs[j]<=q then trss := trss+sqr(W[j,n+i]);
        end;
        writeln('Residual sum of squares =',trss);
      end;
    end;
  end;
end;
{main program}
var
i, j, k, m, n, nRHS : integer; {order of problem, number of right hand sides}
Vtranspose : rmatrix;
rssvec : rvector;  {to hold the residual sum of squares}
svs : rvector; {to hold the singular values}
B : rmatrix; {the matrix of least squares solutions }
A : rmatrix; {needed only for svdtst}
fail, sym : boolean;
t1: real;
begin
banner:='dr04.pas -- run Algorithm 4 problems -- Givens'+chr(39)+
' reduction,';
GivSVD( n, nRHS, B, rssvec, svs, Vtranspose, m);
end. {dr04.pas == Givens' reduction, svd and least squares soln}
