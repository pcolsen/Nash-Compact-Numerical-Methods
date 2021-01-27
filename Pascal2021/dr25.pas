program dr25(input, output);
{dr25.pas == eigensolutions by minimisation of the Rayleigh quotient

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

{$I matrixin.pas}
{$I vectorin.pas}

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


procedure rqmcg( n : integer;
            A, B : rmatrix;
          var X : rvector;
          var ipr : integer;
          var rq : real);

var
  count, i, itn, itlimit : integer;
  avec, bvec, yvec, zvec, g, t : rvector;
  beta, d, eps, g2, gg, oldg2, pa, pn, s2, step : real;
  t2, ta, tabt, tat, tb, tbt, tol, u, v, w, xat, xax, xbt, xbx : real;
  conv, fail : boolean;

begin
  writeln('alg25.pas -- Rayleigh quotient minimisation');
  itlimit := ipr;
  fail := false;
  conv := false;
  ipr := 0;
  eps := calceps;
  tol := n*n*eps*eps;

  pa := big;
  while (ipr<=itlimit) and (not conv) do
  begin
    matmul(n, A, X, avec);
    matmul(n, B, X, bvec);
    ipr := ipr+1;

    xax := 0.0; xbx := 0.0;
    for i := 1 to n do
    begin
      xax := xax+X[i]*avec[i]; xbx := xbx+X[i]*bvec[i];
    end;
    if xbx<=tol then halt;
    rq := xax/xbx;
    write(ipr,' products -- ev approx. =',rq:18);
    if rq<pa then
    begin
      pa := rq;
      gg := 0.0;
      for i := 1 to n do
      begin
        g[i] := 2.0*(avec[i]-rq*bvec[i])/xbx; gg := gg+g[i]*g[i];
      end;
      writeln(' squared gradient norm =',gg:8);
      if gg>tol then

      begin

        for i := 1 to n do t[i] := -g[i];
        itn := 0;
        repeat
          itn := itn+1;
          matmul(n, A, t, yvec);
          matmul(n, B, t, zvec); ipr := ipr+1;
          tat := 0.0; tbt := 0.0; xat := 0.0; xbt := 0.0;
          for i := 1 to n do
          begin
            xat := xat+X[i]*yvec[i]; tat := tat+t[i]*yvec[i];
            xbt := xbt+X[i]*zvec[i]; tbt := tbt+t[i]*zvec[i];
          end;

          u := tat*xbt-xat*tbt; v := tat*xbx-xax*tbt;
          w := xat*xbx-xax*xbt; d := v*v-4.0*u*w;
          if d<0.0 then halt;

          d := sqrt(d);
          if v>0.0 then step := -2.0*w/(v+d) else step := 0.5*(d-v)/u;

          count := 0;
          xax := 0.0; xbx := 0.0;
          for i := 1 to n do
          begin
            avec[i] := avec[i]+step*yvec[i];
            bvec[i] := bvec[i]+step*zvec[i];
            w := X[i]; X[i] := w+step*t[i];
            if (reltest+w)=(reltest+X[i]) then count := count+1;
            xax := xax+X[i]*avec[i]; xbx := xbx+X[i]*bvec[i];
          end;
          if xbx<=tol then halt
                  else pn := xax/xbx;
          if (count<n) and (pn<rq) then
          begin
            rq := pn; gg := 0.0;
            for i := 1 to n do
            begin
              g[i] := 2.0*(avec[i]-pn*bvec[i])/xbx; gg := gg+g[i]*g[i];
            end;
            if gg>tol then
            begin
              xbt := 0.0; for i := 1 to n do xbt := xbt+X[i]*zvec[i];

              tabt := 0.0; beta := 0.0;
              for i := 1 to n do
              begin
                w := yvec[i]-pn*zvec[i]; tabt := tabt+t[i]*w;
                beta := beta+g[i]*(w-g[i]*xbt);
              end;
              beta := beta/tabt;

              for i := 1 to n do t[i] := beta*t[i]-g[i];
            end;
          end

          else
          begin
            if itn=1 then conv := true;
            itn := n+1;
          end;
        until (itn>=n) or (count=n) or (gg<=tol) or conv;
      end
      else conv := true;
    end
    else
    begin
      conv := true;
    end;
    ta := 0.0;
    for i := 1 to n do ta := ta+sqr(X[i]); ta := 1.0/sqrt(ta);
    for i := 1 to n do X[i] := ta*X[i];
  end;
  if ipr>itlimit then ipr := -ipr;
  writeln;
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


var
  A, B : rmatrix;
  X : rvector; {eigenvector}
  Y : rvector; {for residuals}
  avec : smatvec; {for matrixin only}
  sym : boolean; {to tell if matrix symmetric}
  ch : char;
  i, j, n, itcount : integer;
  ev, t, s : real;

begin
  banner:='dr25.pas -- minimise Rayleigh quotient';
  write('Order of problem =');
  readln(n); writeln(n); 
  writeln('Matrix A');
??
  matrixin(n, n, A, avec, sym);
  if not sym then
  begin
    writeln('Matrix not symmetric -- halting');
    halt;
  end;
??
  writeln('Metric matrix B');
  matrixin(n, n, B, avec, sym);
??  if not sym then
  begin
    writeln('Matrix not symmetric -- halting');
    writeln(confile,'Matrix not symmetric -- halting');
    halt;
  end;
??
  writeln('Initial eigenvector approximation');
??  vectorin(n, X);
  itcount:=100*n; {safety setting}
  rqmcg( n, A, B, X, itcount, ev);
  writeln('Solution after ',itcount,' products. Est. eigenvalue =',ev);
  for i:=1 to n do
  begin
    write(X[i]:10:7,' ');
    if (7 * (i div 7) = i) and (i<n) then  writeln;
    t:=0.0;
    for j:=1 to n do t:=t+B[i,j]*X[j];
    Y[i]:=ev*t; {to save eigenvalue * matrix-vector product for residuals}
  end;
  writeln;
  s := resids(n, n, A, Y, X, true);
end. {dr25.pas}

