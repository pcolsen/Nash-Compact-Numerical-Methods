program dr26(input,output);

{dr26.pas == eigensolutions of a complex matrix by Eberlein's
          complex Jacobi procedure

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

procedure comeig( n : integer;
          var itcount: integer;
          var A, Z, T, U : rmatrix);

var
  Rvec : rvector;

  i, itlimit, j, k, k1, m, n1 : integer;
  aki, ami, bv, br, bi : real;
  c, c1i, c1r, c2i, c2r, ca, cb, ch, cos2a, cot2x, cotx, cx : real;
  d, de, di, diag, dr, e, ei, er, eps, eta, g, hi, hj, hr : real;
  isw, max, nc, nd, root1, root2, root : real;
  s, s1i, s1r, s2i, s2r, sa, sb, sh, sig, sin2a, sx : real;
  tanh, tau, te, tee, tem, tep ,tse, zki, zmi : real;
  mark : boolean;

begin

  writeln('alg26.pas -- comeig');
  eps := Calceps;
  mark := false; n1 := n-1;
  for i := 1 to n do
  begin
    for j := 1 to n do
    begin
      T[i,j] := 0.0; U[i,j] := 0.0; if i=j then T[i,i] := 1.0;
    end;
  end;
  itlimit := itcount;
  itcount := 0;
  while (itcount<=itlimit) and (not mark) do
  begin
    itcount := itcount+1;
    tau := 0.0;
    diag := 0.0;
    for k := 1 to n do
    begin
      tem := 0;
      for i := 1 to n do if i<>k then tem := tem+ABS(A[i,k])+ABS(Z[i,k]);
      tau := tau+tem; tep := abs(A[k,k])+abs(Z[k,k]);
      diag := diag+tep;
      Rvec[k] := tem+tep;
    end;
    writeln('TAU=',tau,'  AT ITN ',itcount);
    for k := 1 to n1 do
    begin
      max := Rvec[k]; i := k; k1 := k+1;
      for j := k1 to n do
      begin
        if max<Rvec[j] then
        begin
          max := Rvec[j]; i := j;
        end;
      end;
      if i<>k then
      begin
        Rvec[i] := Rvec[k];
        for j := 1 to n do
        begin
          tep := A[k,j]; A[k,j] := A[i,j]; A[i,j] := tep; tep := Z[k,j];
          Z[k,j] := Z[i,j]; Z[i,j] := tep;
        end;
        for j := 1 to n do
        begin
          tep := A[j,k]; A[j,k] := A[j,i]; A[j,i] := tep; tep := Z[j,k];
          Z[j,k] := Z[j,i]; Z[j,i] := tep; tep := T[j,k]; T[j,k] := T[j,i];
          T[j,i] := tep; tep := U[j,k]; U[j,k] := U[j,i]; U[j,i] := tep;
        end;
      end;
    end;
    if tau>=100.0*eps then
    begin
      mark := true;
      for k := 1 to n1 do
      begin
        k1 := k+1;
        for m := k1 to n do
        begin
          hj := 0.0; hr := 0.0; hi := 0.0; g := 0.0;
          for i := 1 to n do
          begin
            if (i<>k) and (i<>m) then
            begin
              hr := hr+A[k,i]*A[m,i]+Z[k,i]*Z[m,i];
              hr := hr-A[i,k]*A[i,m]-Z[i,k]*Z[i,m];
              hi := hi+Z[k,i]*A[m,i]-A[k,i]*Z[m,i];
              hi := hi-A[i,k]*Z[i,m]+Z[i,k]*A[i,m];
              te := A[i,k]*A[i,k]+Z[i,k]*Z[i,k]+A[m,i]*A[m,i]+Z[m,i]*Z[m,i];
              tee := A[i,m]*A[i,m]+Z[i,m]*Z[i,m]+A[k,i]*A[k,i]+Z[k,i]*Z[k,i];
              g := g+te+tee; hj := hj-te+tee;
            end;
          end;
          br := A[k,m]+A[m,k]; bi := Z[k,m]+Z[m,k]; er := A[k,m]-A[m,k];
          ei := Z[k,m]-Z[m,k]; dr := A[k,k]-A[m,m]; di := Z[k,k]-Z[m,m];
          te := br*br+ei*ei+dr*dr; tee := bi*bi+er*er+di*di;
          if te>=tee then
          begin
            isw := 1.0; c := br; s := ei; d := dr; de := di;
            root2 := sqrt(te);
          end
          else
          begin
            isw := -1.0; c := bi; s := -er; d := di; de := dr;
            root2 := sqrt(tee);
          end;
          root1 := sqrt(s*s+c*c); sig := -1.0; if d>=0.0 then sig := 1.0;
          sa := 0.0; ca := -1.0; if c>=0.0 then ca := 1.0;
          if root1<=eps then
          begin
            sx := 0.0; sa := 0.0; cx := 1.0; ca := 1.0;
            if isw<=0.0 then
            begin
              e := ei; bv := -br;
            end
            else
            begin
              e := er; bv := bi;
            end;
            nd := d*d+de*de;
          end
          else
          begin
            if abs(s)>eps then
            begin
              ca := c/root1; sa := s/root1;
            end;
            cot2x := d/root1; cotx := cot2x+(sig*sqrt(1.0+cot2x*cot2x));
            sx := sig/sqrt(1.0+cotx*cotx); cx := sx*cotx;

            eta := (er*br+ei*bi)/root1; tse := (br*bi-er*ei)/root1;
            te := sig*(tse*d-de*root1)/root2; tee := (d*de+root1*tse)/root2;
            nd := root2*root2+tee*tee; tee := hj*cx*sx; cos2a := ca*ca-sa*sa;
            sin2a := 2.0*ca*sa; tem := hr*cos2a+hi*sin2a;
            tep := hi*cos2a-hr*sin2a; hr := hr*cx*cx-tem*sx*sx-ca*tee;
            hi := hi*cx*cx+tep*sx*sx-sa*tee;
            bv := isw*te*ca+eta*sa; e := ca*eta-isw*te*sa;
          end;

          s := hr-sig*root2*e; c := hi-sig*root2*bv; root := sqrt(c*c+s*s);
          if root<eps then
          begin
            cb := 1.0; ch := 1.0; sb := 0.0; sh := 0.0;
          end
          else
          begin
            cb := -c/root; sb := s/root; tee := cb*bv-e*sb; nc := tee*tee;
            tanh := root/(g+2.0*(nc+nd)); ch := 1.0/sqrt(1.0-tanh*tanh);
            sh := ch*tanh;
          end;
          tem := sx*sh*(sa*cb-sb*ca); c1r := cx*ch-tem; c2r := cx*ch+tem;
          c1i := -sx*sh*(ca*cb+sa*sb); c2i := c1i; tep := sx*ch*ca;
          tem := cx*sh*sb; s1r := tep-tem; s2r := -tep-tem; tep := sx*ch*sa;
          tem := cx*sh*cb; s1i := tep+tem; s2i := tep-tem;
          tem := sqrt(s1r*s1r+s1i*s1i); tep := sqrt(s2r*s2r+s2i*s2i);
          if tep>eps then mark := false;
          if (tep>eps) and (tem>eps) then
          begin
            for i := 1 to n do
            begin
              aki := A[k,i]; ami := A[m,i]; zki := Z[k,i]; zmi := Z[m,i];
              A[k,i] := c1r*aki-c1i*zki+s1r*ami-s1i*zmi;
              Z[k,i] := c1r*zki+c1i*aki+s1r*zmi+s1i*ami;
              A[m,i] := s2r*aki-s2i*zki+c2r*ami-c2i*zmi;
              Z[m,i] := s2r*zki+s2i*aki+c2r*zmi+c2i*ami;
            end;
            for i := 1 to n do
            begin
              aki := A[i,k]; ami := A[i,m]; zki := Z[i,k]; zmi := Z[i,m];
              A[i,k] := c2r*aki-c2i*zki-s2r*ami+s2i*zmi;
              Z[i,k] := c2r*zki+c2i*aki-s2r*zmi-s2i*ami;
              A[i,m] := -s1r*aki+s1i*zki+c1r*ami-c1i*zmi;
              Z[i,m] := -s1r*zki-s1i*aki+c1r*zmi+c1i*ami;
              aki := T[i,k]; ami := T[i,m]; zki := U[i,k]; zmi := U[i,m];
              T[i,k] := c2r*aki-c2i*zki-s2r*ami+s2i*zmi;
              U[i,k] := c2r*zki+c2i*aki-s2r*zmi-s2i*ami;
              T[i,m] := -s1r*aki+s1i*zki+c1r*ami-c1i*zmi;
              U[i,m] := -s1r*zki-s1i*aki+c1r*zmi+c1i*ami;
            end;
          end;
        end;
      end;
    end
    else mark := true;
  end;
  if itcount>itlimit then itcount := -itcount;
end;

procedure stdceigv(n: integer;
                var T, U: rmatrix);

var
  i, k, m : integer;
  b, e, g, s : real;

begin
  writeln('alg11.pas -- standardized eigensolutions');
  for i := 1 to n do
  begin
    g := T[1,i]*T[1,i]+U[1,i]*U[1,i];
    k := 1;
    if n>1 then
    begin
      for m := 2 to n do
      begin
        b := T[m,i]*T[m,i]+U[m,i]*U[m,i];
        if b>g then
        begin
          k := m;
          g := b;
        end;
      end;
    end;
    e := T[k,i]/g;
    s := -U[k,i]/g;
    for k := 1 to n do
    begin
      g := T[k,i]*e-U[k,i]*s; U[k,i] := U[k,i]*e+T[k,i]*s; T[k,i] := g;
    end;
  end;
end;

procedure comres( i, n: integer;
                  A, Z, T, U, Acopy, Zcopy : rmatrix);
var
  j, k: integer;
  g, s, ss : real;

begin
  writeln('alg12.pas -- complex eigensolution residuals');
  ss := 0.0;
  for j := 1 to n do
  begin
    s := -A[i,i]*T[j,i]+Z[i,i]*U[j,i]; g := -Z[i,i]*T[j,i]-A[i,i]*U[j,i];
    for k := 1 to n do
    begin
      s := s+Acopy[j,k]*T[k,i]-Zcopy[j,k]*U[k,i];
      g := g+Acopy[j,k]*U[k,i]+Zcopy[j,k]*T[k,i];
    end;
    writeln('(',s,',',g,')'); 
    ss := ss+s*s+g*g;
  end;
  writeln('Sum of squares = ',ss); 
end;

{Main program}
var
  A, Z, Acopy, Zcopy, T, U : rmatrix;
  i, it, j, k, n : integer;
  sym : boolean;
  avec : smatvec; {for compatibility of Matrixin only}

begin
  banner:='dr26.pas -- Eigensolutions of a general complex matrix';
  write(' Order of matrix = '); readln(n);
  writeln(n);
  writeln('Provide real part of matrix (A)');
  matrixin(n,n,A,avec,sym);
  writeln('Provide imaginary part of matrix (Z)');
  matrixin(n,n,Z,avec,sym);
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
    Acopy[i,j]:=A[i,j]; Zcopy[i,j]:=Z[i,j];
    write('(',A[i,j]:10:5,',',Z[i,j]:10:5,') ');
    if (3 * (j div 3) = j) and (j<n) then writeln;
    end; {copy loop j}
    writeln;
  end; {copy loop i}
  it:=50; {allow a maximum of 50 iterations}
  comeig( n, it, A, Z, T, U);
  if it>0 then writeln('Converged in ',it,' iterations')
    else writeln('Not converged after ',it,' iterations');
  stdceigv(n, T, U); {standardize the eigensolutions -- alg11.pas}
  for i:=1 to n do
  begin
    writeln('EIGENVALUE ',i,'=(',A[i,i],',',Z[i,i],')');
    writeln('VECTOR');
    for k:=1 to n do
    begin
    writeln('(',T[k,i],',',U[k,i],')');
    end; { loop on  k}
    comres( i, n, A, Z, T, U, Acopy, Zcopy); {residuals -- alg12.pas}
  end; {loop on i}
end. {dr26.pas == eigensolutions of a complex matrix}
