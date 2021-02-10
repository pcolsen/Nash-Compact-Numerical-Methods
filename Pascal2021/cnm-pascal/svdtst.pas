Procedure svdtst( A, U, V: rmatrix; Z: rvector;
              nRow, UCol, VCol: integer);
{svdtst.pas
  == This routine tests the results of a singular value
  decomposion calculation.  The matrix A is presumed to contain
  the matrix of which the purported decomposition is

        U  Z  V-transpose

  This routine tests column orthogonality of U and V,
  row orthogonality of V, and the reconstruction suggested
  by the decomposition. It does not carry out the tests of
  the Moore-Penrose inverse A+, which can be computed as

      A+ :=  V  Z  U-transpose.

  FORTRAN codes for the conditions

          A+ A A+  = ? =  A+
          A  A+ A  = ? =  A
          (A+  A)-transpose  = ? =  A+ A
          (A  A+)-transpose  = ? =  A  A+

  are given in Nash, J.C. and Wang, R.L.C. (1986)
}

var
  i,j,k:integer;
  t1: real;
  imax, jmax: integer;
  valmax: real;

begin
  writeln('Column orthogonality of U');
  writeln(confile,'Column orthogonality of U');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to UCol do
  begin
    for j:=i to UCol do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1;
      for k:=1 to nRow do t1:=t1+U[k,i]*U[k,j];
        {writeln('Inner product ',i,',',j,' = ',t1);}
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,' = ',valmax);
  writeln(confile,'Largest inner product is ',imax,',',jmax,' = ',valmax);
  write('Hit [cr] to continue ');readln(infile);
  writeln(confile,'Hit [cr] to continue ');
  writeln('Row orthogonality of U (NOT guaranteed in svd)');
  writeln(confile,'Row orthogonality of U (NOT guaranteed in svd)');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to nRow do
  begin
    for j:=i to nRow do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1;
      for k:=1 to UCol do t1:=t1+U[i,k]*U[j,k];
      {writeln('Inner product ',i,',',j,' = ',t1);}
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,' = ',valmax);
  writeln(confile,'Largest inner product is ',imax,',',jmax,' = ',valmax);
  write('Hit [cr] to continue ');readln(infile);
  writeln('Column orthogonality of V');
  writeln(confile,'Hit [cr] to continue ');
  writeln(confile,'Column orthogonality of V');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to VCol do
  begin
    for j:=i to VCol do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1.0;
      for k:=1 to VCol do t1:=t1+V[k,i]*V[k,j];
      {writeln('Inner product ',i,',',j,' = ',t1); }
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,' = ',valmax);
  write('Hit [cr] to continue ');readln(infile);
  writeln('Row orthogonality of V');
  writeln(confile,'Largest inner product is ',imax,',',jmax,' = ',valmax);
  writeln(confile,'Hit [cr] to continue ');
  writeln(confile,'Row orthogonality of V');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to VCol do
  begin
    for j:=i to VCol do
    begin
      t1:=0.0; {accumulate inner products}
      if i=j then t1:=-1;
      for k:=1 to VCol do t1:=t1+V[i,k]*V[j,k];
      {writeln('Inner product ',i,',',j,' = ',t1);}
      if abs(t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=t1;
      end;
    end;
  end;
  writeln('Largest inner product is ',imax,',',jmax,' = ',valmax);
  write('Hit [cr] to continue ');readln(infile);
  writeln('Reconstruction of initial matrix');
  writeln(confile,'Largest inner product is ',imax,',',jmax,' = ',valmax);
  writeln(confile,'Hit [cr] to continue ');
  writeln(confile,'Reconstruction of initial matrix');
  valmax:=0.0;
  imax:=0;
  jmax:=0;
  for i:=1 to nRow do
  begin
    for j:=1 to VCol do
    begin
      t1:=0;
      for k:=1 to VCol do
        t1:=t1+U[i,k]*Z[k]*V[j,k]; { U * S * V-transpose}
      {writeln('A[',i,',',j,']=',A[i,j],' Recon. =',t1,'  error=',A[i,j]-t1);}
      if abs(A[i,j]-t1)>abs(valmax) then
      begin
        imax:=i; jmax:=j; valmax:=A[i,j]-t1;
      end;
    end;
  end;
  writeln('Largest error is ',imax,',',jmax,' = ',valmax);
  writeln(confile,'Largest error is ',imax,',',jmax,' = ',valmax);
end; {svdtst.pas}
