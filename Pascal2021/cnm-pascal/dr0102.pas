{$I turbo.cnm}
Program runsvd(input,output);
{dr0102.pas == Calculation of Singular values and vectors of an arbitrary
          real matrix, solution of linear least squares approximation
          problem.

  Modifies a method due to Kaiser.  See Nash and Shlien (1987): Simple
  algorithms for the partial singular value decomposition.  Computer
  Journal, vol.30, pp.268-275.

          Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
          Minor fix for zero rows 1993-07-14
}


{$I constype.def} {definitions of constants and types}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}  {compute machine precision}
{$I resids.pas}   {compute residuals -- used in alg02.pas}
{$I matcopy.pas}  {copy matrix}
{$I psvdres.pas}  {PrtSVDResults procedure}
{$I svdtst.pas}   {test svd results -- needs nice summary}
{$I matrixin.pas} {input or generate a matrix of reals}
{$I vectorin.pas} {input or generate a vector of reals}
{$I alg01.pas}    {NashSVD by column rotations}
{$I alg02.pas}    {svd least squares solution}
{$I startup.pas}    {get control input file and console image file}

{main program}
var
  nRow, nCol : integer;
  A, V, U : rmatrix;
  W : wmatrix; {a working matrix which will contain U Zd in the
    upper nRow rows, and V in the bottom nCol rows, where Zd
    is the diagonal matrix of singular values.  That is, W
    becomes

                  (  U  Zd )
                  (        )
                  (    V   )

  }
  Z : rvector; {Z will contain either the squares of singular
          values or the singular values themselves}
  Y : rvector; {Y will contain the 'right hand side' of the
          least squares problem, i.e.  the vector to be
          approximated }
  Bvec : rvector; {the least squares solution }
  fail: boolean;
  sym : boolean;
  inchar : char;
  i,j,k, imax, jmax : integer;
  t1, t2: real;
  tvec : smatvec; {needed only for Matrixin}

begin
  banner:='dr0102.pas -- driver for svd and least squares solution';
  startup; {perform startup operations}
  repeat
    Write('Give dimensions of matrix (no. rows  no. cols.):');
    readln(infile,nRow, nCol); if infname<>'con' then writeln(nRow,' ',nCol);
    Writeln(confile,'Give dimensions of matrix (no. rows  no. cols.):',
              nRow,' ',nCol);
  until (nRow<=20) and (nCol<=20) and (nRow>0) and (nCol>0);
  {Note: the limits here are chosen for the author's convenience.  Pascal
  requires that all quantities and variables be defined before they are
  used.  The reader is advised to exercise caution if the working array
  bounds are changed.}
  if nRow<nCol then
  begin
    writeln('WARNING -- number of rows less than number of columns');
    writeln('     To get a correct svd, form the svd of the transpose');
    writeln('     of your matrix, and transpose its decomposition.');
    writeln(confile,'WARNING -- number of rows less than number of columns');
    writeln(confile,'     To get a correct svd, form the svd of the transpose');
    writeln(confile,'     of your matrix, and transpose its decomposition.');
  end;
  Matrixin(nRow, nCol, A, tvec, sym); {A matrix of size nRow by nCol is
        read into array A.  If it is symmetric, then the matrix is also
        read into the vector tvec in row-wise lower-triangular form.  }
  Vectorin(nRow, Y); {A vector of nRow elements is read into Y}
  Matcopy(nRow,nCol, A, W); {The matrix A is copied into working array W.}
  NashSVD( nRow, nCol, W, Z); {The singular value decomposition is
        computed for matrix A by columnwise orthogonalization of the
        working array W, to which a unit matrix of order nCol is added
        in order to form the matrix V in the bottom nCol rows of W.}
  write('Do you want to print / test the svd results ? ([cr]=no) ');
  readln(infile,inchar); if infname<>'con' then writeln(inchar);
  writeln(confile,
      'Do you want to print / test the svd results ? ([cr]=no) ',inchar);
  if (inchar='Y') or (inchar='y') then
  begin
    for j:=1 to nCol do
    begin
      Z[j]:= sqrt(Z[j]);
      {MOD: 930714}
      for i:=1 to nRow do if Z[j]>0.0 then U[i,j]:=W[i,j]/Z[j]
                                      else U[i,j]:=0.0;
      for i:=1 to nCol do V[i,j]:=W[i+nRow,j];
    end;
    PrtSVDResults( nRow, nCol, U, V,Z);
    write('Do you want to test the svd results ? ([cr]=no) ');
    readln(infile,inchar); if infname<>'con' then writeln(inchar);
    writeln(confile,
        'Do you want to test the svd results ? ([cr]=no) ',inchar);
    if (inchar='Y') or (inchar='y') then
    begin
      svdtst(A,U,V,Z,nRow,nCol,nCol);
      write('Hit [cr] to continue ');readln(infile);
      if infname<>'con' then writeln;
      writeln('Reconstruction of initial matrix from Nash working form');
      writeln(confile,'Hit [cr] to continue ');
      writeln(confile,
          'Reconstruction of initial matrix from Nash working form');
      t2:=0.0; {to store largest error in reconstruction}
      for i:=1 to nRow do
      begin
        for j:=1 to nCol do
        begin
          t1:=0.0;
          for k:=1 to nCol do
            t1:=t1+W[i,k]*W[j+nRow,k]; { U * S * V-transpose}
        {writeln('A[',i,',',j,']=',A[i,j],' Recon. =',t1,' error=',A[i,j]-t1);}
          t1:=A[i,j]-t1; {to compute the residual}
          if abs(t1)>t2 then
          begin
            t2:=abs(t1); imax:=i; jmax:=j; {to save biggest element}
          end;
        end; {loop over columns}
      end; {loop over rows}
      writeln('Largest error is ',imax,',',jmax,' = ',t2);
      writeln(confile,'Largest error is ',imax,',',jmax,' = ',t2);
    end; {test svd results}
  end; {print results}
  svdlss(nRow, nCol, W, Y, Z, A, Bvec); {This computes the least squares
      solution to the problem A Z ~= Bvec, which will over-write the
      vector Z (currently storing the singular values).  The procedure
      also computes the residuals and their sum of squares for the
      problem.}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr0102.pas == svd and least squares solution}
