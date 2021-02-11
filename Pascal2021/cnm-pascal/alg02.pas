procedure svdlss(nRow, nCol: integer;
                 W : wmatrix;
                 Y: rvector;
                 Z : rvector;
                 A : rmatrix;
                 var Bvec: rvector);

var
 i, j, k : integer;
 q, s : real;

begin
  writeln('alg02.pas == svdlss');
  writeln(confile,'alg02.pas == svdlss');
  repeat
    writeln; writeln(confile);
    writeln('Singular values');
    writeln(confile,'Singular values');
    for j := 1 to nCol do
    begin
      write(sqrt(Z[j]):18,' ');
      write(confile,sqrt(Z[j]):18,' ');
      if j = 4 * (j div 4) then writeln;
      if j = 4 * (j div 4) then writeln(confile);
    end;
    writeln;
    write('Enter a tolerance for zero singular value ');
    readln(infile,q);
    if infname<>'con' then writeln(q);
    writeln(confile);
    writeln(confile,'Enter a tolerance for zero singular value ',q);
    if q>=0.0 then
    begin
      q := q*q;
      for i := 1 to nCol do
      begin
        s := 0.0;
        for j := 1 to nCol do
        begin
          for k := 1 to nRow do
          begin
            if Z[j]>q then
              s := s + W[i+nRow,j]*W[k,j]*Y[k]/Z[j];


          end;
        end;
        Bvec[i] := s;
      end;
      writeln('Least squares solution');
      writeln(confile,'Least squares solution');
      for j := 1 to nCol do
      begin
        write(Bvec[j]:12,' ');
        if j = 5 * (j div 5) then writeln;
        write(confile,Bvec[j]:12,' ');
        if j = 5 * (j div 5) then writeln(confile);
      end;
      writeln;
      writeln(confile);
      s := resids(nRow, nCol, A, Y, Bvec, true);
    end;
  until q<0.0;
end ;
