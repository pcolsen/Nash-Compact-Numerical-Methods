procedure svdlss(nRow, nCol: integer;
                 W : wmatrix;
                 Y: rvector;
                 Z : rvector;
                 A : rmatrix;
                 var Bvec: rvector;
                 q : real);

var
 i, j, k : integer;
 s : real;

begin
  writeln('alg02.pas == svdlss');
{  write('Y:');
  for i := 1 to nRow do 
  begin
    write(Y[i],' ');
  end;
  writeln;

  for i := 1 to (nRow+nCol) do
  begin
     write('W row ',i,':');
     for j:= 1 to nCol do
     begin
       write(W[i,j],' ');
     end;
     writeln;
   end;
}
{    writeln('Singular values');
    for j := 1 to nCol do
    begin
      write(Z[j]:18,' ');
      if j=4 * (j div 4) then writeln;
    end;
    writeln;
}
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
                       { V   S+   U'  y }

          end;
        end;
        Bvec[i] := s;
      end;
      writeln('Least squares solution');
      for j := 1 to nCol do
      begin
        write(Bvec[j]:12,' ');
        if j=5 * (j div 5) then writeln;
      end;
      writeln;
      s := resids(nRow, nCol, A, Y, Bvec, true);
    end;
end;

