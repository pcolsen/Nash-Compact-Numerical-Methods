Procedure PrtSVDResults( nRow, nCol:integer;
                U, V: rmatrix; Z: rvector);
{psvdres.pas
  == routine to display svd results and print them to confile
}
var
  i, j : integer;

begin
  writeln(' Singular values and vectors:');
  writeln(confile,' Singular values and vectors:');
  for j := 1 TO nCol do
  begin
    writeln('Singular value (',j,') =', Z[j]);
    writeln('Principal coordinate (U):');
    writeln(confile,'Singular value (',j,') =', Z[j]);
    writeln(confile,'Principal coordinate (U):');
    for i := 1 to nRow do
    begin
      write(U[i,j]:10:7);
      write(confile,U[i,j]:10:7);
      if (7 * (i div 7) = i) and (i<nRow) then writeln;
      if (7 * (i div 7) = i) and (i<nRow) then
      writeln(confile);
    end;
    writeln;
    writeln(confile);
    writeln('Principal component (V):');
    writeln(confile,'Principal component (V):');
    for i:=1 to nCol do
    begin
      write(V[i,j]:10:7);
      write(confile,V[i,j]:10:7);
      if (7 * (i div 7) = i) and (i<nCol) then writeln;
      if (7 * (i div 7) = i) and (i<nCol) then
      writeln(confile);
    end;
    writeln;
    writeln(confile);
  end;
end; {psvdres == print svd results via procedure PrtSVDResults }
