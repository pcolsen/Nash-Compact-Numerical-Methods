Procedure matcopy(nRow ,nCol: integer; A: rmatrix; var B:wmatrix);
{matcopy.pas
  -- copies matrix A, nRow by nCol, into matrix B }
var i,j: integer;
begin
  for i:=1 to nRow do
    for j:=1 to nCol do
      B[i,j]:=A[i,j];
end;{matcopy.pas}
