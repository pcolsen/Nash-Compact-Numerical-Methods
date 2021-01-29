function rayquo(n :integer;
            A,B : rmatrix;
            Y : rvector): real;
{rayquo.pas
  == compute Rayleigh quotient. If the denominator of
      the quotient is zero, set the quotient to a large
      negative number (-big)
}
var
  s,t : real;
  i,j : integer;

begin
  s:=0.0; t:=0.0;
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
      s:=s+Y[i]*A[i,j]*Y[j];
      t:=t+Y[i]*B[i,j]*Y[j];
    end; {loop on j}
  end; {loop on i}
  if t>0.0 then rayquo:=s/t else rayquo:=-big;
  {note the safety value for the result}
end; {rayquo.pas == compute Rayleigh quotient}
