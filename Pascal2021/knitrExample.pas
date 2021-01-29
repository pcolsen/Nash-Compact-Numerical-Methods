program paskio(input,output);

var
  n, i: integer;
  x: real;

begin
  n := 1000;
  for i:=1 to n do
  begin
    x := exp(sin(cos(i)));
  end;
  writeln('x=',x);
end.
