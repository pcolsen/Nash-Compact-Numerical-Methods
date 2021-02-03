procedure cholback(n: integer;
                   a: smatvec;
                   var x: rvector);
var
 i,j,q  : integer;

begin
  if a[1]=0.0 then x[1]:=0.0
              else x[1]:=x[1]/a[1];
  if n>1 then
  begin
    q:=1;
    for i:=2 to n do
    begin
      for j:=1 to (i-1) do
      begin
        q:=q+1; x[i]:=x[i]-a[q]*x[j];
      end;
      q:=q+1;
      if a[q]=0.0 then x[i]:=0.0
                  else x[i]:=x[i]/a[q];
    end;
  end;

  if a[n*(n+1) div 2]=0.0 then x[n]:=0.0
                           else x[n]:=x[n]/a[n*(n+1) div 2];
  if n>1 then
  begin
    for i:=n downto 2 do
    begin
      q:=i*(i-1) div 2;
      for j:=1 to (i-1) do x[j]:=x[j]-x[i]*a[q+j];
      if a[q]=0.0 then x[i-1]:=0.0
                  else x[i-1]:=x[i-1]/a[q];
    end;
  end;
end;

var
  A : rmatrix;
  avector : smatvec;
  i, j, k, nCol, nRow : integer;
  sym : boolean;
  Y, Ycopy : rvector; {to store the right hand side of the equations}
  singmat : boolean; {set true if matrix discovered to be computationally
              singular during alg07.pas}
  s : real; {an accumulator}

begin
  banner:='dr0708 -- Choleski decomposition and back-substitution';
  write('order of problem = ');  
  readln(nRow);
  writeln(nRow);
  nCol:=nRow; {use symmetric matrix in Choleski}
  Matrixin(nRow,nCol,A,avector,sym);
  writeln;
  writeln('returned matrix of order ',nRow);
  if not sym then halt; {must have symmetric matrix}
  begin
    writeln('Symmetric matrix -- Vector form');
    k:=0;
    for i:=1 to nRow do
    begin
      for j:=1 to i do
      begin
         k:=k+1;
        write(avector[k]:10:5,' ');
        if (7 * (j div 7) = j) and (j<i) then writeln;
      end;
      writeln;
    end;
  end;
  writeln('Enter right hand side of equations');
  vectorin(nRow, Y);
  for i:=1 to nRow do Ycopy[i]:=Y[i];
  writeln;
  choldcmp(nRow,avector, singmat); {decompose matrix}
  begin
    writeln('Decomposed matrix -- Vector form');
    k:=0;
    for i:=1 to nRow do
    begin
      for j:=1 to i do
      begin
        k:=k+1;
        write(avector[k]:10:5,' ');
        if (7 * (j div 7) = j) and (j<i) then writeln;
      end;
      writeln;
    end;
  end;
  if not singmat then
  begin
    Cholback(nRow,avector,Y);
    writeln('Solution');
    for i:=1 to nRow do
    begin
      write(Y[i]:10:5,' ');
      if (7 * (i div 7) = i) and (i<nRow) then writeln;
      writeln;
    end;
    s:=resids(nRow,nCol,A,Ycopy,Y,true);
  end {non-singular case}
  else
  begin
    writeln('Matrix computationally singular -- solution not possible');
  end;
end. {dr0708.pas}
