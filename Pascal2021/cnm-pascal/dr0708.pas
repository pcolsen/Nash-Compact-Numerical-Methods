{$I turbo.cnm}
program dr0708(input,output);
{dr0708.pas == driver program to test procedures for Choleski (Alg07)
      and Choleski back-substitution (Alg08)

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I matrixin.pas}
{$I vectorin.pas}
{$I resids.pas}
{$I alg07.pas}
{$I alg08.pas}
{$I startup.pas}

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
  startup;
  write('order of problem = ');  readln(infile,nRow);
  writeln(confile,'order of problem = ',nRow);
  if infname<>'con' then writeln(nRow);
  nCol:=nRow; {use symmetric matrix in Choleski}
  Matrixin(nRow,nCol,A,avector,sym);
  writeln;
  writeln(confile);
  writeln('returned matrix of order ',nRow);
  writeln(confile,'returned matrix of order ',nRow);
  if not sym then halt; {must have symmetric matrix}
  begin
    writeln('Symmetric matrix -- Vector form');
    writeln(confile,'Symmetric matrix -- Vector form');
    k:=0;
    for i:=1 to nRow do
    begin
      for j:=1 to i do
      begin
         k:=k+1;
        write(avector[k]:10:5,' ');
        write(confile,avector[k]:10:5,' ');
        if (7 * (j div 7) = j) and (j<i) then
        begin
          writeln;
          writeln(confile);
        end;
      end;
      writeln;
      writeln(confile);
    end;
  end;
  writeln('Enter right hand side of equations');
  writeln(confile,'Enter right hand side of equations');
  vectorin(nRow, Y);
  for i:=1 to nRow do Ycopy[i]:=Y[i];
  writeln;
  writeln(confile);
  choldcmp(nRow,avector, singmat); {decompose matrix}
  begin
    writeln('Decomposed matrix -- Vector form');
    writeln(confile,'Decomposed matrix -- Vector form');
    k:=0;
    for i:=1 to nRow do
    begin
      for j:=1 to i do
      begin
        k:=k+1;
        write(avector[k]:10:5,' ');
        write(confile,avector[k]:10:5,' ');
        if (7 * (j div 7) = j) and (j<i) then
        begin
          writeln;
          writeln(confile);
        end;
      end;
      writeln;
      writeln(confile);
    end;
  end;
  if not singmat then
  begin
    Cholback(nRow,avector,Y);
    writeln('Solution');
    writeln(confile,'Solution');
    for i:=1 to nRow do
    begin
      write(Y[i]:10:5,' ');
      write(confile,Y[i]:10:5,' ');
      if (7 * (i div 7) = i) and (i<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
      writeln;
      writeln(confile);
    end;
    s:=resids(nRow,nCol,A,Ycopy,Y,true);
  end {non-singular case}
  else
  begin
    writeln('Matrix computationally singular -- solution not possible');
    writeln(confile,
        'Matrix computationally singular -- solution not possible');
  end;
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr0708.pas}
