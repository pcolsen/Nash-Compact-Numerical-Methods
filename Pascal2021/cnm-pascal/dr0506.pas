{$I turbo.cnm}
program dr0506(input,output);
{dr05606.pas == driver for gelim (alg05) and gebacksub (alg06),
    Gauss elimination and linear equations solution

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}
{$I matrixin.pas}
{$I alg05.pas}
{$I alg06.pas}
{$I startup.pas}

var
  i,j,k, nRHS,nRow,nCol : integer;
  rss : real; {to accumulate residual sum of squares}
  s : real;
  tol : real; {tolerance for zero pivot}
  A, Acopy : rmatrix;
  avector : smatvec; {needed only for matrixin}
  sym : boolean;

begin
  banner:='dr0506.pas -- Gauss elimination and linear equations';
  startup;
  begin
    write('order of problem (n) = '); readln(infile,nRow);
    if infname<>'con' then writeln(nRow);
    write('number of right hand sides (nRHS) = '); readln(infile,nRHS);
    if infname<>'con' then writeln(nRHS);
    writeln(confile,'order of problem (n) = ',nRow);
    writeln(confile,'number of right hand sides (nRHS) = ',nRHS);
    nCol:=nRow+nRHS;
    Matrixin(nRow, nCol, A, avector, sym);
  end;
  writeln;
  writeln(confile);
  writeln('Data matrix :',nRow,' by ',nRow+nRHS);
  writeln(confile,'Data matrix :',nRow,' by ',nRow+nRHS);
  for i:=1 to nRow do
  begin
    writeln('Row ',i);
    writeln(confile,'Row ',i);
    for j:=1 to (nRow+nRHS) do
    begin
      write(A[i,j]:10:5,' ');
      write(confile,A[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nRow+nRHS) then
      begin
        writeln;
        writeln(confile);
      end;
      Acopy[i,j]:=A[i,j];
    end;
    writeln;
    writeln(confile);
  end;
  writeln;
  writeln(confile);
  {Note: a negative value for tol will cause the machine epsilon times
  a matrix norm to be used for the actual tolerance.}
  writeln('Enter a tolerance for zero pivot. A negative value will ');
  writeln(' cause machine precision times a matrix norms to be used.');
  write(' Value for tol = '); readln(infile,tol);
  writeln(confile,'Enter a tolerance for zero pivot. A negative value will ');
  writeln(confile,' cause machine precision times a matrix norms to be used.');
  writeln(confile,' Value for tol = ',tol);
  if infname<>'con' then writeln(tol);
  if tol<0.0 then
  begin
    tol:=0.0;
    for i:=1 to nRow do
    begin
      s:=0.0;
      for j:=1 to nRow do s:=s+abs(A[i,j]);
      if s>tol then tol:=s;
    end;
    tol:=tol*calceps; {tol has maximum row norm * EPSILON}
    writeln('Replacing tol with ',tol);
    writeln(confile,'Replacing tol with ',tol);
  end;
  Gelim(nRow,nRHS,A,tol);
  {writeln('after Gauss elimination');
  writeln('returned matrix ',nRow,' by ',nRow+nRHS);
  writeln(confile,'after Gauss elimination');
  writeln(confile,'returned matrix ',nRow,' by ',nRow+nRHS);
  for i:=1 to nRow do
  begin
    writeln('Row ',i);
    writeln(confile,'Row ',i);
    for j:=1 to (nRow+nRHS) do
    begin
      write(A[i,j]:10:5,' ');
      write(confile,A[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<nRow+nRHS) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  writeln;
  writeln(confile);}
  GEBacksub(nRow,nRHS,A);
  writeln;
  writeln(confile);
  for i:=1 to nRHS do
  begin
    writeln('Solution ',i);
    for j:=1 to nRow do
    begin
      write(A[j,nRow+i]:10:5,' ');
      write(confile,A[j,nRow+i]:10:5,' ');
    if (7 * (j div 7) = j) and (j<nRow) then
    begin
      writeln;
      writeln(confile);
    end;
    end;
    writeln;
    writeln(confile);
    writeln('Residuals');
    writeln(confile,'Residuals');
    rss:=0.0; {initialize sum of squared residuals}
    {Rather than use resids.pas, this simple code has been placed in line.}
    for j:=1 to nRow do
    begin
      s:=Acopy[j,nRow+i]; {right hand side}
      for k:=1 to nRow do s:=s-Acopy[j,k]*A[k,nRow+i];
      write(s:10,' ');
      write(confile,s:10,' ');
      rss:=rss+s*s;
      if (7 * (j div 7) = j) and (j<nRow) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
    writeln('Sum of squared residuals = ',rss);
    writeln(confile,'Sum of squared residuals = ',rss);
    writeln;
    writeln(confile);
  end; {loop over solutions}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr0506.pas}
