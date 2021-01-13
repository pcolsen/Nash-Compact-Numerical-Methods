{$I turbo.cnm}
program dr10(input,output);
{dr10.pas == driver to use Gauss elimination for inverse iteration
calculation of matrix eigensolutions

          Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps}
{$I matrixin}
{$I vectorin.pas}
{$I genevres.pas}
{$I rayquo.pas}
{$I alg05.pas}
{$I alg10.pas}
{$I startup.pas}

var
  A, Acopy, B, Bcopy : rmatrix;
  avector : smatvec;
  X, Y : rvector;
  errcode, i, itcount, j, k, m, n, nRHS, nRow, nCol : integer;
  rq, s, shift, ss : real;
  innum : string[20];
  sym, vectorOK : boolean;

begin
  banner:='dr10.pas -- inverse iteration via Gauss elimination';
  startup;
  write('order of problem (n) = '); readln(infile,n);
  if infname<>'con' then writeln(n);
  writeln(confile,'order of problem (n) = ',n);
  nRow:=n; nCol:=2*n; {store matrices in an array n by 2n}
  writeln('Provide the A matrix');
  writeln(confile,'Provide the A matrix');
  Matrixin(nRow, nCol, Acopy, avector, sym);
  writeln('A matrix');
  writeln(confile,'A matrix');
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
      write(Acopy[i,j]:10:5,' ');
      write(confile,Acopy[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<n) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  writeln('Provide the B matrix');
  writeln(confile,'Provide the B matrix');
  Matrixin(nRow, nCol, Bcopy, avector, sym);
  writeln('B matrix');
  writeln(confile,'B matrix');
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
       write(Bcopy[i,j]:10:5,' ');
       write(confile,Bcopy[i,j]:10:5,' ');
      if (7 * (j div 7) = j) and (j<n) then
      begin
        writeln;
        writeln(confile);
      end;
    end;
    writeln;
    writeln(confile);
  end;
  shift:=0.0; {rem initial and safety value for the eigenvalue shift}
  vectorOK:=false; {approximate eigenvector not yet defined}
  repeat
    if vectorOK then
    begin
      write('Do you wish to re-define the trial vector ([cr]=no) ');
      readln(infile,innum);
      writeln(confile,
        'Do you wish to re-define the trial vector ([cr]=no) ',innum);
      if infname<>'con' then writeln(innum);
      if (length(innum)>0) and
        (not ( (copy(innum,1,1)='n') or (copy(innum,1,1)='N') ) )
          then vectorOK:=false;   {to force user to re-define the
            trial eigenvector}
    end; {check if eigenvector approximation to be entered}
    if (not vectorOK) then
    begin
      writeln('Provide a starting vector for inverse iteration');
      writeln(confile,'Provide a starting vector for inverse iteration');
      vectorin(n, Y);
    end;
    vectorOK:=true; {set flag to indicate a vector is now in place}
    writeln('Enter a shift for eigenvalues ([cr] = ',shift,') ');
    write(' A value > 1E+30 halts execution. Entry = ');
    readln(infile,innum);
    writeln(confile,'Enter a shift for eigenvalues ([cr] = ',shift,') ');
    writeln(confile,' A value > 1E+30 halts execution. Entry = ',innum);
    if infname<>'con' then writeln(innum);
    if length(innum)>0 then val(innum, shift, errcode);
    if errcode<>0 then halt; {safety check}
    if shift<=1e30 then
    begin
      for i:=1 to n do  {copy matrices into working matrices}
      begin
        for j:=1 to n do
        begin
          A[i,j]:=Acopy[i,j]; B[i,j]:=Bcopy[i,j];
          A[i,j+n]:=B[i,j]; {to provide work matrix for ALG10}
        end; {loop on j}
      end; {loop on i}
      itcount:=100; {rem fairly liberal bound}
      gii(n,  A , Y, shift, itcount);
      writeln;
      writeln(confile);
      if itcount > 0 then
      begin
        writeln(
          'Converged to eigenvalue =',shift,'  in  ',itcount,' iterations');
        writeln(confile,
          'Converged to eigenvalue =',shift,'  in  ',itcount,' iterations');
      end
      else
      begin
        writeln('Not converged. Approximate eigenvalue=',shift,
                ' after ',-itcount,' iterations');
        writeln(confile,'Not converged. Approximate eigenvalue=',shift,
                ' after ',-itcount,' iterations');
      end; {else not converged}
      writeln('Eigenvector');

      writeln(confile,'Eigenvector');
      for i:=1 to n do
      begin
        write(Y[i]:10:5,' ');
        write(confile,Y[i]:10:5,' ');
        if (7* (i div 7) = i) and (i<n) then
        begin
          writeln;
          writeln(confile);
        end;
      end; {loop on i}
      writeln;
      writeln(confile);
      ss:=genevres(n, Acopy, Bcopy, shift, Y, true);
      rq:=rayquo(n, Acopy, Bcopy, Y);
      writeln('Rayleigh quotient = ',rq);
      writeln(confile,'Rayleigh quotient = ',rq);
    end; {if shift <=1e30}
  until (shift>1e30); {end of loop over possible shifts}
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr10.pas}
