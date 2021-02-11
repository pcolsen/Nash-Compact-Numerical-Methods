{$I turbo.cnm}
program dr24ii(input, output);
{dr24ii.pas == modification of conjugate gradients method to solve
          generalized symmetric matrix eigenvalue problems by inverse
          iteration

  Note, this approach ignores the definiteness of the matrices
  which is supposedly needed for proper operation of the conjugate
  gradients method.

              Copyright 1988 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}
{$I matrixin.pas}
{$I vectorin.pas}
{$I matmul.pas}
{$I alg24.pas}
{$I rayquo.pas}
{$I genevres.pas}
{$I startup.pas}

var
  A, B, H : rmatrix;
  Y : rvector; {RHS in solution of linear equations}
  X : rvector; {solution}
  RHS : rvector; {used for solution of linear equations sub-problem}
  avec : smatvec; {for matrixin only}
  shift : real; {for the eigenvalue shift}
  evalue : real; {for the eigenvalue}
  sym : boolean; {to tell if matrix symmetric}
  ch : char;
  i, j, m, n, itcount : integer;
  iicount: integer; {to count iterations}
  ssmin, t, s, ymax : real;
  tol : real;
  mchoice : integer; {to control flow in program}
  innum : string[20]; {to enter numbers as strings}
  evcon : boolean;
  oldev : real;
  iilimit : integer;
  rq : real; {Rayleigh quotient value}

begin
  banner:='dr24ii.pas -- cg inverse iteration for matrix eigenproblem';
  startup;
  write('Order of problem = '); readln(infile,n);
  if infname<>'con' then writeln(n);
  writeln(confile,'Order of problem = ',n);
  writeln('Matrix A'); writeln(confile,'Matrix A');
  matrixin(n, n, A, avec, sym);
  writeln('Matrix B'); writeln(confile,'Matrix B');
  matrixin(n, n, B, avec, sym);
  if not sym then halt;
  shift:=0.0; {safety setting}
  write('Shift value for eigenvalues ([cr] =',shift,') = ');
  readln(infile,innum);
  writeln(confile,
      'Shift value for eigenvalues ([cr] =',shift,') = ',innum);
  if infname<>'con' then writeln(innum);
  i:=0;
  if length(innum)>0 then val(innum,shift,i);
  if i<>0 then halt; {error checking}
  mchoice:=1; {to indicate we are going to input new vector}
  evalue:=shift; {initial guess}
  while (mchoice<>0) do
  begin
    if mchoice>0 then
    begin
    writeln('Initial guess for eigenvector');
    writeln(confile,'Initial guess for eigenvector');
    vectorin(n, Y);
    end; {new initial vector}
    {This is where inverse iteration loop begins.}
    ssmin:=genevres(n,A,B,evalue,Y,false);
    writeln('Residual sumsquares from trial solution = ',ssmin);
    writeln(confile,'Residual sumsquares from trial solution = ',ssmin);
    iicount:=0; iilimit:=4*n; {loose limit here}
    tol:=calceps; tol:=tol*(abs(shift)+tol); {to be used as a small number}
    oldev:=big; {to ensure no convergence in one iteration}
    writeln('   ev itns':10,'   cg itns':10,'   cgss':14,
            '   ev approx':16,'   RQ':16);
    repeat {inverse iteration loop}
      iicount:=iicount+1;
      for i:=1 to n do {store shifted matrix in H}
        for j:=1 to n do H[i,j]:=A[i,j]-shift*B[i,j];
      for i:=1 to n do
      begin {form RHS:=B * X as rhs of local equations}
        t:=0.0;
        for j:=1 to n do t:=t+B[i,j]*Y[j];
        RHS[i]:=t; X[i]:=Y[i]; {save eigenvector elements for comparison}
      end;
      itcount:=n+1; {safety setting -- short iteration}
      lecg( n, H, RHS, Y, itcount, ssmin);
      {call to alg24.pas -- linear equations by conjugate gradients}
      write(iicount:8,'  ',itcount:8,'    ',ssmin:12,'  ');
      s:=abs(Y[1]); m:=1; {index of current maximum}
      for i:=2 to n do {find largest eigenvector component}
      begin
        if abs(Y[i]) > s then
        begin
          m:=i; s:=abs(Y[i]);
        end;
      end; {loop on i}
      evalue:=shift+X[m]/Y[m];
      write(evalue:14,'  ');  write(confile,evalue:14,'  ');
      {Normalize eigenvector.}
      ymax:=Y[m];
      for i:=1 to n do Y[i]:=Y[i]/ymax; {to normalize the eigenvector}
      rq:=rayquo(n,A,B,Y);    {Compute Rayleigh quotient.}
      writeln(rq:14); writeln(confile,rq:14);
      evcon:= ( (reltest+evalue) = (reltest+oldev) ); {check if evalue converged}
      oldev:=evalue; {to save eigenvalue}
    until ((ssmin<=tol) and (itcount>0) and evcon)
          or (iicount>iilimit) or (ssmin=-big);
    {loop until convergence or failure}
    {Note this is a strict criterion. It could be relaxed somewhat and
    still be satisfactory for many purposes.}
    if iicount>iilimit then
    begin
      writeln('Iteration limit exceeded');
      writeln(confile,'Iteration limit exceeded');
    end;
    ssmin:=genevres(n,A,B,evalue,Y,false);
    writeln('Residual sumsquares of normalized vector = ',ssmin);
    writeln(confile,'Residual sumsquares of normalized vector = ',ssmin);
    writeln('Enter a new shift, or [cr] for ',shift,' ');
    write(' or V for new initial eigenvector, or S to stop ');
    readln(infile,innum); if infname<>'con' then writeln(innum);
    writeln(confile,'Enter a new shift, or [cr] for ',shift,' ');
    writeln(confile,' or V for new initial eigenvector, or S to stop ',
            innum);
    while ((length(innum)>0) and (copy(innum,1,1)=' ')) do
        innum:=copy(innum,2,length(innum)-1);
    if length(innum)=0 then mchoice:=-1
    else
      if (copy(innum,1,1)='S') or (copy(innum,1,1)='s') then mchoice:=0
        else
          if (copy(innum,1,1)='V') or (copy(innum,1,1)='v')
                then mchoice:=1
          else
          begin
            i:=0; {for error code}
            val(innum,shift,i);
            if i<>0 then halt; {error check}
            mchoice:=-1;
          end;
  end; {while mchoice<>=0}
  writeln('Eigensolution for eigenvalue =',evalue);
  writeln('           Rayleigh quotient =',rq);
  writeln(confile,'Eigensolution for eigenvalue =',evalue);
  writeln(confile,'           Rayleigh quotient =',rq);
  for i:=1 to n do
  begin
    write(Y[i]:10:7,' '); write(confile,Y[i]:10:7,' ');
    if (7 * (i div 7) = i) and (i<n) then
    begin
      writeln; writeln(confile);
    end;
  end;
  writeln; writeln(confile);
  ssmin:=genevres(n,A,B,evalue,Y,true);
  flush(confile); close(confile); if infname<>'con' then close(infile);
end. {dr24ii.pas == inverse iteration by conjugate gradients}
