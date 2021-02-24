Procedure matrixin(var m, n: integer; var A: rmatrix;
              var avector: smatvec; var sym :boolean);

{matrixin.pas --

  This procedure generates an m by n real matrix in both
  A or avector.

  A is of type rmatrix, an array[1..nmax, 1..nmax] of real
  where nmax >= n for all possible n to be provided.

  avector is of type rvector, an array[1..nmax*(nmax+1)/2]
  of real, with nmax as above.

  sym is set true if the resulting matrix is symmetric.
}

var
  temp : real;
  i,j,k: integer;
  inchar: char;
  mtype: integer;
  mn : integer;

begin
  if (m=0) or (n=0) then
  begin
    writeln;
    writeln('******* Matrix dimensions zero *********');
    writeln(confile);
    writeln(confile,'******* Matrix dimensions zero *********');
    halt;
  end;
  writeln('Matrixin.pas -- generate or input a real matrix ',m,' by ',n);
  writeln('Possible matrices to generate:');
  writeln('0) Keyboard or console file input');
  writeln('1) Hilbert segment');
  writeln('2) Ding Dong');
  writeln('3) Moler');
  writeln('4) Frank symmetric');
  writeln('5) Bordered symmetric');
  writeln('6) Diagonal');
  writeln('7) Wilkinson W+');
  writeln('8) Wilkinson W-');
  writeln('9) Constant');
  writeln('10) Unit');
  writeln(confile,'Matrixin.pas -- generate or input a real matrix ',m,' by ',n);
  writeln(confile,'Possible matrices to generate:');
  writeln(confile,'0) Keyboard or console file input');
  writeln(confile,'1) Hilbert segment');
  writeln(confile,'2) Ding Dong');
  writeln(confile,'3) Moler');
  writeln(confile,'4) Frank symmetric');
  writeln(confile,'5) Bordered symmetric');
  writeln(confile,'6) Diagonal');
  writeln(confile,'7) Wilkinson W+');
  writeln(confile,'8) Wilkinson W-');
  writeln(confile,'9) Constant');
  writeln(confile,'10) Unit');
{ Note: others could be added.}
  mn:=n;
  if m>mn then mn:=m; {mn is maximum of m and n}
  write('Enter type to generate ');
  readln(infile,mtype);
  if infname<>'con' then writeln(mtype);
  writeln(confile,'Enter type to generate ',mtype);
  case mtype of
    0: begin
      sym:=false;
      if m=n then
      begin
        write('Is matrix symmetric? ');  readln(infile,inchar);
        if infname<>'con' then writeln(inchar);
        writeln(confile,'Is matrix symmetric? ',inchar);
        if (inchar='y') or (inchar='Y') then sym:=true else sym:=false;
      end; {ask if symmetric}
      if sym then
      begin
        for i:=1 to n do
        begin
          writeln('Row ',i,' lower triangle elements');
          writeln(confile,'Row ',i,' lower triangle elements');
          for j:=1 to i do
          begin
            read(infile,A[i,j]);
            if infname<>'con' then write(A[i,j]:10:5,' ') else write(' ');
            write(confile,A[i,j]:10:5,' ');
            A[j,i]:=A[i,j];
            if (7*(j div 7) = j) and (j<i) then
            begin
              writeln;
              writeln(confile);
            end;
          end;
          writeln;
          writeln(confile);
        end;
      end {symmetric matrix}
      else
      begin {not symmetric}
        for i:=1 to m do
        begin
          writeln('Row ',i);
          writeln(confile,'Row ',i);
          for j:=1 to n do
          begin
            read(infile,A[i,j]);
            if infname<>'con' then write(A[i,j]:10:5,' ') else write(' ');
            write(confile,A[i,j]:10:5,' ');
          end; {loop on j}
          writeln;
          writeln(confile);
        end; {loop on i}
      end; {else not symmetric}
    end; {case 0 -- input of matrix}
    1: begin {Hilbert}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=1.0/(i+j-1.0);
      if m=n then sym:=true;
    end;
    2: begin {Ding Dong}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=0.5/(1.5+n-i-j);
      if m=n then sym:=true;
    end;
    3: begin {Moler}
      for i:=1 to mn do
      begin
        for j:=1 to i do
        begin
          A[i,j]:=j-2.0;
          A[j,i]:=j-2.0;
        end;
        A[i,i]:=i;
        if m=n then sym:=true;
      end;
    end;
    4: begin {Frank symmetric}
      for i:=1 to mn do
        for j:=1 to i do
        begin
          A[i,j]:=j;
          A[j,i]:=j;
        end;
        if m=n then sym:=true;
    end;
    5: begin {Bordered}
      temp:=2.0;
      for i:=1 to (mn-1) do
      begin
        temp:=temp/2.0; {2^(1-i)}
        for j:=1 to mn do
          A[i,j]:=0.0;
        A[i,mn]:=temp;
        A[mn,i]:=temp;
        A[i,i]:=1.0;
      end;
      A[mn,mn]:=1.0;
      if m=n then sym:=true;
    end;
    6: begin {Diagonal}
      for i:=1 to mn do
      begin
        for j:=1 to mn do
          A[i,j]:=0.0;
        A[i,i]:=i;
      end;
      if m=n then sym:=true;
    end;
    7: begin {W+}
      k:=mn div 2; {[n/2]}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=0.0;
      if m=n then sym:=true;
      for i:=1 to k do
      begin
        A[i,i]:=k+1-i;
        A[mn-i+1,mn-i+1]:=k+1-i;
      end;
      for i:=1 to mn-1 do
      begin
        A[i,i+1]:=1.0;
        A[i+1,i]:=1.0;
      end;
    end;
    8: begin {W-}
      k:=mn div 2; {[n/2]}
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=0.0;
      if m=n then sym:=true;
      for i:=1 to k do
      begin
        A[i,i]:=k+1-i;
        A[mn-i+1,mn-i+1]:=i-1-k;
      end;
      for i:=1 to mn-1 do
      begin
        A[i,i+1]:=1.0;
        A[i+1,i]:=1.0;
      end;
      if m=n then sym:=true;
    end;
    9: begin {Constant}
      write('Set all elements to a constant value = ');
      readln(infile,temp);
      if infname<>'con' then writeln(temp);
      for i:=1 to mn do
        for j:=1 to mn do
          A[i,j]:=temp;
      if m=n then sym:=true;
    end;
    10: begin {Unit}
      for i:=1 to mn do
      begin
        for j:=1 to mn do A[i,j]:=0.0;
        A[i,i]:=1.0;
      end;
      if m=n then sym:=true;
    end;
    else {case statement else
    begin
      writeln;
      writeln(confile);
      writeln('*** ERROR *** unrecognized option');
      writeln(confile,'*** ERROR *** unrecognized option');
      halt;
    end; {else of case statement}
  end; {case statement}
  if sym then
  begin {convert to vector form}
    k:=0; {index for vector element}
    for i:=1 to n do
    begin
      for j:=1 to i do
      begin
        k:=k+1;
        avector[k]:=A[i,j];
      end;
    end;
  end;
end; {matrixin}
