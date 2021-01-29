procedure vectorin(n: integer; var Y: rvector);
{vectorin.pas
  == enter or generate a vector of n real elements
}
var
  i, j, k, m, nt : integer;
  x : real;

begin
  write('vectorin.pas');
  writeln('   -- enter or generate a real vector of ',n,' elements');
  writeln('Options:');
  writeln('   1) constant');
  writeln('   2) uniform random in [0,user_value) ');
  writeln('   3) user entered from console ');
  writeln('   4) entered from RHS columns in matrix file ');
  write('   Choose option :');
  readln(infile,i);
  if infname<>'con' then writeln(i);
  write(confile,'vectorin.pas');
  writeln(confile,'   -- enter or generate a real vector of ',n,' elements');
  writeln(confile,'Options:');
  writeln(confile,'   1) constant');
  writeln(confile,'   2) uniform random in [0,user_value) ');
  writeln(confile,'   3) user entered ');
  writeln(confile,'   4) entered from RHS columns in matrix file ');
  writeln(confile,'   Choose option :',i);
  Case i of
    1 : begin
        write('Enter constant value ='); readln(infile,x);
        if infname<>'con' then writeln(x);
        writeln(confile,'Enter constant value =',x);
        for j:=1 to n do Y[j]:=x;
    end;
    2 : begin
        write('Enter the upper bound to the generator =');
        readln(infile,x);
        if infname<>'con' then writeln(x);
        writeln(confile,'Enter the upper bound to the generator =',x);
        for j:=1 to n do Y[j]:=Random;
        {According to the Turbo Pascal manual, version 3.0, Random
          returns a number in [0,1). My experience is that most such
          pseudo-random number generators leave a lot to be desired
          in terms of statistical properties. I do NOT recommend it
          for serious use in Monte Carlo calculations or other situations
          where a quality generator is required. For a better generator
          in Pascal, see Wichman B. and Hill, D., (1987)}
    end;
    3 : begin
        writeln('Enter elements of vector one by one');
        writeln(confile,'Enter elements of vector one by one');
        for j:=1 to n do
        begin
          write('Y[',j,']=');
          readln(infile,Y[j]);
          if infname<>'con' then writeln(Y[j]);
          writeln(confile,'Y[',j,']=',Y[j]);
        end;
    end;
    4 : begin  {Get elements of RHS from a matrix + vectors file}
        write('Datafile ');
        readln(infile,dfname);
        if infname<>'con' then writeln(dfname);
        writeln(confile,'Datafile ',dfname);
        write('Which RHS vector should be retrieved? ');
        readln(j);
        if infname<>'con' then writeln(j);
        writeln(confile,'Which RHS vector should be retrieved? ',j);
        {We will use the least squares data file a a form of input}
        if length(dfname)>0 then
        begin
          assign(dfile, dfname);
          reset(dfile);
          read(dfile, nt, i);{reading i rhs elements}
          writeln('Number of columns =',nt);
          writeln(confile,'Number of columns =',nt);
          m:=0; {to count the number of rows}
          while (not eof(dfile)) do
          begin
            m:=m+1;
            for k:=1 to nt do read(dfile,x); {ignore coefficient matrix}
            for k:=1 to j do
            begin
              read(dfile,x); if k=j then Y[m]:=x;
            end;
          end; {while}
          close(dfile);
          writeln('Found ',m,' elements in vector');
          writeln(confile,'Found ',m,' elements in vector');
          if m<>n then
          begin
            writeln('*** ERROR *** not in agreement with procedure call');
            writeln(confile,
                '*** ERROR *** not in agreement with procedure call');
            halt;
          end;
        end {if length(dfname)}
      end; {case 4}
  end {case};
end {vectorin.pas};


