{logistic.pas
  == suite of procedures and functions defining the 3 parameter logistic
    modelling function and residuals
}
procedure fminset(var n:integer;var Bvec: rvector; var Workdata: probdata);
{sets up problem and defines starting values of Bvec}
{setup for 3 parameter logistic function}
var i: integer;
begin
  writeln('Function: 3 parameter logistic growth curve');
  writeln(confile,'Function: 3 parameter logistic growth curve');
  n:=3;
  with Workdata do
  begin
    readln(infile,m);
    writeln(' m =',m:0);
    writeln(confile,' m =',m:0);
    nvar:=2;
    for i:=1 to m do
    begin
      readln(infile, Ydata[i,1]);
      write(Ydata[i,nvar]:10:5,'  ');
      write(confile,Ydata[i,1]:10:5,'  ');
      if (5 * int(i/5) = i ) then
      begin
        writeln; writeln(confile);
      end;
    end; {for i}
    writeln; writeln(confile);
    readln(infile,Bvec[1]);
    readln(infile,Bvec[2]);
    readln(infile,Bvec[3]);
    writeln('Starting values:',Bvec[1],' ',Bvec[2],' ',Bvec[3]);
    writeln(confile,'Starting values:',Bvec[1],' ',Bvec[2],' ',Bvec[3]);
  end; {with Workdata}
end; {fminset from rosen.pas}

function nlres(i, n : integer; Bvec: rvector; var nocomp: boolean;
                                            var Workdata: probdata): real;
{computes residuals for the nonlinear least squares 3 parameter logistic
   function growth curve  logistic.pas}
{NOTE: Workdata has been added to the procedure call. It has been
  declared with the 'var' modifier because intermediate data is
  saved withing the Workdata structure (the residuals).}
var
  temp: real;
begin
  if abs(0.1*i*Bvec[3])>50.0 then nocomp:=true
  else
  begin
   nocomp:=false; {computable function}
   temp:=100.0*Bvec[1]/(1.0+10.0*Bvec[2]*exp(-0.1*Bvec[3]*I))-
        Workdata.Ydata[I,1];
   Workdata.Ydata[i,2]:=temp;
   nlres := temp; {assign residual}
  end; {else}
end; {nlres from rosen.pas}
procedure nljac(i, n: integer; Bvec: rvector; var jacrow: rvector;
                                            var Workdata: probdata);
{computes derivatives of residuals for the nonlinear least squares
  form of the Rosenbrock function from rosen.pas}
{NOTE: Workdata has been added to the procedure call. It has been
  declared with the 'var' modifier in case intermediate data is
  saved withing the Workdata structure.}
var
  t1, t2: real;
begin
  jacrow[1]:=100./(1.+10.*Bvec[2]*EXP(-0.1*Bvec[3]*I));
  jacrow[2]:=-Bvec[1]*0.1*sqr(jacrow[1])*EXP(-0.1*Bvec[3]*I);
  jacrow[3]:=-0.1*jacrow[2]*Bvec[2]*I;
end; {nljac from logistic.pas}
{end of logistic.pas test function code suite}
