{JJACF.PAS
  == suite of procedures and functions defining the Jaffrelot
    first order autocorrelation minimisation function.
    Example 14-1. This is designed to use data in file
    EX1920J.CNM with driver program DR1920.PAS, or file
    EX27J.CNM with driver program DR27.PAS.

    The driver programs MUST have the INCLUDE statements
    modified to invoke the present code (JJACF.PAS).
}
procedure fminset(var n:integer;var Bvec: rvector; var Workdata: probdata);
{sets up problem and defines starting values of Bvec}
var
  i: integer;
  descn : string[80];
{setup for problem from JJACF.PAS}
begin
  with Workdata do
  begin
    writeln('Function: Jaffrelot Minimisation of First Order ACF');
    writeln(confile,'Function: Jaffrelot Minimisation of First Order ACF');
    nvar:=2;
    write('Number of points in data series = ');
    readln(infile,m);
    if (infname<>'con') then writeln(m);
    writeln(confile,'Number of points in data series = ',m);
    readln(infile,descn);
    if (infname<>'con') then writeln(descn);
    writeln(confile,descn);
    for i:=1 to m do
    begin
      read(infile,Ydata[i,1]);
      if infname<>'con' then write(Ydata[i,1]:10:5);
      write(confile,Ydata[i,1]:10:5);
      if (7*(i div 7)=i) and (i<m) then
      begin
        writeln;
        writeln(confile);
      end;
    end; {for i}
    writeln;
    writeln(confile);
  end; {with Workdata}
  n:=2;
  writeln('Enter starting parameters');
  writeln(confile,'Enter starting parameters');
  readln(infile,Bvec[1],Bvec[2]);
  if infname<>'con' then writeln(Bvec[1],' ',Bvec[2]);
  writeln('starting point (',Bvec[1],',',Bvec[2],')');
  writeln(confile,'starting point (',Bvec[1],',',Bvec[2],')');
end; {fminset from JJACF.PAS}

function fminfn(n: integer; var Bvec: rvector; Workdata:probdata;
            var nocomp:boolean):real;
{this is the function from JJACF.PAS}
var
  i, j : integer;
  u, s, t, r1 : real;

begin
  nocomp:=false; {never undefined here}
  with Workdata do
  begin
    for i:=3 to m do
    begin
      Ydata[i,2]:=Ydata[i,1]-Bvec[1]*Ydata[i-1,1]-Bvec[2]*Ydata[i-2,1];
    end;
    u:=0.0;
    for i:=3 to m do u:=u+Ydata[i,2];
    u:=u/(m-2); {mean in Kendall definition}
    s:=0.0; t:=0.0;
    for i:=4 to m do s:=s+(Ydata[i,2]-u)*(Ydata[i-1,2]-u);
    for i:=3 to m do t:=t+sqr(Ydata[i,2]-u);
    if t=0.0 then
    begin
      writeln('Variance is zero -- stopping');
      writeln(confile,'Variance is zero -- stopping');
      halt;
    end;
    r1:=(m-2)*s/((m-3)*t);
  end; {with Workdata}
  fminfn:=sqr(r1);
end; {fminfn from JJACF.PAS}

procedure fmingr(n:integer;Bvec:rvector; Workdata:probdata; var g:rvector);
{computes the gradient at point Bvec from JJACF.PAS}
begin
  halt; {not provided here}
end; {fmingrad from JJACF.PAS}

function nlres(i, n : integer; Bvec: rvector; var nocomp: boolean;
                                            var Workdata: probdata): real;
{computes residuals for the nonlinear least squares form of the
  function from JJACF.PAS}
var
  temp: real;
begin
  halt; {not provided here}
end; {nlres from JJACF.PAS}

procedure nljac(i, n: integer; Bvec: rvector; var jacrow: rvector;
                                              var Workdata:probdata);
{computes derivatives of residuals for the nonlinear least squares
  form of the function from JJACF.PAS}
var
  t1, t2: real;
begin
  halt; {not provided here}
end; {nljac from JJACF.PAS}
{end of JJACF.PAS test function code suite}
