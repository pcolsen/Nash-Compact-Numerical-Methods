{rosen.pas
  == suite of procedures and functions defining the Rosenbrock
    banana shaped valley problem.
}
procedure fminset(var n:integer;var Bvec: rvector; var Workdata: probdata);
{sets up problem and defines starting values of Bvec}
{setup for Rosenbrock problem from rosen.pas}
begin
  writeln('Function: Rosenbrock Banana Valley');
  writeln(confile,'Function: Rosenbrock Banana Valley');
  n:=2;
  Workdata.m:=2; {for nonlinear least squares problems}
  Workdata.nvar:=0;
  Bvec[1]:=-1.2;
  Bvec[2]:=1.0;
  writeln('Classical starting point (-1.2,1)');
  writeln(confile,'Classical starting point (-1.2,1)');
end; {fminset from rosen.pas}

function fminfn(n: integer; var Bvec: rvector; var Workdata:probdata;
            var nocomp:boolean):real; {JN910113 var added to Workdata}
{this is the Rosenbrock banana valley function from rosen.pas}
begin
  nocomp:=false; {never undefined here}
  fminfn:=sqr(Bvec[2]-sqr(Bvec[1]))*100.0+sqr(1.0-Bvec[1]);
end; {fminfn from rosen.pas}
procedure fmingr(n:integer;Bvec:rvector; var Workdata:probdata; 
                           var g:rvector); {JN910113 var added to Workdata}
{computes the gradient of the Rosenbrock banana valley at point Bvec
  from rosen.pas}
begin
  g[1]:=-400.0*Bvec[1]*(Bvec[2]-sqr(Bvec[1]))-2.0*(1.0-Bvec[1]);
  g[2]:=200.0*(Bvec[2]-sqr(Bvec[1]));
end; {fmingrad from rosen.pas}

function nlres(i, n : integer; Bvec: rvector; var nocomp: boolean;
              var Workdata: probdata): real; {JN910113 var added to Workdata}
{computes residuals for the nonlinear least squares form of the
  Rosenbrock function from rosen.pas}
var
  temp: real;
begin
  nocomp:=false; {never set here}
  case i of
    1: begin
      temp:=10.0*(Bvec[2]-sqr(Bvec[1]));
    end;
    2: begin
      temp:=1.0-Bvec[1];
    end;
    else halt; {safety stop}
  end; {case}
  nlres := temp; {assign residual}
end; {nlres from rosen.pas}
procedure nljac(i, n: integer; Bvec: rvector; var jacrow: rvector;
           var Workdata: probdata); {JN910113 var added to Workdata}
{computes derivatives of residuals for the nonlinear least squares
  form of the Rosenbrock function from rosen.pas}
var
  t1, t2: real;
begin
  case i of
    1: begin
      jacrow[1]:=-20.0*Bvec[1];
      jacrow[2]:=10.0;
    end;
    2: begin
      jacrow[1]:=-1.0;
      jacrow[2]:=0.0;
    end;
    else halt; {safety stop}
  end; {case}
end; {nljac from rosen.pas}
{end of rosen.pas test function code suite}

