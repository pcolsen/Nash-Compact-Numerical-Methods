function fn1d(x: real;var nocomp: boolean): real;
  {real valued test function of x for [1D] minimisation}
{htanfn.pas
  == This is the hyperbolic tangent Example 13.2 from
    Compact Numerical Methods.

If this function is called with nocomp TRUE, the root will
be displayed. Otherwise, the function at x will be computed.
}
var
  val,y: real;
  s,t,w,z: real; {constants which define the function}

begin
  s:=1; t:=0.5; w:=0.2; z:=100.0; {fix a particular function}
  {This function is well-behaved.}
  {The following altenative set of parameters gives a function which is
  nearly a step-function.}
{  s:=100; t:=0.5; w:=0.99; z:=100.0;}
  if nocomp then
  begin
    writeln('root is at ',t + ln((1-w)/(1+w))/(2*s)  );
    writeln(confile,'root is at ',t + ln((1-w)/(1+w))/(2*s)  );
  end
  else
  begin
    y:=s*(x - t);
    val:= exp(y);{ partial result for building tanh}
    val:=(val - 1.0/val)/(val + 1.0/val); {tanh}
    val:=z*(val+w);
    nocomp:=false; {never not computable}
    fn1d:=val; {assign function value for return}
  end;
end; {fn1d.pas == htanfn.pas}
