function fn1d(x: real;var nocomp: boolean): real;
  {real valued test function of x for [1D] minimisation and rootfinding}
{quadfn.pas
  == This is a quadratic function (10*x - 5)*x - 8.
  Minimum is at 0.5.
  Roots are at approx. -0.6787088, 1.1787088.
}
var
  val: real;

begin
  nocomp:=false; {always set for this function}
  fn1d:=(10.0*x-5.0)*x-8.0;
end; {fn1d.pas == quadfn.pas}
