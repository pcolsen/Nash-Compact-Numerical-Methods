function fn1d(x: real;var nocomp: boolean): real;
  {real valued test function of x for [1D] minimisation}
{
This is a cubic function x*(x*x - 2.) - 5.
From Forsythe, Malcolm & Moler (1977) page 184.
Minimum is at 0.81650 .
}
var
  val: real;

begin
  nocomp:=false; {always set for this function}
  fn1d:=(x*x-2.0)*x-5.0;
end; {fn1d.pas == cubefn.pas}
