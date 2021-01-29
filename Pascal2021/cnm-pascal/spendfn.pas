function fn1d(x: real;var nocomp: boolean): real;
  {real valued test function of x for [1D] minimisation}
{spendfn.pas
  == This is the expenditure example, illustrated in
  Example 12.5 and Example 13.1 in Compact Numerical Methods.
}
var
  val: real;

begin
  nocomp:=false; {always set for this function}
  {Note: the lack of a power function is a difficulty with Pascal.
    Here we use the log and exp functions. This is NOT recommended
    for general use. See Cody and Waite (1980) for better codes.}
  val:=-250.0*x+22500.0* exp(x * ln(1.00949));
  fn1d:=val;
end; {fn1d.pas == spendfn.pas}
