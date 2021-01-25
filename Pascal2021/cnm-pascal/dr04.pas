{$I turbo.cnm}
Program runsvd(input,output);
{dr04.pas == driver for singular value decomposition and least squares
          solution using row-wise entry of data

          Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}
{$I constype.def}
{$I tdstamp.pas}  {time and date stamp}
{$I calceps.pas}  {compute machine precision}
{$I getobsn.pas}  {allows entry of an observation (row of problem)}
{$I alg04.pas}    {NashSVD by column rotations}
{$I startup.pas}

{main program}
var
  i, j, k, m, n, nRHS : integer; {order of problem, number of right hand sides}
  Vtranspose : rmatrix;
  rssvec : rvector;  {to hold the residual sum of squares}
  svs : rvector; {to hold the singular values}
  B : rmatrix; {the matrix of least squares solutions }
  A : rmatrix; {needed only for svdtst}
  fail, sym : boolean;
  t1: real;

begin
  banner:='dr04.pas -- run Algorithm 4 problems -- Givens'+chr(39)+
                  ' reduction,';
  startup;
  GivSVD( n, nRHS, B, rssvec, svs, Vtranspose, m);
  flush(confile); close(confile);
  if infname<>'con' then close(infile);
end. {dr04.pas == Givens' reduction, svd and least squares soln}
