procedure Getobsn( n: integer;
            nRHS: integer;
            var W : rmatrix;
            k: integer; {row into which data is loaded}
            endflag: real;
            var enddata: boolean); {set true when end of data reached}

{Getobsn.pas
  -- reads one row of data for a least squares problem from console
    into row k of a working array W. The boolean variable enddata is
    set true when an end of date flag is encountered.
    This procedure assumes that the order of the problem, number of
    right hand sides, etc. are already defined, and that the datafile
    is already open.
}
var
  i : integer;

begin
  begin
    for i:=1 to (n+nRHS) do read(infile,W[k,i]);
    if W[k,1]=endflag then enddata:=true else enddata:=false;
  end;
end; {Getobsn.pas}
