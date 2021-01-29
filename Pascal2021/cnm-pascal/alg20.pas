procedure axissrch(n: integer;
               var Bvec: rvector;
               var Fmin: real;
               var lowerfn: boolean;
                   Workdata: probdata);

var
 cradius, eps, f, fplus, step, temp, tilt : real;
 i : integer;
 notcomp : boolean;

begin
  writeln('alg20.pas -- axial search');
  writeln(confile,'alg20.pas -- axial search');
  eps := calceps;
  eps := sqrt(eps);
  writeln('  Axis':6,' Stepsize  ':14,'function +  ':14,
        'function -  ':14,' rad. of curv.':14,'   tilt');
  writeln(confile,'  Axis':6,' Stepsize  ':14,'function +  ':14,
        'function -  ':14,' rad. of curv.':14,'   tilt');
  lowerfn := false;
  for i := 1 to n do
  begin
    if (not lowerfn) then
    begin
      temp := Bvec[i];
      step := eps*(abs(temp)+eps);
      Bvec[i] := temp+step;
      f := fminfn(n, Bvec, Workdata, notcomp);
      if notcomp then f := big;
      write(i:5,' ',step:12,'  ',f:12,'  ');
      write(confile,i:5,' ',step:12,'  ',f:12,'  ');
    end;
    if f<fmin then lowerfn := true;
    if (not lowerfn) then
    begin
      fplus := f;
      Bvec[i] := temp-step;
      f := fminfn(n,Bvec,Workdata,notcomp);
      if notcomp then f := big;
      write(f:12,'  ');  write(confile,f:12,'  ');
    end;
    if f<fmin then lowerfn := true;
    if (not lowerfn) then
    begin
      Bvec[i] := temp;

      temp := 0.5*(fplus-f)/step;
      fplus := 0.5*(fplus+f-2.0*fmin)/(step*step);

      if fplus<>0.0 then
      begin
        cradius := 1.0+temp*temp;
        cradius := cradius*sqrt(cradius)/fplus;
      end
      else
        cradius := big;
      tilt := 45.0*arctan(temp)/arctan(1.0);
      write(cradius:12,'  ',tilt:12);
      write(confile,cradius:12,'  ',tilt:12);
    end;
    writeln; writeln(confile);
  end;
end;
