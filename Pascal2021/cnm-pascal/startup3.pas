procedure startup;
{startup -- startup code modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}
begin
(*
  assign(con,'CON');
  rewrite(con);
  assign(input,''); reset(input);
  assign(output,''); rewrite(output); {needed for I/O redirection}
*)
  writeln(banner);
  {display the program banner}
  tdstamp(CON); {displays a time and date stamp}
  write('File for input of control data ([cr] for keyboard) ');
  readln(infname); if length(infname)=0 then infname:='con';
  {Gets a filename for a file with data for running a problem.}
  assign(infile, infname); reset(infile);
  {Opens this file}
  write('File for console image ([cr] = nul) ');
  readln(infile,confname);
  {Gets a name for a file to which console output is repeated.}
  if length(confname)=0 then confname:='nul';
  {Changes a carriage return to a 'nul' so console copy is omitted.}
  if infname<>'con' then writeln(confname);
  {Writes out the console file name on the screen if it comes from a file.
  Also moves to new line when reading from a file rather than keyboard.}
  assign(confile,confname); rewrite(confile);
  {Opens to console image file.}
  writeln(confile,banner);
  tdstamp(confile); {Repeat banner and time stamp -- it is a little bit
    different from the console timestamp.}
  writeln(confile,'File for input of control data ([cr] for keyboard) ',
                infname);
  writeln(confile,'File for console image ([cr] = nul) ',confname);
  {Repeat the setup information to the console image file.}
  {Now start the program proper.}
end; {startup.pas}
