{********TDSTAMP.PAS***********************************}
procedure leadzero(value:word; var outstr: str2);
{
          Copyright 1990 J.C.Nash
}
begin   {code needed by tdstamp.pas}
  if value<10 then
  begin
    str(value:1,outstr);
    outstr:='0'+outstr;
  end else str(value:2,outstr);
end; {leadzero}
{******************************************************}
procedure tdstamp(var outfile : text);
{tdstamp.pas == writes a time/date stamp to outfile

          Modified for Turbo Pascal 5.0

          Copyright 1988, 1990 J.C.Nash
}
 var
  tdstr: string[20];
  year: word;
  month: word;
  day: word;
  dayofweek: word;
  yrstr: string[4];
  mnstr: str2;
  daystr: str2;
  hrstr : str2;
  minstr: str2;
  secstr: str2;
  hour: word;
  minute: word;
  second: word;
  sec100: word;

begin
  getdate(year, month, day, dayofweek);
  gettime(hour, minute, second, sec100);
  writeln(year,month,day,hour,minute,second);{??}
  if year<yearwrit then
  begin
    writeln('Program senses that clock is not up to date.');
    readln;
    halt;
  end;
  str(year:4,yrstr);
  leadzero(month,mnstr);
  leadzero(day,daystr);
  leadzero(hour,hrstr);
  leadzero(minute,minstr);
  leadzero(second,secstr);
  tdstr:=yrstr+'/'+mnstr+'/'+daystr+' '+hrstr+':'+minstr+':'+secstr;
  writeln(outfile,tdstr);
end; {tdstamp.pas}
