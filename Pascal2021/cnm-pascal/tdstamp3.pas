procedure tdstamp(var outfile : text);
{tdstamp.pas == writes a time/date stamp to outfile

          Copyright 1988 J.C.Nash
}
type
  DateStr = string[10];
  TimeStr = string[10];
  regpack = record
          ax,bx,cx,dx,bp,si,di,ds,es,flags: integer;
        end;

var
  recpack:       regpack;                {record for MsDos call}
  month,day:     string[2];
  year:          string[4];
  dx,cx:         integer;
  hour,min,sec : string[2];
  date : DateStr;
  time : TimeStr;

begin
  with recpack do
  begin
    ax := $2a shl 8;
  end;
  MsDos(recpack);                        { call function }
  with recpack do
  begin
    str(cx,year);                        {convert to string}
    str(dx mod 256,day);                 {convert to string}
    if length(day)<2 then day:='0'+day;
    str(dx shr 8,month);                 {convert to string}
    if length(month)<2 then month:='0'+month;
  end;
  date := year+'/'+month+'/'+day;
  with recpack do
  begin
    ax := $2C shl 8;
  end;
  intr($21,recpack);                     {call interrupt}
  with recpack do
  begin
    str(cx shr 8,hour);
    if length(hour)<2 then hour:='0'+hour;
    str(cx and $FF,min);
    if length(min)<2 then min:='0'+min;
    str(dx shr 8,sec);
    if length(sec)<2 then sec:='0'+sec;
  end;
  time:=hour+':'+min+':'+sec;
  writeln(outfile,date,'   ',time);
end; {tdstamp.pas}
