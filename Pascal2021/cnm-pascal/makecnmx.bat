rem makecnmx.bat
rem   prepare distribution diskette for cnm 2
rem     90-1-18
rem
if %1==a: goto doa
if %1==b: goto dob
goto errors
:dob
rem put a blank diskette in b:
pause
tamrof b:
goto rest
:doa
rem put a blank diskette in a:
pause
tamrof a:/n:9
:rest
label %1nashcnm2
%1
pkxarc c:\cnm\diskcnm readme.cnm 
md \cnm
cd \cnm
pkxarc c:\cnm\diskcnm *.*
c:
cd %1..
goto endit
:errors
rem SOMETHING WRONG
:endit
rem done!!