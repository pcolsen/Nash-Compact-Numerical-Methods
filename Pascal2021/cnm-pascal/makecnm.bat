rem makecnm.bat
rem   prepare distribution diskette for cnm 2
rem     90-1-18
rem
rem put a blank diskette in b:
pause
tamrof b:
label b:nashcnm2
b:
pkxarc c:\cnm\diskcnm readme.cnm 
md \cnm
cd \cnm
pkxarc c:\cnm\diskcnm *.*
c:
cd b:..
