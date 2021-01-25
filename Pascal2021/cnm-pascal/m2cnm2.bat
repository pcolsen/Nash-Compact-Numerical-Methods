rem m2cnm2.bat -- make cnm2 diskette on m2 machine
rem  and put on diskette in drive %1
f:
cd \
cd cnm2
cd cnmfiles
rem PUT EMPTY FORMATTED DISKETTE IN DRIVE %1
pause
md %1cnm
%1
cd cnm
lharc x f:cnmdisk
bac readme*.* ..
bac errata*.* ..
cd ..
f:
optune %1 /op
d %1 /6
d %1\cnm /6
copy %1*.* nul
copy %1\cnm\*.* nul
rem done!
