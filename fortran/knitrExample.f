C  Test a computation
       do 10, i=1,1000
        y=exp(sin(cos(dble(i))))
10     continue
       write(*,100) y
 100   format('0 last value of y = ',1pe16.8)
       end