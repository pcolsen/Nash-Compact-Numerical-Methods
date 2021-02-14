      program timings
      implicit none
      integer, parameter :: ntests = 1
      integer :: n 
      real(kind=8) :: t1, t2, elapsed_time, x, y
      integer(kind=8) :: tclock1, tclock2, clock_rate
      integer :: i,j,k,itest

      call system_clock(tclock1)
      
      n = 54545454
      print *, "Will run exp(sin(cos(dble(i)))) for ",n," loops"
      
      call cpu_time(t1)   ! start cpu timer
      do itest=1,ntests
          y = 0.0
          do j = 1,n
            x=exp(sin(cos(dble(j))))
            y=y+x
          enddo
      enddo
      call cpu_time(t2)   ! end cpu timer
      print *, "CPU time =",(t2-t1)
      call system_clock(tclock2, clock_rate)
      elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
      print *, "Elapsed time = ",elapsed_time
      end
