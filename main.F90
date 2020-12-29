!-----------------------------------------
!!>  this subroutine calculates wave function 
!  rm=mn*mp*amu/(mn+mp)= 469.5662127985447 Mev
!  delta: controlled variable
!  h: step length
!  c: Normalized parameter
!  y: the wave function
!  dyin(dyout):the left(right) Derivative at xm
       program npscattering
        use inputfile 
        use coulfunc
        use numwf
        use comwf

        implicit none
        real :: start, finish
        call cpu_time(start)
        call input()
        call compute()
      !  call output()
        
        WRITE(*,*) 'smoothie version Date: ',VERDATE
        write(*,*) "smoothie git version: ", VERREV
        write(*,*) "smoothie compile Date:", COMPDATE
        call cpu_time(finish)
        print '("Time = ",f6.3," seconds.")',finish-start
       end program
 !------------------------------------------------------------------------
 
 
