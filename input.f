      module inputfile 
!   defines input for potential, grid, and t-matrix
!!>  this subroutine generates some constants 
!  Parameters :
!  xmin(xmax):the lower(upper) limit of integral
!  tol:  error range
!  hbarc: h*c/(2*pi)=197.32698
!  rm : rm: reduced mass 
!       rm=mp*mn/(mp+mn)
!  xm: match point
        implicit none
        integer::n,lmax
        real*8,parameter::hbarc=197.327d0
        real*8,parameter::mp=1.0078d0
        real*8,parameter::mn=1.0086d0
        real*8,parameter::amu=931.49432d0
        real*8::xmin,xmax
        real*8::rm,ecm,h,k,elab
      contains
!***********************************************************************
      subroutine input()
        implicit none
        namelist /global/ n,lmax,xmin,xmax,elab 
      !   namelist /grids/ mp,mn,hbarc
        n=600;lmax=5;xmin=0.d0;xmax=60.0d0;elab=10.0d0
        read (5,nml=global)
        write (6,nml=global)
        
      !   mp=1.0078d0;mn=1.0086d0;hbarc=197.327d0
      !   read (5,nml=grids)
      !   write (6,nml=grids)
        rm=mp*mn*amu/(mp+mn)
        write(*,*)'the reduced mass =',rm
        ecm=elab*mn/(mp+mn)
        write(*,*)'Energy of the central mass system =',ecm
        h=(xmax-xmin)/n
        write(*,*)'the step length =',h
        k=sqrt(2*rm*ecm)/hbarc      

        write (6,*)
        write (6,*) '============  Output  ==================='
        write (6,*)
      end subroutine
***********************************************************************
      end module 