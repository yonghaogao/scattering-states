!     module for potential 
      module comwf
      implicit none 
     
      contains
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
      subroutine compute()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use inputfile 
      use coulfunc
      use numwf
      implicit none
      integer::i,l,ifail
      real*8::zero,eta,rho
      complex*16::b,dyn
      real*8,dimension(0:3000)::x
      real*8,dimension(0:5)::nfc,ngc,nfcp,ngcp
      complex*16,dimension(0:3000)::y,yreal
      complex*16,dimension(0:5)::s,c
      complex*16,dimension(0:5)::hin,hout
      complex*16,dimension(0:5)::dhin,dhout
      b=(0.,1.0)
      eta=0.
      zero=0.
      ifail=0
      rho=xmax*k
      write(*,*)'ARGUMENT FOR COULOMB FUNCTIONS',rho
      write(*,*)'the SOMMERFELD PARAMETER',eta
      write(*,*)'the ifail =',ifail
      call COUL90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
      if (ifail/=0) then
      write(*,*) 'coul90: ifail=',ifail; stop
      endif
      write(*,*)'OUTPUT S-METRIX'
      do l=0,lmax
          call numerov(l,x,y)
          dyn=(y(n-2)-8.*y(n-1)+8.*y(n+1)-y(n+2))/(12.*h)  
          hin(l)=ngc(l)+b*nfc(l)
          hout(l)=ngc(l)-b*nfc(l)
          dhin(l)=ngcp(l)+b*nfcp(l)
          dhout(l)=ngcp(l)-b*nfcp(l)
          s(l)=(k*y(n)*dhout(l)-dyn*hout(l))/(k*y(n)*dhin(l)-dyn*hin(l))
          c(l)=b*(hout(l)-s(l)*hin(l))/(2*y(n))
          yreal=y*c(l)
          write(*,*)real(s(l)),aimag(s(l)),l
          do i=0,n
               write(11,*)x(i),real(yreal(i)),aimag(yreal(i))  
          end do
          write(11,*)'&'
      end do
      end subroutine

      
    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      
      
      
      end module 