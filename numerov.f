!     module for potential 
      module numwf 
      implicit none 
      contains
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 !------------------------------------------------------------------------
 !
 !   Numerov algorithm
 !
 !   INPUT:
 !   e: Test eigenvalues
 !   xmathch: match point
 !   del:adjustable parameter
 !   x(i):ridus
 !   k2:  k2=k**2
 
 !   OUTPUT:
 !   y(i): wave function
      subroutine numerov(l,x,y)
      use inputfile 
      implicit none
      integer::i,l
      real*8,dimension(0:3000)::x
      complex*16,dimension(0:3000)::y
      complex*16,dimension(0:3000)::V,k2
      y(0)=0
      y(1)=h
      call gausspot(h,v(1))
      k2(1)=2*rm*(ecm-v(1))/hbarc**2-l*(l+1)*1./h**2
      y(2)=2.*y(1)-h**2*k2(1)*y(1)
      x(0)=0.d0
      do i=1,n
         x(i)=xmin+i*h
         call gausspot(x(i),v(i))
         k2(i)=2*rm*(ecm-v(i))/hbarc**2-l*(l+1)*1./x(i)**2
      end do
      do i=2,n+2
         y(i+1)=2*(1-5*h**2*k2(i)/12.)*y(i)-(1+h**2*k2(i-1)/12.)*y(i-1)
         y(i+1)=y(i+1)/(1+h**2*k2(i+1)/12.)
      end do
      end subroutine
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>    this subroutine calculates gauss potential.
!     c *** Gaussian Potential
      subroutine gausspot(r,out)
      implicit none
      real*8 r,v0,r0,a
      complex*16::out
      v0=-72.08512d0
      r0=0.
      a=1.484d0
       if (a.gt.1e-6) then
         out=V0*exp(-(r-r0)**2/a**2)
           else
             write(*,*)'a too small in gausspot!'
             stop
       endif
       return
      end subroutine
    
!--------------------------------------------------------------------------   
      
      end module 