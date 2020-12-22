!!>  this subroutine generates some constants 
!  Parameters :
!  xmin(xmax):the lower(upper) limit of integral
!  tol:  error range
!  hbarc: h*c/(2*pi)=197.32698
!  rm : rm: reduced mass 
!       rm=mp*mn/(mp+mn)
!  xm: match point
   module input
    implicit none
    integer,parameter::n=2000
    integer,parameter::lmax=5
    real*8,parameter::xmin=0.d0
    real*8,parameter::xmax=20.d0
    real*8,parameter::hbarc=197.327d0
    real*8,parameter::mp=1.0078d0
    real*8,parameter::mn=1.0086d0
    real*8,parameter::amu=931.49432d0
    real*8,parameter::ecm=10.0d0
    !real*8:: rm,e
 end module
!-----------------------------------------
!!>  this subroutine calculates wave function 
!  rm=mn*mp*c**2/(mn+mp)= 469.5662127985447 Mev
!  delta: controlled variable
!  h: step length
!  c: Normalized parameter
!  y: the wave function
!  dyin(dyout):the left(right) Derivative at xm
program main
 use input
 use coulfunc
 implicit none
 integer::i,l,ifail
 real*8,dimension(0:3000)::x,y,yreal
 real*8,dimension(0:3000)::k2
 complex*8,dimension(1:3)::s
 real*8,dimension(1:3)::nfc,ngc,nfcp,ngcp,c
 real*8::h,zero,eta,rho
 real*8::dyn,rm,e
 complex*8,dimension(1:3)::hin,hout
 complex*8,dimension(1:3)::dyhin,dyhout
 complex*8::b
 real*8,external::gausspot
 rm=mp*mn*amu/(mp+mn)
 e=ecm*(mp+mn)/mn
 h=(xmax-xmin)/n
 b=(0.,1.0)
 eta=0.
 zero=0.
 ifail=0
 do l=0,lmax
    k2(n)=2*rm*(e-gausspot(xmax))/hbarc**2-l*(l+1)*1./xmax**2
    rho=xmax*k2(n)
    call COUL90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
    if (ifail/=0) then
    write(*,*) 'coul90: ifail=',ifail; stop
    endif
    call left(rm,e,h,l,x,y)
    dyn=(y(n-2)-8.*y(n-1)+8.*y(n+1)-y(n+2))/(12.*h)  
    hin(l)=ngc(l)+b*nfc(l)
    hout(l)=ngc(l)-b*nfc(l)
    dyhin(l)=ngcp(l)+b*nfcp(l)
    dyhout(l)=ngcp(l)-b*nfcp(l)
    write(14,*)l,dyn,hin(l),hout(l),dyhin(l),dyhout(l)
    !  compute the s matrix
    s(l)=(y(n)*dyhout(l)-dyn*hout(l))/(y(n)*dyhin(l)-dyn*hin(l))
    !  compute c
    c(l)=b*(hout(l)-s(l)*hin(l))/(2*y(n))
    !  OUTPUT S-METRIX C & WAVE FUNCTION
    write(13,*)s(l),c(l),l
    do i=0,n
      write(11,*)x(i),y(i)  
    end do
    write(11,*)'&'
     !yreal=y*c(lmax)
     !do i=0,n
     ! write(10,*)yreal(i)  
     !end do
 end do
end program
 !------------------------------------------------------------------------
 
 !------------------------------------------------------------------------
 !
 !   Numerov algorithm
 !   LEFT
 !
 !   INPUT:
 !   e: Test eigenvalues
 !   xmathch: match point
 !   del:adjustable parameter
 !   x(i):ridus
 !   k2:  k2=k**2
 
 !   OUTPUT:
 !   y(i): wave function
 
     subroutine left(rm,e,h,l,x,y)
     use input
     implicit none
     integer::i,l
     real*8::h,rm,e
     real*8,dimension(0:3000)::x,y
     real*8,dimension(0:3000)::V,k2
     real*8,external::gausspot
     y(0)=0
     y(1)=h
     k2(1)=2*rm*(e-gausspot(h))/hbarc**2-l*(l+1)*1./h**2
     y(2)=2.*y(1)-h**2*k2(1)*y(1)
     write(12,*)l,y(2)
     write(12,*)'k(i)'
     do i=1,n
         x(i)=xmin+i*h
         V(i)=gausspot(x(i))
         k2(i)=2*rm*(e-v(i))/hbarc**2-l*(l+1)*1./x(i)**2
         write(12,*)k2(i)
     end do
     do i=2,n+2
         y(i+1)=2*(1-5*h**2*k2(i)/12.)*y(i)-(1+h**2*k2(i-1)/12.)*y(i-1)
         y(i+1)=y(i+1)/(1+h**2*k2(i+1)/12.)
     end do
     end subroutine
     
!-------------------------------------------------------------------------
!>    this subroutine calculates gauss potential.
!     c *** Gaussian Potential
     function gausspot(r)
       implicit none
        real*8 r,v0,r0,gausspot,a
        v0=-72.16d0
        r0=0.
        a=1.484d0
          if (a.gt.1e-6) then
            gausspot=V0*exp(-(r-r0)**2/a**2)
              else
                write(*,*)'a too small in gausspot!'
                stop
          endif
          return
     end function
       
!-------------------------------------------------------------------------- 
