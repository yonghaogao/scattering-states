!!>  this subroutine generates some constants 
!  Parameters :
!  xmin(xmax):the lower(upper) limit of integral
!  tol:  error range
!  hbarc: h*c/(2*pi)=197.32698
!  rm : rm: reduced mass 
!       rm=mp*mn/(mp+mn)
!  xm: match point
   module para
      implicit none
      integer,parameter::n=2000
      real*8,parameter::xmin=0.d0
      real*8,parameter::xmax=20.d0
      real*8,parameter::Tol=1.E-5
      real*8,parameter::hbarc=197.327d0
      real*8,parameter::rm=469.566213d0
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
   use para
   use coulfunc
   implicit none
   integer::i,l,lmax,ifail
   real*8,dimension(0:3000)::x,y,yreal
   real*8,dimension(0:3000)::R,V,k2
   complex*8,dimension(1:3)::s,c
   real*8,dimension(1:3)::nfc,ngc,nfcp,ngcp
   real*8::h,zero,eta,rho
   real*8::E
   real*8::dyn
   complex*8,dimension(1:3)::hin,hout
   complex*8,dimension(1:3)::dyhin,dyhout
   complex*8::b
   real*8,external::gausspot
   e=-2.224d0
   h=(xmax-xmin)/n
   b=(0.,1.0)
   eta=0.
   lmax=3
   zero=0.
   ifail=0
   open(10,file='data')
   do i=1,n
     r(i)=xmin+i*h
     V(i)=gausspot(r(i))
     k2(i)=2*rm*(e-v(i))/hbarc**2-l*(l+1)*1./x(i)**2
   end do
   rho=xmax*sqrt(abs(k2(n)))
   call coul90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
   if (ifail/=0) then
    write(*,*) 'coul90: ifail=',ifail; stop
   endif
   do l=1,lmax
      call left(e,h,l,x,y)
      dyn=(y(n-2)-8.*y(n-1)+8.*y(n+1)-y(n+2))/(12.*h)  
      hin(l)=ngc(l)-b*nfc(l)
      hout(l)=ngc(l)+b*nfc(l)
      dyhin(l)=ngcp(l)-b*nfcp(l)
      dyhout(l)=ngcp(l)+b*nfcp(l)
      !  compute the s matrix
      s(l)=(y(n)*dyhin(l)-dyn*hin(l))/(y(n)*dyhout(l)-dyn*hout(l))
      !  compute c
      c(l)=b*(hin(l)-s(l)*hout(l))/(2*y(n))
   end do
   !  OUTPUT S-METRIX C & WAVE FUNCTION
   do l=1,lmax
      write(10,*)s(l),c(l)
   end do
   do i=0,n
    write(10,*)x(i),y(i)  
   end do
   close(10)
   !yreal=y*c(lmax)
   !do i=0,n
   ! write(10,*)yreal(i)  
   !end do
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
   
       subroutine left(e,h,l,x,y)
       use para
       implicit none
       integer::i,l
       real*8::h
       real*8,dimension(0:3000)::x,y
       real*8,dimension(0:3000)::V,k2
       real*8::gausspot
       real*8::E
       y(0)=0
       y(1)=h
       y(2)=h**(l+1)
       do i=1,n
           x(i)=xmin+i*h
           V(i)=gausspot(x(i))
           k2(i)=2*rm*(e-v(i))/hbarc**2-l*(l+1)*1./x(i)**2
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
!!>     this subroutine calculates standard gauss Legendre points
!!>     between x1 and x2 (usually -1.0d0 and  1.0d0)
!!>     N is the number of mesh points required.
!!>     The grid and the weights are stored in the arrays X and W
!!>     @param[in] x1 lower boundary
!!>     @param[in] x2 upper boundary
!!>     @param[in] N number of grid points
!!>     @param[out] X grid points
!!>     @param[out] W integration weights
       SUBROUTINE gauleg(N,x1,x2,X,W)
         IMPLICIT NONE
         INTEGER N
         REAL*8 x1,x2,X(N),W(N)
         REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
         INTEGER m,i,j
 
         pi=acos(-1.0)
         tol=1.E-12
 
         m=(n+1)/2
         xm=0.5*(x2+x1)
         xl=0.5*(x2-x1)
 
         DO 10 i=1,m
          z=cos(pi*(i-0.25)/(N+0.5))
 
  20      CONTINUE
          p1=1.0E0
          p2=0.0E0
          DO 30 j=1,N
           p3=p2
           p2=p1
           p1=((2*j-1)*z*p2-(j-1)*p3)/j
  30      CONTINUE
          pp=N*(z*p1-p2)/(z*z-1.0E0)
          z1=z
          z=z1-p1/pp
          IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende
 
          X(i) = xm - xl*z
          X(n+1-i) = xm + xl*z
          W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
          W(n+1-i) = W(i)
  10     CONTINUE
        END SUBROUTINE gauleg
!--------------------------------------------------------------------------
     FUNCTION FFR4(Y,F,N) 
         IMPLICIT REAL*8(A-H,O-Z)
         REAL*8 F(N),P,P1,P2,Q,X,FFR4
         REAL*8 Y
         PARAMETER(X=.16666666666667)
         P=Y
         I=P
         IF(I.LE.0) GO TO 2
         IF(I.GE.N-2) GO TO 4
       1 P=P-I
         P1=P-1.
         P2=P-2.
         Q=P+1.
         FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
         RETURN
       2 IF(I.LT.0) GO TO 3
         I=1
         GO TO 1
       3 FFR4=F(1)
         RETURN
       4 IF(I.GT.N-2) GO TO 5
         I=N-3
         GO TO 1
       5 FFR4=F(N)
         RETURN
     END function

!--------------------------------------------------------------------------     
   