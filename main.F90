!!>  this subroutine generates some points of  
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
      real*8,parameter::xmin=0.d0,xmax=20.d0
      real*8,parameter::Tol=1.E-5
      real*8,parameter::hbarc=197.327d0
      real*8,parameter::rm=469.566213d0
   end module
!-----------------------------------------
!!>  this subroutine calculates standard gauss Legendre points 
!  rm=mn*mp*c**2/(mn+mp)= 469.5662127985447 Mev
!  delta: controlled variable
!  h: step length
!  c: Normalized parameter
!  yl:wave function (left)
!  yr:wave function (right)
!  y: the whole wave function
!  sum: the integral of the wave function
!  dyin(dyout):the left(right) Derivative at xm
program main
   use para
   implicit none
   integer::i,j,l
   real*8,dimension(1:3000)::xl,yl
   real*8,dimension(1:3000)::y
   real*8,dimension(1:3000)::X,W
   real*8,dimension(1:3000)::r,V
   real*8::xmathch,h
   real*8::c,E,n0
   real*8::delta,ddelta
   real*8::sum
   real*8::dyin,dyout
   real*8,external::FFR4,gausspot
   real*8,parameter::v0=-72.25d0
   e=-2.224d0
   xmathch=2.0d0
   delta=1.
   h=(xmax-xmin)/n
   do l=1,3
      call left(e,h,xmathch,delta,l,x,y)

     
      !  The normalization 
      c=0.
      do i=1,2*n
         r(i)=xmin+(i-1)*h
      end do   
      do i=1,2*n
         c=c+(r(i+1)-r(i))/2.*(y(i)**2+y(i+1)**2)
      end do   
      do i=1,2*n+1
          r(i)=xmin+(i-1)*h
          y(i)=y(i)/sqrt(c)
          write(10,*) r(i),y(i)
      end do
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
   
       subroutine left(e,h,xmathch,delta,l,x,y)
       use para
       implicit none
       integer::i,l
       real*8::delta,h,xmathch
       real*8,dimension(0:2000)::x,y
       real*8,dimension(0:2000)::V,k2
       real*8::gausspot
       real*8::E
       y(0)=0
       y(1)=h
       y(2)=h**(l+1)
       do i=0,n
           x(i)=xmin+i*h
           V(i)=gausspot(x(i))
           k2(i)=2*rm*(l*(l+1)*1./x(i)**2+e-v(i))/hbarc**2
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
          v0=-150d0
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
   