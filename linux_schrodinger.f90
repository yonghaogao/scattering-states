!!>  this subroutine generates some points of  
!  Parameters :
!  xmin(xmax):the lower(upper) limit of integral
!  tol:  error range
!  hbarc: h*c/(2*pi)=197.32698
!  rm : rm: reduced mass 
!       rm=mp*mn/(mp+mn)
!  E: e:  Energy eigenvalue
!  xm: match point
   module para
      implicit none
      integer::n
      real*8::xmin,xmax,Tol
      real*8::hbarc,rm
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
   program schrodinger
      use para
      implicit none
      integer::i,j
      real*8,dimension(1:3000)::xl,xr,yl,yr
      real*8,dimension(1:3000)::y,y1
      real*8,dimension(1:3000)::X,W
      real*8,dimension(1:3000)::r,V
      real*8::xmathch,h
      real*8::c,E,n0
      real*8::delta,ddelta
      real*8::sum
      integer::num1,num2
      real*8::dyin,dyout
      real*8,external::FFR4,gausspot
      real*8,parameter::v0=-150d0
      n=1000
      xmin=0.
      xmax=20.
      tol=1.E-5
      e=-2.224d0
      hbarc=197.327d0
      rm=469.566213d0
      xmathch=2.0d0
      delta=1.
      h=(xmax-xmin)/(2*n)
      num1=200
      num2=1800
      call left(e,xmathch,delta,xl,yl)
      call right(e,xmathch,delta,xr,yr)
!     whole wave function
      do i=1,num1
         y(i)=yl(i)
      end do
      !n0=yl(n+1)/yr(1),for
      n0=yl(num1+1)/yr(1)
      yr=n0*yr
      do i=num1,2*n+1
         y(i)=yr(i+1-num1)
      end do
      do i=1,2*n+1
         r(i)=xmin+(i-1)*h
         V(i)=delta*gausspot(r(i))
      end do  
!     integeral & interpolation
      call  gauleg(2*n,xmin,xmax,X,W)
      sum=0.
      do i=1,2*n+1
         y1(i)=FFR4(x(i)/h,y(i),2*n)
         sum=sum+w(i)*v(i)*y1(i)**2
      end do   
!     Derivative at xm------------------------
      dyin=(yl(num1-1)-8.*yl(num1)+8.*yl(num1+2)-yl(num1+3))/(12.*h)
      dyout=(-11./6.*yr(1)+3.*yr(2)-3./2.*yr(3)+1./3.*yr(4))/h
      ddelta=y(num1+1)*(dyout-dyin)/sum
      !--------------------------------
      do while(abs(ddelta)>tol)
         ddelta=y(num1+1)*(dyout-dyin)/sum
         delta=delta+ddelta
         call left(e,xmathch,delta,xl,yl)
         call right(e,xmathch,delta,xr,yr)
         do i=1,num1+1
            y(i)=yl(i)
         end do
         !n0=yl(n+1)/yr(1),for
         n0=yl(num1+1)/yr(1)
         yr=n0*yr
         do i=num1,2*n+1
            y(i)=yr(i+1-num1)
         end do
         do i=1,2*n+1
            r(i)=xmin+(i-1)*h
            V(i)=delta*gausspot(r(i))
         end do  
!      integeral & interpolation
         sum=0.
         do i=1,2*n+1
         y1(i)=FFR4(x(i)/h,y(i),2*n)
         sum=sum+w(i)*v(i)*y1(i)**2
         end do   
         dyin=(yl(num1-1)-8.*yl(num1)+8.*yl(num1+2)-yl(num1+3))/(12.*h)
         dyout=(-11./6.*yr(1)+3.*yr(2)-3./2.*yr(3)+1./3.*yr(4))/h
      end do
      write(*,*)delta*V0
!     The normalization 
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
   
       subroutine left(e,xmathch,delta,x,y)
       use para
       implicit none
       integer::i,num1
       real*8::delta,h,xmathch
       real*8,dimension(1:2000)::x,y
       real*8,dimension(1:2000)::V,k2
       real*8::gausspot
       real*8::E
       num1=200
       h=(xmathch-xmin)/num1
       y(1)=0
       y(2)=h
       do i=1,num1+3
           x(i)=xmin+(i-1)*h
           V(i)=delta*gausspot(x(i))
           k2(i)=2*rm*abs(e-v(i))/hbarc**2
       end do
       do i=2,num1+2
           y(i+1)=2*(1-5*h**2*k2(i)/12.)*y(i)-(1+h**2*k2(i-1)/12.)*y(i-1)
           y(i+1)=y(i+1)/(1+h**2*k2(i+1)/12.)
       end do
       end subroutine
       
   !--------------------------------------------------------------------------
   !   RIGHT	
   !
   !   INPUT:
   !   e: Test eigenvalues
   !   xm: match point
   !   delta: adjustable parameter
   
   !   HETA : Sommerfeld parameter
   !   hbar: h*c/(2*pi)
   !   ec: e=-1.60218E-19
   !   OUTPUT:
   !   x(i):ridus
   !   y(i): wave function
       subroutine right(e,xmathch,delta,x,y)
       use para
       implicit none
       integer::i,num2
       integer::zero
       real*8::delta,h,xmathch
       real*8::HETA
       real*8,dimension(1:2000)::x,y,dy
       real*8,dimension(1:2000)::V,k2
       real*8::gausspot
       real*8::temp,kn,knn
       real*8::E
       heta=0.
       zero=0
       num2=1800
       h=(xmax-xmathch)/num2
       temp=xmax-h
       !heta=z1*z2*e2*niu*/(hbar*sqrt(w*niu*Energy))
       do i=1,num2+1
           x(i)=xmathch+(i-1)*h
           V(i)=delta*gausspot(x(i))
           k2(i)=2*rm*(e-v(i))/hbarc**2
       end do
       kn=sqrt(2*rm*abs(e)/hbarc**2)
       call WHIT(heta,xmax,kn,e,zero,y(num2+1),dy(num2),zero)
       call WHIT(heta,temp,kn,e,zero,y(num2),dy(num2-1),zero)
       do i=num2,2,-1
         y(i-1)=2*(1-5./12*k2(i)*h**2)*y(i)-(1+1./12*k2(i+1)*h**2)*y(i+1)
         y(i-1)=y(i-1)/(1+1./12*k2(i-1)*h**2)
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
    SUBROUTINE WHIT(HETA,R,XK,E,LL,F,FD,IE)
      !
      !     CALCULATES  WHITTAKER  FUNCTION  WL(K,R)  WITH
      !     ASYMPTOTIC  FORM  EXP(-(KR + ETA(LOG(2KR)))
      !     E  IS  NEGATIVE
      !     If IE = 0, allowed to return result e**IE larger than Whittaker,
      !                for the IE value returned.
      !     If IE > 0, must scale results by that amount.
      !
      !   input : 
      !           HETA : Sommerfeld parameter
      !           R : radius 
      !           XK: module of wavenumber in fm^{-1}
      !           E :  C.M. energy in MeV 
      !           LL :  partial wave 
      !           IE :  normally set to 0 
      !   output:
      !           F(LL+1) : WHITTAKER  FUNCTION
      !           FD(LL+1) : derivative WHITTAKER  FUNCTION
      
      
      !	use drier !  AMoro
            IMPLICIT REAL*8 (A-H,O-Z)
            DIMENSION F(LL+1),FD(LL+1) ,T(12),S(7)
      !! AMoro: to replace drier module
            REAL*8 FPMAX
      !	acc8 = epsilon(acc8);
            fpmax = huge(acc8)**0.8d0 
      !! ------------------------------
      
            L = LL+1
      !              NOW L = NO. OF VALUES TO FIND
            EE=-1.0
            AK=XK
            ETA=HETA
            LP1=L+1
            RHO=AK*R
         S(:) = 0
            IF(L-50)1,1,2
          1 LM=60
            GO TO 3
          2 LM=L+10
          3 LMP1=LM+1
            IS=7
            PJE=30.0*RHO+1.0
            H=max(INT(PJE),4)
            H=RHO/H
            RHOA=10.0*(ETA+1.0)
            IF(RHOA-RHO)13,13,14
         13 IFEQL=1
            RHOA=RHO
            GO TO 15
         14 IFEQL=0
         15 PJE=RHOA/H+0.5
            RHOA=H*INT(PJE)
            IF(IFEQL)16,16,18
         16 IF(RHOA-RHO-1.5*H)17,18,18
         17 RHOA=RHO+2.0*H
         18 IF(EE)55,55,19
         19 STOP 'WHIT'
         27 A=2.0-10.0/12.0*H*H*EE
            B=1.0/6.0*H*ETA
            C=1.0+1.0/12.0*H*H*EE
            M1=INT(RHOA/H-0.5)
            M2=INT(RHO/H-1.5)
            T(2)=B/FLOAT(M1+1)
            T(3)=B/FLOAT(M1)
            JS=M1
            DO 29 IS=M2,M1
            DO 28 I=1,6
            S(I)=S(I+1)
         28 CONTINUE
            T(1)=T(2)
            T(2)=T(3)
            T(3)=B/FLOAT(JS-1)
            S(7)=((A+10.0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))
            JS=JS-1
            IF(ABS(S(7)).LE.FPMAX) GO TO 29
             DO 285 I=2,7
        285   S(I) = S(I) / FPMAX
         29 CONTINUE
            T(1)=S(4)
            T(2)=(1.0/60.0*(S(1)-S(7))+0.15*(S(6)-S(2))+0.75*(S(3)-S(5)))/H
            GO TO 60
         55 C=1.0/RHOA
            A=1.0
            B=1.0-C*ETA
            F(1)=A
            FD(1)=B
            DO 56 M=1,26
            D=0.5*(ETA+FLOAT(M-1))*(ETA+FLOAT(M))*C/FLOAT(M)
            A=-A*D
            B=-B*D-A*C
            F(1)=F(1)+A
            FD(1)=FD(1)+B
         56 CONTINUE
            A=-ETA*LOG(2.0*RHOA)-RHOA
            FPMINL = -LOG(FPMAX)
            if(IE.eq.0.and.A.LT.FPMINL) IE = INT(FPMINL-A)
            A=EXP(A+IE)
            F(1)=A*F(1)
      !      FD(1)=A*FD(1)
            FD(1)=A*FD(1) * (-1d0 - 2*ETA/(RHOA))
            IF(IFEQL)57,57,61
         57 S(IS)=F(1)
            IF(IS-7)27,58,27
         58 IS=6
            RHOA=RHOA+H
            GO TO 55
         60 F(1)=T(1)
            FD(1)=T(2)
         61 C=1.0/RHO
            DO 63 M=1,L-1
            A=ETA/FLOAT(M)
            B=A+C*FLOAT(M)
            F(M+1)=(B*F(M)-FD(M))/(A+1.0)
            FD(M+1)=(A-1.0)*F(M)-B*F(M+1)
         63 CONTINUE
            DO 65 M=1,L
            FD(M)=AK*FD(M)
         65 CONTINUE
            RETURN
            END SUBROUTINE
         
   