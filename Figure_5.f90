!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!	program kuramoto3DHDFZloop.f90
!
!   includes up to four body interactions and external force
!   W is in the z-direction and F=(F0*cos(sigma*t),F0*sin(sigma*t),0)
!   Equations in the rotating reference frame
!
!   Loops over F and sigma
!
!	Marcus A.M. de Aguiar - 26/02/2025
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! module defining global variables
MODULE globals
INTEGER(4), SAVE :: nosc,n
REAL*8, SAVE :: r(3),anosc,k1,k2,k3
REAL*8, SAVE :: force,sigma
REAL*8, ALLOCATABLE, SAVE :: omega(:)
END MODULE

PROGRAM kuramoto3D
USE globals
IMPLICIT REAL*8(A-H,O-Z)
INTEGER :: iseed(33)
REAL*8 :: avj
REAL*8, ALLOCATABLE :: dy(:),y(:),yscal(:)
REAL*8 :: rav(3),rsqd(3),rmodav,rmodsqd
REAL*8 :: ravold(3),rsqdold(3),rmodavold,rmodsqdold

! constants
pi = dacos(-1.0D0)
pi2 = 2.0D0*pi

! interactions
k3 = 0.0D0
force = 0.2D0
sigma = 1.0D0

nosc = 5000
delta = 0.1D0
xf = 300.0D0
xcut = 0.75D0*xf
ifq = 2
omegamod = 0.0D0
inicond = 0

! Precision and initial integration step
eps=1.D-08

! Allocate vectors
n = 3*nosc
anosc = dfloat(nosc)
ALLOCATE (dy(n),y(n),yscal(n))
ALLOCATE (omega(n))

! initialize random number generator
OPEN(UNIT=50,FILE='seed.in',STATUS='OLD')  
	READ(50,*) iseed(1),iseed(2)
CLOSE(50)
CALL RANDOM_SEED(put=iseed)
CALL RANDOM_NUMBER(aux)

! Frequency distributions
omega = 0.0D0
IF(ifq == 2) THEN
    ! Entries are Gaussian distributed around omegamod in z direction
    avj = 0.D0
    DO i=1,nosc
        ii=3*(i-1)
        CALL RANDOM_NUMBER(aux)
        CALL gaussian(delta,aux,deltaomega)
        omega(ii+3) = deltaomega
        avj = avj + deltaomega
    END DO
    ! Make sure average is exactly omegamod
    avj = avj/anosc
    DO i=1,nosc
        ii=3*(i-1)
        omega(ii+3) = omega(ii+3) - avj + omegamod
    END DO
ELSE
    ! Identical oscillators
    DO i=1,nosc
        j=3*(i-1)
        omega(j+1) = 0.D0
        omega(j+2) = 0.D0
        omega(j+3) = omegamod
    END DO
END IF

! Output files
OPEN(UNIT=17,FILE='Figure_5_A.dat',STATUS='UNKNOWN')

!loops over k1 and k2

DO l1=1,81
    k1 = 0.0D0 + (dfloat(l1)-1.0D0)*0.025D0
    ! initial conditions on the sphere
    CALL RANDOM_NUMBER(aux)
    aux0 = aux
    CALL RANDOM_NUMBER(aux)
    aux1 = aux
    DO i=1,nosc
        CALL RANDOM_NUMBER(auxa)
        CALL RANDOM_NUMBER(auxb)
        IF(inicond == 0) THEN
            phi = auxa*pi2
            ct = 2.0D0*auxb - 1.D0    ! cos(theta)
            st = dsqrt(1.D0 - ct**2)  ! sin(theta)
        ELSE
            phi = aux0 + auxa*0.01
            ct = 2.0D0*aux1 - 1.D0    ! cos(theta)
            st = dsqrt(1.D0 - ct**2)  ! sin(theta)
        END IF
        j=3*(i-1)
        y(j+1) = st*dcos(phi)
        y(j+2) = st*dsin(phi)
        y(j+3) = ct
    END DO
    htry=0.001D0
    DO l2=1,401
        k2 = 10.0D0 - 0.025D0*(dfloat(l2)-1.0D0)
        ravold = 0.D0
        rsqdold = 0.D0
        rmodavold = 0.D0
        rmodsqdold = 0.D0
        icount = 0
        aicount = 0.D0
        ! Time evolution
        x = 0.0D0
        DO WHILE (x < xf)
            DO i=1,n
                yscal(i)=y(i)
                IF(yscal(i) == 0.0) yscal(i)=0.01
            END DO
            CALL DERIV(X,Y,DY)
            CALL RKQS(y,dy,n,x,htry,eps,yscal,hdid,hnext,DERIV)
            IF(X+HNEXT.EQ.X) THEN
                print*, 'Stepsize not significant in RKDUMB.'
            END IF
            HTRY=HNEXT
            ! ensure vectors remain on the sphere
            DO i=1,nosc
                j=3*(i-1)
                aux = dsqrt(y(j+1)**2 + y(j+2)**2 + y(j+3)**2)
                y(j+1) = y(j+1)/aux
                y(j+2) = y(j+2)/aux
                y(j+3) = y(j+3)/aux
            END DO
            rmod = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
            ! compute averages for x > xcut
            if(x > xcut) then
                rav = (aicount*ravold + r)/(aicount+1.D0)
                rmodav = (aicount*rmodavold + rmod)/(aicount+1.D0)
                if(icount > 0) then
                    rsqd = (rsqdold*(aicount-1.D0) + (r-ravold)*(r-rav))/aicount
                    rmodsqd = (rmodsqdold*(aicount-1.D0) + (rmod-rmodavold)*(rmod-rmodav))/aicount
                end if
                icount = icount + 1
                aicount = dfloat(icount)
                ravold = rav
                rmodavold = rmodav
                rsqdold = rsqd
                rmodsqdold = rmodsqd
            end if
        END DO
        rsqd = dsqrt(rsqd)
        rmodsqd = dsqrt(rmodsqd)
        drmod = rmodsqd/rmod
        dmax = maxval(rsqd)
        WRITE(17,'(f7.3,f7.3,f10.7,f10.7,f10.7)')k1,k2,rmodav,rmodsqd,dmax
    END DO
END DO
CLOSE(10)

! write seed
CALL RANDOM_SEED(get=iseed)
OPEN(UNIT=50,FILE='seed.in',STATUS='OLD', POSITION='REWIND')
    WRITE (50,*) iseed(1),iseed(2)
close(50)

END PROGRAM kuramoto3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations of motion
SUBROUTINE DERIV(x,y,dy)
USE globals
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 x,y(n),dy(n),r2(3),r3(3),rf(3)
REAL*8 :: rp

r = 0.0D0
DO i=1,nosc
    j = 3*(i-1)
    r(1) = r(1) + y(j+1)
    r(2) = r(2) + y(j+2)
    r(3) = r(3) + y(j+3)
END DO
r = r/anosc
rf(1) = k1*r(1) + force
rf(2) = k1*r(2)
rf(3) = k1*r(3)
rmod2 = r(1)**2 + r(2)**2 + r(3)**2

r2 = 0.D0
DO i=1,nosc
    j = 3*(i-1)
    aux = (r(1)*y(j+1) + r(2)*y(j+2) + r(3)*y(j+3))
    r2(1) = r2(1) + aux*y(j+1)
    r2(2) = r2(2) + aux*y(j+2)
    r2(3) = r2(3) + aux*y(j+3)
END DO
r2(1) = 2.0D0*r2(1)/anosc - r(1)
r2(2) = 2.0D0*r2(2)/anosc - r(2)
r2(3) = 2.0D0*r2(3)/anosc - r(3)

r3 = rmod2*r

DO i=1,nosc
    j = 3*(i-1)
    rs = rf(1)*y(j+1) + rf(2)*y(j+2) + rf(3)*y(j+3)
    rs2 = r2(1)*y(j+1) + r2(2)*y(j+2) + r2(3)*y(j+3)
    rs3 = r3(1)*y(j+1) + r3(2)*y(j+2) + r3(3)*y(j+3)
    t1 = (rf(1) - rs*y(j+1)) + k2*(r2(1)-rs2*y(j+1)) + k3*(r3(1)-rs3*y(j+1))
    t2 = (rf(2) - rs*y(j+2)) + k2*(r2(2)-rs2*y(j+2)) + k3*(r3(2)-rs3*y(j+2))
    t3 = (rf(3) - rs*y(j+3)) + k2*(r2(3)-rs2*y(j+3)) + k3*(r3(3)-rs3*y(j+3))
    dy(j+1) = -(omega(j+3)-sigma)*y(j+2) + t1
    dy(j+2) = (omega(j+3)-sigma)*y(j+1) + t2
    dy(j+3) = t3
END DO

RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrotines from numerical recipies
SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
IMPLICIT REAL*8(A-H,O-Z)
dimension dydx(n),y(n),yscal(n)
EXTERNAL derivs
PARAMETER (NMAX=60000)
dimension yerr(NMAX),ytemp(NMAX)
PARAMETER (SAFETY=0.9D0,PGROW=-.2D0,PSHRNK=-.25D0,ERRCON=1.89D-4)
h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
errmax=0.D0
do i=1,n
    errmax=max(errmax,abs(yerr(i)/yscal(i)))
end do
errmax=errmax/eps
if(errmax.gt.1.D0)then
    h=SAFETY*h*(errmax**PSHRNK)
    if(h.lt.0.1D0*h)then
        h=.1D0*h
    endif
    xnew=x+h
    if(xnew.eq.x) then
        write(6,*) 'stepsize underflow in rkqs'
        stop
    end if
    goto 1
else
    if(errmax.gt.ERRCON)then
        hnext=SAFETY*h*(errmax**PGROW)
    else
        hnext=5.D0*h
    endif
    hdid=h
    x=x+h
    do i=1,n
        y(i)=ytemp(i)
    end do
20  return
endif
END

! Subrotine from numerical recipies
SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
IMPLICIT REAL*8(A-H,O-Z)
dimension dydx(n),y(n),yerr(n),yout(n)
EXTERNAL derivs
PARAMETER (NMAX=60000)
dimension ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),ytemp(NMAX)
PARAMETER (A2=.2D0,A3=.3D0,A4=.6D0,A5=1.D0,A6=.875D0,&
     & B21=.2D0,B31=3.D0/40.D0,B32=9.D0/40.D0,B41=.3D0,B42=-.9D0,&
     & B43=1.2D0,B51=-11.D0/54.D0,B52=2.5D0,B53=-70.D0/27.D0,&
     & B54=35.D0/27.D0,B61=1631.D0/55296.D0,B62=175.D0/512.D0,&
     & B63=575.D0/13824.D0,B64=44275.D0/110592.D0,B65=253.D0/4096.D0,&
     & C1=37.D0/378.D0,C3=250.D0/621.D0,C4=125.D0/594.D0,&
     & C6=512.D0/1771.D0,DC1=C1-2825.D0/27648.D0,&
     & DC3=C3-18575.D0/48384.D0,DC4=C4-13525.D0/55296.D0,&
     & DC5=-277.D0/14336.D0,DC6=C6-.25D0)
do i=1,n
    ytemp(i)=y(i)+B21*h*dydx(i)
end do
call derivs(x+A2*h,ytemp,ak2)
do i=1,n
    ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
end do
call derivs(x+A3*h,ytemp,ak3)
do i=1,n
    ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
end do
call derivs(x+A4*h,ytemp,ak4)
do i=1,n
    ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
end do
call derivs(x+A5*h,ytemp,ak5)
do i=1,n
    ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+ B65*ak5(i))
end do
call derivs(x+A6*h,ytemp,ak6)
do i=1,n
    yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
end do
do i=1,n
    yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
end do
return
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gaussian(a,y,x)
IMPLICIT REAL*8(A-H,O-Z)
REAL*8 a,y,x,error
REAL*8 x0,delta_x,error_x
REAL*8 F,rho,anorm

error=0.001D0
pi = dacos(-1.0D0)
b = 1.0D0/(2.0D0*a*a)
anorm = dsqrt(b/pi)
x0 = 0.0D0
F = 0.5D0*(1+erf(x0*dsqrt(b))) - y
rho = anorm*dexp(-b*x0*x0)

DO WHILE(abs(F) > error)
    delta_x = - F/rho
    x = x0 + delta_x
    F = 0.5*(1+erf(x*sqrt(b))) - y
    rho = anorm*exp(-b*x*x)
    x0 = x
END DO

RETURN
END
