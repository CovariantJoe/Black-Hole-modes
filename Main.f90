MODULE QNMFinder
    ! ========================================================================================
    !
    !       Filename:  Main.f90
    !
    !    Description:  This module contains two algorithms to obtain the Quasi Normal Modes
    !                  (QNM) of a Schwarzschild's black hole by Leaver's and WKB methods.
    !                  Behavior is controlled by the file Configs.nml
    !    
    !        Version:  1.0 Changes:
    !        Created:  06/25/2024
    !       Revision:  08/08/2024
    !       Compiler:  gfortran
    !         Author:  Covariant Joe
    !        Company: 
    !
    ! ========================================================================================
    !

    IMPLICIT NONE
    SAVE   
    CONTAINS

    RECURSIVE FUNCTION Leaver(N,s,l,i) RESULT(ans)
    COMPLEX, INTENT(IN) :: s
    REAL, OPTIONAL, INTENT(IN) :: i
    REAL, INTENT(IN) :: l,N
    REAL :: index
    COMPLEX :: ans
    
    IF (.NOT. PRESENT(i)) THEN
        index = 0
    ELSE
        index = i
    END IF

    IF (index >= N) THEN
        ans = beta(index,s,l) - alpha(index,s)*gamma(index+1.0,s)
    ELSE
        ans = beta(index,s,l) - alpha(index,s)*gamma(index+1.0,s)/Leaver(N,s,l,index+1.0)
    END IF
    
END FUNCTION
    SUBROUTINE PossibleRoots()
    REAL :: L, Resolution, Terms
    REAL, DIMENSION(2) :: Interval
    REAL, ALLOCATABLE, DIMENSION(:) :: x,y
    INTEGER :: NRoots,err, Nx, Ny, i,j,k, r
    CHARACTER (Len = 64) :: filename
    CHARACTER (Len = 80) :: errMsg
    COMPLEX :: Val, zk
    COMPLEX, ALLOCATABLE, DIMENSION(:) :: roots
    NAMELIST / LeaverMethod / L, Interval, Resolution, NRoots, filename
    Interval = [-2.0,2.0]
    OPEN(unit=10, file='Configs.nml', IOSTAT = err, status = "old",IOMSG = errMsg)
    IF (err /= 0) THEN
        WRITE(*,*) "Error opening Configs.nml"
        WRITE(*,*) errMsg
        STOP
    END IF
    READ(10,NML = LeaverMethod, IOSTAT = err, IOMSG = errMsg)
    CLOSE(10)
    IF (err /= 0) THEN
        WRITE(*,*) errMsg
        STOP
    END IF
    Terms = 50.0*(MAX(abs(Interval(1)),abs(Interval(2)))) + 30.0
    Nx = CEILING((Interval(2)-Interval(1))/resolution)
    Ny = CEILING((Interval(2)-0)/resolution)
    ALLOCATE(x(Nx))
    ALLOCATE(y(Ny))
    ALLOCATE(roots(2*NRoots))
    roots = CMPLX(Interval(2)) + CMPLX(1,1) 

    DO i = 1,Nx
        x(i) = Interval(1) + (i-1)*resolution
    END DO
    DO i = 1,Ny
        y(i) = (i-1)*resolution
    END DO

    k = 0
    DO i = 1,Nx
        DO j = 1,Ny
            Val = Leaver(Terms,CMPLX(x(i),y(j)),L)
            IF (CABS(Val) < 0.1) THEN
                zk = CMPLX(x(i),y(j))
                DO r = 1, 2*NRoots
                    IF(abs(CABS(roots(r))-CABS(zk)) < 0.01) THEN
                        EXIT
                    ELSE IF (r == 2*NRoots) THEN
                        k = k+1
                        roots(k) = zk
                    END IF 
                END DO
            END IF
        END DO
    END DO
    DEALLOCATE(x)
    DEALLOCATE(y)
    IF (k == 0) THEN
        WRITE(*,*) "Found no modes within the square domain"
        WRITE(*,*) "Try adjusting the grid resolution"
        STOP
    END IF
    CALL MullerRoots(Nroots,Terms,l,roots,k,filename)
    END SUBROUTINE

    SUBROUTINE MullerRoots(Nroots,Terms,l,posroot,k,filename)
    INTEGER, INTENT(IN) :: Nroots,k
    INTEGER :: MaxIter, i,j,r = 0
    INTEGER, DIMENSION(1) :: indx
    REAL, INTENT(IN) :: l,Terms
    COMPLEX, INTENT(INOUT), DIMENSION(:) :: posroot
    COMPLEX, DIMENSION(2) :: Temp
    COMPLEX :: z1,z2,z3, h1,h2, d1,d2, a,b,c, z4,denom
    CHARACTER(Len = 64),INTENT(IN) :: filename

    MaxIter = 30
    DO i = 1,k
        z1 = posroot(i); z2 = z1 + CMPLX(0,5e-3); z3 = z1 + CMPLX(5e-3,0)
        DO j = 1, MaxIter
            h1 = z2-z1; h2 = z3-z2
            d1 = (Leaver(Terms,z2,L) - Leaver(Terms,z1,L))/h1; d2 = (Leaver(Terms,z3,L) - Leaver(Terms,z2,L))/h2  
            a = (d2-d1)/(h2+h1); b = d2 + h2*a; c = Leaver(Terms,z3,L)
            
            Temp = [b + (b*b-4*a*c)**0.5,b - (b*b-4*a*c)**0.5]
            indx =  MAXLOC(CABS(Temp))
            denom = Temp(indx(1))
            z4 = z3 - 2.0*c/denom

            IF(CABS(z4-z3) < 1e-4 .OR. CABS(Leaver(Terms,z4,L)) < 1e-4) THEN
                WRITE(*,*) z4
                r = r+1
                posroot(r) = z4
                EXIT
            ELSE
                z1 = z2; z2 = z3; z3 = z4
            END IF
        END DO
        IF (r == Nroots) THEN
            EXIT
        END IF
    END DO
    OPEN(unit=30, file= filename, status='replace', form='unformatted',access = "stream", action='write')
                WRITE(*,*) "Saving ",r ,"found modes to ",  filename
                WRITE(30) posroot(:r)
                CLOSE(30)
    END SUBROUTINE
 
    SUBROUTINE WKBMethod()
    INTEGER :: i,j,Nx,N,err, Flag
    REAL :: L,M,xn,step = 1e-2
    CHARACTER (Len = 64) :: filename
    CHARACTER (Len = 80) :: errMsg
    COMPLEX, ALLOCATABLE, DIMENSION(:) :: wn
    NAMELIST / WKB / L,M, N, filename
    Nx = CEILING(9/step)

    OPEN(unit=10, file='Configs.nml', IOSTAT = err, status = "old",IOMSG = errMsg)
    IF (err /= 0) THEN
        WRITE(*,*) "Error opening Configs.nml"
        WRITE(*,*) errMsg
        STOP
    END IF
    READ(10,NML = WKB, IOSTAT = err, IOMSG = errMsg)
    CLOSE(10)
    IF (err /= 0) THEN
        WRITE(*,*) errMsg
        STOP
    END IF
    ALLOCATE(wn(N+1))
    DO  i = 1,Nx
        IF (abs(V((i-1)*step,l,M,1)) < 0.1) THEN
            CALL NewtonRaphson(0.0,(i-1.0)*step,l,M,Flag,xn)
            IF (Flag == 1 .AND. V(xn,l,M,2) > 0.0) THEN
                DO j = 0,n
                wn(j+1) = V(xn,l,M)/M**2 + CMPLX(0,1)*SQRT(2.0*V(xn,l,M,2))*(j+0.5)/M**2
                WRITE(*,*) CSQRT(wn(j+1))
                END DO
                OPEN(unit=30, file= filename, status='replace', form='unformatted',access = "stream", action='write')
                WRITE(*,*) "Saving modes to ",  filename
                WRITE(30) wn
                CLOSE(30)
                RETURN
            END IF
        END IF
    END DO
    WRITE(*,*) "Could not find maximum value of potential within r = 0->9"

    END SUBROUTINE

    SUBROUTINE NewtonRaphson(a,b,l,M,Flag,xn)
    INTEGER :: MaxIter,j
    INTEGER, INTENT(OUT) :: Flag
    REAL, INTENT(IN) :: a,b,l,M
    REAL :: x0
    REAL, INTENT(OUT) :: xn
    Flag = 0
    MaxIter = 30
    x0 = b
    DO j = 1, MaxIter
        xn = x0 - V(x0,l,M,1)/V(x0,l,M,2)*(1.0-2.0*M/x0)
        IF (abs(V(xn,l,M,1)) < 1e-5 .AND. V(xn,l,M,2)*(1.0-2.0*M/xn) > 0.0) THEN
            Flag = 1
            EXIT
        END IF
        x0 = xn
    END DO
    END SUBROUTINE

    FUNCTION alpha(n,s)
        REAL, INTENT(IN) :: n
        COMPLEX, INTENT(IN) :: s
        COMPLEX :: alpha
        alpha = n**2 +2.0*n*(s+1.0) +2.0*s +1.0
    END FUNCTION

    FUNCTION beta(n,s,l)
        REAL, INTENT(IN) :: n,l
        REAL :: m = 2.0
        COMPLEX, INTENT(IN) :: s
        COMPLEX :: beta
        beta = -2.0*n**2 -2.0*n*(4.0*s+1.0) -8.0*s*s -4.0*s - l*(l+1.0) + m**2 -1.0
    END FUNCTION

    FUNCTION gamma(n,s)
        REAL, INTENT(IN) :: n
        COMPLEX, INTENT(IN) :: s
        REAL :: m = 2.0
        COMPLEX :: gamma
        gamma = n**2 + 4.0*s*n +4.0*s*s -m**2
    END FUNCTION
    
    FUNCTION F(r,M,k)
    INTEGER, OPTIONAL :: k
    REAL :: r,M
    REAL :: F

    IF(.NOT. PRESENT(k)) THEN
        F = 1.0-2.0*M/r
    ELSE IF (k == 1) THEN
        F = 2.0*M/r**2
    ELSE IF (k == 2) THEN
        F = -4.0*M/r**3
    ELSE IF (k == 3) THEN
        F = 12.0*M/r**4
    END IF
    END FUNCTION

    FUNCTION Lambda(n,r0,l,M)
    INTEGER, INTENT(IN) :: n
    REAL :: a
    REAL, INTENT(IN) :: r0,l,M
    COMPLEX :: Lambda
    a = REAL(n) + 0.5
    Lambda = 1.0/csqrt(CMPLX(2.0*V(r0,l,m,2)))*(1.0/8.0*(V(r0,l,m,4)/V(r0,l,m,2))*(0.25 + a**2) &
    - 1.0/288.0*(7.0+60.0*a**2)*(V(r0,l,m,3)/V(r0,l,m,2))**2 )

    END FUNCTION

    FUNCTION Omega(n,r0,l,M)
        INTEGER, INTENT(IN) :: n
        REAL :: a
        REAL, INTENT(IN) :: r0,l,M
        COMPLEX :: Omega
        a = REAL(n) + 0.5
        Omega = a/(2.0*V(r0,l,m,2)) *(5.0/6912.0*(77.0 + 188.0*a**2) &
        *(V(r0,l,m,3)/V(r0,l,m,2))**4 - 1.0/384.0*(51.0+100.0*a**2) &
        *(V(r0,l,m,4)*V(r0,l,m,3)**2/V(r0,l,m,2)**3) + 1.0/2304.0*(67.0+68.0*a**2) &
        *(V(r0,l,m,4)/V(r0,l,m,2))**2 + 1.0/288.0*( (19.0+28.0*a**2) &
        *(V(r0,l,m,3)*V(r0,l,m,5)/V(r0,l,m,2)**2) - (5.0 +4.0*a**2)*(V(r0,l,m,6)/V(r0,l,m,2))))
    
        END FUNCTION

    FUNCTION V(r,l,M,k)
    REAL:: l,M,lam
    INTEGER, OPTIONAL :: k
    REAL :: r, V
    lam = l*(l+1.0)
    IF(.NOT. PRESENT(k)) THEN
        V = F(r,M)*(l*(l+1.0)/r**2 - 6.0*M/r**3)
    ELSE IF (k == 1) THEN
        V = -6.0*F(r,M)*(-lam*r**2*1.0/3.0 + M*r*(lam+3.0) -8.0*M**2)/r**5
    ELSE IF (k == 2) THEN
        V = -60.0*F(r,M)*(0.1*lam*r**3 -2.0*M*r**2*(lam + 9.0/5.0)/3.0 +(lam+7.0)*r*M**2 - 0.2*48.0*M**3)/r**7
    ELSE IF (k == 3) THEN
        V = ( (24.0*l**2 + 24.0*l)*r**2 - 120.0*M*(l**2 + l +3.0)*r + 1440.0*M**2)/r**7
    ELSE IF (k == 4) THEN
        V = ( (-120.0*l**2 - 120.0*l)*r**2 + 720.0*M*(l**2 + l +3.0)*r - 10080.0*M**2)/r**8
    ELSE IF (k == 5) THEN
        V = ( (720.0*l**2 + 720.0*l)*r**2 - 5040.0*M*(l**2 + l +3.0)*r + 80640.0*M**2)/r**9
    ELSE IF (k == 6) THEN
        V = ( (-5040.0*l**2 - 5040.0*l)*r**2 + 40320.0*M*(l**2 + l +3.0)*r - 725760.0*M**2)/r**10
    END IF
    END FUNCTION
END MODULE 

PROGRAM MainProgram
    USE QNMFinder
    USE ZerilliSolver
    IMPLICIT NONE
    CHARACTER :: Mode
    REAL :: t0, tf
    WRITE(*,*) "Call Leaver's method (l), WKB (w) or the solver for the zerilli equation (z)?"
    READ(*,*) Mode
    CALL cpu_time(t0)
    SELECT CASE(Mode)
        CASE('L','l')
            CALL PossibleRoots
        CASE('W','w')
            CALL WKBMethod
        CASE('Z','z')
            Call Solver(V)
        END SELECT
    CALL cpu_time(tf)
    WRITE(*,*) "tiempo "
    WRITE(*,*) tf-t0
END PROGRAM
