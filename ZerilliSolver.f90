MODULE ZerilliSolver
    ! =====================================================================================
    !
    !       Filename:  ZerilliSolver.f90
    !
    !    Description:  This module contains the solver for the Zerilli equation
    !                  u_rr + (w**2 - V(r)) u = 0 by Runge Kuta 4th order.
    !                  Behaviour is controlled by the file Configs.nml
    !    
    !        Version:  1.0 Changes:
    !        Created:  06/05/2024
    !       Revision:  06/28/2024
    !       Compiler:  gfortran
    !         Author:  Covariant Joe
    !        Company:
    !
    ! =====================================================================================
    !

    IMPLICIT NONE
    SAVE
    
    CONTAINS
    SUBROUTINE Solver(V)
    DOUBLE PRECISION  :: xf,x0
    DOUBLE PRECISION :: dx,u_0,z_0
    REAL, EXTERNAL :: V
    REAL :: L,M
    INTEGER :: Ln,err,i
    DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:) :: x
    COMPLEX(kind = 8) , ALLOCATABLE, DIMENSION(:) :: Phi
    COMPLEX(kind = 8) :: w
    CHARACTER (Len = 64) :: filename
    CHARACTER (Len = 80) :: errMsg
    NAMELIST / ODE / L,M,x0,xf,dx,u_0,z_0,w,Filename
    OPEN(unit=10, file='Configs.nml', IOSTAT = err, status = "old",IOMSG = errMsg)
    IF (err /= 0) THEN
        WRITE(*,*) "Error openning Configs.nml"
        WRITE(*,*) errMsg
        STOP
    END IF
    READ(10,NML = ODE, IOSTAT = err, IOMSG = errMsg)
    CLOSE(10)
    IF (err /= 0) THEN
        WRITE(*,*) errMsg
        STOP
    END IF
    Ln = CEILING((xf-x0)/dx)
    ALLOCATE(x(Ln))
    ALLOCATE(Phi(Ln))
    DO i = 1,Ln
        x(i) = x0 + (i-1)*dx
    END DO
        Phi(1) = u_0;

        CALL RK45(x,dx,w,Phi,z_0,V,L,M)
        
    OPEN(unit=30, file= filename, status='replace', form='unformatted',access = "stream", action='write')
    WRITE(*,*) "Saving solution to ",  filename
    WRITE(30) Phi
    CLOSE(30)
    END SUBROUTINE


    SUBROUTINE RK45(x,dx,w,u,z_0,V,L,M)
        REAL, EXTERNAL :: V
        REAL :: L,M
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: x
        COMPLEX(kind = 8) , INTENT(INOUT), DIMENSION(:) :: u
        COMPLEX(kind = 8), ALLOCATABLE, DIMENSION(:) :: z
        DOUBLE PRECISION, INTENT(IN) :: dx,z_0
        COMPLEX(kind = 8) :: w
        INTEGER :: k
        COMPLEX(kind = 8) :: M1,M2,M3,M4
        ALLOCATE(z(SIZE(x)))
        z(1) = z_0; u(1) = 1

        DO k = 2,SIZE(x)
            M1 = (-V(x(k-1),L,M) + w**2)*u(k-1)
            M2 = (-V(x(k-1)+dx*0.5,L,M) + w**2)*(u(k-1) + 0.5*M1*dx)
            M3 = (-V(x(k-1)+dx*0.5,L,M) + w**2)*(u(k-1) + 0.5*M2*dx)
            M4 = (-V(x(k-1)+dx,L,M) + w**2)*(u(k-1) + M3*dx)
            z(k) = z(k-1) + dx*(M1 +2.0*M2 + 2.0*M3 + M4)/6.0

            M1 = z(k-1)
            M2 = z(k-1) + 0.5*M1*dx
            M3 = z(k-1) + 0.5*M2*dx
            M4 = z(k-1) + M3*dx
            u(k) = u(k-1) + dx*(M1 +2.0*M2 + 2.0*M3 + M4)/6.0
        END DO
        DEALLOCATE(z)
    END SUBROUTINE
END MODULE
