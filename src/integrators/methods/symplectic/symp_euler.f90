! ************************************************************
!! Symplectic Euler Method
!
!> Objectives
!  Application of the Symplectic Euler Method.
!
!> Modified
!  2026.06.15
!
!> Created
!  2026.06.15
!
!> Author
!  oap
!
MODULE symp_euler_integrator_mod
    USE types_mod
    USE integrator_mod
    IMPLICIT NONE
    PUBLIC symp_euler_integrator_type
    PRIVATE

    TYPE, EXTENDS(integrator_type) :: symp_euler_integrator_type
        CONTAINS
            PROCEDURE :: method
    END TYPE

CONTAINS

SUBROUTINE method (self, qs, ps, fs)
    CLASS(symp_euler_integrator_type), INTENT(INOUT) :: self
    REAL(pf), DIMENSION(:,:), INTENT(INOUT) :: qs, ps, fs

    ! apply the Symplectic Euler Method
    ps = ps + self % dt * fs
    qs = qs + self % dt * ps * self % inv_masses
    
    ! update the forces
    fs = self % forces_method % eval(self % N, self % masses, qs)
END SUBROUTINE

END MODULE