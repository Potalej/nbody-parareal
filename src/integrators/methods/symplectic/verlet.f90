! ************************************************************
!! (Stormer/velocity) Verlet Method
!
!> Objectives
!  Application of the (symplectic) Stormer/velocity-Verlet Method.
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
MODULE verlet_integrator_mod
    USE types_mod
    USE integrator_mod
    IMPLICIT NONE
    PUBLIC verlet_integrator_type
    PRIVATE

    TYPE, EXTENDS(integrator_type) :: verlet_integrator_type
        CONTAINS
            PROCEDURE :: method
    END TYPE

CONTAINS

SUBROUTINE method (self, qs, ps, fs)
    CLASS(verlet_integrator_type), INTENT(INOUT) :: self
    REAL(pf), DIMENSION(:,:), INTENT(INOUT) :: qs, ps, fs

    ! apply the Velocity-Verlet Method
    ps = ps + 0.5_pf * self % dt * fs
    qs = qs + self % dt * ps * self % inv_masses
    
    ! update the forces
    fs = self % forces_method % eval(self % N, self % masses, qs)

    ps = ps + 0.5_pf * self % dt * fs
END SUBROUTINE

END MODULE