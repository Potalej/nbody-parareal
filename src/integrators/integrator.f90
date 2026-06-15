! ************************************************************
!! Integrator
!
!> Objectives
!  This module contains the integrator_type, the base for all
!  integrators in the program.
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
MODULE integrator_mod
    USE configs_mod
    USE types_mod
    USE forces_mod
    IMPLICIT NONE
    PUBLIC integrator_type
    PRIVATE

    TYPE :: integrator_type
        CHARACTER(30) :: name
        REAL(pf) :: dt
        INTEGER  :: N
        REAL(pf), ALLOCATABLE :: masses(:), inv_masses(:,:)
        
        CHARACTER(10) :: forces_method_NAME
        REAL(pf)      :: forces_angle ! parameter 'theta' of Barnes-Hut
        LOGICAL       :: forces_parallel
        CLASS(forces_type), POINTER :: forces_method => NULL()

        CONTAINS
            PROCEDURE :: init, apply, method
    END TYPE

CONTAINS

SUBROUTINE init (self, m, dt, fm, fa, fp, G, softening)
    CLASS(integrator_type), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: m(:) ! masses
    REAL(pf), INTENT(IN) :: dt ! timestep
    REAL(pf), INTENT(IN) :: fa ! forces angle
    LOGICAL,  INTENT(IN) :: fp ! forces parallel
    CHARACTER(LEN=*), INTENT(IN) :: fm ! forces method
    REAL(pf), INTENT(IN) :: G, softening ! parameters of the Newtoninan potential
    LOGICAL :: use_bh
    INTEGER :: a

    CALL msg("[integrator] initializing", 2)

    self % N = SIZE(m)
    self % dt = dt
    self % forces_method_name = fm
    self % forces_angle = fa
    self % forces_parallel = fp
    
    ! masses
    ALLOCATE(self % masses(self % N))
    ALLOCATE(self % inv_masses(self % N,3))
    self % masses = m
    DO a = 1, self % N
        self % inv_masses(a,:) = 1.0_pf / m(a)
    END DO

    ! initialize the forces
    use_bh = .FALSE.
    IF (TRIM(fm) == "bh") use_bh = .TRUE.
    CALL msg("[integrator] initializing forces", 2)
    ALLOCATE(self % forces_method)
    CALL self % forces_method % init(fp, use_bh, G, softening, fa)
END SUBROUTINE

SUBROUTINE method (self, qs, ps, fs)
    CLASS(integrator_type), INTENT(INOUT) :: self
    REAL(pf), DIMENSION(:,:), INTENT(INOUT) :: qs, ps, fs

    WRITE (*,*) 'OPS'
END SUBROUTINE

SUBROUTINE apply (self, qs, ps, fs, eval_forces)
    CLASS(integrator_type), INTENT(INOUT) :: self
    REAL(pf), DIMENSION(self % N, 3), INTENT(INOUT) :: qs, ps, fs
    LOGICAL, INTENT(INOUT), OPTIONAL :: eval_forces

    IF (PRESENT(eval_forces)) THEN
        IF (eval_forces) THEN
            fs = self % forces_method % eval(self % N, self % masses, qs)
            eval_forces = .FALSE.
        ENDIF
    ENDIF

    CALL self % method(qs, ps, fs)
END SUBROUTINE

END MODULE