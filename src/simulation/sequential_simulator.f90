! ************************************************************
!! Sequential in time simulator
!
!> Objectives
!  This module contains a sequential in time simulator, ie, a
!  classical N-body simulator.
!
!> Modified
!  2026.06.15
!
!> Author
!  oap
!
MODULE sequential_simulator_mod
    USE configs_mod
    USE types_mod
    USE integrator_mod
    USE integrators_mod
    USE utils_mod
    USE omp_lib
    IMPLICIT NONE
    PUBLIC sequential_simulation_type
    PRIVATE

    TYPE :: sequential_simulation_type
        INTEGER :: N, checkpoints
        REAL(pf) :: G, softening, t, dt, tf
        CLASS(integrator_type), POINTER :: integrator
        REAL(pf), ALLOCATABLE, DIMENSION(:)   :: masses
        REAL(pf), ALLOCATABLE, DIMENSION(:,:) :: qs0, ps0
        REAL(pf), ALLOCATABLE, DIMENSION(:,:) :: qs, ps, fs

        CONTAINS
            PROCEDURE :: init, run
    END TYPE
CONTAINS

SUBROUTINE init (self, N, m, qs0, ps0, dt, tf, method, checkpoints, G, fsoft, fangle, fpar, fmethod)
    CLASS(sequential_simulation_type), INTENT(INOUT) :: self
    ! simulation parameters
    INTEGER,      INTENT(IN) :: N, checkpoints
    REAL(pf),     INTENT(IN) :: dt, tf
    CHARACTER(*), INTENT(IN) :: method
    ! forces parameters
    REAL(pf),     INTENT(IN) :: G, fsoft, fangle
    LOGICAL,      INTENT(IN) :: fpar
    CHARACTER(*), INTENT(IN) :: fmethod
    ! initial values
    REAL(pf), DIMENSION(N),   INTENT(IN) :: m
    REAL(pf), DIMENSION(N,3), INTENT(IN) :: qs0, ps0

    self % N = N
    self % dt = dt
    self % tf = tf
    self % t = 0
    self % G = G
    self % softening = fsoft
    self % checkpoints = checkpoints

    ALLOCATE(self % masses(N))
    self % masses = m

    ALLOCATE(self % qs0(N, 3))
    ALLOCATE(self % ps0(N, 3))
    ALLOCATE(self % qs(N, 3))
    ALLOCATE(self % ps(N, 3))
    self % qs0 = qs0
    self % ps0 = ps0
    self % qs = qs0
    self % ps = ps0

    ALLOCATE(self % fs(N, 3))

    ! integrator
    CALL msg('[seq. sim.] defining the method', 2)
    CALL define_method(self % integrator, method)
    CALL msg('[seq. sim.] initializing the integrator', 2)
    CALL self % integrator % init(m, dt, fmethod, fangle, fpar, G, fsoft)
END SUBROUTINE

SUBROUTINE run (self, output_file_par)
    CLASS(sequential_simulation_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN), OPTIONAL :: output_file_par
    
    INTEGER :: output_file
    LOGICAL :: eval_forces

    REAL(pf) :: init_E, init_J(3), init_P(3), error_E, error_J, error_P
    REAL(pf) :: timer_0

    INTEGER :: total_steps, number_substeps, substep

    output_file = -1
    IF (PRESENT(output_file_par)) output_file = output_file_par

    total_steps = NINT(ABS(self % tf / self % dt))
    number_substeps = NINT(1.0_pf * total_steps / self % checkpoints)

    init_E = total_energy(self%masses, self%qs0, self%ps0, self % G, self % softening)
    init_J = total_linear_momentum(self%ps0)
    init_P = total_angular_momentum(self%qs0, self%ps0)

    self % t = 0
    self % qs = self % qs0
    self % ps = self % ps0
    
    eval_forces = .TRUE.
    substep = 0
    timer_0 = omp_get_wtime()
    DO WHILE (self % t < self % tf)
        self % t = self % t + self % dt

        CALL self % integrator % apply(self % qs, self % ps, self % fs, eval_forces)
        IF (output_file > 0) WRITE(output_file,*) self % qs, self % ps

        substep = substep + 1
        IF (substep == number_substeps) THEN
            substep = 0
            error_E = ABS(total_energy(self%masses, self%qs, self%ps, self % G, self % softening) - init_E)
            error_J = NORM2(total_linear_momentum(self%ps) - init_J)
            error_P = NORM2(total_angular_momentum(self%qs, self%ps) - init_P)

            WRITE(*,"(A4,E12.5, A17,E12.3)") "t = ", self % t, " / time elapsed: ", omp_get_wtime() - timer_0
            WRITE(*,'(A15," = ",ES12.4)') '||E - E_0||', error_E
            WRITE(*,'(A15," = ",ES12.4)') '||J - J_0||', error_J
            WRITE(*,'(A15," = ",ES12.4)') '||P - P_0||', error_P
            WRITE(*,*)
        ENDIF
    END DO

    WRITE(*,"(A,E12.3)") "Simulation complete! Total time:", omp_get_wtime() - timer_0

    error_E = ABS(total_energy(self%masses, self%qs, self%ps, self % G, self % softening) - init_E)
    error_J = NORM2(total_linear_momentum(self%ps) - init_J)
    error_P = NORM2(total_angular_momentum(self%qs, self%ps) - init_P)
    WRITE(*,'(A15," = ",ES12.4)') '||E - E_0||', error_E
    WRITE(*,'(A15," = ",ES12.4)') '||J - J_0||', error_J
    WRITE(*,'(A15," = ",ES12.4)') '||P - P_0||', error_P
END SUBROUTINE

END MODULE