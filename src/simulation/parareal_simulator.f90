! ************************************************************
!! Parareal
!
!> Objectives
!  This module contains the Parareal simulator.
!
!> Modified
!  2026.06.18
!
!> Created
!  2026.06.15
!
!> Author
!  oap
!
MODULE parareal_simulator_mod
    USE configs_mod
    USE types_mod
    USE integrator_mod
    USE integrators_mod
    USE utils_mod
    USE omp_lib
    IMPLICIT NONE
    PUBLIC parareal_simulation_type
    PRIVATE

    TYPE :: parareal_simulation_type
        INTEGER  :: N
        REAL(pf) :: G, softening, tf
        REAL(pf), ALLOCATABLE, DIMENSION(:)   :: masses
        REAL(pf), ALLOCATABLE, DIMENSION(:,:) :: qs0, ps0
        REAL(pf), ALLOCATABLE, DIMENSION(:,:) :: qs, ps

        ! Parareal configs
        INTEGER  :: number_of_windows, N_itermax
        REAL(pf) :: tolerance
        LOGICAL  :: fine_parallel

        CLASS(integrator_type), POINTER :: coarse_method => NULL()
        INTEGER :: coarse_dt_mult
        
        CLASS(integrator_type), POINTER :: fine_method => NULL()
        INTEGER :: fine_dt_mult

        CONTAINS
            PROCEDURE :: init, run
            PROCEDURE :: apply_coarse, apply_fine
    END TYPE

CONTAINS

SUBROUTINE init (self, N, m, qs0, ps0, & ! initial values
                tf, G, & ! simulation configs
                num_of_windows, N_itermax, tolerance, fine_parallel, & ! parareal configs
                co_method, co_dt_mult, co_fm, co_fa, co_fs, co_fp, & ! coarse
                fi_method, fi_dt_mult, fi_fm, fi_fa, fi_fs, fi_fp &  ! fine
                )
    CLASS(parareal_simulation_type), INTENT(INOUT) :: self
    ! simulation parameters
    INTEGER,      INTENT(IN) :: N
    REAL(pf),     INTENT(IN) :: tf, G
    ! initial values
    REAL(pf), DIMENSION(N),   INTENT(IN) :: m
    REAL(pf), DIMENSION(N,3), INTENT(IN) :: qs0, ps0
    ! parareal configs
    INTEGER,  INTENT(IN) :: num_of_windows, N_itermax
    REAL(pf), INTENT(IN) :: tolerance
    LOGICAL,  INTENT(IN) :: fine_parallel ! fine parallel
    ! coarse method parameters
    CHARACTER(*), INTENT(IN) :: co_method, co_fm ! coarse method, coarse forces method
    INTEGER,      INTENT(IN) :: co_dt_mult ! coarse how many steps
    REAL(pf),     INTENT(IN) :: co_fa ! coarse forces angle
    REAL(pf),     INTENT(IN) :: co_fs ! coarse forces softening
    LOGICAL,      INTENT(IN) :: co_fp ! coarse forces parallel
    ! fine method parameters
    CHARACTER(*), INTENT(IN) :: fi_method, fi_fm ! fine method, fine forces method
    INTEGER,      INTENT(IN) :: fi_dt_mult ! fine how many steps
    REAL(pf),     INTENT(IN) :: fi_fa ! fine forces angle
    REAL(pf),     INTENT(IN) :: fi_fs ! fine forces softening
    LOGICAL,      INTENT(IN) :: fi_fp ! fine forces parallel

    ! local variables
    REAL(pf) :: coarse_dt, fine_dt

    self % N = N
    self % tf = tf
    
    ! initial values
    ALLOCATE(self % masses(N))
    self % masses = m
    ALLOCATE(self % qs0(N,3))
    ALLOCATE(self % qs(N,3))
    self % qs0 = qs0
    self % qs = qs0
    ALLOCATE(self % ps0(N,3))
    ALLOCATE(self % ps(N,3))
    self % ps0 = ps0
    self % ps = ps0
    self % G = G
    self % softening = co_fs

    ! parareal configs
    self % number_of_windows = num_of_windows
    self % N_itermax = N_itermax
    self % tolerance = tolerance
    self % fine_parallel = fine_parallel

    ! coarse method configs
    CALL msg('[parareal] defining the coarse method', 2)
    ALLOCATE(self % coarse_method)
    CALL define_method(self % coarse_method, co_method)
    CALL msg('[parareal] initializing the coarse integrator', 2)
    self % coarse_dt_mult = co_dt_mult
    coarse_dt = tf / (co_dt_mult * num_of_windows)
    CALL self % coarse_method % init(m, coarse_dt, co_fm, co_fa, co_fp, G, co_fs)

    ! fine method configs
    CALL msg('[parareal] defining the fine method', 2)
    ALLOCATE(self % fine_method)
    CALL define_method(self % fine_method, fi_method)
    CALL msg('[parareal] initializing the fine integrator', 2)
    self % fine_dt_mult = fi_dt_mult
    fine_dt = tf / (fi_dt_mult * num_of_windows)
    CALL self % fine_method % init(m, fine_dt, fi_fm, fi_fa, fi_fp, G, fi_fs)

    CALL msg('[parareal] initialized', 2)
END SUBROUTINE

SUBROUTINE run (self, output_file_par)
    CLASS(parareal_simulation_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN), OPTIONAL :: output_file_par
    INTEGER :: output_file
    REAL(pf), ALLOCATABLE :: U0(:,:)
    REAL(pf), ALLOCATABLE :: U(:,:,:,:)
    REAL(pf), ALLOCATABLE :: U_tilde(:,:,:)
    REAL(pf), ALLOCATABLE :: U_comp(:,:,:)
    INTEGER  :: i, iterations
    REAL(pf) :: error
    INTEGER :: N

    REAL(pf) :: timer_0, timer_1, timer_total
    REAL(pf) :: init_E, init_J(3), init_P(3), error_E, error_J, error_P
    REAL(pf), ALLOCATABLE :: q(:,:), p(:,:)

    N = self % N
    ALLOCATE(q(N,3))
    ALLOCATE(p(N,3))

    init_E = total_energy(self%masses, self%qs0, self%ps0, 1.0_pf, SQRT(self%coarse_method%forces_method%softening2))
    init_J = total_angular_momentum(self%qs0, self%ps0)
    init_P = total_linear_momentum(self%ps0)

    output_file = -1
    IF (PRESENT(output_file_par)) output_file = output_file_par

    CALL msg('[parareal.run] running parareal', 2)

    ALLOCATE(U0(2*N,3))
    U0(1:N, :) = self % qs0
    U0(N+1:,:) = self % ps0

    ALLOCATE(U(self % N_itermax+1, self % number_of_windows, 2*N, 3))
    U = 0.0_pf

    ALLOCATE(U_tilde(self % number_of_windows, 2*N, 3))
    ALLOCATE(U_comp(self % number_of_windows, 2*N, 3))

    ! coarse application (zeroth iteration)
    U(1,1,:,:) = U0
    DO i = 1, self % number_of_windows - 1
        U(1,i+1,:,:) = self % apply_coarse(U(1,i,:,:))
    END DO

    timer_total = 0
    DO iterations = 1, self % N_itermax
        timer_0 = omp_get_wtime()
        U(iterations + 1, 1, :, :) = U0
        U_tilde(1,:,:) = U0

        WRITE(*,'(A)') "---------------------------------------------------"

        IF (self % fine_parallel) THEN
            !$OMP PARALLEL PRIVATE(i)
            !$OMP DO
            DO i = 1, self % number_of_windows - 1
                U_tilde(i+1,:,:) = self % apply_fine(U(iterations,i,:,:))
            END DO
            !$OMP END DO
            !$OMP END PARALLEL
        ELSE
            DO i = 1, self % number_of_windows - 1
                U_tilde(i+1,:,:) = self % apply_fine(U(iterations,i,:,:))
            END DO
        ENDIF

        DO i = 1, self % number_of_windows - 1
            U(iterations+1,i+1,:,:) = self % apply_coarse(U(iterations+1,i,:,:)) + U_tilde(i+1,:,:)
            U(iterations+1,i+1,:,:) = U(iterations+1,i+1,:,:) - self % apply_coarse(U(iterations,i,:,:))
            ! PRINT *, NORM2(U(iterations+1,i+1,:,:) - U(iterations,i+1,:,:))
        END DO

        timer_1 = omp_get_wtime()
        timer_total = timer_total + timer_1 - timer_0

        error = NORM2(U(iterations+1,self%number_of_windows,:,:)-U(iterations, self%number_of_windows,:,:))/(6*N)

        ! CALL msg('error:'//TRIM(error_str)//" / time:",, 1)
        WRITE(*,'(A6,I4,A9,ES12.4,A8,ES12.4)') "Iter.:", iterations, " / Error:", error, " / Time:", timer_1 - timer_0
        
        q = u(iterations+1,self%number_of_windows,1:N,:)
        p = u(iterations+1,self%number_of_windows,N+1:,:)
        error_E = ABS(total_energy(self%masses, q, p, 1.0_pf,SQRT(self%coarse_method%forces_method%softening2)) - init_E)
        error_J = NORM2(total_angular_momentum(q, p) - init_J)
        error_P = NORM2(total_linear_momentum(p) - init_P)
        WRITE(*,'(A15," = ",ES12.4)') '||E - E_0||', error_E
        WRITE(*,'(A15," = ",ES12.4)') '||J - J_0||', error_J
        WRITE(*,'(A15," = ",ES12.4)') '||P - P_0||', error_P

        IF (error <= self % tolerance) EXIT
    END DO

    WRITE(*,"(A)") "Simulation complete!"
    WRITE(*,"(A,E12.3)") "Total time:", omp_get_wtime() - timer_0
    WRITE(*,"(A,I4)") "Iterations:", iterations
    WRITE(*,"(A,E12.3)") "Final rel. error:", error
    WRITE(*,'(A15," = ",ES12.4)') '||E - E_0||', error_E
    WRITE(*,'(A15," = ",ES12.4)') '||J - J_0||', error_J
    WRITE(*,'(A15," = ",ES12.4)') '||P - P_0||', error_P

    IF (output_file > 0) THEN
        DO i = 1, self % number_of_windows
            self % qs = U(iterations+1, i, 1:N, :)
            self % ps = U(iterations+1, i, N+1:, :)
            WRITE (output_file, *) self % qs, self % ps
        END DO
    ENDIF

    DEALLOCATE(U, U0, U_tilde)
END SUBROUTINE

FUNCTION apply_coarse (self, U) RESULT(U_out)
    CLASS(parareal_simulation_type), INTENT(INOUT) :: self
    REAL(pf), DIMENSION(2*self % N, 3), INTENT(IN) :: U
    REAL(pf), DIMENSION(2*self % N, 3) :: U_out
    REAL(pf), DIMENSION(self % N, 3) :: qs, ps, fs
    INTEGER :: step
    LOGICAL :: eval_forces
    
    qs = U(1:self%N,:)
    ps = U(self%N+1:,:)

    eval_forces = .TRUE.
    DO step = 1, self % coarse_dt_mult
        CALL self % coarse_method % apply(qs, ps, fs, eval_forces)
    END DO

    U_out(1:self%N,:) = qs
    U_out(self%N+1:,:) = ps
END FUNCTION

FUNCTION apply_fine (self, U) RESULT(U_out)
    CLASS(parareal_simulation_type), INTENT(INOUT) :: self
    REAL(pf), DIMENSION(2*self % N, 3), INTENT(IN) :: U
    REAL(pf), DIMENSION(2*self % N, 3) :: U_out
    REAL(pf), DIMENSION(self % N, 3) :: qs, ps, fs
    INTEGER :: step
    LOGICAL :: eval_forces

    qs = U(1:self%N,:)
    ps = U(self%N+1:,:)
    eval_forces = .TRUE.

    DO step = 1, self % fine_dt_mult
        CALL self % fine_method % apply(qs, ps, fs, eval_forces)
    END DO
    
    U_out(1:self%N,:) = qs
    U_out(self%N+1:,:) = ps
END FUNCTION

END MODULE