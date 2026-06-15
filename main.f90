! ************************************************************
!! NBODY-PARAREAL
!
!> Objectives
!  This file calls all other routines and do the management
!  of input parameters.
!
!> Modified
!  2026.06.15
!
!> Author
!  oap
!
PROGRAM main
    USE configs_mod
    USE types_mod
    USE utils_mod
    USE json_utils_mod
    USE integrators_mod, ONLY : list_available_integrators
    USE sequential_simulator_mod
    USE parareal_simulator_mod
    IMPLICIT NONE

    CHARACTER(256) :: action

    CALL header
    WRITE (*,"(A)") 'v'//version_string//'_'//precision//' ('//build_date//' '//build_time//')'
    WRITE (*,*)

    IF (command_argument_count() == 0) THEN
        STOP "No arguments were provided! For more information, use -h or --help"
    ENDIF

    ! Get the action argument
    CALL get_command_argument(1, action)
    SELECT CASE(TRIM(action))
        CASE ('-h', '--help')
            CALL help

        CASE ('-r', '--run')
            CALL run

        CASE DEFAULT
            STOP 'Action "'//TRIM(action)//'" not found. For more information, use -h or --help'
    END SELECT
CONTAINS

SUBROUTINE help ()
    WRITE (*,"(A)") "USE:"
    WRITE (*,"(A)") "    ./nbody_parareal [--help|-h]"
    WRITE (*,"(A)") "    ./nbody_parareal [--run|-r] <input_filename> <output_filename:optional>"
    WRITE (*,*)
    WRITE (*,"(A)") "AVAILABLE NUMERICAL INTEGRATORS:"
    WRITE (*,"(A)",ADVANCE='NO') "    "
    CALL list_available_integrators
    WRITE (*,*)
    WRITE (*,"(A)") "EXAMPLES:"
    WRITE (*,"(A)") "    ./nbody_parareal -r examples/sequential_random.json"
    WRITE (*,"(A)") "    ./nbody_parareal -r examples/sequential_iv.json"
    WRITE (*,"(A)") "    ./nbody_parareal -r examples/parareal_random.json"
    WRITE (*,"(A)") "    ./nbody_parareal -r examples/parareal_iv.json"
    WRITE (*,*)
END SUBROUTINE

SUBROUTINE run ()
    TYPE(json_value), POINTER :: json_data
    CHARACTER(LEN=20) :: file_type, integration_type
    INTEGER :: N
    INTEGER :: file_out_int
    REAL(pf), ALLOCATABLE :: m(:), qs(:,:), ps(:,:)

    IF (command_argument_count() < 1 + 1) THEN
        STOP "You need to provide an input file. For more information, use -h or --help"
    ENDIF

    !> manage the input and output files
    CALL manage_files_parameters(json_data, file_out_int)

    !> initial values type case
    CALL msg('[main.run] selecting file type case...', 2)
    file_type = json_get_string(json_data, "type")
    
    IF (TRIM(file_type) == "random") THEN
        CALL msg('[main.run] case random selected, generating random initial values', 2)
        CALL random_generate(json_data, m, qs, ps)
        CALL msg('[main.run] initial values generated!', 2)

    ELSE IF (TRIM(file_type) == "initial_values") THEN
        CALL msg('[main.run] case initial values selected, reading values', 2)
        N = json_get_int(json_data, 'N')
        ALLOCATE(m(N))
        ALLOCATE(qs(N,3))
        ALLOCATE(ps(N,3))
        m = json_get_float_vec(json_data, 'initial_values.masses')
        qs = json_get_float_matrix(json_data, 'initial_values.positions', N, 3)
        ps = json_get_float_matrix(json_data, 'initial_values.momenta', N, 3)
        CALL msg('[main.run] initial values readed', 2)

    ELSE
        STOP 'File type "'//TRIM(file_type)//'" not found!'
    ENDIF

    !> integration type case
    CALL msg('[main.run] selecting integration type case...', 2)
    integration_type = json_get_string(json_data, "integration.type")

    IF (TRIM(integration_type) == "sequential") THEN
        CALL msg('[main.run] case sequential selected, running...', 2)
        CALL run_sequential(json_data, file_out_int, m, qs, ps)

    ELSE IF (TRIM(integration_type) == "parareal") THEN
        CALL msg('[main.run] case parareal selected, running...', 2)
        CALL run_parareal(json_data, file_out_int, m, qs, ps)

    ELSE
        STOP 'Integration type "'//TRIM(integration_type)//'" not found!'
    ENDIF
END SUBROUTINE

SUBROUTINE manage_files_parameters (json_data, file_out_int)
    TYPE(json_value), POINTER :: json_data
    INTEGER, INTENT(INOUT) :: file_out_int
    CHARACTER(LEN=256) :: file_in, file_out
    LOGICAL :: file_in_exists

    !> input file
    CALL get_command_argument(2, file_in)
    CALL msg('[main.run] reading the input file "' // TRIM(file_in) // '"', 2)
    INQUIRE(FILE=TRIM(file_in), EXIST=file_in_exists)
    IF (.NOT. file_in_exists) THEN
        STOP 'The input file "'//TRIM(file_in)//'" was not found!'
    ENDIF
    CALL read_json(TRIM(file_in), json_data)
    CALL msg('[main.run] json readed', 2)

    !> output file
    IF (command_argument_count() < 3) THEN
        ! if no out file provided
        file_out_int = 0
    ELSE
        file_out_int = 23
        CALL get_command_argument(3, file_out)
        CALL msg('[main.run] creating the output file "' // TRIM(file_out) // '"', 2)
        OPEN(file_out_int, file = TRIM(file_out), status="replace")
    ENDIF
END SUBROUTINE

! to generate random initial values
SUBROUTINE random_generate (json_data, m, qs, ps)
    TYPE(json_value), POINTER, INTENT(IN) :: json_data
    REAL(pf), ALLOCATABLE, INTENT(INOUT)  :: m(:), qs(:,:), ps(:,:)
    INTEGER :: N, masses_value
    REAL(pf) :: qs_range, ps_range
    LOGICAL :: nm, two_dimensional

    N = json_get_int(json_data, "N")

    qs_range = json_get_float(json_data, "random_parameters.qs_range")
    ps_range = json_get_float(json_data, "random_parameters.ps_range")
    nm = json_get_logical(json_data, "random_parameters.normalized_masses")
    two_dimensional = json_get_logical(json_data, "random_parameters.2d")

    ALLOCATE(m(N))
    ALLOCATE(qs(N,3))
    ALLOCATE(ps(N,3))

    CALL generate_initial_values(N, m, qs, ps, qs_range, ps_range, nm, two_dimensional)

    IF (masses_value > 0) m = masses_value
END SUBROUTINE

! to simulate sequentially (in time)
SUBROUTINE run_sequential (json_data, output_file, m, qs, ps)
    TYPE(json_value), POINTER :: json_data
    INTEGER,  INTENT(IN)    :: output_file
    REAL(pf), INTENT(INOUT) :: m(:), qs(:,:), ps(:,:)

    INTEGER  :: N, checkpoints
    REAL(pf) :: dt, tf
    CLASS(sequential_simulation_type), POINTER :: simulation
    REAL(pf)      :: G, fsoft, fangle
    LOGICAL       :: fpar
    CHARACTER(30) :: fmethod
    CHARACTER(30) :: method

    CALL msg('[main.run_sequential] getting information about the simulation', 2)
    method = json_get_string(json_data, "integration.method")
    N = json_get_int(json_data, 'N')
    dt = json_get_float(json_data, 'integration.timestep')
    tf = json_get_float(json_data, 'integration.tf')
    G = json_get_float(json_data, 'G')
    checkpoints = json_get_int(json_data, 'integration.checkpoints')
    fsoft = json_get_float(json_data, 'integration.softening')
    fangle = json_get_float(json_data, 'integration.forces_angle')
    fpar = json_get_logical(json_data, 'integration.forces_parallel')
    fmethod = json_get_string(json_data, 'integration.forces_method')

    CALL msg('[main.run_sequential] initializing the sequential simulator', 2)
    ALLOCATE(simulation)
    CALL simulation % init(N, m, qs, ps, & ! initial values
        dt, tf, method, checkpoints, & ! method parameters
        G, fsoft, fangle, fpar, fmethod & ! forces parameters
    )
    CALL simulation % run(output_file)
END SUBROUTINE

! to simulate parallel (in time) via parareal
SUBROUTINE run_parareal (json_data, output_file, m, qs, ps)
    TYPE(json_value), POINTER :: json_data
    INTEGER,  INTENT(IN)    :: output_file
    REAL(pf), INTENT(INOUT) :: m(:), qs(:,:), ps(:,:)
    
    ! simulation
    INTEGER  :: N
    REAL(pf) :: tf, G, softening
    
    ! parareal config
    CLASS(parareal_simulation_type), POINTER :: parareal
    INTEGER  :: number_of_windows, N_itermax
    REAL(pf) :: tolerance
    LOGICAL  :: fine_parallel

    ! coarse method
    CHARACTER(30) :: co_method, co_fmethod
    INTEGER       :: co_dt_mult
    REAL(pf)      :: co_fa, co_fs
    LOGICAL       :: co_fp

    ! fine method
    CHARACTER(30) :: fi_method, fi_fmethod
    INTEGER       :: fi_dt_mult
    REAL(pf)      :: fi_fa, fi_fs
    LOGICAL       :: fi_fp
    
    CALL msg('[main.run_parareal] getting information about the simulation', 2)
    N = json_get_int(json_data, 'N')
    G = json_get_float(json_data, 'G')
    tf = json_get_float(json_data, 'integration.tf')
    softening = json_get_float(json_data, 'integration.softening')
    
    ! parareal
    number_of_windows = json_get_int(json_data, 'integration.number_of_windows')
    N_itermax = json_get_int(json_data, 'integration.N_itermax')
    tolerance = json_get_float(json_data, 'integration.tolerance')
    fine_parallel = json_get_logical(json_Data, 'integration.parareal_parallel')

    ! coarse
    co_method = json_get_string(json_data, 'integration.coarse.method')
    co_dt_mult = json_get_int(json_data, 'integration.coarse.dt_mult')
    co_fmethod = json_get_string(json_data, 'integration.coarse.forces_method')
    co_fa = json_get_float(json_data, 'integration.coarse.forces_angle')
    co_fp = json_get_logical(json_data, 'integration.coarse.forces_parallel')
    co_fs = softening

    ! fine
    fi_method = json_get_string(json_data, 'integration.fine.method')
    fi_dt_mult = json_get_int(json_data, 'integration.fine.dt_mult')
    fi_fmethod = json_get_string(json_data, 'integration.fine.forces_method')
    fi_fa = json_get_float(json_data, 'integration.fine.forces_angle')
    fi_fp = json_get_logical(json_data, 'integration.fine.forces_parallel')
    fi_fs = softening

    CALL msg('[main.run_parareal] initializing the parareal simulator', 2)
    ALLOCATE(parareal)
    CALL parareal % init(N, m, qs, ps, & ! initial values
        tf, G, & ! simulation configs
        number_of_windows, N_itermax, tolerance, fine_parallel, & ! parareal configs
        co_method, co_dt_mult, co_fmethod, co_fa, co_fs, co_fp, & ! coarse
        fi_method, fi_dt_mult, fi_fmethod, fi_fa, fi_fs, fi_fp &  ! fine
    )
    CALL parareal % run(output_file)
END SUBROUTINE

END PROGRAM