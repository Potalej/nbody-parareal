! ************************************************************
!! Forces
!
!> Objectives
!  This module contains direct-forces routines and using the
!  Barnes-Hut method. It also defines the forces function of
!  the integrator_type.
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
MODULE forces_mod
    USE configs_mod
    USE octree_mod
    USE types_mod
    IMPLICIT NONE
    PUBLIC forces_type

    ABSTRACT INTERFACE
        FUNCTION forces_function (N, ms, qs, G, eps2, angle2)
            IMPORT :: pf
            IMPLICIT NONE
            INTEGER,                  INTENT(IN) :: N
            REAL(pf), DIMENSION(N),   INTENT(IN) :: ms
            REAL(pf), DIMENSION(N,3), INTENT(IN) :: qs
            REAL(pf),                 INTENT(IN) :: G, eps2, angle2
            REAL(pf), DIMENSION(N,3) :: forces_function
        END FUNCTION
    END INTERFACE

    TYPE :: forces_type
        REAL(pf) :: softening2
        REAL(pf) :: G
        REAL(pf) :: angle2
        PROCEDURE(forces_function), POINTER, NOPASS :: forces_function => NULL()
        CONTAINS
            PROCEDURE :: init => init_forces
            PROCEDURE :: eval
    END TYPE

CONTAINS

SUBROUTINE init_forces (self, use_parallel, use_bh, G, softening, angle)
    CLASS(forces_type), INTENT(INOUT) :: self
    LOGICAL,  INTENT(IN) :: use_parallel, use_bh
    REAL(pf), INTENT(IN) :: G, softening, angle

    CALL msg("[forces] initializing", 2)

    IF (use_parallel) THEN
        IF (use_bh) THEN
            CALL msg("[forces] allocating bh par", 2)
            self % forces_function => compute_forces_bh_par
        ELSE
            CALL msg("[forces] allocating direct par", 2)
            self % forces_function => compute_forces_direct_par
        ENDIF
    ELSE
        IF (use_bh) THEN
            CALL msg("[forces] allocating bh", 2)
            self % forces_function => compute_forces_bh
        ELSE
            CALL msg("[forces] allocating direct", 2)
            self % forces_function => compute_forces_direct
        ENDIF
    ENDIF

    ! parameters of the Newtonian potential
    self % G = G
    self % softening2 = softening * softening

    ! parameter of the Barnes-Hut method
    self % angle2 = angle * angle
END SUBROUTINE

FUNCTION eval (self, N, ms, qs) RESULT(forces)
    CLASS(forces_type), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: ms(N), qs(N,3)
    REAL(pf) :: forces(N,3)

    forces = self % forces_function(N, ms, qs, self % G, self % softening2, self % angle2)
END FUNCTION

FUNCTION compute_forces_direct (N, ms, qs, G, eps2, angle2) RESULT(forces)
    INTEGER,  INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: ms(N), qs(N,3), G, eps2, angle2
    REAL(pf) :: forces(N,3), dist2, dist, f
    REAL(pf) :: dx, dy, dz
    INTEGER  :: a, b

    forces = 0.0_pf

    DO a = 2, N
        DO b = 1, a - 1
            dx = qs(b,1) - qs(a,1)
            dy = qs(b,2) - qs(a,2)
            dz = qs(b,3) - qs(a,3)
            
            dist2 = dx*dx + dy*dy + dz*dz + eps2
            dist = SQRT(dist2)
            f = G * ms(a) * ms(b) / (dist * dist2)

            ! forces over 'a'
            forces(a,1) = forces(a,1) + f * dx
            forces(a,2) = forces(a,2) + f * dy
            forces(a,3) = forces(a,3) + f * dz

            ! forces over 'b'
            forces(b,1) = forces(b,1) - f * dx
            forces(b,2) = forces(b,2) - f * dy
            forces(b,3) = forces(b,3) - f * dz
        END DO
    END DO
END FUNCTION

FUNCTION compute_forces_direct_par (N, ms, qs, G, eps2, angle2) RESULT(forces)
    INTEGER,  INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: ms(N), qs(N,3), G, eps2, angle2
    REAL(pf) :: forces(N,3), dist2, dist, f
    REAL(pf) :: dx, dy, dz
    INTEGER  :: a, b

    forces = 0.0_pf

    !$OMP PARALLEL SHARED(forces) PRIVATE(a, b, dx, dy, dz, dist2, dist, f)
    !$OMP DO
    DO a = 1, N
        DO b = 1, N
            IF (a == b) CYCLE

            dx = qs(b,1) - qs(a,1)
            dy = qs(b,2) - qs(a,2)
            dz = qs(b,3) - qs(a,3)
            
            dist2 = dx*dx + dy*dy + dz*dz + eps2
            dist = SQRT(dist2)
            f = G * ms(a) * ms(b) / (dist * dist2)

            ! forces over 'a'
            forces(a,1) = forces(a,1) + f * dx
            forces(a,2) = forces(a,2) + f * dy
            forces(a,3) = forces(a,3) + f * dz
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END FUNCTION

! compute forces using the barnes-hut method (sequential)
FUNCTION compute_forces_bh (N, ms, qs, G, eps2, angle2) RESULT (forces)
    INTEGER,  INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: ms(N), qs(N,3)
    REAL(pf), INTENT(IN) :: G, eps2, angle2

    TYPE(OctreeType) :: tree
    INTEGER :: p
    REAL(pf) :: forces(N, 3)

    CALL tree % init(ms, qs(:,1), qs(:,2), qs(:,3))
    CALL tree % evaluate_quad()

    DO p = 1, N
        forces(p,:) = tree % forces(p, angle2, G, eps2)
    END DO
END FUNCTION

! compute forces using the barnes-hut method (parallel)
FUNCTION compute_forces_bh_par (N, ms, qs, G, eps2, angle2) RESULT (forces)
    INTEGER,  INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: ms(N), qs(N,3)
    REAL(pf), INTENT(IN) :: G, eps2, angle2

    TYPE(OctreeType) :: tree
    INTEGER :: p
    REAL(pf) :: forces(N, 3)

    CALL tree % init(ms, qs(:,1), qs(:,2), qs(:,3))
    CALL tree % evaluate_quad()

    !$OMP PARALLEL SHARED(forces) PRIVATE(p) 
    !$OMP DO
    DO p = 1, N
        forces(p,:) = tree % forces(p, angle2, G, eps2)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END FUNCTION

END MODULE