! ************************************************************
!! Utilities
!
!> Objectives
!  This module contains some utilities, like observables and
!  the generation of random initial values.
!
!> Modified
!  2026.06.15
!
!> Author
!  oap
!
MODULE utils_mod
    USE types_mod
    IMPLICIT NONE
    PUBLIC
CONTAINS

SUBROUTINE generate_initial_values (N, m, qs, ps, size_qs, size_ps, nm, two_dim)
    INTEGER,  INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: size_qs, size_ps
    REAL(pf), INTENT(INOUT) :: m(N), qs(N,3), ps(N,3)
    LOGICAL, INTENT(IN) :: nm, two_dim

    CALL RANDOM_SEED()

    ! normalized masses
    IF (nm) THEN
        m = 1.0_pf / N
    ELSE
        CALL random_number(m)
        DO WHILE (MINVAL(m) == 0.0_pf)
            CALL random_number(m)
        END DO
    ENDIF

    CALL RANDOM_NUMBER(qs)
    qs = size_qs * (2.0_pf * qs - 1.0_pf)

    CALL RANDOM_NUMBER(ps)
    ps = size_ps * (2.0_pf * ps - 1.0_pf)

    IF (two_dim) THEN
        qs(:,3) = 0.0_pf
        ps(:,3) = 0.0_pf
    ENDIF
END SUBROUTINE

FUNCTION total_energy (m, qs, ps, G, eps)
    REAL(pf), INTENT(IN) :: m(:), qs(:,:), ps(:,:), G, eps
    REAL(pf) :: total_energy, V, rab
    INTEGER :: a, b

    total_energy = 0.0_pf
    DO a = 1, SIZE(m)
        total_energy = total_energy + 0.5_pf*DOT_PRODUCT(ps(a,:),ps(a,:))/m(a)
        DO b = 1, a - 1
            rab = SQRT(NORM2(qs(b,:)-qs(a,:))**2 + eps*eps)
            V = - G * m(a) * m(b) / rab
            total_energy = total_energy + V
        END DO
    END DO 
END FUNCTION

FUNCTION total_linear_momentum (ps) RESULT(P)
    REAL(pf), INTENT(IN) :: ps(:,:)
    REAL(pf) :: P(3)
    INTEGER :: a

    P = 0.0_pf
    DO a = 1, SIZE(ps,1)
        P = P + ps(a,:)
    END DO
END FUNCTION

FUNCTION cross_product (u, v)
  REAL(pf), DIMENSION(:), INTENT(IN) :: u, v
  REAL(pf), DIMENSION(3)             :: cross_product

  cross_product(1) =  u(2)*v(3)-v(2)*u(3)
  cross_product(2) = -u(1)*v(3)+v(1)*u(3)
  cross_product(3) =  u(1)*v(2)-v(1)*u(2)
END FUNCTION

FUNCTION total_angular_momentum (qs, ps) RESULT(J)
    REAL(pf), INTENT(IN) :: qs(:,:), ps(:,:)
    REAL(pf) :: J(3)
    INTEGER :: a

    J = 0.0_pf
    DO a = 1, SIZE(ps,1)
        J = J + cross_product(qs(a,:), ps(a,:))
    END DO
END FUNCTION

END MODULE