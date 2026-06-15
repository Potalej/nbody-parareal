! ************************************************************
!! Types
!
!> Objectives
!  This module contains the precision types used in the program.
!
!> Modified
!  2026.06.15
!
!> Author
!  oap
!
MODULE types_mod
  IMPLICIT NONE
  INTEGER, PARAMETER, PUBLIC :: pf32  = SELECTED_REAL_KIND(6, 37)    ! simple precision (~32 bits)
  INTEGER, PARAMETER, PUBLIC :: pf64  = SELECTED_REAL_KIND(15, 307)  ! double precision (~64 bits)
  INTEGER, PARAMETER, PUBLIC :: pf128 = SELECTED_REAL_KIND(33, 4931) ! quadruple precision (~128 bits)

#ifdef REAL32
  INTEGER, PARAMETER, PUBLIC :: pf = pf32   !! Standard real type [4 bytes]
#elif REAL64
  INTEGER, PARAMETER, PUBLIC :: pf = pf64   !! Standard real type [8 bytes]
#elif REAL128
  INTEGER, PARAMETER, PUBLIC :: pf = pf128  !! Standard real type [16 bytes]
#else
  INTEGER, PARAMETER, PUBLIC :: pf = pf64   !! Standard real type if not specified [8 bytes]
#endif
END MODULE