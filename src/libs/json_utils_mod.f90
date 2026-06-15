! ************************************************************
!! JSON Utils
!
!> Objectives
!  This module contains some utilities for the JSON-Fortran
!  API. Its essencially an API.
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
MODULE json_utils_mod
    USE json_module, only: json_core, json_value, json_ck
    USE types_mod
    IMPLICIT NONE

    TYPE(json_core) :: json
    PUBLIC
CONTAINS

SUBROUTINE read_json (file, data)
    CHARACTER(LEN=*), INTENT(IN) :: file
    TYPE(json_value), POINTER, INTENT(INOUT) :: data
    CALL json % parse(file=file, p=data)
END SUBROUTINE

FUNCTION json_get_string (dados, chave) RESULT(texto)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    CHARACTER(LEN=:, KIND=json_ck), ALLOCATABLE :: texto_ck
    CHARACTER(LEN=:), ALLOCATABLE :: texto
    INTEGER :: i

    CALL json % get(dados, chave, texto_ck)

    ALLOCATE(CHARACTER(LEN=LEN(texto_ck)) :: texto)
    DO i = 1, LEN(texto)
        texto(i:i) = ACHAR(IACHAR(texto_ck(i:i)))
    END DO
END FUNCTION json_get_string

FUNCTION json_get_int (dados, chave) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    INTEGER :: valor
    CALL json % get(dados, chave, valor)
END FUNCTION

FUNCTION json_get_logical (dados, chave) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    LOGICAL :: valor
    CALL json % get(dados, chave, valor)
END FUNCTION

FUNCTION json_get_float (dados, chave) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    REAL(pf)   :: valor
    CALL json % get(dados, chave, valor)
END FUNCTION json_get_float

FUNCTION json_get_float_vec (dados, chave) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    REAL(pf), ALLOCATABLE   :: valor(:)
    CALL json % get(dados, chave, valor)
END FUNCTION json_get_float_vec

FUNCTION json_get_float_matrix (dados, chave, lines, columns) RESULT(value)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    INTEGER, INTENT(IN) :: lines, columns
    INTEGER       :: line
    CHARACTER(32) :: line_str
    REAL(pf), ALLOCATABLE :: value_line(:), value(:,:)

    ALLOCATE(value(lines, columns))
    ALLOCATE(value_line(columns))
    DO line = 1, lines
        WRITE(line_str, *) line
        value_line = json_get_float_vec(dados, chave//"["//line_str//"]")
        value(line,:) = value_line
    END DO
    DEALLOCATE(value_line)
END FUNCTION json_get_float_matrix

SUBROUTINE json_clone (entrada, saida)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: entrada
    TYPE(json_value), POINTER, INTENT(OUT) :: saida
    call json % clone(entrada, saida)
END SUBROUTINE json_clone

END MODULE json_utils_mod