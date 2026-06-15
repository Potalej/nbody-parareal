# ============================================================================
# CMake File for the integrators
# ----------------------------------------------------------------------------
# This file is called by the root CMakeLists.txt and contains the set of
# available numerical integrators.
# When an integrator is added, one just need to change here and recompile the
# program to the new integrator became available.
# ============================================================================

set(INTEGRATORS_SOURCE
  # Basic symplectic methods
  "./src/integrators/methods/symplectic/symp_euler.f90" # (Symplectic Euler)
  "./src/integrators/methods/symplectic/verlet.f90" # (Verlet)
CACHE STRING "List of integrators")

# Formated lists to Fortran
set(INTEGRATORS_USE "")
set(INTEGRATORS_CLASSES "")
set(INTEGRATORS_LIST "")
set(INTEGRATORS_SELECT "")

foreach(src ${INTEGRATORS_SOURCE})
  # get the file
  get_filename_component(fname ${src} NAME_WE)
  string(TOUPPER ${fname} FNAME_UP)

  # 'USE' in the script
  string(APPEND INTEGRATORS_USE "  USE ${fname}_integrator_mod\n")

  # TYPE in the script
  string(APPEND INTEGRATORS_CLASSES "  TYPE(${fname}_integrator_type), TARGET, SAVE :: INT_${fname}\n")

  # CASE of the method in the script
  string(APPEND INTEGRATORS_SELECT "    CASE ('${fname}'); integrator => INT_${fname}\n")
endforeach()

SET(COUNTER 0)
list(LENGTH INTEGRATORS_SOURCE NUMBER_OF_INTEGRATORS)
foreach(src ${INTEGRATORS_SOURCE})
  # get the file
  get_filename_component(fname ${src} NAME_WE)
  string(TOUPPER ${fname} FNAME_UP)

  # add to the counter
  math(EXPR COUNTER "${COUNTER} + 1")

  # add to the string
  if (COUNTER EQUAL 1)
    string(APPEND INTEGRATORS_LIST "  WRITE(*,'(A)',ADVANCE='NO') '${fname}, '\n")
  elseif(COUNTER EQUAL NUMBER_OF_INTEGRATORS)
    string(APPEND INTEGRATORS_LIST "  WRITE(*,'(A)') '${fname}.'\n")
  else()
    string(APPEND INTEGRATORS_LIST "  WRITE(*,'(A)',ADVANCE='NO') '${fname}, '\n")
  endif()
endforeach()

# generate the file from the template
configure_file(
  ${CMAKE_SOURCE_DIR}/src/integrators/integrators.f90.in
  ${CMAKE_BINARY_DIR}/src/integrators/integrators.f90
  @ONLY
)