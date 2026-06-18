cmake_minimum_required(VERSION 3.15)
project("nbody_parareal" VERSION 0.0.4
        DESCRIPTION "Nbody solver using the Parareal method!"
        LANGUAGES Fortran)
enable_language(Fortran)

# Compilation time
string(TIMESTAMP BUILD_DATE "%Y-%m-%d")
string(TIMESTAMP BUILD_TIME "%H:%M:%S")

# ============================================================================
# Compilation
# ----------------------------------------------------------------------------
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -fopenmp")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fcheck=all -Wall -fbacktrace -O3 -march=native -fopenmp") # Debug
message("\n******************************************")
message("#          nbody-parareal  v${CMAKE_PROJECT_VERSION}        #")
message("******************************************")
# ============================================================================

# ============================================================================
# GProf
# ----------------------------------------------------------------------------
option(GPROF "gprof profile" OFF)

if (GPROF)
  message(STATUS "[debug] GNU Profiler (gprof) enabled")
  add_compile_options(-pg)
  add_link_options(-pg)
endif()
# ============================================================================

# ============================================================================
# Precision
# ----------------------------------------------------------------------------
set(PRECISION "64" CACHE STRING "Choose the precision: 32, 64 or 128")
set_property(CACHE PRECISION PROPERTY STRINGS 32 64 128)

# Define pre-processor flags for each precision
if(PRECISION STREQUAL "32")
  add_compile_definitions(REAL32)
elseif(PRECISION STREQUAL "64")
  add_compile_definitions(REAL64)
elseif(PRECISION STREQUAL "128")
  add_compile_definitions(REAL128)
else()
  message(FATAL_ERROR "Invalid precision: ${PRECISION}. Enter 32, 64 or 128")
endif()
# ============================================================================

# ============================================================================
# Configs
# ----------------------------------------------------------------------------
set(PRINT_LEVEL "1" CACHE STRING "Choose the print level: 0, 1 or 2")
set_property(CACHE PRINT_LEVEL PROPERTY STRINGS 0 1 2)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configs.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/configs.f90
  @ONLY
)
add_library(configs
  ${CMAKE_CURRENT_BINARY_DIR}/configs.f90
)
# ============================================================================

# ============================================================================
# Executables
# ----------------------------------------------------------------------------
add_executable(nbody_parareal
  "main.f90"
)
# target_link_libraries(version)
target_include_directories(nbody_parareal PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# ============================================================================

# ============================================================================
# Inclui o JSONFortran
# ----------------------------------------------------------------------------
message("\n> External lib: JSON-Fortran")
option(FETCHJSONFORTRAN "Compile json-fortran locally" OFF)

if (FETCHJSONFORTRAN OR PRECISION EQUAL "32" OR PRECISION EQUAL "128" OR NOT DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Baixando do repositório para compilar localmente. Isso pode demorar um pouco...")
  
  include(FetchContent)

  FetchContent_Declare(
    jsonfortran
    GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
    GIT_TAG 9.2.0
    GIT_SHALLOW TRUE
  )

  set(JSON_REAL_KIND "REAL${PRECISION}" CACHE INTERNAL "Tipo da variavel real do JSON-Fortran")
  
  FetchContent_MakeAvailable(jsonfortran)

  target_include_directories(jsonfortran
    INTERFACE
      ${jsonfortran_BINARY_DIR}
  )

elseif (DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Ambiente conda ativo, verificando instalacao...")
  set(JSONFORTRAN_LIB "$ENV{CONDA_PREFIX}/lib/libjsonfortran.so")
  set(JSONFORTRAN_INC "$ENV{CONDA_PREFIX}/include")
  
  if (EXISTS "${JSONFORTRAN_LIB}")
    add_library(jsonfortran::jsonfortran SHARED IMPORTED)
    set_target_properties(jsonfortran::jsonfortran PROPERTIES
      IMPORTED_LOCATION "${JSONFORTRAN_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${JSONFORTRAN_INC}"
    )
    add_library(jsonfortran ALIAS jsonfortran::jsonfortran)
    message(STATUS "JSON-Fortran habilitado via conda!")
  else()
    message(FATAL_ERROR "JSON-Fortran nao encontrado no conda! Instale no conda ou compile localmente atraves da flag -DFETCHJSONFORTRAN=ON")
  endif()
endif()
# ============================================================================

# ============================================================================
# =-=-=-=-=-=-=-=-=-=-=-=-=-=- Modules =-=-=-=-=-=-=-=-=-=-=-=-=-=-
# ============================================================================

# ============================================================================
#> Types
# ----------------------------------------------------------------------------
add_library(types
"./src/types.F90"
)
# ============================================================================

# ============================================================================
#> json_utils_mod
# ----------------------------------------------------------------------------
add_library(json_utils_mod 
  "./src/libs/json_utils_mod.f90"
)
target_link_libraries(json_utils_mod PRIVATE types)
target_link_libraries(json_utils_mod PRIVATE jsonfortran)
# ============================================================================

# ============================================================================
#> Octree
# ----------------------------------------------------------------------------
add_library(tree
  "./src/barneshut/octree.f90"
)
target_link_libraries(tree PRIVATE types)
# ============================================================================

# ============================================================================
#> Forces
# ----------------------------------------------------------------------------
add_library(forces
  "./src/forces/forces.f90"
)
target_link_libraries(forces PRIVATE types)
target_link_libraries(forces PRIVATE tree)
# ============================================================================

# ============================================================================
#> Numerical integration
# ----------------------------------------------------------------------------
add_subdirectory(src/integrators) # Metodos de integracao
add_library(integration
  "./src/integrators/integrator.f90"
  ${CMAKE_CURRENT_BINARY_DIR}/src/integrators/integrators.f90
  ${INTEGRATORS_SOURCE}
  
)
target_link_libraries(integration PRIVATE types)
target_link_libraries(integration PRIVATE forces)

# ============================================================================

# ============================================================================
#> utils
# ----------------------------------------------------------------------------
add_library(utils
  "./src/utils.f90"
)
target_link_libraries(utils PRIVATE types)

# ============================================================================

# ============================================================================
#> Simulation
# ----------------------------------------------------------------------------
add_library(simulation
  "./src/simulation/sequential_simulator.f90"
  "./src/simulation/parareal_simulator.f90"
)
target_link_libraries(simulation PRIVATE types)
target_link_libraries(simulation PRIVATE integration)
target_link_libraries(simulation PRIVATE utils)

# ============================================================================

# ============================================================================
# Link the libs to the executable
# ----------------------------------------------------------------------------
target_link_libraries(nbody_parareal
  configs
  types
  json_utils_mod
  simulation
  utils
)
# ============================================================================

# ============================================================================
# After the compilation, copy the executable to the root
# ----------------------------------------------------------------------------
add_custom_command(TARGET nbody_parareal POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
     "$<TARGET_FILE:nbody_parareal>"
     "../$<CONFIGURATION>/$<TARGET_FILE_NAME:nbody_parareal>" 
  COMMENT "Copying the executable to the root.")
# ----------------------------------------------------------------------------
