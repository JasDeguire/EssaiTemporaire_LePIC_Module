cmake_minimum_required(VERSION 3.23)
project(LePIC_3D LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

get_filename_component(FORTRAN_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

set(USING_MPI_WRAPPER FALSE)
if(FORTRAN_COMPILER_NAME MATCHES "^(mpiifx|mpiifort|mpifort|mpif90|mpifc)$")
  set(USING_MPI_WRAPPER TRUE)
endif()

message(STATUS "Fortran compiler          : ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler ID       : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Fortran compiler basename : ${FORTRAN_COMPILER_NAME}")
message(STATUS "Using MPI wrapper         : ${USING_MPI_WRAPPER}")

if(USING_MPI_WRAPPER)
  message(STATUS "MPI provided by wrapper: ${FORTRAN_COMPILER_NAME}")
else()
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

# ------------------------------------------------------------
# Compiler flags
# ------------------------------------------------------------
set(OMP_COMPILE_FLAG "")
set(OMP_LINK_FLAG "")

if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")

  set(OMP_COMPILE_FLAG "-qopenmp")
  set(OMP_LINK_FLAG    "-qopenmp")

  set(INTEL_RELEASE_FLAGS
    -O3
    -ipo
    -no-prec-div
    -m64
    -xHost
    -qopenmp
    -auto
    -fp-model=fast=2
    -unroll
  )

  set(INTEL_DEBUG_FLAGS
    -O0
    -g
    -traceback
    -check
    all
    -fpe0
    -warn
    all
    -qopenmp
    -auto
    -fno-inline
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(${INTEL_DEBUG_FLAGS})
  else()
    add_compile_options(${INTEL_RELEASE_FLAGS})
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  set(OMP_COMPILE_FLAG "-fopenmp")
  set(OMP_LINK_FLAG    "-fopenmp")

  set(GNU_RELEASE_FLAGS
    -O3
    -march=native
    -funroll-loops
    -fopenmp
  )

  set(GNU_DEBUG_FLAGS
    -O0
    -g
    -fbacktrace
    -fcheck=all
    -ffpe-trap=invalid,zero,overflow
    -Wall
    -Wextra
    -fopenmp
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(${GNU_DEBUG_FLAGS})
  else()
    add_compile_options(${GNU_RELEASE_FLAGS})
  endif()

else()
  message(WARNING "Unknown compiler: ${CMAKE_Fortran_COMPILER_ID}")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g)
  else()
    add_compile_options(-O3)
  endif()
endif()

add_library(build_config INTERFACE)

target_include_directories(build_config INTERFACE
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_options(build_config INTERFACE
      -O0
      -g
      -traceback
      -qopenmp
    )
  else()
    target_link_options(build_config INTERFACE
      -O3
      -ipo
      -no-prec-div
      -m64
      -xHost
      -qopenmp
      -auto
      -fp-model=fast=2
      -unroll
    )
  endif()
elseif(OMP_LINK_FLAG)
  target_link_options(build_config INTERFACE ${OMP_LINK_FLAG})
endif()

if(NOT USING_MPI_WRAPPER)
  target_link_libraries(build_config INTERFACE MPI::MPI_Fortran)
  target_include_directories(build_config INTERFACE ${MPI_Fortran_INCLUDE_DIRS})
endif()

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
include_directories(${MODULES_DIR}/legacy)

add_library(core_deps OBJECT
  ${MODULES_DIR}/core/mod_constants.f90
  ${MODULES_DIR}/core/mod_utils.f90
  ${MODULES_DIR}/core/mod_RNG.f90
)
target_link_libraries(core_deps PUBLIC build_config)
target_include_directories(core_deps PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${MODULES_DIR}/legacy
)

add_library(config_base OBJECT
  ${MODULES_DIR}/config/mod_config.f90
)
target_link_libraries(config_base PUBLIC build_config core_deps)
target_include_directories(config_base PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

add_library(geometry_base OBJECT
  ${MODULES_DIR}/geometry/mod_domain.f90
)
target_link_libraries(geometry_base PUBLIC
  build_config
  core_deps
  config_base
)
target_include_directories(geometry_base PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(geometry_base config_base)

add_library(fields_base OBJECT
  ${MODULES_DIR}/fields/mod_fields.f90
  ${MODULES_DIR}/fields/mod_density.f90
)
target_link_libraries(fields_base PUBLIC
  build_config
  core_deps
  config_base
  geometry_base
)
target_include_directories(fields_base PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(fields_base geometry_base)

add_library(pic_core_obj OBJECT
  ${MODULES_DIR}/config/mod_readConditions.f90

  ${MODULES_DIR}/core/mod_intro.f90
  ${MODULES_DIR}/core/mod_part_info.f90
  ${MODULES_DIR}/core/mod_debug_checks.f90

  ${MODULES_DIR}/fields/mod_PoissonDecomposition.f90
  ${MODULES_DIR}/fields/mod_PoissonSolver_legacy.f90
  ${MODULES_DIR}/fields/mod_charge_weights.f90
  ${MODULES_DIR}/fields/mod_chargeDeposition.f90
  ${MODULES_DIR}/fields/mod_magneticField.f90
  ${MODULES_DIR}/fields/mod_electricField.f90

  ${MODULES_DIR}/geometry/mod_generateBoundary.f90
  ${MODULES_DIR}/geometry/mod_boundary.f90

  ${MODULES_DIR}/particles/mod_chemistryState.f90
  ${MODULES_DIR}/particles/mod_particles.f90

  ${MODULES_DIR}/collisions/mod_collisionProducts.f90
  ${MODULES_DIR}/collisions/mod_collisionUtilities.f90
  ${MODULES_DIR}/collisions/mod_collisionDiagnostics.f90
  ${MODULES_DIR}/collisions/mod_MCCcollisions.f90
  ${MODULES_DIR}/collisions/mod_BMCCcollisions.f90
  ${MODULES_DIR}/collisions/mod_DSMCcollisions.f90
  ${MODULES_DIR}/collisions/mod_collisionsGwenael.f90
  ${MODULES_DIR}/collisions/mod_collisions.f90
  ${MODULES_DIR}/collisions/mod_CoulombCollisions.f90

  ${MODULES_DIR}/particles/mod_particle_sorting.f90
  ${MODULES_DIR}/particles/mod_particleMover.f90
  ${MODULES_DIR}/particles/mod_particleBC.f90
  ${MODULES_DIR}/particles/mod_particle_loader.f90
  ${MODULES_DIR}/particles/mod_reactions.f90
  ${MODULES_DIR}/particles/mod_reactionsDB.f90
  ${MODULES_DIR}/particles/mod_heating.f90
  ${MODULES_DIR}/particles/mod_injection.f90
  ${MODULES_DIR}/particles/mod_planeMoments.f90

  ${MODULES_DIR}/simulation/mod_simParams.f90
  ${MODULES_DIR}/io/mod_output_2d.f90
  ${MODULES_DIR}/io/mod_restart.f90

  ${MODULES_DIR}/core/mod_state.f90
  ${MODULES_DIR}/simulation/mod_simulation.f90
)
target_link_libraries(pic_core_obj PUBLIC
  build_config
  core_deps
  config_base
  geometry_base
  fields_base
)
target_include_directories(pic_core_obj PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(pic_core_obj
  core_deps
  config_base
  geometry_base
  fields_base
)

add_library(legacy_mg OBJECT
  ${MODULES_DIR}/legacy/pdesolver.f90
  ${MODULES_DIR}/legacy/mg.f90
  ${MODULES_DIR}/legacy/sors.f90
)
target_link_libraries(legacy_mg PUBLIC build_config)
target_include_directories(legacy_mg PUBLIC
  ${MODULES_DIR}/legacy
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(legacy_mg pic_core_obj)

add_executable(run_min
  ${MODULES_DIR}/main/main.f90
  $<TARGET_OBJECTS:core_deps>
  $<TARGET_OBJECTS:config_base>
  $<TARGET_OBJECTS:geometry_base>
  $<TARGET_OBJECTS:fields_base>
  $<TARGET_OBJECTS:pic_core_obj>
  $<TARGET_OBJECTS:legacy_mg>
)

target_link_libraries(run_min PRIVATE
  build_config
)

target_include_directories(run_min PRIVATE
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${MODULES_DIR}/legacy
)

set_target_properties(run_min PROPERTIES
  OUTPUT_NAME "run_min"
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

option(ENABLE_SAFE_IPO "Try CMake-managed IPO/LTO" OFF)

if(ENABLE_SAFE_IPO AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error LANGUAGES Fortran)

  if(ipo_supported)
    message(STATUS "IPO supported: enabling IPO on run_min only")
    set_property(TARGET run_min PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO not supported: ${ipo_error}")
  endif()
endif()

message(STATUS "Executable will be: ${CMAKE_BINARY_DIR}/run_min")
