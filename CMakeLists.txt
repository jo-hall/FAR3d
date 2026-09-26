cmake_minimum_required(VERSION 3.10)

project(far3d LANGUAGES Fortran C)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(DEBUG_FLAGS "-g -fcheck=all -fbacktrace")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(DEBUG_FLAGS "-g -check all -traceback")
else()
    message(WARNING "Unknown Fortran compiler, using default debug flags")
    set(DEBUG_FLAGS "-g")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${DEBUG_FLAGS}")

# Enable the Fortran preprocessor so files can gate GPU-only code (e.g. OpenACC
# device init in far3d.f90) behind #ifdef _OPENACC, which compilers only define
# when an OpenACC-enabling flag (-fopenacc / -acc=...) is actually passed.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# -------------------------
# Find MPI
# -------------------------
find_package(MPI REQUIRED)

# -------------------------
# Find OpenMP (optional, but often used with MPI)
# -------------------------
find_package(OpenMP REQUIRED)

# -------------------------
# Source files
# -------------------------
# Automatically include all Fortran source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.f90 src/*.F90 src/*.f)

# Everything except the main program goes into a static library, so the unit
# tests in tests/ can link against the same compiled modules as far3d.x.
set(MAIN_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/far3d.f90)
list(REMOVE_ITEM SOURCES ${MAIN_SOURCE})

add_library(far3d_core STATIC ${SOURCES})
set_target_properties(far3d_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_include_directories(far3d_core PUBLIC ${CMAKE_BINARY_DIR}/modules)
target_link_libraries(far3d_core PUBLIC
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
)

# -------------------------
# Add executable
# -------------------------
add_executable(far3d.x ${MAIN_SOURCE})
target_link_libraries(far3d.x PRIVATE far3d_core)

# -------------------------
# Tests (ctest); disable with -DBUILD_TESTING=OFF. See tests/README.md.
# -------------------------
option(BUILD_TESTING "Build the FAR3d unit and regression tests" ON)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
# -------------------------
# Optional: Compiler flags
# -------------------------
# You can use CMake flags set from the outside (e.g., via Makefile), or set defaults here.
# target_compile_options(far3d.x PRIVATE -g -O3 -fopenmp)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE = ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG = ${CMAKE_Fortran_FLAGS_DEBUG}")
