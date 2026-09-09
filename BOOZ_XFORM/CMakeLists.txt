cmake_minimum_required (VERSION 3.29)

project (booz_xform Fortran)

# Find MPI
find_package(MPI REQUIRED)
if(MPI_FOUND)
    include_directories(SYSTEM ${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS} -cpp -DLINUX -DMPI_OPT")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
else()
    message(FATAL_ERROR "MPI not found")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Manually specify NetCDF paths
set(NETCDF_INC_PATH "/opt/local/include")
set(NETCDF_LIB_PATH "/opt/local/lib")
set(NetCDF_Fortran_LIBRARY "/opt/local/lib/libnetcdff.dylib")
message(STATUS "Using netcdf")
add_definitions(-DNETCDF)
include_directories(${NETCDF_INC_PATH})
link_directories(${NETCDF_LIB_PATH})

add_executable (xbooz_xform)

target_link_libraries (xbooz_xform PUBLIC ${MPI_Fortran_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran ${NetCDF_Fortran_LIBRARY})

add_subdirectory (Sources)

install (TARGETS xbooz_xform)
