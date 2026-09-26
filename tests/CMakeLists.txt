# FAR3d test suite. Configured from the top-level CMakeLists.txt; run with
#   ctest --output-on-failure            (from the build directory)
#   ctest -L unit | -L regression        (one group)
# See tests/README.md.

# ---------------------------------------------------------------------------
# Fortran unit tests: small programs linked against far3d_core. Each prints
# PASS/FAIL lines and exits nonzero on failure.
# ---------------------------------------------------------------------------
add_library(far3d_testing STATIC unit/testing_mod.f90)
set_target_properties(far3d_testing PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_link_libraries(far3d_testing PUBLIC far3d_core)

foreach(t test_block_tridiag test_radial_grid test_spectral_ops test_numerics_tools)
    add_executable(${t} unit/${t}.f90)
    target_link_libraries(${t} PRIVATE far3d_testing)
    set_target_properties(${t} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
    add_test(NAME unit.${t} COMMAND ${t})
    set_tests_properties(unit.${t} PROPERTIES LABELS unit ENVIRONMENT "OMP_NUM_THREADS=1")
endforeach()

# ---------------------------------------------------------------------------
# Regression tests: Python (stdlib unittest + numpy) scripts that run far3d.x
# on small cases and compare against tests/regression/reference/.
# ---------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter)
set(FAR3D_TEST_PYTHON "${Python3_EXECUTABLE}" CACHE FILEPATH
    "Python interpreter (with numpy) used for the regression tests")

if(FAR3D_TEST_PYTHON)
    execute_process(COMMAND ${FAR3D_TEST_PYTHON} -c "import numpy"
                    RESULT_VARIABLE _numpy_missing OUTPUT_QUIET ERROR_QUIET)
endif()

if(NOT FAR3D_TEST_PYTHON OR _numpy_missing)
    message(WARNING "No Python with numpy found: regression tests disabled. "
                    "Set -DFAR3D_TEST_PYTHON=/path/to/python3 to enable them.")
else()
    set(_env
        "FAR3D_EXE=$<TARGET_FILE:far3d.x>"
        "FAR3D_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "FAR3D_REPO=${PROJECT_SOURCE_DIR}"
        "FAR3D_WORKDIR=${CMAKE_BINARY_DIR}/tests/runs"
        "OMP_NUM_THREADS=1")
    foreach(t test_io test_equilibrium test_linear test_nonlinear)
        add_test(NAME regression.${t}
                 COMMAND ${FAR3D_TEST_PYTHON} -m unittest -v ${t}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/regression)
        set_tests_properties(regression.${t} PROPERTIES
            LABELS regression
            ENVIRONMENT "${_env}"
            TIMEOUT 900)
    endforeach()
endif()
