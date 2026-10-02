MAKE_CMD := make
ifdef VERBOSE
  MAKE_CMD += VERBOSE=1
endif

# Build directories
BUILD_RELEASE := ../build-release
BUILD_DEBUG   := ../build-debug
BUILD_GPU     := ../build-gpu

# GPU build flags (NVIDIA HPC SDK). GPU_FC defaults to the NVHPC-bundled MPI
# wrapper; override with e.g. `make gpu GPU_FC=ftn GPU_CC=cc` on Cray systems.
GPU_FC      ?= mpif90
GPU_CC      ?= mpicc
GPU_FFLAGS  := -O3 -mp -acc=gpu -gpu=managed -Minstrument -traceback -cuda
# Linked only into Fortran executables: libnvhpcwrapnvtx needs the Fortran
# runtime, so passing it via CMAKE_EXE_LINKER_FLAGS breaks CMake's C compiler
# check. (-cuda reaches the link line through GPU_FFLAGS.)
GPU_FLIBS   := -lnvhpcwrapnvtx

.PHONY: all debug gpu test perf clean

# Python (with numpy) for the regression tests: make test PYTHON=/path/to/python3
PYTHON_FLAG := $(if $(PYTHON),-DFAR3D_TEST_PYTHON=$(PYTHON))

all: $(BUILD_RELEASE)
	@cd $(BUILD_RELEASE) && \
	cmake -DCMAKE_BUILD_TYPE=Release $(PYTHON_FLAG) \
	  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
	  -DCMAKE_Fortran_COMPILER=mpif90 .. && \
	$(MAKE_CMD)

debug: $(BUILD_DEBUG)
	@cd $(BUILD_DEBUG) && \
	cmake -DCMAKE_BUILD_TYPE=Debug .. && \
	$(MAKE_CMD)

gpu: $(BUILD_GPU)
	@cd $(BUILD_GPU) && \
	export OMPI_MCA_opal_cuda_support=0 && \
	cmake -DCMAKE_BUILD_TYPE=Release -DFAR3D_GPU=ON \
	  -DCMAKE_C_COMPILER=$(GPU_CC) \
	  -DCMAKE_Fortran_COMPILER=$(GPU_FC) \
	  -DCMAKE_Fortran_FLAGS="$(GPU_FFLAGS)" \
	  -DCMAKE_Fortran_STANDARD_LIBRARIES="$(GPU_FLIBS)" \
	  .. && \
	$(MAKE_CMD)

# Build the release configuration, then run the unit and regression tests
# (tests/README.md). Pass CTEST_ARGS, e.g. make test CTEST_ARGS="-L unit".
test: all
	@cd $(BUILD_RELEASE) && ctest --output-on-failure $(CTEST_ARGS)

# Build the release configuration, then run the performance benchmarks
# (tests/performance/README.md) on it. Pass options in PERF_ARGS, e.g.
#   make perf PYTHON=/opt/homebrew/bin/python3.11 PERF_ARGS="--np 1 2 4 --label mymachine"
# Benchmark the GPU build with run_benchmarks.py --exe ../build-gpu/far3d.x directly.
perf: all
	@cd ../tests/performance && $(or $(PYTHON),python3) run_benchmarks.py \
	  --exe $(abspath $(BUILD_RELEASE))/far3d.x $(PERF_ARGS)

# Directory creation rules
$(BUILD_RELEASE) $(BUILD_DEBUG) $(BUILD_GPU):
	@mkdir -p $@

clean:
	@echo "Cleaning build folders"
	@rm -rf $(BUILD_RELEASE) $(BUILD_DEBUG) $(BUILD_GPU)
	@rm -f far3d.x

