SHELL = /bin/sh

SHIELD_BUILD = /SHiELD/SHiELD_build
include $(SHIELD_BUILD)/site/gnu.mk

FC ?= mpif90
CC ?= mpicc
CXX ?= mpicxx
LD ?= mpif90

ARTIFACTS=$(SHIELD_BUILD)/Build/exec/shield_wrapper_nh.prod.64bit.gnu
FMS_64BIT=$(SHIELD_BUILD)/Build/libFMS/gnu/64bit
FFLAGS += -I$(ARTIFACTS) -I$(FMS_64BIT)
MPI_FLAGS=`$(FC) -show | sed 's/\(-I[^ ]*\|^[^ ]*\)//g'`
NCEPLIBS = -L$(SHIELD_BUILD)/Build/nceplibs/gnu -lbacio  -lsp_d  -lw3emc_d  -lw3nco_d

TEMPLATES_RELATIVE = dynamics_data.F90 flagstruct_data.F90 physics_data.F90 _wrapper.pyx
TEMPLATES_IN = $(addprefix templates/,$(TEMPLATES_RELATIVE))
TEMPLATES_OUT = $(addprefix lib/,$(TEMPLATES_RELATIVE))
WRAPPER_OBJECTS=$(TEMPLATES_OUT:.F90=.o) lib/coupler_lib.o

PROPERTIES_FILES = \
	shield/wrapper/dynamics_properties.json \
	shield/wrapper/flagstruct_properties.json \
	shield/wrapper/physics_properties.json

SETUP_PY_LIBS = -L$(ARTIFACTS) \
                -L$(FMS_64BIT) \
                -lfv3 \
                -lFMS \
                -lgfs \
                $(NCEPLIBS) \
                $(LDFLAGS) \
                $(MPI_FLAGS) \
                -lmvec \
                -lm

$(TEMPLATES_OUT): $(TEMPLATES_IN) $(PROPERTIES_FILES)
	# fill_templates.py needs to be run in the same directory as the outputs
	cd lib && python3 ../fill_templates.py $(TEMPLATES_RELATIVE)

build:$(WRAPPER_OBJECTS)
	SETUP_PY_LIBS="$(SETUP_PY_LIBS)" python3 setup.py build_ext --inplace
	SETUP_PY_LIBS="$(SETUP_PY_LIBS)" python3 setup.py bdist_wheel

install:
	SETUP_PY_LIBS="$(SETUP_PY_LIBS)" python3 setup.py install

# Specify dependency structure of the wrapper object files
lib/physics_data.o: lib/physics_data.F90 lib/dynamics_data.o 
lib/coupler_lib.o: lib/coupler_lib.F90 lib/dynamics_data.o
lib/flagstruct_data.o: lib/flagstruct_data.F90
lib/dynamics_data.o: lib/dynamics_data.F90

test:
	pytest -vv tests/test_all_mpi_requiring.py
	pytest -vv tests/pytest

test_regtest_reset:
	pytest -vv tests/test_all_mpi_requiring.py
	pytest -vv tests/pytest --regtest-reset
