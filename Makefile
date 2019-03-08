# some system dependent settings...

ifndef SYSTEM
SYSTEM = linux
endif

ifeq ($(SYSTEM),linux)
DBG=
OPENMP=
FFLAGS=$(DBG) -O3 $(OPENMP) -fPIC
BLAS = -lblas
LAPACK = -llapack
PREO = 
POSTO = $(BLAS) $(LAPACK)
LDFLAGS = -shared
SHELL = /bin/sh
OBJSUF=o
MODSUF=mod
FC = gfortran
MEX=mex
MWRAP=mwrap
MWFLAGS=-c99complex
endif

ifeq ($(SYSTEM),mac)
DBG=
OPENMP=
FFLAGS=$(DBG) -O3 $(OPENMP) -fPIC
PREO = -framework accelerate
POSTO =
LDFLAGS = -shared
SHELL = /bin/sh
OBJSUF=o
MODSUF=mod
FC = gfortran
MEX=mex
MWRAP=mwrap
MWFLAGS=-c99complex
endif

# locations

SRC_DIR = src
BIN_DIR = bin
TEST_DIR = test
TMP_DIR = tmp
MDIR = matlab
MWRAP_DIR = mwrap
RGG_DIR = external/rgg_tools
VPATH = $(SRC_DIR):$(BIN_DIR):$(RGG_DIR)

# names

LIBBASE = chunkie
LIBNAME  = lib$(LIBBASE).so
LIBLINK = -l$(LIBBASE)
#FSOURCES = $(shell echo $(SRC_DIR)/*.f)
#F90SOURCES = $(shell echo $(SRC_DIR)/*.f90)
#OBJS = $(patsubst $(SRC_DIR)/%.f,%.o,$(FSOURCES)) $(patsubst $(SRC_DIR)/%.f90,%.o,$(F90SOURCES))

GATEWAY = $(LIBBASE)gateway
MWRAPFILE = $(LIBBASE)

MODS =
OBJS = 	chunks_ders.o \
	chunks.o \
	legeexps.o \
	chunks_quads.o \
	corners.o \
	qerrfun.o \
	gammanew_eval.o \
	prini.o adapgaus.o pplot.o

LOBJS = $(patsubst %.o,../$(BIN_DIR)/%.o,$(OBJS))
# targets

all: lib

.PHONY : all lib profile release \
  install install-strip uninstall clean distclean setup_dir\
	deepclean

setup_dir: 
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(TMP_DIR)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

lib: setup_dir $(MODS) $(OBJS) 
	cd $(BIN_DIR); $(FC) $(LDFLAGS) -o $(LIBNAME) $(OBJS)

clean:
	cd $(BIN_DIR); rm -f *
	cd $(MWRAP_DIR); rm -f $(GATEWAY).c *.mex*


distclean: clean
	cd $(BIN_DIR); rm -f $(LIBNAME)

printflags: 
	@echo $(DBG) $(OPENMP)

## tests

# note that, without install, we need to point to the library
# at run-time as well... (see the -Wl,-rpath part)
test%: setup_dir lib
	cd $(TMP_DIR); $(FC) $(FFLAGS) -o $@ ../$(TEST_DIR)/$@.f $(PREO) ../$(BIN_DIR)/$(LIBNAME) $(POSTO) -Wl,-rpath=../$(BIN_DIR)
	cd $(TMP_DIR); ./$@

## matlab

$(GATEWAY).c: $(MWRAP_DIR)/$(MWRAPFILE).mw Makefile
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw


# all mex-ing

# easier just to link against object files than library,
# if it's not being installed (I think)
mexfile: $(GATEWAY).c $(OBJS) Makefile
	cd $(MWRAP_DIR); $(MEX) $(GATEWAY).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm
