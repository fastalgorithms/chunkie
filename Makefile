# some system dependent settings...

ifndef SYSTEM
SYSTEM = linux
endif

ifeq ($(SYSTEM),linux)
DBG=
OPENMP= -fopenmp
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

# locations

SRC_DIR = src
BIN_DIR = bin
TEST_DIR = test
TMP_DIR = tmp
MDIR = chunkie
DEV_DIR = dev
MWRAP_DIR = mwrap
RGG_DIR = external/rgg_tools
VPATH = $(SRC_DIR):$(BIN_DIR):$(RGG_DIR)

MEX_DEST = $(DEV_DIR)

# names

BASE = utils
GATEWAY = $(BASE)gateway
MWRAPFILE = $(BASE)

BASE2 = fasthank
GATEWAY2 = $(BASE2)gateway
MWRAPFILE2 = $(BASE2)



MODS =
OBJS = 	chunks_ders.o \
	chunks.o \
	legeexps.o \
	chunks_quads.o \
	corners.o \
	qerrfun.o \
	gammanew_eval.o \
	prini.o adapgaus.o pplot.o \
	alpert_quads.o \
	hank103.o

LOBJS = $(patsubst %.o,../$(BIN_DIR)/%.o,$(OBJS))
# targets

all: mexfile

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

clean:
	cd $(BIN_DIR); rm -f *
	cd $(MWRAP_DIR); rm -f $(GATEWAY).c *.mex* *.m


distclean: clean
	cd $(BIN_DIR); rm -f $(LIBNAME)

printflags: 
	@echo $(DBG) $(OPENMP)

## matlab

$(GATEWAY).c: $(MWRAP_DIR)/$(MWRAPFILE).mw Makefile
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw

$(GATEWAY2).c: $(MWRAP_DIR)/$(MWRAPFILE2).mw Makefile
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb $(MWRAPFILE2).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw


# all mex-ing

# easier just to link against object files than library,
# if it's not being installed (I think)
mexfile: $(GATEWAY).c $(OBJS) Makefile
	cd $(MWRAP_DIR); $(MEX) $(GATEWAY).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm
	mv $(MWRAP_DIR)/*.m $(MEX_DEST)/
	mv $(MWRAP_DIR)/*.mex* $(MEX_DEST)/


NAMESPACESTR = fasthank
SEDSTRING = s/$(GATEWAY2)/$(NAMESPACESTR).$(GATEWAY2)/g
mexfile2: $(GATEWAY2).c $(OBJS) Makefile
	mkdir -p $(MDIR)/+$(NAMESPACESTR)
	cd $(MWRAP_DIR); $(MEX) $(GATEWAY2).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm -D_OPENMP
	cd $(MWRAP_DIR); sed -i '$(SEDSTRING)' *.m
	mv $(MWRAP_DIR)/*.m $(MDIR)/+$(NAMESPACESTR)/
	mv $(MWRAP_DIR)/*.mex* $(MDIR)/+$(NAMESPACESTR)/
