PYTHON_VERSION = python2.7
PYTHON_LANG_LEVEL = 2
PYTHON_INCLUDE_PATH = /Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PYTHON_LD_PATH = /Library/Frameworks/Python.framework/Versions/2.7/lib
CYTHONFLAGS = -DNPY_NO_DEPRECATED_API
OUTEXT= o
LIBEXT= so

NUMPY_INCLUDE_PATH = /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/

ifeq ($(OS),Windows_NT)
REMOVE = del
MOVE = move
else
REMOVE = rm -rf
MOVE = mv
endif

# cython compiler PYX->C
CYTHON = cython -$(PYTHON_LANG_LEVEL)
# gnu c compiler C->O
#CC=gcc
CC = /usr/local/opt/llvm/bin/clang
#CCFLAGS = -O -fopenmp -DAdd_ $(CYTHONFLAGS) -fPIC
CCFLAGS = -I/usr/local/opt/llvm/include -fopenmp
# gnu fortran compiler (linker) O->SO
FL=gfortran
#FLFLAGS = -O -fopenmp -shared -fPIC
FLFLAGS = -L/usr/local/opt/llvm/lib
FLLIBS = -ldmumps -lmumps_common -lpord -llapack -lmpiseq -lblas -lpthread

MUMPS_INCLUDE_PATH = /usr/local/opt/brewsci-mumps/lib
MUMPS_LIBRARY_PATH = /usr/local/opt/brewsci-mumps/include
MUMPS_LIBRSEQ_PATH = /usr/local/opt/brewsci-mumps/libexec/include
LAPACK_LIBRARY_PATH = /usr/local/opt/lapack/lib

all: build_mumps

build_mumps:
	@echo "Building dMUMPS Inteface module"
	$(CYTHON) _dMUMPS_.pyx
	$(CC) $(CCFLAGS) -c _dMUMPS_.c -o _dMUMPS_.$(OUTEXT) -I. -I$(PYTHON_INCLUDE_PATH) -I$(NUMPY_INCLUDE_PATH) -L$(PYTHON_LD_PATH) -l$(PYTHON_VERSION) -I$(MUMPS_INCLUDE_PATH) -I$(MUMPS_LIBRSEQ_PATH)
	$(FL) -o _dMUMPS_.$(LIBEXT) $(FLFLAGS) _dMUMPS_.$(OUTEXT) -L$(MUMPS_LIBRARY_PATH) -L$(LAPACK_LIBRARY_PATH) -L$(MUMPS_LIBRSEQ_PATH) $(FLLIBS)
	$(MOVE) _dMUMPS_.$(LIBEXT) ../

clean:
	$(REMOVE) ../_dMUMPS_.$(LIBEXT) _dMUMPS_.$(OUTEXT) _dMUMPS_.c

source_clean:
	$(REMOVE) _dMUMPS_.$(OUTEXT) _dMUMPS_.c

