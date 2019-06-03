
# Check OS enviroment
ifeq ($(OS),Windows_NT)
	# Check USER environment
	ifeq ($(USERNAME),feruanos)
		INC_ARMADILLO_DIR = "C:/ScientificLibraries/cpp/armadillo940/include"
		INC_HDF5_DIR = "C:/Program Files/HDF5-1.10.5-win64/include"
		IDIRS = -I$(INC_ARMADILLO_DIR) -I$(INC_HDF5_DIR)
		
		LIB_LAPACK_DIR = "C:/ScientificLibraries/fortran90/lapack380"
		LIB_BLAS_DIR = "C:/ScientificLibraries/fortran90/blas380"
		LIB_HDF5_DIR = "C:/Program Files/HDF5-1.10.5-win64/lib"
		LDIRS = -L$(LIB_LAPACK_DIR) -L$(LIB_BLAS_DIR) -L$(LIB_HDF5_DIR)
		
		LIBS = -llapack -lblas -lgfortran "C:/Program Files/HDF5-1.10.5-win64/lib/hdf5.lib"
		
	else ifeq ($(USERNAME),rodriguezlua)
		LIB_LAPACK_DIR = "C:/ScientificLibraries/fortran90/lapack380"
		LIB_BLAS_DIR = "C:/ScientificLibraries/fortran90/blas380"
	endif
else
	VAR_ECHOED="Cluster environment"
	USER_NAME=$(USER)
endif

ODIR=obj
SDIR=src
BDIR=bin

CC=g++
CFLAGS=$(IDIRS) -std=c++14 -O2 -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5

_DEPS=Lines/Lines.hpp Bodies/Bodies.hpp Spring/Spring.hpp Hydro/Hydro.hpp BCPs/BCPs.hpp BCPs/Winchies.hpp BCPs/WinchiesController.hpp ODE_solvers/ODE_solvers.hpp SEM_math/quadrule.hpp
DEPS=$(patsubst %,$(SDIR)/%,$(_DEPS))
_OBJS=main.o SEM_math/quadrule.o BCPs/BCPs.o BCPs/Winchies.o BCPs/WinchiesController.o Lines/Lines_Dyn.o Lines/Lines_QS.o Bodies/Bodies.o Spring/Spring.o Hydro/Hydro.o ODE_solvers/ODE_solvers.o
OBJS=$(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
all: oasis
	
oasis: $(OBJS)
	@echo $(LDIRS)
	$(CC) -o $(BDIR)/$@.exe $^ $(LDIRS) $(LIBS)

.PHONY: clean

clean:
ifeq ($(OS),Windows_NT)
	del /S/F *.o
endif
	
