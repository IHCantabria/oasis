# Check OS enviroment
ifeq ($(OS),Windows_NT)
# Check USER environment
ifeq ($(USERNAME),feruanos)
INC_ARMADILLO_DIR = "C:/ScientificLibraries/cpp/armadillo940/include"
INC_HDF5_DIR = "C:/Program Files/HDF5-1.10.5-win64/include"
IDIRS = -I$(INC_ARMADILLO_DIR) -I$(INC_HDF5_DIR)
		
LIB_LAPACK_DIR = "C:/ScientificLibraries/fortran90/lapack380"
LIB_BLAS_DIR = "C:/ScientificLibraries/fortran90/blas380"
LIB_ARMADILLO_DIR = "C:/ScientificLibraries/cpp/armadillo940/lib"
LIB_HDF5_DIR = "C:/Program Files/HDF5-1.10.5-win64/lib"
LDIRS = -L$(LIB_ARMADILLO_DIR) -L$(LIB_HDF5_DIR)
		
LIBS = -larmadillo -lgfortran "C:/Program Files/HDF5-1.10.5-win64/lib/hdf5.lib"
		
else ifeq ($(USERNAME),rodriguezlua)
LIB_LAPACK_DIR = "C:/ScientificLibraries/fortran90/lapack380"
LIB_BLAS_DIR = "C:/ScientificLibraries/fortran90/blas380"
endif
else ifeq ($(OS),centos)
INC_ARMADILLO_DIR="/home/projects/energia/ArmadilloIH/armadillo-9.100.5/include"
INC_HDF5_DIR=$(HDF5_DIR)/include
IDIRS=-I$(INC_ARMADILLO_DIR) -I$(INC_HDF5_DIR)
	
LIB_OPENBLAS_DIR=-L$(EBROOTOPENBLAS)/lib
LDIRS=$(LIB_OPENBLAS_DIR)
	
LIBS=-lopenblas -lhdf5
	
else
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
ifeq ($(OS),Windows_NT)
	$(CC) -o $(BDIR)/$@.exe $^ $(LDIRS) $(LIBS)
else ifeq ($(OS),centos)
	$(CC) -o $(BDIR)/$@ $^ $(LDIRS) $(LIBS)
endif

.PHONY: clean

clean:
ifeq ($(OS),Windows_NT)
	del /S/F *.o
	del $(BDIR)\oasis.exe
else ifeq ($(OS),centos)
	find . -name "*.o" -type f -delete
	rm $(BDIR)/oasis
endif