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
		
SCI_LIBS = -larmadillo -lgfortran -lquadmath "C:/Program Files/HDF5-1.10.5-win64/lib/hdf5.lib"
		
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
	
SCI_LIBS=-lopenblas -lhdf5
	
else ifeq ($(OS),ubuntu)

ifeq ($(USERNAME),feruanos)

LIBS=-lopenblas -larmadillo -lhdf5
		
else ifeq ($(USERNAME),rodriguezlua)
# The following lines are needed to choose the serial version of HDF5
INC_HDF5_DIR="/usr/include/hdf5/serial"
IDIRS=-I$(INC_HDF5_DIR)
LIB_HDF5_DIR="/usr/lib/x86_64-linux-gnu/hdf5/serial"
LDIRS=-L$(LIB_HDF5_DIR)

LIBS=-larmadillo -lhdf5 -lopenblas
endif

endif

ODIR=obj
SDIR=src
BDIR=bin

CC=g++
CFLAGS=$(IDIRS) -std=c++14 -O2 -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5
LIBS=$(SCI_LIBS) -lstdc++fs

_DEPS=Lines/Lines.hpp Bodies/Bodies.hpp Spring/Spring.hpp Hydro/HydroDatabase.hpp BCPs/BCPs.hpp BCPs/Winchies.hpp BCPs/WinchiesController.hpp ODE_solvers/ODE_solvers.hpp SEM_math/quadrule.hpp os_tools.hpp Exceptions/Exception.hpp Simulations/Simulation.hpp MathTools.hpp CommonTools.hpp
DEPS=$(patsubst %,$(SDIR)/%,$(_DEPS))
_OBJS=main.o SEM_math/quadrule.o BCPs/BCPs.o BCPs/Winchies.o BCPs/WinchiesController.o Lines/Lines_Dyn.o Lines/Lines_QS.o Bodies/Bodies.o Spring/Spring.o Hydro/HydroDatabase.o ODE_solvers/ODE_solvers.o os_tools.o Exceptions/Exception.o Simulations/Simulation.o MathTools.o CommonTools.o
OBJS=$(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
all: oasis
	
oasis: $(OBJS)
ifeq ($(OS),Windows_NT)
	$(CC) -static -o $(BDIR)/$@.exe $^ $(LDIRS) $(LIBS)
else ifeq ($(OS),centos)
<<<<<<< HEAD
	$(CC) -static -o $(BDIR)/$@ $^ $(LDIRS) $(LIBS)
=======
	$(CC) -o $(BDIR)/$@ $^ $(LDIRS) $(LIBS)
else ifeq ($(OS),ubuntu)
	$(CC) -o $(BDIR)/$@ $^ $(LDIRS) $(LIBS)
>>>>>>> master
endif

.PHONY: clean

clean:
ifeq ($(OS),Windows_NT)
	del /S/F *.o
	del $(BDIR)\oasis.exe
else ifeq ($(OS),centos)
	find . -name "*.o" -type f -delete
	rm $(BDIR)/oasis
else ifeq ($(OS),ubuntu)
	find . -name "*.o" -type f -delete
	rm $(BDIR)/oasis
endif