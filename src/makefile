CC = g++ # Compiler macro
CFLAGS = -std=c++11 -O2 -I /home/users/rodriguezlua/ArmadilloIH/armadillo-9.100.5/bin/include -DARMA_DONT_USE_WRAPPER -lopenblas # List of flags to pass to the compilation command
DEPS = Lines/Lines.hpp BCPs/BCPs.hpp ODE_solvers/ODE_solvers.hpp SEM_math/quadrule.hpp  # Set of .h files on which the .cpp files depend macro
OBJ = main.o SEM_math/quadrule.o BCPs/BCPs.o Lines/Lines_Dyn.o Lines/Lines_QS.o ODE_solvers/ODE_solvers.o  # Object files macro

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

SEM_Lines.exe: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)


# Si el proyecto crece tanto como para organizarlo en varios directorios, mirar el tutorial:
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
