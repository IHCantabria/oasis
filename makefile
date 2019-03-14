CC = g++                # Compiler macro
CFLAGS = -std=c++11 -O2 -I /home/users/rodriguezlua/ArmadilloIH/armadillo-9.100.5/bin/include -DARMA_DONT_USE_WRAPPER -lopenblas # List of flags to pass to the compilation command macro (. -> directorio)
DEPS = classes.h quadrule.hpp math_lib.h        # Set of .h files on which the .cpp files depend macro
OBJ = main.o quadrule.o math_lib.o SetLine_lib.o QS_lib.o BCP_lib.o  # Object files macro

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Moorings.exe: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)


# Si el proyecto crece tanto como para organizarlo en varios directorios, mirar el tutorial:
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/