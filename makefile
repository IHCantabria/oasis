CC = g++               # Compiler macro
CFLAGS = -I.           # List of flags to pass to the compilation command macro (. -> directorio)
DEPS = Catline.h       # Set of .h files on which the .cpp files depend macro
OBJ = Moorings.o setLine.o qs.o  # Object files macro

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Moorings.exe: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)


# Si el proyecto crece tanto como para organizarlo en varios directorios, mirar el tutorial:
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/