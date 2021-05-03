# Remove reference to mingw64/lib if compiling using Linux system.

LIBS = -L/mingw64/lib -llapack -lblas

all:
	g++ ctf.cpp -o ctf $(LIBS)
