all:
	g++ ctf.cpp -o ctf -L/mingw64/lib -llapack -lblas
	# g++ ctf.cpp -o ctf -larmadillo -llapack -lblas -pthread
