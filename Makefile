CXX = g++
CPPFLAGS=-std=c++11 -O3
CPPFLAGSMAIN=-std=c++11 -O3 -fopenmp
CFLAGS= -O3  -fopenmp
INCLFLAGS=-I ./eigen-eigen-5a0156e40feb/

dieta: main.o parser.o hamiltonian.o lanczos.o greens.o vmatrix.o cpt.o gauss_legendre.o density.o meanfield.o potential.o bzIntegrator.o
	$(CXX) $(INCLFLAGS) $(CPPFLAGSMAIN) -o dieta main.o parser.o hamiltonian.o lanczos.o greens.o vmatrix.o cpt.o gauss_legendre.o density.o meanfield.o potential.o bzIntegrator.o
	
main.o: main.cxx parser.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c main.cxx
	
parser.o: parser.cxx parser.h
	$(CXX) $(INCLFLAGS) $(CPPFLAGS) -c parser.cxx
	
hamiltonian.o: hamiltonian.cxx hamiltonian.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGSMAIN) -c hamiltonian.cxx
	
lanczos.o: lanczos.cxx lanczos.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c lanczos.cxx
	
greens.o: greens.cxx greens.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c greens.cxx
	
vmatrix.o: vmatrix.cxx vmatrix.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c vmatrix.cxx
	
cpt.o: cpt.cxx cpt.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c cpt.cxx
	
gauss_legendre.o: gauss_legendre/gauss_legendre.c gauss_legendre/gauss_legendre.h
	$(CXX) $(INCLFLAGS) $(CFLAGS) -c gauss_legendre/gauss_legendre.c
	
density.o: density.cxx density.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c density.cxx
	
potential.o: potential.cxx potential.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c potential.cxx
	
meanfield.o: meanfield.cxx meanfield.h
	$(CXX) $(INCLFLAGS)  $(CPPFLAGS) -c meanfield.cxx
	
bzIntegrator.o: bzIntegrator.cxx bzIntegrator.h
	$(CXX) $(INCLFLAGS) $(CPPFLAGS) -c bzIntegrator.cxx
	
clean: 
	rm -f *.o dieta
	echo Clean done
