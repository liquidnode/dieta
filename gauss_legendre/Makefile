CXX = g++
CPPFLAGS=-O3
INCLFLAGS=-I ../eigen-eigen-5a0156e40feb/

main: example.o gauss_legendre.o
	$(CXX) $(INCLFLAGS) $(CPPFLAGS) -o main example.o gauss_legendre.o
	
example.o: example.c
	$(CXX) $(INCLFLAGS) $(CPPFLAGS) -c example.c
	
gauss_legendre.o: gauss_legendre.c gauss_legendre.h
	$(CXX) $(INCLFLAGS) $(CPPFLAGS) -c gauss_legendre.c
	
clean: 
	rm -f *.o main
	echo Clean done
