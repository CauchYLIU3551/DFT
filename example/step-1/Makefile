CPPFLAGS := -I../../include -I$(HOME)/local/include -I/usr/local/include

CXXFLAGS := -O2 -g -fopenmp

LDFLAGS := -L$(HOME)/local/lib -L/usr/local/lib -L/usr/local/dealii-8.1.0/lib

LDLIBS := ../../lib/libDFT.so -lAFEPack -ldeal_II -ldl -llapack -leigsolver

main : main.o
	g++ -o main main.o $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -rf main.o main *~
