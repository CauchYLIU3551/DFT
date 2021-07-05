src := $(wildcard *.cpp)
objects := $(patsubst %.cpp, %.o, $(wildcard *.cpp))

CXX:=mpic++
#CXX:=g++

CPPFLAGS := -I$(HOME)/local/include -I$(HOME)/local/include/deal.II 
LDFLAGS := -L$(HOME)/local/lib
LIBS := -lAFEPack -ldeal_II -ldl  -lCG
CXXFLAGS := -O2 -fPIC -g -fopenmp

all : $(objects)
	g++ -shared -fPIC -o ../lib/libDFT.so $(objects)

clean:
	rm -rf *~ *.o *.so
