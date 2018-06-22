OS := $(shell uname)
ifeq ($(OS),Darwin)
CC = g++-6 -DEIGEN_DONT_PARALLELIZE -fdiagnostics-color=always  -O2  #-ftree-vectorize  -ftree-vectorizer-verbose=7 -fopt-info-vec-missed
else
CC = g++ -fopenmp -std=c++11
endif
CINCLUDE = -I/opt/local/lib/ 
CINCLUDE += -ISrc -I/usr/local/hdf5/include -I/usr/local/Cellar/hdf5/1.10.1_2/include/ -I/usr/include/hdf5/serial #HDF5
CINCLUDE += -I$(HOME)/include/eigen3 -I/opt/local/include/eigen3 -I/usr/include/eigen3/ -I/usr/local/include/eigen3 #EIGEN3 
CLIBS = -L/usr/local/hdf5/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/local/Cellar/hdf5/1.10.1_2/lib/ 
CLIBS += -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5  
VPATH = .
OBJS = *.o
DEBUG=1
VERBOSE=1

all:    clean
	cd Src; $(CC) $(CFLAGS) $(CINCLUDE) -c *.cpp -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}
	@echo "linking..."
	cd Src; $(CC) $(OBJS) $(CLIBS) $(CFLAGS) -o ../KITE-tools
	rm -f Src/*.o
	cd ..

clean:
	rm -f Src/*.o Src/*~ core *~
