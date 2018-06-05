CC = g++
CFLAGS =   -fopenmp -std=gnu++11 -march=native  -O3 -ffast-math -Wall -DEIGEN_DONT_PARALLELIZE 
CINCLUDE = -I/opt/local/include -I/usr/local/include -ISrc 
CINCLUDE+= -I/usr/local/hdf5/include -I/usr/include/hdf5/serial 
CINCLUDE+= -I$(HOME)/include/eigen3 -I/opt/local/include/eigen3 
CLIBS = -L/usr/local/hdf5/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5  
VPATH = .
OBJS = *.o
DEBUG=0
VERBOSE=1

all:    clean
	cd Src; $(CC) $(CFLAGS) $(CINCLUDE) -c *.cpp -DDEBUG=${DEBUG} -DVERBOSE=${VERBOSE}
	@echo "linking..."
	cd Src; $(CC) $(OBJS) $(CLIBS) $(CFLAGS) -o ../KITEx
	rm -f Src/*.o
	cd ..



clean:
	rm -f Src/*.o Src/*~ core *~
