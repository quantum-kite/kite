
        ##############################################################################      
        #                        KITE | Pre-Release version 0.1                      #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018                      #      
        #                                                                            #      
        ##############################################################################      


OS := $(shell uname)

ifeq ($(OS),Darwin)
CC = g++-6  -DEIGEN_DONT_PARALLELIZE -fdiagnostics-color=always  -O2  #-ftree-vectorize  -ftree-vectorizer-verbose=7 -fopt-info-vec-missed
else
CC = g++ -march=native  -O2 
endif	
# This makefile has reduntant paths so that is compiles the code in both ubuntu 16.04 and Mac OSX with homebrew. 
#If you know what you are doing, feel free to edit and remove the unnecessary paths

CFLAGS = -DEIGEN_DONT_PARALLELIZE -fopenmp -std=gnu++11 -Wall -ffast-math 
CINCLUDE = -I/opt/local/lib/ 
CINCLUDE += -ISrc -I/usr/local/hdf5/include -I/usr/local/Cellar/hdf5/1.10.1_2/include/ -I/usr/include/hdf5/serial #HDF5
CINCLUDE += -I$(HOME)/include/eigen3 -I/opt/local/include/eigen3 -I/usr/include/eigen3/ -I/usr/local/include/eigen3 #EIGEN3 
CLIBS = -L/usr/local/hdf5/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/local/Cellar/hdf5/1.10.1_2/lib/ 
CLIBS += -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5  
VPATH = .
OBJS = *.o

compile_main=1
verbose=1
debug=0
estimate_time=1

all:    clean
	cd Src; $(CC) $(CFLAGS) $(CINCLUDE) -DDEBUG=$(debug) -DCOMPILE_MAIN=$(compile_main) -DVERBOSE=$(verbose) -DESTIMATE_TIME=$(estimate_time)  -c *.cpp  
	@echo "linking..."
	cd Src; $(CC) $(OBJS) $(CLIBS) $(CFLAGS) -o ../KITEx
	rm -f Src/*.o
	cd ..

debug:  clean  
	cd Src; $(CC) $(CFLAGS) $(CINCLUDE) $(CDEFS)  -DMEM1=$(MEM1) -DMEM2=$(MEM2) -c *.cpp   -g
	@echo "linking..."
	cd Src; $(CC) $(CLIBS) $(CFLAGS) $(CINCLUDE) $(CDEFS) $(OBJS)  -o ../KITEx -g
	cd ..

clean:
	rm -f Src/*.o Src/*~ core *~
