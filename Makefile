# MODEL=ISING|SG L=8|16... make -e
# ENTROPY=CFP|BERG|EXACT DINAMICS=SSF|NFOLD|CLUSTER ./comando
# -std=c99 

OS := $(shell uname)

ifeq ($(OS),Darwin)
CC = g++-mp-6  -DEIGEN_DONT_PARALLELIZE -Wall  -fdiagnostics-color=always -march=native  -O2 #-ffast-math #-ftree-vectorize  -ftree-vectorizer-verbose=7 -fopt-info-vec-missed
else
CC = g++ -Wall -DEIGEN_DONT_PARALLELIZE -g -march=native  -O3 -ffast-math 

endif	

CFLAGS =   -fopenmp -std=gnu++11
CINCLUDE = -I/opt/local/lib/ -ISrc -I/usr/local/hdf5/include -I/homes/jvlopes/include/eigen3 -I/opt/local/include/eigen3  
CLIBS = -L/usr/local/hdf5/lib/ -lhdf5 -lhdf5_cpp
VPATH = .
OBJS = *.o
EFERMI = 0.0
REDUCTION  = 1


all:    clean
	cd Src; $(CC) $(CFLAGS) $(CINCLUDE) -c *.cpp  
	@echo "linking..."
	cd Src; $(CC) $(CLIBS) $(CFLAGS) $(CINCLUDE) $(OBJS)  -o ../pp
	rm -f Src/*.o
	cd ..



debug:  clean  
	cd Src; $(CC) $(CFLAGS) $(CINCLUDE) $(CDEFS)  -DMEM1=$(MEM1) -DMEM2=$(MEM2) -c *.cpp   -g
	@echo "linking..."
	cd Src; $(CC) $(CLIBS) $(CFLAGS) $(CINCLUDE) $(CDEFS) $(OBJS)  -o ../pp -g
	cd ..

clean:
	rm -f Src/*.o Src/*~ core *~


tratar: clean
	$(CC) -Wall $(CLIBS) $(CFLAGS) $(CINCLUDE) $(CDEFS)  -DMEM1=$(MEM1) -DMEM2=$(MEM2) -DREDUCTION=$(REDUCTION) -DEFERMI=$(EFERMI) tratar.cpp  -o tratar
