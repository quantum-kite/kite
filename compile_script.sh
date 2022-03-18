#!/bin/bash
echo "Compiling KITEx and KITE-tools. To have the full functionality, you need to have at least version 8 of gcc."
echo "If you do not, you will not be able to run guassian_wavepacket. To enable compiling with this feature, please edit this file and set WAVEPACKET=1" 
echo "By default, the flag WAVEPACKET is set to 0"
WAVEPACKET=0

libraries="-L$HOME/lib -L/opt/local/lib  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5"
headers="-I/opt/local/include/eigen3 -I$HOME/include -I$HOME/include/eigen3 -I/opt/local/include"

# make the directory structure
mkdir -p build
mkdir -p tools/build

echo "Compiling KITEx"
cd build
for i in ../Src/*.cpp; do
  echo "Compiling $i"
  g++ -DCOMPILE_WAVEPACKET=$WAVEPACKET -std=gnu++11 -O2 $i  -I../Src/ $headers  -fopenmp -c
  done 
echo "Linking"
g++ -std=gnu++11 -O2 *.o $libraries -fopenmp -o ../KITEx
echo "Done."
cd ..


echo "Compiling KITE-tools"
cd tools/build
for i in ../src/*.cpp; do
  echo "Compiling $i"
  g++ -std=gnu++11 -O2 $i -I../src/ -I/../src/cond_2order $headers -fopenmp -c
  done 

for i in ../src/cond_2order/*.cpp; do
  echo "Compiling $i"
  g++ -std=gnu++11 -O2 $i -I../src/ -I/../src/cond_2order $headers -fopenmp -c
  done 

echo "Linking"
g++ -std=gnu++11 -O2 *.o $libraries -fopenmp -o ../../KITE-tools
echo "Done."

cd ..
