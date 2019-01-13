echo "Compiling KITEx and KITE-tools. To have the full functionality, you need to have at least version 8 of gcc. If you do not, then please edit this file and set the flag WAVEPACKET to 0"
WAVEPACKET=0

mkdir -p build
mkdir -p tools/build

echo "Compiling KITEx"
cd build
for i in ../Src/*.cpp
do
  echo "Compiling $i"
  g++ -DCOMPILE_WAVEPACKET=$WAVEPACKET -std=gnu++11 $i  -I../Src/ -I/opt/local/include/eigen3 -I$HOME/include -I$HOME/include/eigen3 -I/opt/local/include  -fopenmp -c
  done 
echo "Linking"
g++ -std=gnu++11 *.o   -L$HOME/lib -L/opt/local/lib  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -fopenmp -o Kite-tools
echo "Done."
cd ..


echo "Compiling KITE-tools"
cd tools/build
for i in ../src/*.cpp
do
  echo "Compiling $i"
  g++ -std=gnu++11 $i -I../src/ -I/../src/cond_2order -I/opt/local/include/eigen3 -I$HOME/include -I$HOME/include/eigen3 -I/opt/local/include  -fopenmp -c
  done 

for i in ../src/cond_2order/*.cpp
do
  echo "Compiling $i"
  g++ -std=gnu++11 $i -I../src/ -I/../src/cond_2order -I/opt/local/include/eigen3 -I$HOME/include -I$HOME/include/eigen3 -I/opt/local/include  -fopenmp -c
  done 

echo "Linking"
g++ -std=gnu++11 *.o   -L$HOME/lib -L/opt/local/lib  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -fopenmp -o Kite-tools
echo "Done."






cd ..
