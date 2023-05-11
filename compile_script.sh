echo "Compiling KITEx and KITE-tools. To have the full functionality, you need to have at least version 8 of gcc. If you do not, you will not be able to run guassian_wavepacket. To enable compiling with this feature, please edit this file and set WAVEPACKET=1. By default, the flag WAVEPACKET is set to 0"
WAVEPACKET=0

# create the directory structure
mkdir -p build
mkdir -p tools/build

# Common locations for the eigen headers. Edit this line if eigen is not in any of these places
eigen="-I/opt/local/include/eigen3 -I$HOME/include/eigen3 -I/usr/include/eigen3"

# Library and header locations
libraries="-L$HOME/lib -L/opt/local/lib  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5"
headers="-I$HOME/include -I/opt/local/include $eigen" 

# Compile KITEx
include="-I../Src -I../Src/Hamiltonian -I../Src/Tools -I../Src/Vector -I../Src/Simulation -I ../Src/Lattice"
sources="../Src/*.cpp ../Src/*/*.cpp"

echo "Compiling KITEx"
cd build
for i in $sources; do
  echo "Compiling $i"
  g++ -DCOMPILE_WAVEPACKET=$WAVEPACKET -std=gnu++11 -O2 $i  $include $headers  -fopenmp -c
done 
echo "Linking"
g++ -std=gnu++11 -O2 *.o $libraries -fopenmp -o ../KITEx
echo "Done."
cd ..

# Compile KITE-tools
include="-I../Src -I../Src/OptCond_1order -I../Src/OptCond_2order -I../Src/Spectral -I../Src/Tools -I ../Src/CondDC"
sources="../Src/*.cpp ../Src/*/*.cpp"

echo "Compiling KITE-tools"
cd tools/build
for i in $sources; do
  echo "Compiling $i"
  g++ -std=gnu++11 -O2 $i $include $headers -fopenmp -c
done 
echo "Linking"
g++ -std=gnu++11 -O2 *.o $libraries -fopenmp -o ../../KITE-tools
echo "Done."
cd ../..
