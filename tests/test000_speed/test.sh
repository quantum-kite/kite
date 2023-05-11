/usr/bin/time -v ../KITEx config.h5 &> out.dat

cat out.dat | grep Elapsed | awk '{print $8}'

# A normal time is 1.5s for 512x512 with 512 polys
