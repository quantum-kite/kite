file1=

# Search for the list of energies in one of the files. All of the files have the same list of energies
awk '/Energies/,/ARPES/' arpes0.dat | tail -n +2 | head -n -1 > energies.dat

# Search for the list of k-vectors and the arpes matrix
# remove the last column of each file so they stitch together nicely
for i in {0..17}; do
  awk '/k-vectors/,/Energies/' arpes${i}.dat | head -n -2 | tail -n +2 >> file_k.dat
  awk '/ARPES/,/G/' arpes${i}.dat | tail -n +2 | rev | cut -d" " -f2- | rev > file${i}.dat
done

paste file{0..17}.dat > joined.dat

final=arpes_joined.dat
echo "k-vectors:" > $final
cat file_k.dat >> $final
echo "Energies:" >> $final
cat energies.dat >> $final
echo "ARPES:" >> $final
cat joined.dat >> $final

sed -i 's/\t/ /g' $final
