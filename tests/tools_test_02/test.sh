#!/bin/bash

tools=../KITE-tools

#echo Testing DoS for the same Hamiltonian with and without energy shift.
$tools config_sym.h5 --DOS -N dos_sym.dat -X > log
$tools config.h5 --DOS -N dos.dat -X >> log

head -n 900 dos.dat > dos1.dat
tail -n 900 dos_sym.dat > dos2.dat

result=$(paste dos1.dat dos2.dat | awk '
BEGIN {max=0; ORS=" "} 
{sum += ($2-$5)^2} 
sqrt(($2-$5)^2)>max {max = sqrt(($2-$5)^2)} 
END { print max; print sqrt(sum); ORS="\n"; print ""}')

printable_result=$(echo $result | awk '{print $1}')
rm dos1.dat dos2.dat

acceptable=$(echo $result | awk '{
if($1<0.02)
  print 1;
else
  print 0;
}')

if [[ "$acceptable" == "1" ]]; then 
  echo -e "\e[1m\e[32mOK\e[0m"
else 
  echo -e "\e[1m\e[91m$printable_result\e[0m"
fi
