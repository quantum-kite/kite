#!/bin/bash
 
# This script will compare the .h5 files of all the tests
printall=0


echo -e "Testing \e[1mKITEx\e[0m functionalities. \e[1m\e[32mOK\e[0m means the test passed perfectly."
for i in test0*; do
    cd $i
    echo -n "  "
    cat description | tr '\n' ' '
    ./test.sh quick
    cd ..
done
echo -e "\n"

