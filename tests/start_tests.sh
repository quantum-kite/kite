#!/bin/bash
 
# This script will compare the .h5 files of all the tests
printall=0


echo -e "Testing \e[1mKITEx\e[0m functionalities. \e[1m\e[32mOK\e[0m means the test passed perfectly."
for i in test*; do
    cd $i
    echo -n "  "
    cat description | tr '\n' ' '
    result=$(./test.sh script | awk 'BEGIN {sum=0} {sum+=$5} END {print sum}')
    if [[ "$result" == "0" ]]; then 
      echo -e "\e[1m\e[32mOK\e[0m"
    else 
      echo -e "\e[1m\e[91m$result\e[0m"
    fi
    cd ..
done
echo -e "\n"

echo -e "Testing \e[1mKITE-tools\e[0m functionalities. \e[1m\e[32mOK\e[0m means the test passed."
for i in tools_test*; do
    cd $i
    echo -n "  "
    cat description | tr '\n' ' '
    ./test.sh script 
    cd ..
done
