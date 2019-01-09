#!/bin/bash
 
# This script will compare the .h5 files of all the tests


echo "Objects:    dif_norm    max_dif      norm1      norm2       %"
for i in test*; do
    cat $i/description
    if [[ "$1" == "redo" ]]; then
        ./KITEx $i/config.h5 > /dev/null
    fi

    while read line; do 
        #echo $line; 
        file=$line
        a=$(basename $line)
        printf "%-12s" "$a"


        python compare.py $i/configREF.h5 $file $i/config.h5 $file
    done < $i/function

done
