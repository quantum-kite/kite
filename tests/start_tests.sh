#!/bin/bash
 
# This script will compare the .h5 files of all the tests


S=3  # the seed that will be used by KITE



echo "Objects:    dif_norm    max_dif      norm1      norm2       %"
for i in test*; do
    cd $i
    cat description
    if [[ "$1" == "redo" ]]; then
        cp configORIG.h5 config.h5
        chmod 755 config.h5
    else 
        if [[ "$1" == "script" ]]; then
            python config.py > /dev/null
        fi
    fi
    cd ..
    SEED=$S ../build/KITEx $i/config.h5 > /dev/null

    while read line; do 
        file=$line
        a=$(basename $line)
        printf "%-12s" "$a"


        python compare.py $i/configREF.h5 $file $i/config.h5 $file
    done < $i/function

done
