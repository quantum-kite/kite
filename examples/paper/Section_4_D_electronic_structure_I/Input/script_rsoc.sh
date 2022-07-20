for i in {0..17}; do
    python config.py 1024 4096 0.005 10.0 5.0 78.0 16.0 1.0 0.01 31.0 9.0 $i 0.2
    ./KITEx config.h5 
    ./KITE-tools config.h5 --ARPES -E -10 10 16384 -F 100
    mv arpes.dat arpes${i}.dat
done
./join_files.sh
python process_arpes1.py arpes_joined.dat
