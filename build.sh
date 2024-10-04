#! /bin/bash

nthreads=4
avgload=3
nmaxfailures=2 # 0 is infinity

# old option was -l1 -k2 works immediately, but very slow
ninja -C build -l$avgload -k$nmaxfailures
exit_code=$?

echo exit_code=$exit_code

while [[ $exit_code -ne 0 && $nthreads -ne 0 ]]
do
    echo "Build failed, trying with less threads"
    [ $avgload -ge 2 ]  && avgload=$(($avgload-1))
    nthreads=$(($nthreads-1))
    echo "Running: ninja -C build -l$avgload -k$nmaxfailures -j$nthreads"
    ninja -C build -l$avgload -k$nmaxfailures -j$nthreads
    exit_code=$?
done

[[ $exit_code -ne 0 ]] && echo "Building failed" || echo "Build succeded"