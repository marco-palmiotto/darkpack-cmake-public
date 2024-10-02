#! /bin/bash

# old option was -l1 -k2 works immediately, but very slow
ninja -C build -l1 
exit_code=$?

if [[ !$exit_code -ne 0 ]]
then
    echo "First build failed"
    nthreads=4
    avgload=2
    nmaxfailures=0 # 0 is infinity
fi

while [[ $exit_code -ne 0 && $nthreads -ne 0 ]]
do
    echo "Build failed, trying with less threads"
    echo "Running: ninja -C build -l$avgload -k$nmaxfailures -j$nthreads"
    ninja -C build -l$avgload -k$nmaxfailures -j$nthreads
    exit_code=$?

    [ $avgload -ge 2 ]  && avgload=$(($avgload-1))
    nthreads=$(($nthreads-1))
done

[[ !$exit_code ]] && echo "Building failed" || echo "Build succeded"