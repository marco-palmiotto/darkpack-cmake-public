#!/bin/bash

distros=( "fedora36" "debian" "ubuntu2204" "fedora41" "fedora41gcc13")
suffixes=( "dev" "scan" "deploy") #deploy must be the last one

for distro in ${distros[@]}
do
    for suffix in ${suffixes[@]}
    do
        name_output=darkpack-$suffix-$distro
        filename=Dockerfile-${distro}-$suffix
        echo Working on $filename

        [[ ! -f $filename ]] && echo "File $filename not found" && continue

        podman build -t $name_output:latest -f $filename

        if [ $suffix == "deploy" -o $suffix == "scan" ]
        then
            podman tag localhost/${name_output} docker.io/marcopalmiotto/darkpack-$suffix:$distro
            podman push docker.io/marcopalmiotto/darkpack-$suffix:$distro &
        fi
    done
done

wait
