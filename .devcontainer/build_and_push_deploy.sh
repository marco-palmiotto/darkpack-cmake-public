#!/bin/bash

distros=( "fedora36" "debian" )
suffix=( "dev" "scan" "deploy") #deploy must be the last one

for distro in ${distros[@]}
do
    for suffix in ${suffixes[@]}
    do
        name_output=darkpack-$suffix-$distro
        filename=Dockerfile-${distro}-$suffix
        echo Working on $filename

        [ -f $filename ] && podman build -t $name_output:latest -f $filename
    done

    # [ -f $filename ] && podman tag localhost/${name_output} docker.io/marcopalmiotto/darkpack-deploy:$distro

    # [ -f $filename ] && podman push docker.io/marcopalmiotto/darkpack-deploy:$distro

    # To save the image in a sif file
    # podman save --format oci-archive -o $name_output.tar localhost/$name_output:latest
    # sudo apptainer build $name_output.sif $name_output.tar
done



