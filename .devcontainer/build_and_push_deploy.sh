#!/bin/bash

distros=( "fedora36" "debian" )
suffix=deploy

for distro in ${distros[@]}
do
    echo Working on $distro
    name_output=darkpack-$suffix-$distro

    podman build -t $name_output:latest -f Dockerfile-${distro}-$suffix

    podman tag localhost/${name_output} docker.io/marcopalmiotto/darkpack-deploy:$distro

    podman push docker.io/marcopalmiotto/darkpack-deploy:$distro

    # To save the image in a sif file
    # podman save --format oci-archive -o $name_output.tar localhost/$name_output:latest
    # sudo apptainer build $name_output.sif $name_output.tar
done



