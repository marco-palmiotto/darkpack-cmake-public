#!/bin/bash

distros=( "fedora36" "debian" )

for distro in ${distros[@]}
do
    echo Working on $distro
    name_output=darkpack-launch-$distro

    podman build -t $name_output:latest -f Dockerfile-${distro}-launch

    podman tag localhost/${name_output} docker.io/marcopalmiotto/darkpack-deploy:$distro

    podman push docker.io/marcopalmiotto/darkpack-deploy:$distro
done


# podman save --format oci-archive -o $name_output.tar localhost/$name_output:latest
# sudo apptainer build $name_output.sif $name_output.tar
