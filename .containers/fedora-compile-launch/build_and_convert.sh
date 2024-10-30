#!/bin/bash

podman build -t fedora36:dp .
podman save --format oci-archive -o fedora36.tar localhost/fedora36:dp
sudo apptainer build fedora36.sif oci-archive://fedora36.tar
scp fedora36.sif mpalmiot@cca.in2p3.fr:/pbs/home/m/mpalmiot/code/
