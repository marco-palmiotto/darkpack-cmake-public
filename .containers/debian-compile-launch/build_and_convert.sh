#!/bin/bash

podman build -t darkpack-launch:latest .
podman save --format oci-archive -o darkpack-launch.tar localhost/darkpack-launch:latest
sudo apptainer build darkpack-launch.sif darkpack-launch.tar
scp darkpack-launch.sif mpalmiot@cca.in2p3.fr:/pbs/home/m/mpalmiot/code/
