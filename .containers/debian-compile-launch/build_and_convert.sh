#!/bin/bash

name_output=darkpack-launch-debian

podman build -t $name_output:latest ./Dockerfile-debian-launch && \
podman save --format oci-archive -o $name_output.tar localhost/$name_output:latest
# sudo apptainer build $name_output.sif $name_output.tar
