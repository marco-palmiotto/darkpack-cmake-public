#!/bin/bash

podman build -t darkpack-launch:latest .
# podman save -o darkpack-launch.tar darkpack-launch:latest
#&& sudo apptainer build darkpack-launch.sif docker-archive://darkpack-launch.tar
