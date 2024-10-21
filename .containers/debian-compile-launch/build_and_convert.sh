#!/bin/bash

docker build -t darkpack-launch:latest . && \
docker save -o darkpack-launch.tar darkpack-launch:latest && \
sudo apptainer build darkpack-launch.sif docker-archive://darkpack-launch.tar
