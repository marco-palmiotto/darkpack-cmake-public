docker run --rm --user $(id -u):$(id -g)  --mount type=bind,src=.,target=/home/darkpack  -w /home/darkpack  docker.io/marcopalmiotto/darkpack-deploy:fedora36  /bin/bash -c "./run_CI.sh"     
