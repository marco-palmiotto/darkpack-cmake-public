# # Limita la memoria virtuale a 8GB
# ulimit -v 4000000
#
# # Limita il numero di processi a 100
# ulimit -u 10
#
# # Limita l'utilizzo della CPU a 10 minuti
# ulimit -t 600

sudo apptainer build darkpack-launch.sif docker-archive://darkpack-launch.tar
scp darkpack-launch.sif mpalmiot@cca.in2p3.fr:/pbs/home/m/mpalmiot/code/
