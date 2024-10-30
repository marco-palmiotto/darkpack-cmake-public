sudo yum install gcc -y || exit 1
sudo yum install gfortran -y || exit 1
sudo yum install g++ -y || exit 1
sudo yum install gsl -y || exit 1
sudo yum install gsl-devel -y || exit 1
sudo yum install texlive -y || exit 1
sudo yum install ImageMagick -y || exit 1
sudo yum install dvipng -y || exit 1
sudo yum install coreutils -y || exit 1
sudo yum install qtchooser -y || exit 1
sudo yum install cmake -y || exit 1

# dependencies for darkpack
sudo yum install ninja-build -y || exit 1
# sudo yum install libgomp.x86_64 libomp-devel.x86_64 -y || exit 1

# dependencies for cma
sudo yum install eigen3-devel -y || exit 1

