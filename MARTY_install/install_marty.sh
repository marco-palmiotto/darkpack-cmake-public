#!/bin/bash

nthreads=1

[ $# -eq 1 ] && nthreads=$1

work_dir=$(pwd)
sourcepath=$(pwd)/marty-public
installpath=$sourcepath/build/.install 

echo Checking gcc version

# Get the version of GCC
gcc_version=$(gcc --version | grep -oP '\d+\.\d+\.\d+' | head -n1 | cut -f 1 -d '.')

echo You have installed $gcc_version

if [[ $gcc_version -le 9 || $gcc_version -ge 12 ]]
then
  echo "The required version must be between gcc-9 and gcc-12"
  exit 1
fi

echo Cleaning possible past installations in marty-public

[ -d $sourcepath ] && rm -rf $sourcepath

download_marty()
{
  git clone https://github.com/docbrown1955/marty-public.git 
  cd $sourcepath
  if [ $? -ne 0 ] ; then
    echo "Impossible to change directory in $sourcepath, aborting" ; exit 5
  fi
  git pull
  git checkout plugin/boltzmann_solver

  if [ $? -ne 0 ] ; then
    echo "Error in downloading MARTY from the repository" ; cd .. ; exit 1
  fi
}

cd $work_dir
download_marty 

mkdir -p $installpath/bin 
if [ $? -ne 0 ] ; then
  echo "Error in creating $installpath/bin" ; exit 2
fi

mkdir -p $installpath/include 
if [ $? -ne 0 ] ; then
  echo "Error in creating $installpath/include" ; exit 3
fi
mkdir -p $installpath/lib 
if [ $? -ne 0 ] ; then
  echo "Error in creating $installpath/lib" ; exit 4
fi

cd $sourcepath

cmake -DCMAKE_INSTALL_PREFIX=.install -S . -B build
if [ $? -ne 0 ] ; then
 echo "Error in calling cmake (line 16), aborting" ; exit 8
fi

echo "Starting compilation of MARTY with $nthreads threads"

cd build && make -j$nthreads 
if [ $? -ne 0 ] ; then
 echo "Error in calling make (line 17), aborting" ; exit 9
fi

make install 
if [ $? -ne 0 ] ; then
 echo "Error in calling make install (line 18), aborting" ; exit 10
fi

echo 'Installation ended succesfully'
echo \n
echo 'Generating configuration file'

# If $INSTALLMARTYPATH is not defined, we define it in the ~/.bashrc
echo "# MARTY envirnmental variables" > $work_dir/.martyconfig
echo "export INSTALLMARTYPATH=$installpath # path where MARTY is installed">> $work_dir/.martyconfig
echo 'export PATH=$INSTALLMARTYPATH/bin:$PATH # Find the command grafed'>> $work_dir/.martyconfig
echo 'export LD_LIBRARY_PATH=$INSTALLMARTYPATH/lib:$LD_LIBRARY_PATH # Find libraries at run-time (Linux)'>>  $work_dir/.martyconfig
echo 'export LIBRARY_PATH=$INSTALLMARTYPATH/lib:$LIBRARY_PATH # Find libraries at compile-time' >> $work_dir/.martyconfig
echo 'export C_INCLUDE_PATH=$INSTALLMARTYPATH/include:$C_INCLUDE_PATH # Find headers at compile time' >> $work_dir/.martyconfig
echo 'export CPATH=$INSTALLMARTYPATH/include:$CPATH # Find headers at compile time' >> $work_dir/.martyconfig

echo \n

echo "IMPORTANT: before proceeding, copy the content of $work_dir/.martyconfig in your .bashrc file for the enviromental variables!"

echo \n

if [ -f compile_with_marty.sh ]
then 
  marty_script_path="$(pwd)/compile_with_marty.sh"
elif [ -f "$work_dir/compile_with_marty.sh" ]
then
  marty_script_path="$work_dir/compile_with_marty.sh"
else
  marty_script_path=""
fi

[ $marty_script_path -ne "" ] && sed "s|export INSTALLMARTYPATH=*|export INSTALLMARTYPATH=$installpath|g" -i compile_with_marty.sh && echo "Alternatively, you can try compiling your MARTY .cpp sources by using the script $marty_script_path, adding in your .bashrc only the path of the shared libraries, needed at runtime"


