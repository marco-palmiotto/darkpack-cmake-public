#!/bin/bash

if [ $# -ne 1 ]
then 
  echo "This script requires one argument:"
  echo " 1. The name of the numerical library"
  exit 1
fi

[ ! -d $1 ] && echo "The library does not exist" && exit 2


barebone_lib_destination_path=../$1
if [ ! -d $barebone_lib_destination_path ]
then
    echo Creating the new library
    cp -r $1 $barebone_lib_destination_path
else
    echo A library already exists in /src, skipping
fi


darkpack_lib_destination_path=../dp_$1
if [ ! -d $darkpack_lib_destination_path ]
then 
    echo Creting the new DarkPACK library
    mkdir $darkpack_lib_destination_path
fi

cp -r auxiliary_library/$1/ $darkpack_lib_destination_path
cp auxiliary_library/blank_files/leshouchesfrommarty.hpp $darkpack_lib_destination_path/leshouchesfrommarty.hpp

