#!/bin/bash -x

mkdir -p lib
cd lib
if [ $? -ne 0 ]
then
    echo "lib directory not found"
    exit
fi

cd mimalloc
if [ $? -ne 0 ]
then
    echo "lib/mimalloc directory not found"
    git clone -b v1.2.0 https://github.com/microsoft/mimalloc.git
    if [ $? -ne 0 ]
    then
        echo "Unable to perform git clone https://github.com/microsoft/mimalloc.git"
        exit
    fi
    cd mimalloc
fi

echo "CMake is required for installing mimalloc"
unset LD_PRELOAD
mkdir -p out/release
cd out/release
if [ $? -ne 0 ]
then
    echo "out/release directory not found"
    exit
fi

cmake ../..
if [ $? -ne 0 ]
then
    echo "cmake not found"
    exit
fi
make

echo "################################################################################"
echo "Mimalloc installed."
echo "Ensure that you export LD_PRELOAD=`pwd`/libmimalloc.so"
echo "################################################################################"
