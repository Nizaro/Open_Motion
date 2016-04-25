#!/bin/bash

#get current directory
CURRENT_DIR=`pwd`
NDK_BUILD_dir=/home/thomas/Android/Sdk/ndk-bundle/

#build project and install it

mkdir build
cd build/
cmake ..

sudo make
sudo make install


