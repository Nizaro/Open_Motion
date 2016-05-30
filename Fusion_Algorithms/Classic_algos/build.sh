#!/bin/bash

#get current directory
CURRENT_DIR=`pwd`

#build project and install it
mkdir build
cd build/
cmake .. 
sudo make
sudo make install

cd $CURRENT_DIR

#create android-native libraries

cp $CURRENT_DIR/src/*.c $CURRENT_DIR/android-native/jni
cp $CURRENT_DIR/src/*.h $CURRENT_DIR/android-native/jni/include


