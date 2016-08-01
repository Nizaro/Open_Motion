# OpenMotion : Inertial fusion Library

Synopsis
========

OpenMotion is an open source library written in C providing data fusion and attitude estimation algorithms. OpenMotion was designed for computational efficiency and with a strong focus on real-time applications.
A set of functions are provided in order to help one implementing its own algortithms and to compare its results.
The implemented algorithms use the output of an IMU composed by a 3D gyroscope, a 3D accelerometer and a 3D magnetometer and can be used in real-time applications on embedded systems.


Motivation
==========

This opensource project aims to provide a simple and usefull way for scientists ans users to implement data fusion algorithms and to test them.
An other feature is to provide a tool for benchmarking as far as data fusion algorithms are concerned. Thus, scientists and industrials have now the possibility to achieve a scientificaly reliable performance comparison between diferent algorithms. 
On the contary of must 3D libraries, OpenMotion library has been written in C language in order to enable better portability and easier bindings like with MatLab for instance in order to run simulations.

Installation
============

Source code for OpenMotion is kept on GitHub : https://github.com/Nizaro/Open_Motion

Linux
------

Clone the GitHub OpenMotion repository :

	git clone https://github.com/Nizaro/Open_Motion

Get into the corresponding folder and create a build repository :
	
	cd OpenMotion/Fusion_Algorithms/Classic_algos

Then execute the "build.sh"	script which will make use of cmake in odrer to create the OpenMotion library.

	./build.sh
	
If cmake is not installed or is deprecated, you need to run the following lines before :
	
	sudo apt-get install cmake

Or :
	
	sudo apt-get upgrade
	
In order to use OpenMotion as a linux library use this command after the build is complete
	
	sudo ldconfig -v
	

Windows
-------

First you need to install cmake : https://cmake.org/download/

Then you need to clone GitHub OpenMotion repository: https://github.com/Nizaro/Open_Motion

Create a "build" repository under OpenMotion repository.

On CMake (cmake-gui), the source code is under OpenMotion repository, the binaries are to be built under "build" repository and the generator should be MinGW Makefiles.

Then press "Generate".

Mac OS
------

The library has not been tested in a Mac OS environment yet

API Reference
=============

**Link towards API reference docs live** :

Code Example
============

An sample of usage is available here in the folder ./OpenMotion/Fusion_Algorithms/samples/.
A makefile is provided. But you can compile by yourself with the command :
	
	gcc -o om_sample om_sample.c -std=c99 -lopenmotion -lm -D_BSD_SOURCE
	
And execute with 

	./om_sample -i imu_data.csv	

The file imu_data.csv contains simulated and noisy values of a IMU such as we have access to the ground truth too.
However, you can use the library in a real situation.

Android Example
============

An example of android application is available too. In this application, OpenMotion is used as a native library. 
However, in order to use the code properly, several modifications must be done before. Otherwise the application will not compile and work in a optimal manner.

- A compilation with ndk-build is required of the Library. Go to android-native and launch the command.
- In the file OpenMotionAndroidExample/app/build.gradle, please replace the path to the OpenMotion library at line 19,39 and 41. Also in files Android.mk and Application.mk
- The Library works only if all the component of the IMU (gyroscope, accelerometer and magnetometer) belong to the same reference frame. 
Then, please check if the IMU of your smartphone or tablet verify this constraint. Do not hesite to change the axis if require.
- The Library will transform the magnetic field (output of the magnetometer) into the geographic north. 
This transformation is based on the declinaison and the inclinaison of the magnetic fiel. This values depends on your location.
Go see http://www.ngdc.noaa.gov/geomag-web/#igrfwmm for more information
- The Library need an estimate of the sensors biases and variances, this values provided here are not optimal for every systems

This application has been tested on an Samsung galaxy tab S II with Android 5.0.2 API 21.

Additional functionality
========================

C functions implementing static calibration for the sensors are also provided in order to achieve the best accuracy possible for the given algorithms.

Future features and improvements
================================

An android application linked to the calibration native C functions is also provided in order to calibrate the 3D sensors (Accelerometer, Gyroscope, Magnetometer) .
This application can also be used "as is" in others applications or projects.

Contributors
============

OpenMotion library is an opensource project. 

All people who want to improve, discuss or implement new functionalities are welcomed.

Here are some usefull links fot those who want to study or understand deeper the implemented algorithms :

- 3D geometry : http://sedris.org/wg8home/Documents/WG80485.pdf
- Kalman filters : https://cse.sc.edu/~terejanu/files/tutorialKF.pdf
- Kalman extended filters : https://homes.cs.washington.edu/~todorov/courses/cseP590/readings/tutorialEKF.pdf
- CGO : http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4608934
- USQUE : http://www.acsu.buffalo.edu/~johnc/uf_att.pdf
- Particle filter : http://arc.aiaa.org/doi/pdf/10.2514/1.47236
- Nonlinear Attitude Estimation Filtering : http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20020060647.pdf


License
=======

OpenMotion is released under a BSD license and hence it’s free for both academic and commercial use :

Copyright (c) <2016>, <OUARTI Nizar, BRAUD Thomas, BILLAUD Vivien>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
