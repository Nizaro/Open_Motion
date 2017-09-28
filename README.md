# OpenMotion : Inertial Fusion Library

Synopsis
========

OpenMotion is an open source library written in C providing data fusion and attitude estimation algorithms. OpenMotion was designed for computational efficiency and with a strong focus on real-time applications.
A set of functions are provided in order to facilitate the implementation of original algorithms and to compare the results with state-of-the-art algorithms already implemented.
The implemented algorithms combines the output of an IMU composed by a 3D gyroscope, a 3D accelerometer and a 3D magnetometer and can be used in real-time applications for embedded systems.


Motivation
==========

This open source project aims to provide a simple and useful way for scientists and users to test and implement data fusion algorithms. It can be also useful for non-experts to learn how to implement these kinds of algorithms.
This software can also be a benchmarking tool for attitude estimation. Thus, scientists and industrials have now the possibility to achieve a scientifically reliable performance comparison between different algorithms. 
OpenMotion library has been written in C language in order to enable better portability and to ease the bindings with MatLab or Android.

Installation
============

Source code for OpenMotion is kept on GitHub : https://github.com/Nizaro/Open_Motion

Linux
------

Clone the GitHub OpenMotion repository :

	git clone https://github.com/Nizaro/Open_Motion

Get into the corresponding folder and create a build repository :
	
	cd Open_Motion/Fusion_Algorithms/Classic_algos
    mkdir build

Then execute the "build.sh"	script which will make use of cmake in order to create the OpenMotion library.

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

First you need to install homebrew :

	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Then use homebrew to install cmake:

	brew install cmake
	
Once cmake is installed, go into the corresponding folder and create a build repository :

	cd Open_Motion/Fusion_Algorithms/Classic_algos
    mkdir build

Then execute the "build.sh"	script which will make use of cmake in order to create the OpenMotion library.
    
    cd Open_Motion/Fusion_Algorithms/C_code
	./build.sh


API Reference
=============

**Link towards API reference docs live** :

Code Example
============

An sample of usage is available here in the folder ./Open_Motion/Fusion_Algorithms/samples/.
A makefile is provided. But you can compile by yourself with the command :
	
	gcc -o om_sample om_sample.c -std=c99 -lopenmotion -lm -D_BSD_SOURCE
	
And execute with 

	./om_sample -i imu_data.csv	

The file imu_data.csv contains noisy simulated values of an IMU. Because the data are here generated, the ground truth is accessible.
However, you can use the library with real data by making your fork of om_sample.c.

Matlab Example
============
A Matlab version is under development.

Android Example
============

An example of android application is available too. In this application, OpenMotion is used as a native library. 
However, in order to use the code properly, several modifications must be done before. Otherwise the application will not compile and work in a optimal manner.

- A compilation with ndk-build is required of the Library. Go to android-native and launch the command.
- In the file OpenMotionAndroidExample/app/build.gradle, please replace the path to the OpenMotion library at line 19,39 and 41. Also in files Android.mk and Application.mk
- The Library works only if all the component of the IMU (gyroscope, accelerometer and magnetometer) belong to the same reference frame. 
Then, please check if the IMU of your smartphone or tablet verify this constraint. Do not hesitate to change the axis frame if required.
- The Library will transform the magnetic field (output of the magnetometer) into the geographic north. 
This transformation is based on the declination and the inclinaison of the magnetic fiel. This values depends on your location.
Go see http://www.ngdc.noaa.gov/geomag-web/#igrfwmm for more information
- The Library need an estimate of the sensors biases and variances, this values provided here are not optimal for every systems

This application has been tested on an Samsung galaxy tab S II with Android 5.0.2 API 21.

Comparative Framework
============

This tool aims to compare the performance of the algorithm on several aspect. Not only on the accuracy, but also on the computation time requirement.

Preliminary
------

In order to use the framework, few installation is required. The framework had been developed under Linux environment. Thus, the following description concern only Linux user. Windows and Mac OS adaptation will come in future version.
Here the required package:

- csvkit 1.0.2(or newer) http://csvkit.readthedocs.io/en/1.0.2/index.html
- R package 'limma' https://bioconductor.org/packages/release/bioc/html/limma.html

How to use it
------

Once everything is installed the framework is quite easy to use. All script files are in the folder /Comparative_framework/scripts/

To launch the simulation of a IMU, use this commande
	
	./generate_data -n number_of_simulation

This will generate csv file in the folder /Comparative_framework/data/csv/ and the associated png showing the data in /Comparative_framework/data/png/

Then, to generate all the performance results, use this commande

	./generate_results 

All the results will be located in /Comparative_framework/results/current-date_current-hour/ 
This folder contains other folder:

- perfo/hist/ : the histogram showing the different score 
- perfo/table/  : the table showing the different score
- plots_angle : contains all plot showing the true and the estimate angle for all methods
- plots_error : contains all plot showing the error for all methods

Do not take into account the other folder

Additional functionality
========================

C functions implementing static calibration for the sensors are also provided in order to achieve the best accuracy possible for the given algorithms.

Future features and improvements
================================

An android application linked to the calibration native C functions is also provided in order to calibrate the 3D sensors (Accelerometer, Gyroscope, Magnetometer) .
This application can also be used "as is" in others applications or projects.

Contributors
============

OpenMotion library is an open source project. 

All people who want to improve, discuss or implement new functionalities are welcomed.

Here are some useful links for those who want to study or understand deeper the implemented algorithms :

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

Copyright (c) <2017>, <OUARTI Nizar, BRAUD Thomas>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




