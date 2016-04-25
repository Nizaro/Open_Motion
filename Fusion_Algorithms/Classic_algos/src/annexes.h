/*
 * annexes.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#ifndef ANNEXES_H_
#define ANNEXES_H_

#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <typeinfo>

#if defined(_OPENMP)
#include <omp.h>
#endif


#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

using namespace std;

namespace om{

///////////////////////////////////////////////////////
/////            Global constants                 /////
///////////////////////////////////////////////////////


#define PI 3.14159265359
#define DELTA_T 0.001
#define G 9.81
#define EPSILON 0.000000001
#define RAD_TO_DEG 180.0/PI
#define DEG_TO_RAD PI/180

///////////////////////////////////////////////////////
/////          Tools display console              /////
///////////////////////////////////////////////////////

std::string prd(const double x, const int decDigits, const int width);
std::string prd(const double x, const int decDigits);
std::string center(const string s, const int w);


///////////////////////////////////////////////////////
/////                 Others                      /////
///////////////////////////////////////////////////////

double signeOf(double x);




}





#endif /* ANNEXES_H_ */
