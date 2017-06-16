/*
 * calibration.h
 *
 *  Created on: 30 May 2017
 *      Author: thomas
 */

#ifndef SRC_CALIBRATION_H_
#define SRC_CALIBRATION_H_


#include "algebra.h"

/////////////////////////////////////////////
///////   Adaptive least square         /////
/////////////////////////////////////////////

void mat_to_vecs(struct omMatrix* A,struct omVector* vecs_A);

double tensor_T(int k,int i,int l,omMatrix* X,double var);

int tensor_R(int p,int q,int i,omMatrix* M);

double function_n_als(int p,int q,omMatrix* M,omMatrix* X);

void om_calibration_adjusted_least_square(omMatrix* X,omMatrix* Q,omVector* b,double* d,double Hm_square);




/////////////////////////////////////////////
///////   Variance estimation           /////
/////////////////////////////////////////////


void om_calibration_noise_estimation(omMatrix* X,omVector* v);



#endif /* SRC_CALIBRATION_H_ */
