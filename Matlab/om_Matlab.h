//
//  om_Matlab.h
//  
//
//  Created by Nizar Ouarti on 13/06/17.
//
//

#ifndef ____om_Matlab__
#define ____om_Matlab__

#include <stdio.h>
#include "utils_.h"
#include "openmotion/om.h"
#include "utils_.h"

double *data;
Init_Params I_P;
char *method;

//double *DELTA_T_;

void buildStruct(double DELTA_T,double var_accel, double var_gyro, double var_magn,double *local_magn,double *bias_accel,double *bias_gyro,double *bias_magn,double *init_quat, Init_Params *I_P);

void Compute(double *data,Init_Params *I_P, char *method, int cols,double *quats);

#endif /* defined(____om_Matlab__) */
