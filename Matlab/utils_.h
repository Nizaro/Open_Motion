//
//  utils_.h
//  
//
//  Created by Nizar Ouarti on 12/06/17.
//
//

#ifndef ____utils___
#define ____utils___

#include <stdio.h>
#include "openmotion/om.h"
#include <string.h>

/*
 * File reader
 */
struct line_reader {
    
    // All members are private.
    FILE	*f;
    char	*buf;
    size_t	 siz;
    
};

// This to move to om.h

/*
 * The different data fusion algorithm present in the library
 */
typedef enum MethodType{
    MEKF,REQUEST,QUEST,USQUE,CGO,PF,GDOF,CFA,CSP
}MethodType;



/*
 * The ground truth composed by the true position in a 3D scene and the true orientation represented by a quaternion
 */
struct Groundtruth{
    
    omVector position;
    omQuaternion q_true;
    
};


/*Initialisation structure : different parameters that allow to initialise the the fusion manager*/
typedef struct Init_Params_{
    
    double DELTA_T;
    
    double var_accel;
    double var_magn;
    double var_gyro;
    
    double local_magn[3];
    
    double bias_accel[3];
    double bias_magn[3];
    double bias_gyro[3];

    double init_quat[4];

} Init_Params;




void lr_init(struct line_reader *lr, FILE *f);
char* next_line(struct line_reader *lr, size_t *len);
void lr_free(struct line_reader *lr);
char** str_split(char* a_str, const char a_delim);
void displayLoadingBar(int index);
double calculErrorOrientation(omQuaternion *q_real,omQuaternion *q_est);
void init_ned(double magn[3]);
void init_manager(omSensorFusionManager *manager,void* filter, MethodType type,Init_Params *I_P);


#endif /* defined(____utils___) */
