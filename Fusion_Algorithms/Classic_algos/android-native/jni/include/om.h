/*
 * om.h
 *
 *  Created on: 25 May, 2016
 *      Author: thomas
 */

#ifndef OM_H_
#define OM_H_


#include "geometry.h"


static omVector ned_magnetic_field;
static omVector ned_gravity;
static omVector ned_geographic_north;

void init_ned_frame();


///////////////////////////////////////////////////////
/////             Structures                      /////
///////////////////////////////////////////////////////



typedef struct omIMUParams{

	double variance_accelerometer;
	double variance_gyroscope;
	double variance_magnetometer;

	omVector bias_accelerometer;
	omVector bias_gyroscope;
	omVector bias_magnetometer;

}omIMUParams;

typedef struct omIMUData{

	omVector data_accelerometer;
	omVector data_gyroscope;
	omVector data_magnetometer;

}omIMUData;


typedef enum omAttitudeReprensatationType{
	Quarternion,Matrix,AxisAngle,EulerAngle
}omAttitudeReprensatationType;


union omAttitudeReprensatation{

	omQuaternion quaternion;
	omEulerAngle euler;
	omAxisAngle axis_angle;
	omMatrix matrix;
};



typedef struct omSensorFusionManager{

	enum omAttitudeReprensatationType type;
	union omAttitudeReprensatation output;
	omIMUData imu_data;
	omIMUParams imu_params;

	void (*process)(struct omSensorFusionManager *manager,void *filter);
	void (*initialization)(struct omSensorFusionManager *manager,void *filter);

}omSensorFusionManager;




///////////////////////////////////////////////////////
/////           NonLinearFilter CGO               /////
///////////////////////////////////////////////////////


typedef struct omNonLinearFilter_CGO{

	double _seed;
	double _k_mag;
	double _k_acc;
	double _k_I;
	double _k_P;

	omQuaternion _q_est;
	omQuaternion _q_pred;

	omVector _bias_est;
	omVector _bias_pred;


}omNonLinearFilter_CGO;

void om_cgo_initialization(struct omSensorFusionManager *manager,void *filter);
void om_cgo_process(struct omSensorFusionManager *manager,void *filter);


///////////////////////////////////////////////////////
/////           NonLinearFilter MEKF              /////
///////////////////////////////////////////////////////

typedef struct omNonLinearFilter_USQUE{

	double _seed;
	double _f;
	double _a;
	double _lambda;

	omQuaternion _q_est;
	omQuaternion _q_pred;

	omVector _x_k;
	omMatrix _P_k;

	omMatrix _Q;
	omMatrix _R;


}omNonLinearFilter_USQUE;

void om_initialization_usque(struct omSensorFusionManager *manager,void *filter);
void om_process_usque(struct omSensorFusionManager *manager,void *filter);
void om_prediction_usque(struct omSensorFusionManager *manager,void *filter);
void om_update_usque(struct omSensorFusionManager *manager,void *filter);



#endif /* OM_H_ */






