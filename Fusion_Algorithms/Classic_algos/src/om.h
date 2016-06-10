/*
 * om.h
 *
 *  Created on: 25 May, 2016
 *      Author: thomas
 */

#ifndef OM_H_
#define OM_H_


#include "geometry.h"


extern omVector ned_magnetic_field;
extern omVector ned_gravity;
extern omVector ned_geographic_north;

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
/////           NonLinearFilter USQUE             /////
///////////////////////////////////////////////////////

typedef struct omNonLinearFilter_USQUE{

	double _seed;
	double _f;
	double _a;
	double _lambda;

	omQuaternion _q_est;
	omQuaternion _q_pred;

	omVector _x_k;
	omVector _x_k_pred;
	omMatrix _P_k;
	omMatrix _P_k_pred;

	omMatrix _Q;
	omMatrix _R;


}omNonLinearFilter_USQUE;

void om_usque_initialization(struct omSensorFusionManager *manager,void *filter);
void om_usque_process(struct omSensorFusionManager *manager,void *filter);
void om_usque_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter);
void om_usque_update(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter,omVector *sigma_points,omQuaternion* sigma_quaternion);


///////////////////////////////////////////////////////
/////           NonLinearFilter REQUEST           /////
///////////////////////////////////////////////////////


typedef struct omNonLinearFilter_REQUEST{

	double _seed;
	double _lambda_m_k;
	double _mu_k;
	double _m_k;

	omQuaternion _q_est;

	omMatrix _P_k;
	omMatrix _P_k_pred;

	omMatrix _K_k;
	omMatrix _K_k_pred;
	omMatrix _d_K_k;

	omMatrix _d_B_k;
	omMatrix _d_S_k;
	omVector _d_z_k;
	double _d_m_k;
	double _d_sigma_k;


	omMatrix _R;
	omMatrix _Q;


	double* _a;
	omVector* _r;
	omVector* _b;


}omNonLinearFilter_REQUEST;


void om_request_initialization(struct omSensorFusionManager *manager,void *filter);
void om_request_process(struct omSensorFusionManager *manager,void *filter);
void om_request_preprocess(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter);
void om_request_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter);
void om_request_update(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter);

void om_request_computeR(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter);
void om_request_computeQ(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter);


#endif /* OM_H_ */






