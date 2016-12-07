/**
 * \file om.h
 * \author Thomas BRAUD, Nizar OUARTI
 * \date 10 june 2016
 * \brief File containing linear algebra methods
 *
 */


#ifndef OM_H_
#define OM_H_


#include "geometry.h"


extern omVector ned_magnetic_field;
extern omVector ned_gravity;




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


typedef enum omAttitudeRepresentationType{

	Quarternion,Matrix,AxisAngle,EulerAngle

}omAttitudeRepresentationType;


union omAttitudeRepresentation{

	omQuaternion quaternion;
	omEulerAngle euler;
	omAxisAngle axis_angle;
	omMatrix matrix;

};



typedef struct omSensorFusionManager{

	enum omAttitudeRepresentationType type;
	union omAttitudeRepresentation output;
	omIMUData imu_data;
	omIMUParams imu_params;

	void (*process_filter)(struct omSensorFusionManager *manager,void *filter);
	void (*initialization_filter)(struct omSensorFusionManager *manager,void *filter);
	void (*free_filter)(void *filter);

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
void om_cgo_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CGO *filter);
void om_cgo_update(struct omSensorFusionManager *manager,omNonLinearFilter_CGO *filter);

void om_cgo_free(void *filter);

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

	omVector *sigma_points;
	omQuaternion* sigma_quaternion;

}omNonLinearFilter_USQUE;

void om_usque_initialization(struct omSensorFusionManager *manager,void *filter);
void om_usque_process(struct omSensorFusionManager *manager,void *filter);
void om_usque_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter);
void om_usque_update(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter);
void om_usque_free(void *filter);


///////////////////////////////////////////////////////
/////           NonLinearFilter USQUE             /////
///////////////////////////////////////////////////////

typedef struct omNonLinearFilter_CSP{

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

	double _k_mag;
	double _k_acc;
	double _k_I;
	double _k_P;

	omVector *sigma_points;
	omQuaternion* sigma_quaternion;

}omNonLinearFilter_CSP;

void om_csp_initialization(struct omSensorFusionManager *manager,void *filter);
void om_csp_process(struct omSensorFusionManager *manager,void *filter);
void om_csp_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CSP *filter);
void om_csp_update(struct omSensorFusionManager *manager,omNonLinearFilter_CSP *filter);
void om_csp_free(void *filter);



///////////////////////////////////////////////////////
/////           NonLinearFilter MEKF              /////
///////////////////////////////////////////////////////


typedef struct omNonLinearFilter_MEKF{

	double _seed;
	double _f;
	double _a;
	double _h;
	double _lambda;

	omQuaternion _q_est;
	omQuaternion _q_pred;

	omVector _x_k;
	omVector _x_k_pred;
	omMatrix _P_k;
	omMatrix _P_k_pred;

	omVector _w_k;
	omVector _v_k;

	omMatrix _F;
	omMatrix _H;

	omMatrix _Q;
	omMatrix _R;
	omMatrix _Q_cho;
	omMatrix _R_cho;

}omNonLinearFilter_MEKF;



void om_mekf_initialization(struct omSensorFusionManager *manager,void *filter);
void om_mekf_process(struct omSensorFusionManager *manager,void *filter);
void om_mekf_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter);
void om_mekf_update(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter);
void om_mekf_f_function(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter,omVector *x,omVector *f_x);
void om_mekf_h_function(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter,omVector *x,omVector *h_x);
void om_mekf_f_jacobian(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter);
void om_mekf_f_jacobian(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter);
void om_mekf_free(void *filter);

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

void om_request_free(void *filter);

///////////////////////////////////////////////////////
/////           NonLinearFilter PF                /////
///////////////////////////////////////////////////////

typedef struct omNonLinearFilter_PF{

	double _seed;
	omQuaternion _q_est;
	omQuaternion _q_pred;

	double _sum_w;
	double _sum_n_eff;

	double _n;
	double _threeshold;
	double _f;
	double _h;

	int _resample;

	omMatrix _R;
	omMatrix _L;

	omVector _x_k;


	double* _particle_w;
	omVector* _particle_wn;
	omVector* _particle_x;
	omQuaternion* _particle_q;


}omNonLinearFilter_PF;

void om_pf_initialization(struct omSensorFusionManager *manager,void *filter);
void om_pf_process(struct omSensorFusionManager *manager,void *filter);
void om_pf_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_PF *filter);
void om_pf_update(struct omSensorFusionManager *manager,omNonLinearFilter_PF *filter);
void om_pf_resampling(omNonLinearFilter_PF *filter);
void om_pf_swap(omNonLinearFilter_PF *filter,int i,int j);
void om_pf_quicksort(omNonLinearFilter_PF *filter,int left, int right);

void om_pf_free(void *filter);


///////////////////////////////////////////////////////
/////           NonLinearFilter GDOF              /////
///////////////////////////////////////////////////////

typedef struct omNonLinearFilter_GDOF{

	double _eta;
	double _beta;

	omQuaternion _q_est;
	omVector _e_m;
	omVector _e_b;

	omVector _f_a;
	omVector _f_b;

	omMatrix _J_a;
	omMatrix _J_b;

	omVector _bias_est;



}omNonLinearFilter_GDOF;

void om_gdof_initialization(struct omSensorFusionManager *manager,void *filter);
void om_gdof_process(struct omSensorFusionManager *manager,void *filter);
void om_gdof_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter);
void om_gdof_update(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter);

void om_gdof_free(void *filter);


///////////////////////////////////////////////////////
/////           NonLinearFilter CFA               /////
///////////////////////////////////////////////////////


typedef struct omNonLinearFilter_CFA{

	double _lambda;
	double _beta;

	omQuaternion _q_est;
	omQuaternion _q_pred;

	omVector _v_acc_pred;
	omVector _v_mag_pred;

}omNonLinearFilter_CFA;


void om_cfa_initialization(struct omSensorFusionManager *manager,void *filter);
void om_cfa_process(struct omSensorFusionManager *manager,void *filter);
void om_cfa_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CFA *filter);
void om_cfa_update(struct omSensorFusionManager *manager,omNonLinearFilter_CFA *filter);

void om_cfa_free(void *filter);



///////////////////////////////////////////////////////
/////           NonLinearFilter QUEST             /////
///////////////////////////////////////////////////////


typedef struct omNonLinearFilter_QUEST{

	omQuaternion _q_est;

	double* _a;
	omVector* _r;
	omVector* _b;

}omNonLinearFilter_QUEST;


void om_quest_initialization(struct omSensorFusionManager *manager,void *filter);
void om_quest_process(struct omSensorFusionManager *manager,void *filter);
void om_quest_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_QUEST *filter);
void om_quest_update(struct omSensorFusionManager *manager,omNonLinearFilter_QUEST *filter);

void om_quest_free(void *filter);



#endif /* OM_H_ */






