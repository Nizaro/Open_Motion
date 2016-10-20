/*
 * om.c
 *
 *  Created on: 25 May, 2016
 *      Author: thomas
 */

#include "om.h"

// global variables
omVector ned_gravity;
omVector ned_geographic_north;


/* initialization of North East Down frame vector */
void init_ned_frame(){

	om_vector_create(&ned_geographic_north,3,1.0,0.0,0.0);
	om_vector_create(&ned_gravity,3,0.0,0.0,1.0);

}


///////////////////////////////////////////////////////
/////           NonLinearFilter CGO               /////
///////////////////////////////////////////////////////

/* initialization of all component used by the nonlinear filter CGO */
void om_cgo_initialization(struct omSensorFusionManager *manager,void *filter){

	// constant initialization
	(*(omNonLinearFilter_CGO*)filter)._seed = 0;
	(*(omNonLinearFilter_CGO*)filter)._k_I = 0.3;
	(*(omNonLinearFilter_CGO*)filter)._k_P = 1.0;
	(*(omNonLinearFilter_CGO*)filter)._k_acc = 10.0;
	(*(omNonLinearFilter_CGO*)filter)._k_mag = 10.0;

	// attitude initialization
	om_quat_create(&(*(omNonLinearFilter_CGO*)filter)._q_est,1.0,0.0,0.0,0.0);
	switch(manager->type){

	// convertion according to user choices
	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_CGO*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_CGO*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_CGO*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_CGO*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	// sensors' biases representation
	om_vector_create(&(*(omNonLinearFilter_CGO*)filter)._bias_est,3);
	om_vector_create(&(*(omNonLinearFilter_CGO*)filter)._bias_pred,3);

	om_vector_clone(&manager->imu_params.bias_gyroscope,&(*(omNonLinearFilter_CGO*)filter)._bias_est);
	om_vector_clone(&manager->imu_params.bias_gyroscope,&(*(omNonLinearFilter_CGO*)filter)._bias_pred);

	// initialization of North East Down frame vector
	init_ned_frame();

}

/* process function of the nonlinear filter CGO */
void om_cgo_process(struct omSensorFusionManager *manager,void *filter){

	om_cgo_prediction(manager,(omNonLinearFilter_CGO*)filter);
	om_cgo_update(manager,(omNonLinearFilter_CGO*)filter);


}



void om_cgo_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CGO *filter){

	//////
	// Prediction step
	/////

	// compute b_(k+1)^- = b_(k)^+
	for(int i=0;i<3;++i)
		om_vector_setValue(&filter->_bias_pred,i,filter->_bias_est._values[i]);

	// compute angular_velocity = y_gyro - b_pred
	omVector angular_velocity;
	om_vector_create(&angular_velocity,3);
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&filter->_bias_pred,&angular_velocity);

	// compute q_(k+1)^- = Omega(Y_gyro - b_pred)q_(k)^+
	om_kinematics_quaternion(&filter->_q_est,&angular_velocity,&filter->_q_pred);

	// free memory
	om_vector_free(&angular_velocity);

}


void om_cgo_update(struct omSensorFusionManager *manager,omNonLinearFilter_CGO *filter){


	//////
	// Update step
	/////

	// variable
	omVector z_acc;
	omVector z_mag;
	omVector v_acc_pred;
	omVector v_mag_pred;
	omVector cross_acc;
	omVector cross_mag;
	omVector omega;
	omVector omega_q;


	// allocation
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);
	om_vector_create(&v_acc_pred,3);
	om_vector_create(&v_mag_pred,3);
	om_vector_create(&omega,3);
	om_vector_create(&omega_q,3);
	om_vector_create(&cross_acc,3);
	om_vector_create(&cross_mag,3);


	// compute v_acc = R(q_(k+1)^-)g^a
	om_rotate_vector_quaternion(&filter->_q_pred,&ned_gravity,&v_acc_pred);

	// compute v_mag = R(q_(k+1)^-)m^a
	om_rotate_vector_quaternion(&filter->_q_pred,&ned_geographic_north,&v_mag_pred);

	// removing
	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	// normalization
	om_vector_normalize(&z_acc);
	om_vector_normalize(&z_mag);
	om_vector_normalize(&v_acc_pred);
	om_vector_normalize(&v_mag_pred);

	// compute cross product
	om_vector_crossProduct(&z_acc,&v_acc_pred,&cross_acc);
	om_vector_crossProduct(&z_mag,&v_mag_pred,&cross_mag);

	om_operator_vector_scal_mul(&cross_acc,filter->_k_acc/2.0,&cross_acc);
	om_operator_vector_scal_mul(&cross_mag,filter->_k_mag/2.0,&cross_mag);

	// compute gain omega
	om_operator_vector_add(&cross_acc,&cross_mag,&omega);

	// compute bias_est =  (_omega*(_k_I)*_delta_t*(-1.0));
	om_operator_vector_scal_mul(&omega,filter->_k_I*DELTA_T*(-1.0),&filter->_bias_est);

	// compute angular_velocity = y_gyro - b_est + omega*k_P
	omVector angular_velocity;
	om_vector_create(&angular_velocity,3);
	om_operator_vector_scal_mul(&omega,filter->_k_P,&omega_q);
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&filter->_bias_est,&angular_velocity);
	om_operator_vector_add(&angular_velocity,&omega_q,&angular_velocity);

	// compute q_(k+1)^+ = Omega(Y_gyro - b_pred + omega*k_P)q_(k)^+
	om_kinematics_quaternion(&filter->_q_est,&angular_velocity,&filter->_q_est);

	// set output
	double qw = filter->_q_est._qw;
	double qx = filter->_q_est._qx;
	double qy = filter->_q_est._qy;
	double qz = filter->_q_est._qz;

	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,qw,qx,qy,qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&filter->_q_est,1.0,0.0,0.0,0.0);
		break;

	}

	//free memory
	om_vector_free(&z_acc);
	om_vector_free(&z_mag);
	om_vector_free(&cross_acc);
	om_vector_free(&cross_mag);
	om_vector_free(&omega);
	om_vector_free(&omega_q);
	om_vector_free(&v_mag_pred);
	om_vector_free(&v_acc_pred);



}



/* release all component used for the nonlinear filter CGO */
void om_cgo_free(void *filter){

	//free memory
	om_vector_free(&(*(omNonLinearFilter_CGO*)filter)._bias_est);
	om_vector_free(&(*(omNonLinearFilter_CGO*)filter)._bias_pred);

}


///////////////////////////////////////////////////////
/////           NonLinearFilter USQUE             /////
///////////////////////////////////////////////////////

/* initialization of all component used by the nonlinear filter USQUE */
void om_usque_initialization(struct omSensorFusionManager *manager,void *filter){

	// constant
	(*(omNonLinearFilter_USQUE*)filter)._seed = 0;
	(*(omNonLinearFilter_USQUE*)filter)._a = 1.0;
	(*(omNonLinearFilter_USQUE*)filter)._f = 4.0;
	(*(omNonLinearFilter_USQUE*)filter)._lambda = 1.0;

	om_quat_create(&(*(omNonLinearFilter_USQUE*)filter)._q_pred,1.0,0.0,0.0,0.0);

	switch(manager->type){

	case Quarternion:

		om_quat_create(&(*(omNonLinearFilter_USQUE*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_USQUE*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_USQUE*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_USQUE*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}


	double var_u = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double var_v = manager->imu_params.variance_gyroscope;

	om_matrix_create(&(*(omNonLinearFilter_USQUE*)filter)._Q,6,6);
	for(int i=0;i<3;++i){

		om_matrix_setValue(&(*(omNonLinearFilter_USQUE*)filter)._Q,i,i,(var_v*DELTA_T + 0.33333*var_u*(pow(DELTA_T,3.0))));
		om_matrix_setValue(&(*(omNonLinearFilter_USQUE*)filter)._Q,i+3,i,-(0.5*var_u*DELTA_T*DELTA_T));
		om_matrix_setValue(&(*(omNonLinearFilter_USQUE*)filter)._Q,i,i+3,-(0.5*var_u*DELTA_T*DELTA_T));
		om_matrix_setValue(&(*(omNonLinearFilter_USQUE*)filter)._Q,i+3,i+3,(var_u*DELTA_T));

	}

	om_matrix_create(&(*(omNonLinearFilter_USQUE*)filter)._R,6,6);
	for(int i=0;i<3;++i){
		om_matrix_setValue(&(*(omNonLinearFilter_USQUE*)filter)._R,i,i,manager->imu_params.variance_accelerometer);
		om_matrix_setValue(&(*(omNonLinearFilter_USQUE*)filter)._R,i+3,i+3,manager->imu_params.variance_magnetometer);
	}


	double b_x = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double b_y = om_vector_getValue(&manager->imu_params.bias_gyroscope,1);
	double b_z = om_vector_getValue(&manager->imu_params.bias_gyroscope,2);

	om_matrix_create(&(*(omNonLinearFilter_USQUE*)filter)._P_k,6,6);
	om_matrix_create(&(*(omNonLinearFilter_USQUE*)filter)._P_k_pred,6,6);
	om_vector_create(&(*(omNonLinearFilter_USQUE*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);
	om_vector_create(&(*(omNonLinearFilter_USQUE*)filter)._x_k_pred,6);

	init_ned_frame();

}


void om_usque_free(void *filter){


	om_vector_free(&(*(omNonLinearFilter_USQUE*)filter)._x_k);
	om_vector_free(&(*(omNonLinearFilter_USQUE*)filter)._x_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_USQUE*)filter)._P_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_USQUE*)filter)._P_k);
	om_matrix_free(&(*(omNonLinearFilter_USQUE*)filter)._Q);
	om_matrix_free(&(*(omNonLinearFilter_USQUE*)filter)._R);

}



void om_usque_process(struct omSensorFusionManager *manager,void *filter){

	om_usque_prediction(manager,(omNonLinearFilter_USQUE*)filter);

}


void om_usque_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter){

	double a = filter->_a;
	double f = filter->_f;
	double lambda = filter->_lambda;

	omMatrix S;
	omMatrix T;

	om_matrix_create(&S,6,6);
	om_matrix_create(&T,6,6);

	// compute S = sqrt( (P + Q)*(n + lambda )
	om_operator_matrix_add(&filter->_P_k,&filter->_Q,&T);
	om_operator_matrix_scal_mul(&T,lambda + 6.0,&T);
	om_matrix_choleskyDecomposition(&T,&S);

	// generate sigma point at time k
	omVector* sigma_points = (omVector*)malloc(13 * sizeof(omVector));
	om_vector_create(&sigma_points[0],6);
	om_vector_clone(&filter->_x_k,&sigma_points[0]);

	for(int i=1;i<=6;++i){

		omVector S_row;
		om_vector_create(&S_row,6);

		om_vector_create(&sigma_points[i],6);
		om_vector_create(&sigma_points[i+6],6);

		om_matrix_getRow(&S,(i-1),&S_row);
		om_operator_vector_add(&filter->_x_k,&S_row,&sigma_points[i]);
		om_operator_vector_sub(&filter->_x_k,&S_row,&sigma_points[i+6]);

		om_vector_free(&S_row);

	}



	// generate sigma quaternion at time k
	omQuaternion* sigma_quaternion = (omQuaternion*)malloc(13 * sizeof(omQuaternion));
	for(unsigned int i=0;i<13;++i){

		omVector dp_i;
		om_vector_create(&dp_i,3,om_vector_getValue(&sigma_points[i],0),om_vector_getValue(&sigma_points[i],1),om_vector_getValue(&sigma_points[i],2));

		if(i==0){

			double qw = filter->_q_est._qw;
			double qx = filter->_q_est._qx;
			double qy = filter->_q_est._qy;
			double qz = filter->_q_est._qz;
			om_quat_create(&sigma_quaternion[i],qw,qx,qy,qz);

		}else{

			// compute error quaternion i
			double dp_i_norm = om_vector_norm(&dp_i);
			double dp_i_norm_sq = dp_i_norm*dp_i_norm;
			double dqw_i = ( ( a*(-1.0)*dp_i_norm_sq ) + (f*sqrt( f*f + (1.0 - a*a)*dp_i_norm_sq  ))) / ( (f*f) + dp_i_norm_sq );

			omVector dqv_i;
			om_vector_create(&dqv_i,3);
			om_operator_vector_scal_mul(&dp_i,( (1.0/f) * (a + dqw_i )),&dqv_i);
			omQuaternion dq_i;
			om_quat_create(&dq_i,dqw_i,om_vector_getValue(&dqv_i,0),om_vector_getValue(&dqv_i,1),om_vector_getValue(&dqv_i,2) );

			// compute sigma quaternion i
			om_operator_quat_mul(&dq_i,&filter->_q_est,&sigma_quaternion[i]);
			om_vector_free(&dqv_i);

		}

		om_vector_free(&dp_i);


	}

	// propagation of sigma quaternion and sigma points
	omQuaternion s_q_0_inv;
	for(int i=0;i<13;++i){

		omVector b_i;
		omVector angular_velocity_i;

		om_vector_create(&b_i,3,om_vector_getValue(&sigma_points[i],3),om_vector_getValue(&sigma_points[i],4),om_vector_getValue(&sigma_points[i],5));
		om_vector_create(&angular_velocity_i,3);

		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);

		// propagation of sigma quaternion
		omQuaternion q_pred_i;
		om_kinematics_quaternion(&sigma_quaternion[i],&angular_velocity_i,&q_pred_i);

		if(i==0){

			om_vector_setValues(&sigma_points[0],6,0.0,0.0,0.0,om_vector_getValue(&b_i,0),om_vector_getValue(&b_i,1),om_vector_getValue(&b_i,2));
			om_quat_inverse(&q_pred_i,&s_q_0_inv);

		}else{

			omQuaternion dq_i;

			om_operator_quat_mul(&q_pred_i,&s_q_0_inv,&dq_i);

			double tmp = (f/ (a + dq_i._qw));
			double dp_i_x = dq_i._qx * tmp;
			double dp_i_y = dq_i._qy * tmp;
			double dp_i_z = dq_i._qz * tmp;

			om_vector_setValues(&sigma_points[i],6,dp_i_x,dp_i_y,dp_i_z,
					om_vector_getValue(&b_i,0),om_vector_getValue(&b_i,1),om_vector_getValue(&b_i,2));

		}


		om_quat_create(&sigma_quaternion[i],q_pred_i._qw,q_pred_i._qx,q_pred_i._qy,q_pred_i._qz);

		om_vector_free(&b_i);
		om_vector_free(&angular_velocity_i);
	}

	// compute mean of x_k
	omVector sum_x;
	omVector sum_a;
	omVector sum_b;

	om_vector_create(&sum_x,6,0.0,0.0,0.0,0.0,0.0,0.0);
	om_vector_create(&sum_a,6);
	om_vector_create(&sum_b,6);

	for(int i=1;i<13;++i){
		om_operator_vector_add(&sum_x,&sigma_points[i],&sum_x);
	}

	om_operator_vector_scal_mul(&sigma_points[0],lambda/(6.0+lambda),&sum_a);
	om_operator_vector_scal_mul(&sum_x,1.0/(2.0*(6.0+lambda)),&sum_b);
	om_operator_vector_add(&sum_a,&sum_b,&filter->_x_k_pred);

	// compute covariance P_k
	omMatrix sum_P;

	om_matrix_create(&sum_P,6,6);

	for(unsigned int i=1;i<13;++i){

		omVector var_x_i;
		omMatrix tmp_i;
		omMatrix tmp_i_t;
		omMatrix P_i;

		//allocate
		om_matrix_create(&P_i,6,6);
		om_vector_create(&var_x_i,6);
		om_matrix_create(&tmp_i,6,1);
		om_matrix_create(&tmp_i_t,1,6);

		//compute  X_k(i) - x_k^-
		om_operator_vector_sub(&sigma_points[i],&filter->_x_k_pred,&var_x_i);
		for(int i=0;i<6;i++)
			om_matrix_setValue(&tmp_i,i,0,om_vector_getValue(&var_x_i,i));

		//compute P_k^- +=  (x_k^- - X_k(i))(x_k^- - X_k(i))^T
		om_matrix_transpose(&tmp_i,&tmp_i_t);
		om_operator_matrix_mul(&tmp_i,&tmp_i_t,&P_i);
		om_operator_matrix_add(&sum_P,&P_i,&sum_P);

		//free memory
		om_vector_free(&var_x_i);
		om_matrix_free(&tmp_i);
		om_matrix_free(&tmp_i_t);
		om_matrix_free(&P_i);

	}


	omVector var_x_0;
	omMatrix tmp_0;
	omMatrix tmp_0_t;
	omMatrix P_0;

	om_vector_create(&var_x_0,6);
	om_matrix_create(&tmp_0,6,1);
	om_matrix_create(&tmp_0_t,1,6);
	om_matrix_create(&P_0,6,6);

	//compute  X_k(0) - x_k^-
	om_operator_vector_sub(&sigma_points[0],&filter->_x_k_pred,&var_x_0);

	for(int i=0;i<6;i++)
		om_matrix_setValue(&tmp_0,i,0,om_vector_getValue(&var_x_0,i));

	om_matrix_transpose(&tmp_0,&tmp_0_t);
	om_operator_matrix_mul(&tmp_0,&tmp_0_t,&P_0);

	omMatrix sum_P_a;
	omMatrix sum_P_b;
	om_matrix_create(&sum_P_a,6,6);
	om_matrix_create(&sum_P_b,6,6);

	//compute P_k^- = (lambda/(6.0+lambda))* (X_k(0) - x_k^-)*(X_k(0) - x_k^-)^T + (1.0/(2.0*(6.0+lambda)))*sum_{i=0}^n( (X_k(i) - x_k^-)(X_k(i) - x_k^-)^Ts  )
	om_operator_matrix_scal_mul(&P_0,lambda/(6.0+lambda),&sum_P_a);
	om_operator_matrix_scal_mul(&sum_P,1.0/(2.0*(6.0+lambda)),&sum_P_b);
	om_operator_matrix_add(&sum_P_a,&sum_P_b,&filter->_P_k_pred);
	om_operator_matrix_add(&filter->_P_k_pred,&filter->_Q,&filter->_P_k_pred);

	omVector b_pred;
	omVector omega;

	om_vector_create(&omega,3);
	om_vector_create(&b_pred,3,om_vector_getValue(&filter->_x_k,3),
							   om_vector_getValue(&filter->_x_k,4),
							   om_vector_getValue(&filter->_x_k,5));

	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_pred,&omega);
	om_kinematics_quaternion(&filter->_q_est,&omega,&filter->_q_pred);



	//correction step
	om_usque_update(manager,filter,sigma_points,sigma_quaternion);


	//free memory
	om_vector_free(&var_x_0);
	om_matrix_free(&tmp_0);
	om_matrix_free(&tmp_0_t);
	om_matrix_free(&P_0);
	om_matrix_free(&sum_P);
	om_matrix_free(&sum_P_a);
	om_matrix_free(&sum_P_b);
	om_matrix_free(&S);
	om_matrix_free(&T);
	om_vector_free(&sum_a);
	om_vector_free(&sum_b);
	om_vector_free(&sum_x);


	for(int i=0;i<13;++i){
		om_vector_free(&sigma_points[i]);
	}

	free(sigma_points);
	sigma_points = 0;
	free(sigma_quaternion);
	sigma_quaternion = 0;

}

void om_usque_update(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter,omVector *sigma_points,omQuaternion *sigma_quaternion){



	//measurement vector;
	omVector z_acc;
	omVector z_mag;
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);

	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_mag);


	omVector z_k;
	om_vector_create(&z_k,6,om_vector_getValue(&z_acc,0),om_vector_getValue(&z_acc,1),om_vector_getValue(&z_acc,2),
							om_vector_getValue(&z_mag,0),om_vector_getValue(&z_mag,1),om_vector_getValue(&z_mag,2));



	double a = filter->_a;
	double f = filter->_f;
	double lambda = filter->_lambda;


	//compute sigma z
	omVector* sigma_z = (omVector*)malloc(13*sizeof(omVector));

	for(unsigned int i=0;i<13;++i){

		omVector tmp_res_m;
		omVector tmp_res_a;

		om_vector_create(&tmp_res_m,3);
		om_vector_create(&tmp_res_a,3);

		om_rotate_vector_quaternion(&sigma_quaternion[i],&ned_gravity,&tmp_res_a);
		om_rotate_vector_quaternion(&sigma_quaternion[i],&ned_geographic_north,&tmp_res_m);

		om_operator_vector_scal_mul(&tmp_res_a,G,&tmp_res_a);


		om_vector_create(&sigma_z[i],6,om_vector_getValue(&tmp_res_a,0),om_vector_getValue(&tmp_res_a,1),om_vector_getValue(&tmp_res_a,2),
									   om_vector_getValue(&tmp_res_m,0),om_vector_getValue(&tmp_res_m,1),om_vector_getValue(&tmp_res_m,2));

        om_vector_free(&tmp_res_m);
        om_vector_free(&tmp_res_a);

	}



	// compute mean of z_k
	omVector mean_z;
	omVector sum_z;
	omVector sum_a;
	omVector sum_b;

	om_vector_create(&sum_z,6,0.0,0.0,0.0,0.0,0.0,0.0);
	om_vector_create(&sum_a,6);
	om_vector_create(&sum_b,6);
	om_vector_create(&mean_z,6);

	for(unsigned int i=1;i<13;++i){
		om_operator_vector_add(&sum_z,&sigma_z[i],&sum_z);
	}

	om_operator_vector_scal_mul(&sigma_z[0],lambda/(6.0+lambda),&sum_a);
	om_operator_vector_scal_mul(&sum_z,1.0/(2.0*(6.0+lambda)),&sum_b);
	om_operator_vector_add(&sum_a,&sum_b,&mean_z);

	om_vector_free(&sum_a);
	om_vector_free(&sum_b);
	om_vector_free(&sum_z);


	// compute covariance and cross-covariance
	omMatrix P_zz;
	omMatrix P_xz;
	omMatrix sum_P_zz;
	omMatrix sum_P_xz;

	omVector var_x_0;
	omMatrix tmp_x_0;
	omVector var_z_0;
	omMatrix tmp_z_0;
	omMatrix tmp_z_0_t;

	omMatrix P_zz_0;
	omMatrix P_xz_0;

	om_matrix_create(&sum_P_zz,6,6);
	om_matrix_create(&sum_P_xz,6,6);
	om_matrix_create(&P_zz,6,6);
	om_matrix_create(&P_xz,6,6);

	om_matrix_create(&P_zz_0,6,6);
	om_matrix_create(&P_xz_0,6,6);
	om_vector_create(&var_x_0,6);
	om_matrix_create(&tmp_x_0,6,1);
	om_vector_create(&var_z_0,6);
	om_matrix_create(&tmp_z_0,6,1);
	om_matrix_create(&tmp_z_0_t,1,6);

	//compute (Z_k(0) - z_k^-)*(Z_k(0) - z_k^-)^T
	// and (X_k(0) - x_k^-)*(Z_k(0) - z_k^-)^T
	om_operator_vector_sub(&sigma_points[0],&filter->_x_k_pred,&var_x_0);
	om_operator_vector_sub(&sigma_z[0],&mean_z,&var_z_0);

	for(int i=0;i<6;i++){
		om_matrix_setValue(&tmp_x_0,i,0,om_vector_getValue(&var_x_0,i));
		om_matrix_setValue(&tmp_z_0,i,0,om_vector_getValue(&var_z_0,i));
	}

	om_matrix_transpose(&tmp_z_0,&tmp_z_0_t);
	om_operator_matrix_mul(&tmp_z_0,&tmp_z_0_t,&P_zz_0);
	om_operator_matrix_mul(&tmp_x_0,&tmp_z_0_t,&P_xz_0);

	for(unsigned int i=1;i<13;++i){

		omVector var_x_i;
		omMatrix tmp_x_i;
		omVector var_z_i;
		omMatrix tmp_z_i;
		omMatrix tmp_z_i_t;

		omMatrix P_zz_i;
		omMatrix P_xz_i;

		om_matrix_create(&P_zz_i,6,6);
		om_matrix_create(&P_xz_i,6,6);

		//allocate
		om_vector_create(&var_x_i,6);
		om_matrix_create(&tmp_x_i,6,1);
		om_vector_create(&var_z_i,6);
		om_matrix_create(&tmp_z_i,6,1);
		om_matrix_create(&tmp_z_i_t,1,6);

		//compute (Z_k(i) - z_k^-)*(Z_k(i) - z_k^-)^T
		// and (X_k(i) - x_k^-)*(Z_k(i) - z_k^-)^T
		om_operator_vector_sub(&sigma_points[i],&filter->_x_k_pred,&var_x_i);
		om_operator_vector_sub(&sigma_z[i],&mean_z,&var_z_i);

		for(int i=0;i<6;i++){
			om_matrix_setValue(&tmp_x_i,i,0,om_vector_getValue(&var_x_i,i));
			om_matrix_setValue(&tmp_z_i,i,0,om_vector_getValue(&var_z_i,i));
		}

		om_matrix_transpose(&tmp_z_i,&tmp_z_i_t);
		om_operator_matrix_mul(&tmp_z_i,&tmp_z_i_t,&P_zz_i);
		om_operator_matrix_mul(&tmp_x_i,&tmp_z_i_t,&P_xz_i);

		om_operator_matrix_add(&sum_P_xz,&P_xz_i,&sum_P_xz);
		om_operator_matrix_add(&sum_P_zz,&P_zz_i,&sum_P_zz);

		//free memory
		om_vector_free(&var_x_i);
		om_vector_free(&var_z_i);

		om_matrix_free(&tmp_z_i);
		om_matrix_free(&tmp_x_i);
		om_matrix_free(&tmp_z_i_t);
		om_matrix_free(&P_xz_i);
		om_matrix_free(&P_zz_i);



	}


	omMatrix sum_P_zz_a;
	omMatrix sum_P_zz_b;
	omMatrix sum_P_xz_a;
	omMatrix sum_P_xz_b;

	om_matrix_create(&sum_P_xz_a,6,6);
	om_matrix_create(&sum_P_xz_b,6,6);
	om_matrix_create(&sum_P_zz_a,6,6);
	om_matrix_create(&sum_P_zz_b,6,6);

	//compute P_k^- = (lambda/(13.0+lambda))* (X_k(0) - x_k^-)*(X_k(0) - x_k^-)^T + (1.0/(2.0*(13.0+lambda)))*sum_{i=0}^n( (X_k(i) - x_k^-)(X_k(i) - x_k^-)^Ts  )
	om_operator_matrix_scal_mul(&P_xz_0,lambda/(6.0+lambda),&sum_P_xz_a);
	om_operator_matrix_scal_mul(&sum_P_xz,1.0/(2.0*(6.0+lambda)),&sum_P_xz_b);
	om_operator_matrix_scal_mul(&P_zz_0,lambda/(6.0+lambda),&sum_P_zz_a);
	om_operator_matrix_scal_mul(&sum_P_zz,1.0/(2.0*(6.0+lambda)),&sum_P_zz_b);

	om_operator_matrix_add(&sum_P_xz_a,&sum_P_xz_b,&P_xz);
	om_operator_matrix_add(&sum_P_zz_a,&sum_P_zz_b,&P_zz);
	om_operator_matrix_add(&filter->_R,&P_zz,&P_zz);


	//Kalman gain
	omMatrix K;
	omMatrix P_zz_inv;

	om_matrix_create(&K,6,6);
	om_matrix_create(&P_zz_inv,6,6);
	om_matrix_inverse(&P_zz,&P_zz_inv);

	om_operator_matrix_mul(&P_xz,&P_zz_inv,&K);

	//correction x_k
	omVector s_z_k;
	omVector corr_x;

	om_vector_create(&s_z_k,6);
	om_vector_create(&corr_x,6);

	//compute x_k_^+ = x_k_^- + K*(z_k - z_k^-)
	om_operator_vector_sub(&z_k,&mean_z,&s_z_k);
	om_operator_matrix_vector_mul(&K,&s_z_k,&corr_x);
	om_operator_vector_add(&filter->_x_k_pred,&corr_x,&filter->_x_k);


	//correction P_k
	omMatrix K_t;
	omMatrix S_tmp;
	omMatrix S;

	om_matrix_create(&S,6,6);
	om_matrix_create(&S_tmp,6,6);
	om_matrix_create(&K_t,6,6);
	om_matrix_transpose(&K,&K_t);

	//compute _P_k^+ = _P_k^- - (K_k*_cov_Z*K_k.transpose());
	om_operator_matrix_mul(&P_zz,&K_t,&S_tmp);
	om_operator_matrix_mul(&K,&S_tmp,&S);
	om_operator_matrix_sub(&filter->_P_k_pred,&S,&filter->_P_k);

	// update of quaternion
	omVector dp_est;
	om_vector_create(&dp_est,3,om_vector_getValue(&filter->_x_k,0),om_vector_getValue(&filter->_x_k,1),om_vector_getValue(&filter->_x_k,2));

	double dp_i_norm = om_vector_norm(&dp_est);
	double dp_i_norm_sq = dp_i_norm*dp_i_norm;
	double dqw_i = ( ( a*(-1.0)*dp_i_norm_sq ) + (f*sqrt( f*f + (1.0 - a*a)*dp_i_norm_sq  ))) / ( (f*f) + dp_i_norm_sq );

	omQuaternion dq;
	om_quat_create(&dq,dqw_i,om_vector_getValue(&dp_est,0)*( (1.0/f) * (a + dqw_i )),
							 om_vector_getValue(&dp_est,1)*( (1.0/f) * (a + dqw_i )),
							 om_vector_getValue(&dp_est,2)*( (1.0/f) * (a + dqw_i )) );

	om_operator_quat_mul(&dq,&filter->_q_pred,&filter->_q_est);

	//set output
	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}


	//free memory
	om_vector_free(&z_k);
	om_vector_free(&dp_est);
	om_vector_free(&mean_z);
	om_vector_free(&s_z_k);
	om_vector_free(&corr_x);
	om_vector_free(&var_x_0);
	om_vector_free(&var_z_0);


	om_matrix_free(&tmp_x_0);
	om_matrix_free(&tmp_z_0);
	om_matrix_free(&tmp_z_0_t);
	om_matrix_free(&K);
	om_matrix_free(&K_t);
	om_matrix_free(&S);
	om_matrix_free(&S_tmp);
	om_matrix_free(&sum_P_xz);
	om_matrix_free(&sum_P_zz);
	om_matrix_free(&P_zz);
	om_matrix_free(&P_zz_inv);
	om_matrix_free(&P_xz);
	om_matrix_free(&sum_P_xz_a);
	om_matrix_free(&sum_P_xz_b);
	om_matrix_free(&sum_P_zz_a);
	om_matrix_free(&sum_P_zz_b);

	for(int i=0;i<13;++i){
		om_vector_free(&sigma_z[i]);
	}

	free(sigma_z);
	sigma_z = 0;

}

///////////////////////////////////////////////////////
/////           NonLinearFilter MEKF              /////
///////////////////////////////////////////////////////


void om_mekf_initialization(struct omSensorFusionManager *manager,void *filter){

	switch(manager->type){

	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_MEKF*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_MEKF*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_MEKF*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_MEKF*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	(*(omNonLinearFilter_MEKF*)filter)._a = 1.0;
	(*(omNonLinearFilter_MEKF*)filter)._f = 4.0;
	(*(omNonLinearFilter_MEKF*)filter)._h = 0.01;
	(*(omNonLinearFilter_MEKF*)filter)._seed = 0.0;

	double var_u = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double var_v = manager->imu_params.variance_gyroscope;

	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._Q,6,6);

	for(int i=0;i<3;++i){

		om_matrix_setValue(&(*(omNonLinearFilter_MEKF*)filter)._Q,i,i,var_v);
		om_matrix_setValue(&(*(omNonLinearFilter_MEKF*)filter)._Q,i+3,i+3,var_u);

	}

	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._R,6,6);
	for(int i=0;i<3;++i){
		om_matrix_setValue(&(*(omNonLinearFilter_MEKF*)filter)._R,i,i,manager->imu_params.variance_accelerometer);
		om_matrix_setValue(&(*(omNonLinearFilter_MEKF*)filter)._R,i+3,i+3,manager->imu_params.variance_magnetometer);
	}

	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._Q_cho,6,6);
	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._R_cho,6,6);

	om_matrix_choleskyDecomposition(&(*(omNonLinearFilter_MEKF*)filter)._Q,&(*(omNonLinearFilter_MEKF*)filter)._Q_cho);
	om_matrix_choleskyDecomposition(&(*(omNonLinearFilter_MEKF*)filter)._R,&(*(omNonLinearFilter_MEKF*)filter)._R_cho);

	double b_x = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double b_y = om_vector_getValue(&manager->imu_params.bias_gyroscope,1);
	double b_z = om_vector_getValue(&manager->imu_params.bias_gyroscope,2);

	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._P_k,6,6);
	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._P_k_pred,6,6);
	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._H,6,6);
	om_matrix_create(&(*(omNonLinearFilter_MEKF*)filter)._F,6,6);
	om_vector_create(&(*(omNonLinearFilter_MEKF*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);
	om_vector_create(&(*(omNonLinearFilter_MEKF*)filter)._x_k_pred,6);
	om_vector_create(&(*(omNonLinearFilter_MEKF*)filter)._v_k,6);
	om_vector_create(&(*(omNonLinearFilter_MEKF*)filter)._w_k,6);



	init_ned_frame();


}


void om_mekf_free(void *filter){

	om_vector_free(&(*(omNonLinearFilter_MEKF*)filter)._w_k);
	om_vector_free(&(*(omNonLinearFilter_MEKF*)filter)._v_k);
	om_vector_free(&(*(omNonLinearFilter_MEKF*)filter)._x_k);
	om_vector_free(&(*(omNonLinearFilter_MEKF*)filter)._x_k_pred);

	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._P_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._P_k);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._Q);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._R_cho);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._Q_cho);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._R);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._F);
	om_matrix_free(&(*(omNonLinearFilter_MEKF*)filter)._H);

}




void om_mekf_f_function(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter,omVector *x,omVector *f_x){

	omVector dp;
	omVector dp_var;
	omVector b;
	omVector omega;

	om_vector_create(f_x,6);

	om_vector_create(&omega,3);
	om_vector_create(&dp_var,3);

	om_vector_create(&dp,3,om_vector_getValue(x,0),om_vector_getValue(x,1),om_vector_getValue(x,2));
	om_vector_create(&b,3,om_vector_getValue(x,3),om_vector_getValue(x,4),om_vector_getValue(x,5));

	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b,&omega);
	om_vector_crossProduct(&omega,&dp,&dp_var);
	om_operator_vector_scal_mul(&dp_var,-1.0,&dp_var);

	for(int i=0;i<3;i++){
		om_vector_setValue(f_x,i,om_vector_getValue(&dp_var,i));
		om_vector_setValue(f_x,i+3,om_vector_getValue(&b,i));
	}

	om_vector_free(&dp);
	om_vector_free(&dp_var);
	om_vector_free(&b);
	om_vector_free(&omega);

}

void om_mekf_h_function(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter,omVector *x,omVector *h_x){

	omQuaternion dq;
	omQuaternion q_var;

	omVector dp;
	omVector tmp_res_a;
	omVector tmp_res_m;

	om_vector_create(h_x,6);
	om_vector_create(&tmp_res_a,3);
	om_vector_create(&tmp_res_m,3);
	om_vector_create(&dp,3,om_vector_getValue(x,0),om_vector_getValue(x,1),om_vector_getValue(x,2));

	double a = filter->_a;
	double f = filter->_f;
	double dp_norm = om_vector_norm(&dp);
	double dp_norm_square = dp_norm*dp_norm;

	double dq_w = ( (a*(-1.0)*dp_norm_square) + (f*sqrt( (f*f) + (1.0 - (a*a))*dp_norm_square))  ) / ( (f*f) + dp_norm_square);
	double dq_x = om_vector_getValue(&dp,0) * ( (1.0/f) * (a + dq_w) );
	double dq_y = om_vector_getValue(&dp,1) * ( (1.0/f) * (a + dq_w) );
	double dq_z = om_vector_getValue(&dp,2) * ( (1.0/f) * (a + dq_w) );

	om_quat_create(&dq,dq_w,dq_x,dq_y,dq_z);
	om_operator_quat_mul(&dq,&filter->_q_pred,&q_var);

	om_rotate_vector_quaternion(&q_var,&ned_geographic_north,&tmp_res_m);
	om_rotate_vector_quaternion(&q_var,&ned_gravity,&tmp_res_a);
	om_operator_vector_scal_mul(&tmp_res_a,9.81,&tmp_res_a);


	for(int i=0;i<3;i++){
		om_vector_setValue(h_x,i,om_vector_getValue(&tmp_res_a,i) + om_vector_getValue(&filter->_v_k,i));
		om_vector_setValue(h_x,i+3,om_vector_getValue(&tmp_res_m,i) + om_vector_getValue(&filter->_v_k,i+3) );
	}

	om_vector_free(&dp);
	om_vector_free(&tmp_res_a);
	om_vector_free(&tmp_res_m);


}

void om_mekf_f_jacobian(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter){

	double h = 0.001;

	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			om_matrix_setValue(&filter->_F,i,j,0.0);


	for (int j = 0; j < 6; j++){

		omVector xhp;
		omVector xhm;
		omVector f_xhm;
		omVector f_xhp;
		omVector value;

		om_vector_create(&value,6);
		om_vector_create(&xhp,6);
		om_vector_create(&xhm,6);

		om_vector_clone(&filter->_x_k,&xhp);
		om_vector_clone(&filter->_x_k,&xhm);

		om_vector_setValue(&xhp,j,om_vector_getValue(&xhp,j)+h);
		om_vector_setValue(&xhm,j,om_vector_getValue(&xhm,j)-h);

		om_mekf_f_function(manager,filter,&xhp,&f_xhp);
		om_mekf_f_function(manager,filter,&xhm,&f_xhm);

		om_operator_vector_sub(&f_xhp,&f_xhm,&value);
		om_operator_vector_scal_div(&value,2.0*h,&value);

		for (int i = 0; i < 6; i++)
			om_matrix_setValue(&filter->_F,i,j,om_vector_getValue(&value,i));

		om_vector_free(&xhp);
		om_vector_free(&xhm);
		om_vector_free(&f_xhp);
		om_vector_free(&f_xhm);
		om_vector_free(&value);


	}

}

void om_mekf_h_jacobian(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter){


	double h = 0.001;

	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			om_matrix_setValue(&filter->_H,i,j,0.0);


	for (int j = 0; j < 6; j++){

		omVector xhp;
		omVector xhm;
		omVector h_xhm;
		omVector h_xhp;
		omVector value;

	    om_vector_create(&value,6);
	    om_vector_create(&xhp,6);
	    om_vector_create(&xhm,6);

		om_vector_clone(&filter->_x_k_pred,&xhp);
		om_vector_clone(&filter->_x_k_pred,&xhm);

		om_vector_setValue(&xhp,j,om_vector_getValue(&xhp,j)+h);
		om_vector_setValue(&xhm,j,om_vector_getValue(&xhm,j)-h);

		om_mekf_h_function(manager,filter,&xhp,&h_xhp);
		om_mekf_h_function(manager,filter,&xhm,&h_xhm);

		om_operator_vector_sub(&h_xhp,&h_xhm,&value);
		om_operator_vector_scal_div(&value,2.0*h,&value);

		for (int i = 0; i < 6; i++)
			om_matrix_setValue(&filter->_H,i,j,om_vector_getValue(&value,i));

		om_vector_free(&xhp);
		om_vector_free(&xhm);
		om_vector_free(&h_xhp);
		om_vector_free(&h_xhm);
		om_vector_free(&value);


	}


}


void om_mekf_process(struct omSensorFusionManager *manager,void *filter){

	om_mekf_prediction(manager,(omNonLinearFilter_MEKF*)filter);
	om_mekf_update(manager,(omNonLinearFilter_MEKF*)filter);
}

void om_mekf_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter){


	om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&filter->_Q_cho,filter->_seed,&filter->_w_k);
	filter->_seed += 6.0;
	om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&filter->_R_cho,filter->_seed,&filter->_v_k);
	filter->_seed += 6.0;

	omVector b;
	omVector omega;

	// propagation of quaternion
	om_vector_create(&b,3,om_vector_getValue(&filter->_x_k,3),om_vector_getValue(&filter->_x_k,4),om_vector_getValue(&filter->_x_k,5));
    om_vector_create(&omega,3);
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b,&omega);
	om_kinematics_quaternion(&filter->_q_est,&omega,&filter->_q_pred);

	// compute the jacobian matrix of function f
	om_mekf_f_jacobian(manager,filter);

	// propagation of state vector

	om_mekf_f_function(manager,filter,&filter->_x_k,&filter->_x_k_pred);

	// propagation of state conversion
	omMatrix F_t;
	omMatrix T;
	om_matrix_create(&T,6,6);
	om_matrix_create(&F_t,6,6);

	om_matrix_transpose(&filter->_F,&F_t);
	om_operator_matrix_mul(&filter->_F,&filter->_P_k,&T);
	om_operator_matrix_mul(&T,&F_t,&filter->_P_k_pred);
	om_operator_matrix_add(&filter->_P_k_pred,&filter->_Q,&filter->_P_k_pred);

	//free memory
	om_vector_free(&b);
	om_vector_free(&omega);
	om_matrix_free(&F_t);
	om_matrix_free(&T);



}

void om_mekf_update(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter){


	// compute the jacobian matrix of function h
	om_mekf_h_jacobian(manager,filter);

	// get the measurement vector
	omVector z;
	om_vector_create(&z,6);
	for(int i=0;i<3;++i){
		om_vector_setValue(&z,i,om_vector_getValue(&manager->imu_data.data_accelerometer,i) - om_vector_getValue(&manager->imu_params.bias_accelerometer,i));
		om_vector_setValue(&z,i+3,om_vector_getValue(&manager->imu_data.data_magnetometer,i));
	}

	// estimate the measurement vector
	omVector z_pred;
	om_mekf_h_function(manager,filter,&filter->_x_k_pred,&z_pred);

	//compute Kalman gain
	omMatrix K;
	omMatrix Ka;
	omMatrix Kb;
	omMatrix Kb_inv;
	omMatrix H_t;

	om_matrix_create(&K,6,6);
	om_matrix_create(&Ka,6,6);
	om_matrix_create(&Kb,6,6);
	om_matrix_create(&Kb_inv,6,6);
	om_matrix_create(&H_t,6,6);
	om_matrix_transpose(&filter->_H,&H_t);

	om_operator_matrix_mul(&filter->_P_k_pred,&H_t,&Ka);
	om_operator_matrix_mul(&filter->_H,&Ka,&Kb);
	om_operator_matrix_add(&Kb,&filter->_R,&Kb);
	om_matrix_inverse(&Kb,&Kb_inv);
	om_operator_matrix_mul(&Ka,&Kb_inv,&K);

	// update state vector
	omVector s_x;
	omVector delta_z;

	om_vector_create(&s_x,6);
	om_vector_create(&delta_z,6);
	om_operator_vector_sub(&z,&z_pred,&delta_z);
	om_operator_matrix_vector_mul(&K,&delta_z,&s_x);

	om_operator_vector_add(&filter->_x_k_pred,&s_x,&filter->_x_k);

	//compute Innovation
	omMatrix I;
	omMatrix S;

	om_matrix_create(&S,6,6);
	om_matrix_createIdentity(&I,6);
	om_operator_matrix_mul(&K,&filter->_H,&S);
	om_operator_matrix_sub(&I,&S,&S);

	// update state covariance
	omMatrix S_t;
	omMatrix K_t;
	omMatrix Pa;
	omMatrix Pb;
	omMatrix T;

	om_matrix_create(&T,6,6);
	om_matrix_create(&Pa,6,6);
	om_matrix_create(&Pb,6,6);
	om_matrix_create(&S_t,6,6);
	om_matrix_create(&K_t,6,6);
	om_matrix_transpose(&S,&S_t);
	om_matrix_transpose(&K,&K_t);

	om_operator_matrix_mul(&S,&filter->_P_k_pred,&T);
	om_operator_matrix_mul(&T,&S_t,&Pa);

	om_operator_matrix_mul(&K,&filter->_R,&T);
	om_operator_matrix_mul(&T,&K_t,&Pb);

	om_operator_matrix_add(&Pa,&Pb,&filter->_P_k);

	// update attitude quaternion
	omQuaternion dq;
	omVector dp_est;
	om_vector_create(&dp_est,3,om_vector_getValue(&filter->_x_k,0),om_vector_getValue(&filter->_x_k,1),om_vector_getValue(&filter->_x_k,2));

	double a = filter->_a;
	double f = filter->_f;
	double dp_norm = om_vector_norm(&dp_est);
	double dp_norm_square = dp_norm*dp_norm;

	double dq_w = ( (a*-1.0*dp_norm_square) + (f*sqrt( (f*f) + (1.0 - (a*a))*dp_norm_square))  ) / ( (f*f) + dp_norm_square);
	double dq_x = om_vector_getValue(&filter->_x_k,0) * ( (1.0/f) * (a + dq_w) );
	double dq_y = om_vector_getValue(&filter->_x_k,1) * ( (1.0/f) * (a + dq_w) );
	double dq_z = om_vector_getValue(&filter->_x_k,2) * ( (1.0/f) * (a + dq_w) );

	om_quat_create(&dq,dq_w,dq_x,dq_y,dq_z);
	om_operator_quat_mul(&dq,&filter->_q_pred,&filter->_q_est);

	// reset error dp to null vector
	for(int i=0;i<3;++i)
		om_vector_setValue(&filter->_x_k,i,0.0);

	//set output
	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}


	//free memory
	om_matrix_free(&S_t);
	om_matrix_free(&S);
	om_matrix_free(&K_t);
	om_matrix_free(&K);
	om_matrix_free(&H_t);
	om_matrix_free(&Pa);
	om_matrix_free(&Pb);
	om_matrix_free(&Ka);
	om_matrix_free(&Kb);
	om_matrix_free(&Kb_inv);
	om_matrix_free(&I);
	om_matrix_free(&T);
	om_vector_free(&z_pred);
	om_vector_free(&z);
	om_vector_free(&s_x);
	om_vector_free(&delta_z);
	om_vector_free(&filter->_x_k_pred);

}



///////////////////////////////////////////////////////
/////           NonLinearFilter REQUEST           /////
///////////////////////////////////////////////////////


void om_request_initialization(struct omSensorFusionManager *manager,void *filter){


	// init q_0
	switch(manager->type){

	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_REQUEST*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_REQUEST*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_REQUEST*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_REQUEST*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	(*(omNonLinearFilter_REQUEST*)filter)._mu_k = 0.01;

	// init weights
	(*(omNonLinearFilter_REQUEST*)filter)._a = (double*)malloc(3*sizeof(double));
	(*(omNonLinearFilter_REQUEST*)filter)._a[0]=0.45;
	(*(omNonLinearFilter_REQUEST*)filter)._a[1]=0.1;
	(*(omNonLinearFilter_REQUEST*)filter)._a[2]=0.45;


	// init measuerment vector & reference vector
	(*(omNonLinearFilter_REQUEST*)filter)._b = (omVector*)malloc(3*sizeof(omVector));
	(*(omNonLinearFilter_REQUEST*)filter)._r = (omVector*)malloc(3*sizeof(omVector));

	for(int i=0;i<3;i++){

		om_vector_create(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],3);
		om_vector_create(&(*(omNonLinearFilter_REQUEST*)filter)._r[i],3,0.0,0.0,0.0);

		om_vector_setValue(&(*(omNonLinearFilter_REQUEST*)filter)._r[i],i,1.0);
		om_rotate_vector_quaternion(&(*(omNonLinearFilter_REQUEST*)filter)._q_est,&(*(omNonLinearFilter_REQUEST*)filter)._r[i],&(*(omNonLinearFilter_REQUEST*)filter)._b[i]);

	}

	// init m_0, sigma_0, B_0 and z_0
	(*(omNonLinearFilter_REQUEST*)filter)._d_m_k = 0.0;
	(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k = 0.0;
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,3,3);
	om_vector_create(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,3,0.0,0.0,0.0);

	for(int i=0;i<3;++i){

		omVector cross;
		omMatrix m_bi;
		omMatrix m_ri_t;
		omMatrix M;

		om_vector_create(&cross,3);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&M,3,3);

		(*(omNonLinearFilter_REQUEST*)filter)._d_m_k += (*(omNonLinearFilter_REQUEST*)filter)._a[i];
		(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k += om_vector_dotProduct(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],&(*(omNonLinearFilter_REQUEST*)filter)._r[i])*(*(omNonLinearFilter_REQUEST*)filter)._a[i];

		// compute _d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		om_vector_crossProduct(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],&(*(omNonLinearFilter_REQUEST*)filter)._r[i],&cross);
		om_operator_vector_scal_mul(&cross,(*(omNonLinearFilter_REQUEST*)filter)._a[i],&cross);
		om_operator_vector_add(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,&cross,&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k);


		//compute _d_B_k = _d_B_k + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],j));
			om_matrix_setValue(&m_ri_t,0,j,om_vector_getValue(&(*(omNonLinearFilter_REQUEST*)filter)._r[i],j));
		}

		om_operator_matrix_mul(&m_bi,&m_ri_t,&M);
		om_operator_matrix_scal_mul(&M,(*(omNonLinearFilter_REQUEST*)filter)._a[i],&M);
		om_operator_matrix_add(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,&M,&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k);

		//free memory
		om_vector_free(&cross);
		om_matrix_free(&m_bi);
		om_matrix_free(&m_ri_t);
		om_matrix_free(&M);


	}

	(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k /= (*(omNonLinearFilter_REQUEST*)filter)._d_m_k;
	om_operator_matrix_scal_div(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,(*(omNonLinearFilter_REQUEST*)filter)._d_m_k,&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k);
	om_operator_vector_scal_div(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,(*(omNonLinearFilter_REQUEST*)filter)._d_m_k,&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k);


	// init S_0
	omMatrix d_B_k_t;
	om_matrix_create(&d_B_k_t,3,3);
	om_matrix_transpose(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,&d_B_k_t);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._d_S_k,3,3);


	//compute _d_S_k =_d_B_k + _d_B_k.transpose();
	om_operator_matrix_add(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,&d_B_k_t,&(*(omNonLinearFilter_REQUEST*)filter)._d_S_k);

	om_matrix_free(&d_B_k_t);

	// init K_est_0
	omMatrix I;
	omMatrix T;

	om_matrix_create(&T,3,3);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._Q,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._R,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._K_k,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._d_K_k,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._K_k_pred,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._P_k,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._P_k_pred,4,4);

	om_matrix_createIdentity(&I,3);
	om_operator_matrix_scal_mul(&I,(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k,&I);
	om_operator_matrix_sub(&(*(omNonLinearFilter_REQUEST*)filter)._d_S_k,&I,&T);


	//set K_est_0 values
	for(int i=0;i<3;i++){

		om_matrix_setValue(&(*(omNonLinearFilter_REQUEST*)filter)._K_k,i,3,om_vector_getValue(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,i));
		om_matrix_setValue(&(*(omNonLinearFilter_REQUEST*)filter)._K_k,3,i,om_vector_getValue(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,i));

		for(int j=0;j<3;j++)
			om_matrix_setValue(&(*(omNonLinearFilter_REQUEST*)filter)._K_k,i,j,om_matrix_getValue(&T,i,j));

	}


	om_matrix_setValue(&(*(omNonLinearFilter_REQUEST*)filter)._K_k,3,3,(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k);

	om_request_computeR(manager,(omNonLinearFilter_REQUEST*)filter);

	om_matrix_clone(&(*(omNonLinearFilter_REQUEST*)filter)._R,&(*(omNonLinearFilter_REQUEST*)filter)._P_k);
	(*(omNonLinearFilter_REQUEST*)filter)._m_k = (*(omNonLinearFilter_REQUEST*)filter)._d_m_k;


	om_matrix_free(&T);
	om_matrix_free(&I);

}

void om_request_free(void *filter){

	om_vector_free(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k);

	for(int i=0;i<3;i++){
		om_vector_free(&(*(omNonLinearFilter_REQUEST*)filter)._r[i]);
		om_vector_free(&(*(omNonLinearFilter_REQUEST*)filter)._b[i]);
	}

	free((*(omNonLinearFilter_REQUEST*)filter)._r);
	free((*(omNonLinearFilter_REQUEST*)filter)._b);
	free((*(omNonLinearFilter_REQUEST*)filter)._a);

	(*(omNonLinearFilter_REQUEST*)filter)._r = 0;
	(*(omNonLinearFilter_REQUEST*)filter)._b = 0;
	(*(omNonLinearFilter_REQUEST*)filter)._a = 0;

	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._K_k);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._K_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._P_k);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._P_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._Q);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._R);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._d_K_k);
	om_matrix_free(&(*(omNonLinearFilter_REQUEST*)filter)._d_S_k);


}






void om_request_process(struct omSensorFusionManager *manager,void *filter){

	om_request_preprocess(manager,(omNonLinearFilter_REQUEST*)filter);
	om_request_prediction(manager,(omNonLinearFilter_REQUEST*)filter);
	om_request_update(manager,(omNonLinearFilter_REQUEST*)filter);

}

void om_request_preprocess(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){

	/* update measurement vector */
	omVector z_acc;
	omVector z_mag;
	omVector y_axis;
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);
	om_vector_create(&y_axis,3);

	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_acc);
	om_vector_normalize(&z_mag);

	om_vector_crossProduct(&z_acc,&z_mag,&y_axis);

	om_vector_normalize(&y_axis);

	om_vector_setValues(&filter->_b[0],3,om_vector_getValue(&z_mag,0),om_vector_getValue(&z_mag,1),om_vector_getValue(&z_mag,2));
	om_vector_setValues(&filter->_b[1],3,om_vector_getValue(&y_axis,0),om_vector_getValue(&y_axis,1),om_vector_getValue(&y_axis,2));
	om_vector_setValues(&filter->_b[2],3,om_vector_getValue(&z_acc,0),om_vector_getValue(&z_acc,1),om_vector_getValue(&z_acc,2));


	om_vector_free(&z_acc);
	om_vector_free(&z_mag);
	om_vector_free(&y_axis);

}

void om_request_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){

	// compute d_m_k, d_sigma_k, d_B_k, d_z_k
	filter->_d_m_k = 0.0;
	filter->_d_sigma_k = 0.0;

	om_vector_setValues(&filter->_d_z_k,3,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j){
			om_matrix_setValue(&filter->_d_B_k,i,j,0.0);
			om_matrix_setValue(&filter->_d_S_k,i,j,0.0);
		}



	for(int i=0;i<3;++i){

		omVector cross;
		omMatrix m_bi;
		omMatrix m_ri_t;
		omMatrix M;

		om_vector_create(&cross,3,0.0,0.0,0.0);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&M,3,3);

		filter->_d_m_k += filter->_a[i];
		filter->_d_sigma_k += om_vector_dotProduct(&filter->_b[i],&filter->_r[i])*filter->_a[i];

		// compute _d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		om_vector_crossProduct(&filter->_b[i],&filter->_r[i],&cross);
		om_operator_vector_scal_mul(&cross,filter->_a[i],&cross);
		om_operator_vector_add(&filter->_d_z_k,&cross,&filter->_d_z_k);


		//compute _d_B_k = _d_B_k + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);

		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&filter->_b[i],j));
			om_matrix_setValue(&m_ri_t,0,j,om_vector_getValue(&filter->_r[i],j));
		}

		om_operator_matrix_mul(&m_bi,&m_ri_t,&M);
		om_operator_matrix_scal_mul(&M,filter->_a[i],&M);
		om_operator_matrix_add(&filter->_d_B_k,&M,&filter->_d_B_k);

		//free memory
		om_vector_free(&cross);
		om_matrix_free(&m_bi);
		om_matrix_free(&m_ri_t);
		om_matrix_free(&M);


	}



	filter->_d_sigma_k /= filter->_d_m_k;
	om_operator_matrix_scal_div(&filter->_d_B_k,filter->_d_m_k,&filter->_d_B_k);
	om_operator_vector_scal_div(&filter->_d_z_k,filter->_d_m_k,&filter->_d_z_k);


	//compute _d_S_k =_d_B_k + _d_B_k.transpose();
	omMatrix d_B_k_t;
	om_matrix_create(&d_B_k_t,3,3);
	om_matrix_transpose(&filter->_d_B_k,&d_B_k_t);
	om_operator_matrix_add(&filter->_d_B_k,&d_B_k_t,&filter->_d_S_k);
	om_matrix_free(&d_B_k_t);

	// compute d_K
	omMatrix I;
	omMatrix T;

	om_matrix_create(&T,3,3);

	for(int i=0;i<4;++i)
		for(int j=0;j<4;++j)
			om_matrix_setValue(&filter->_d_K_k,i,j,0.0);

	om_matrix_createIdentity(&I,3);
	om_operator_matrix_scal_mul(&I,filter->_d_sigma_k,&I);
	om_operator_matrix_sub(&filter->_d_S_k,&I,&T);


	//set d_K_k values
	for(int i=0;i<3;i++){

		om_matrix_setValue(&filter->_d_K_k,i,3,om_vector_getValue(&filter->_d_z_k,i));
		om_matrix_setValue(&filter->_d_K_k,3,i,om_vector_getValue(&filter->_d_z_k,i));

		for(int j=0;j<3;j++)
			om_matrix_setValue(&filter->_d_K_k,i,j,om_matrix_getValue(&T,i,j));

	}

	om_matrix_setValue(&filter->_d_K_k,3,3,filter->_d_sigma_k);

	om_matrix_free(&T);
	om_matrix_free(&I);

	//propagate Q
	om_request_computeQ(manager,filter);

	//propagate R
	om_request_computeR(manager,filter);

	// compute state transition matrix Phy
	omMatrix Phi;
	omMatrix Phi_t;
	omMatrix tmp_K;
	omMatrix tmp_P;
	omVector angular_velocity;

	om_matrix_create(&tmp_K,4,4);
	om_matrix_create(&tmp_P,4,4);
	om_matrix_create(&Phi_t,4,4);
	om_matrix_create(&Phi,4,4);

	om_vector_create(&angular_velocity,3);
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&manager->imu_params.bias_gyroscope,&angular_velocity);

	om_operator_omega_kinematics(&angular_velocity,&Phi);
	om_matrix_transpose(&Phi,&Phi_t);

	// propagate K_k
	om_operator_matrix_mul(&Phi,&filter->_K_k,&tmp_K);
	om_operator_matrix_mul(&tmp_K,&Phi_t,&filter->_K_k_pred);

	// propagate P_k
	om_operator_matrix_mul(&Phi,&filter->_P_k,&tmp_P);
	om_operator_matrix_mul(&tmp_P,&Phi_t,&filter->_P_k_pred);
	om_operator_matrix_add(&filter->_P_k_pred,&filter->_Q,&filter->_P_k_pred);

	// free memory
	om_matrix_free(&Phi);
	om_matrix_free(&Phi_t);
	om_matrix_free(&tmp_K);
	om_matrix_free(&tmp_P);
	om_vector_free(&angular_velocity);



}

void om_request_update(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){

	// compute optimal gain
	double trace_P = om_matrix_trace(&filter->_P_k_pred);
	double trace_R = om_matrix_trace(&filter->_R);
	double a = (filter->_m_k*filter->_m_k)*trace_P;
	double rho_k_opt = ( a )/( a + ((filter->_d_m_k*filter->_d_m_k)*trace_R)  );
	rho_k_opt=0.005;


	// propagate m_k
	double m_kp1 = ((1.0 - rho_k_opt)*filter->_m_k) + (rho_k_opt*filter->_d_m_k);

	// some variable
	double tmp_a = (1.0 - rho_k_opt)*(filter->_m_k/m_kp1);
	double tmp_b = (rho_k_opt)*(filter->_d_m_k/m_kp1);

	// update K
	om_operator_matrix_scal_mul(&filter->_K_k_pred,tmp_a,&filter->_K_k_pred);
	om_operator_matrix_scal_mul(&filter->_d_K_k,tmp_b,&filter->_d_K_k);
	om_operator_matrix_add(&filter->_K_k_pred,&filter->_d_K_k,&filter->_K_k);

	// update P
	om_operator_matrix_scal_mul(&filter->_P_k_pred,tmp_a*tmp_a,&filter->_P_k_pred);
	om_operator_matrix_scal_mul(&filter->_R,tmp_b*tmp_b,&filter->_R);
	om_operator_matrix_add(&filter->_P_k_pred,&filter->_R,&filter->_P_k);

	// get eigen values of K
	omVector* eigen_vector;
	double* eigen_values;
	om_matrix_getEingenValues(&filter->_K_k,&eigen_vector,&eigen_values,50);


	// calcul of lambda_max
	double lambda=0.0;
	for(int i=0;i<4;i++){

		//get the max values
		if(lambda < eigen_values[i]){

			lambda = eigen_values[i];

			double qw = om_vector_getValue(&eigen_vector[i],3);
			double qx = om_vector_getValue(&eigen_vector[i],0);
			double qy = om_vector_getValue(&eigen_vector[i],1);
			double qz = om_vector_getValue(&eigen_vector[i],2);

			// update q_est
			om_quat_create(&filter->_q_est,qw,qx,qy,qz);

		}

		// free memory
		om_vector_free(&eigen_vector[i]);

	}


	//set output
	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}



	// free memory
	free(eigen_values);
	free(eigen_vector);
	eigen_values = 0;
	eigen_vector = 0;



}


void om_request_computeR(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){


	double n_k = 3.0;

	for(int i=0;i<4;++i)
		for(int j=0;j<4;++j)
			om_matrix_setValue(&filter->_R,i,j,0.0);



	/* computation of R22 */
	double R22 = (2.0*filter->_mu_k)/n_k;

	/* computation of R11 */
	omMatrix R11;
	om_matrix_create(&R11,3,3);


	for(int i=0;i<3;i++){

		double dot = om_vector_dotProduct(&filter->_r[i],&filter->_b[i]);

		omMatrix m_bi;
		omMatrix m_ri;
		omMatrix S_ri;
		omMatrix S_ri_t;
		omMatrix I;
		omMatrix tmp_bi_ri;
		omMatrix tmp_ri_bi;
		omMatrix tmp_bi_bi;
		omMatrix m_bi_t;
		omMatrix m_ri_t;
		omMatrix R11_i_a;
		omMatrix R11_i_b;
		omMatrix tmp;

		om_matrix_create(&S_ri,3,1);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri,3,1);

		om_matrix_create(&m_bi_t,1,3);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&S_ri_t,3,3);

		//vectorToMatrix
		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&filter->_b[i],j));
			om_matrix_setValue(&m_ri,j,0,om_vector_getValue(&filter->_r[i],j));
		}

		// compute S_ri = skewSymetricMatrix(_r[i]);
		om_matrix_skewSymmetricMatrix(&filter->_r[i],&S_ri);

		// compute I = identity(3)*(3.0 -(dotProduct(_r[i],_b[i])*dotProduct(_r[i],_b[i])));
		om_matrix_createIdentity(&I,3);
		om_operator_matrix_scal_mul(&I,3.0 - (dot*dot),&I);

		// transpose
		om_matrix_transpose(&m_bi,&m_bi_t);
		om_matrix_transpose(&m_ri,&m_ri_t);
		om_matrix_transpose(&S_ri,&S_ri_t);

		//allocation
		om_matrix_create(&tmp_bi_bi,3,3);
		om_matrix_create(&tmp_ri_bi,3,3);
		om_matrix_create(&tmp_bi_ri,3,3);
		om_matrix_create(&R11_i_a,3,3);
		om_matrix_create(&R11_i_b,3,3);
		om_matrix_create(&tmp,3,3);

		// compute tmp_bi_ri = m_bi*m_ri.transpose()
		om_operator_matrix_mul(&m_bi,&m_ri_t,&tmp_bi_ri);

		// compute tmp_ri_bi = m_ri*m_bi.transpose()
		om_operator_matrix_mul(&m_ri,&m_bi_t,&tmp_ri_bi);

		// compute tmp_bi_bi = m_bi*m_bi.transpose()
		om_operator_matrix_mul(&m_bi,&m_bi_t,&tmp_bi_bi);

		// compute ( I + ( ( (m_bi*m_ri.transpose()) +(m_ri*m_bi.transpose())   )*dotProduct(_b[i],_r[i]) )
		om_operator_matrix_add(&tmp_bi_ri,&tmp_ri_bi,&R11_i_a);
		om_operator_matrix_scal_mul(&R11_i_a,dot,&R11_i_a);
		om_operator_matrix_add(&R11_i_a,&I,&R11_i_a);


		// compute (S_ri*(m_bi*m_bi.transpose())*S_ri.transpose())  )
		om_operator_matrix_mul(&S_ri,&tmp_bi_bi,&tmp);
		om_operator_matrix_mul(&tmp,&S_ri_t,&R11_i_b);

		//compute R11 = R11 + ( I + ( ( (m_bi*m_ri.transpose()) +(m_ri*m_bi.transpose())   )*dotProduct(_b[i],_r[i]) )  + (S_ri*(m_bi*m_bi.transpose())*S_ri.transpose())  );
		om_operator_matrix_add(&R11_i_a,&R11_i_b,&tmp);
		om_operator_matrix_add(&R11,&tmp,&R11);



		//free memory
		om_matrix_free(&m_bi);
		om_matrix_free(&m_ri);
		om_matrix_free(&S_ri);
		om_matrix_free(&S_ri_t);
		om_matrix_free(&I);
		om_matrix_free(&tmp_bi_bi);
		om_matrix_free(&tmp_bi_ri);
		om_matrix_free(&tmp_ri_bi);
		om_matrix_free(&m_bi_t);
		om_matrix_free(&m_ri_t);
		om_matrix_free(&R11_i_a);
		om_matrix_free(&R11_i_b);
		om_matrix_free(&tmp);


	}

	//compute R11 = R11*(_mu_k/n_k);
	om_operator_matrix_scal_mul(&R11,(filter->_mu_k/n_k),&R11);


	//set R values
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			om_matrix_setValue(&filter->_R,i,j,om_matrix_getValue(&R11,i,j));

	om_matrix_setValue(&filter->_R,3,3,R22);

	om_matrix_free(&R11);



}

void om_request_computeQ(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){


	double eta_k = manager->imu_params.variance_gyroscope;

	for(int i=0;i<4;++i)
		for(int j=0;j<4;++j)
			om_matrix_setValue(&filter->_Q,i,j,0.0);



	omMatrix d_B_k_t;
	omMatrix tmp_d_B_k;
	omMatrix tmp_d_B_k_2;
	omMatrix tmp_d_B_k_square;
	omMatrix tmp_d_B_k_t_square;

	om_matrix_create(&d_B_k_t,3,3);
	om_matrix_transpose(&filter->_d_B_k,&d_B_k_t);

	//compute (_d_B_k*_d_B_k.transpose())
	om_matrix_create(&tmp_d_B_k,3,3);
	om_operator_matrix_mul(&filter->_d_B_k,&d_B_k_t,&tmp_d_B_k);

	//compute (_d_B_k.transpose()*_d_B_k)
	om_matrix_create(&tmp_d_B_k_2,3,3);
	om_operator_matrix_mul(&d_B_k_t,&filter->_d_B_k,&tmp_d_B_k_2);

	//compute (_d_B_k*_d_B_k)
	om_matrix_create(&tmp_d_B_k_square,3,3);
	om_operator_matrix_mul(&filter->_d_B_k,&filter->_d_B_k,&tmp_d_B_k_square);

	//compute (_d_B_k.transpose()*_d_B_k.transpose())
	om_matrix_create(&tmp_d_B_k_t_square,3,3);
	om_operator_matrix_mul(&d_B_k_t,&d_B_k_t,&tmp_d_B_k_t_square);



	double norm_d_z_k = om_vector_norm(&filter->_d_z_k);
	double trace = om_matrix_trace(&tmp_d_B_k);

	// computation of Q22
	double Q22 = eta_k *( trace + (filter->_d_sigma_k*filter->_d_sigma_k) + (norm_d_z_k*norm_d_z_k));


	double tmp = (filter->_d_sigma_k*filter->_d_sigma_k) + (norm_d_z_k*norm_d_z_k) - trace;

	omMatrix I;
	omMatrix T;
	omMatrix Q11;

	om_matrix_create(&T,3,3);
	om_matrix_create(&Q11,3,3);

	om_matrix_createIdentity(&I,3);
	om_operator_matrix_scal_mul(&I,tmp,&I);

	om_operator_matrix_sub(&tmp_d_B_k_2,&tmp_d_B_k_square,&T);
	om_operator_matrix_sub(&T,&tmp_d_B_k_t_square,&T);
	om_operator_matrix_scal_mul(&T,2.0,&T);

	// computation of Q11
	om_operator_matrix_add(&I,&T,&Q11);

	/* computation of Q12 et Q21 */
	omMatrix M;
	omMatrix M_t;
	omMatrix D;
	omMatrix I2;
	omMatrix A;
	omVector y;
	omVector x;
	omVector Q12;

	om_vector_create(&x,3);
	om_vector_create(&Q12,3);
	om_matrix_create(&A,3,3);
	om_matrix_create(&M,3,3);
	om_matrix_create(&M_t,3,3);
	om_matrix_create(&D,3,3);
	om_matrix_createIdentity(&I2,3);
	om_operator_matrix_scal_mul(&I2,filter->_d_sigma_k,&I2);

	//compute M = _d_B_k*(_d_B_k - (identity(3)*_d_sigma_k))
	om_operator_matrix_sub(&filter->_d_B_k,&I2,&A);
	om_operator_matrix_mul(&filter->_d_B_k,&A,&M);

	//compute y = operatorVex(M.transpose() - M);
	om_matrix_transpose(&M,&M_t);
	om_operator_matrix_sub(&M_t,&M,&D);
	om_operator_vex(&D,&y);

	//compute  Q12 = (y + (_d_B_k.transpose()*_d_z_k))*(eta_k*(-1));
	om_operator_matrix_vector_mul(&d_B_k_t,&filter->_d_z_k,&x);
	om_operator_vector_add(&y,&x,&Q12);
	om_operator_vector_scal_mul(&Q12,eta_k*(-1.0),&Q12);

	//set Q values
	for(int i=0;i<3;i++){

		om_matrix_setValue(&filter->_Q,i,3,om_vector_getValue(&Q12,i));
		om_matrix_setValue(&filter->_Q,3,i,om_vector_getValue(&Q12,i));

		for(int j=0;j<3;j++)
			om_matrix_setValue(&filter->_Q,i,j,om_matrix_getValue(&Q11,i,j));

	}

	om_matrix_setValue(&filter->_Q,3,3,Q22);
	om_operator_matrix_scal_mul(&filter->_Q,DELTA_T*DELTA_T,&filter->_Q);


	//free memory
	om_matrix_free(&d_B_k_t);
	om_matrix_free(&tmp_d_B_k);
	om_matrix_free(&tmp_d_B_k_2);
	om_matrix_free(&tmp_d_B_k_square);
	om_matrix_free(&tmp_d_B_k_t_square);
	om_matrix_free(&I);
	om_matrix_free(&T);
	om_matrix_free(&Q11);
	om_matrix_free(&M);
	om_matrix_free(&M_t);
	om_matrix_free(&D);
	om_matrix_free(&I2);
	om_matrix_free(&A);

	om_vector_free(&y);
	om_vector_free(&x);
	om_vector_free(&Q12);

}




///////////////////////////////////////////////////////
/////           NonLinearFilter PF                /////
///////////////////////////////////////////////////////



void om_pf_initialization(struct omSensorFusionManager *manager,void *filter){

	// init q_0
	switch(manager->type){

	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_PF*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_PF*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_PF*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_PF*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	(*(omNonLinearFilter_PF*)filter)._seed = 0.0;
	(*(omNonLinearFilter_PF*)filter)._h = 0.1;
	(*(omNonLinearFilter_PF*)filter)._f = 4.0;
	(*(omNonLinearFilter_PF*)filter)._n = 500;
	(*(omNonLinearFilter_PF*)filter)._threeshold = 4.0*( (double)(*(omNonLinearFilter_PF*)filter)._n)/7.0 ;

	// init state
	double b_x = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double b_y = om_vector_getValue(&manager->imu_params.bias_gyroscope,1);
	double b_z = om_vector_getValue(&manager->imu_params.bias_gyroscope,2);

	om_vector_create(&(*(omNonLinearFilter_PF*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);

	// compute Cholesky Decomposition of matrix Q in order to generate noise
	//double var_u = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	//double var_v = manager->imu_params.variance_gyroscope;

	double var_u = 0.000031623;
	double var_v = 0.0031623;


	om_matrix_create(&(*(omNonLinearFilter_PF*)filter)._L,6,6);
	double rate = DELTA_T * 0.01;
	for(int i=0;i<3;i++){
		double tmp = ((var_v*var_v)/rate) + ((var_u*var_u*rate)/12.0);

		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._L,i,i,sqrt(tmp));
		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._L,i+3,i+3, var_u*sqrt(rate));
		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._L,i+3,i,(-0.5)*var_u*sqrt(rate));
	}

	// compute covariance matrix R
	om_matrix_create(&(*(omNonLinearFilter_PF*)filter)._R,6,6);
	for(int i=0;i<3;++i){
		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._R,i,i,manager->imu_params.variance_accelerometer);
		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._R,i+3,i+3,manager->imu_params.variance_magnetometer);
	}


	(*(omNonLinearFilter_PF*)filter)._particle_q = (omQuaternion*)malloc((*(omNonLinearFilter_PF*)filter)._n*sizeof(omQuaternion));
	(*(omNonLinearFilter_PF*)filter)._particle_x = (omVector*)malloc((*(omNonLinearFilter_PF*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_PF*)filter)._particle_wn = (omVector*)malloc((*(omNonLinearFilter_PF*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_PF*)filter)._particle_w = (double*)malloc((*(omNonLinearFilter_PF*)filter)._n*sizeof(double));



	double w_0 = 1.0/ ((double)(*(omNonLinearFilter_PF*)filter)._n);
	double f = (*(omNonLinearFilter_PF*)filter)._f;

	// generation of N particles a time t=0
	for(int i=0;i<(*(omNonLinearFilter_PF*)filter)._n;i++){

		// generate normal random number
		om_vector_create(&(*(omNonLinearFilter_PF*)filter)._particle_wn[i],6,0.0,0.0,0.0,0.0,0.0,0.0);
		om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&(*(omNonLinearFilter_PF*)filter)._L,(*(omNonLinearFilter_PF*)filter)._seed,&(*(omNonLinearFilter_PF*)filter)._particle_wn[i]);
		(*(omNonLinearFilter_PF*)filter)._seed += 6.0;

		//om_vector_display(&(*(omNonLinearFilter_PF*)filter)._particle_wn[i]);

		// perturbation of x_0
		om_vector_create(&(*(omNonLinearFilter_PF*)filter)._particle_x[i],6,0.0,0.0,0.0,0.0,0.0,0.0);
		om_operator_vector_add(&(*(omNonLinearFilter_PF*)filter)._x_k,&(*(omNonLinearFilter_PF*)filter)._particle_wn[i], &(*(omNonLinearFilter_PF*)filter)._particle_x[i]);

		// get d_p_0_i
		omVector dp_0_i;
		om_vector_create(&dp_0_i,3);
		for(int l=0;l<3;l++)
			om_vector_setValue(&dp_0_i,l,om_vector_getValue(&(*(omNonLinearFilter_PF*)filter)._particle_x[i],l));

		// compute d_q_0_i
		omQuaternion dq_0_i;

		double dp_0_i_norm = om_vector_norm(&dp_0_i);
		double dp_0_i_norm_square = dp_0_i_norm*dp_0_i_norm;

		double dq_0_i_w = ( (f*f) - dp_0_i_norm_square)/ ( (f*f) + dp_0_i_norm_square);
		double dq_0_i_x = om_vector_getValue(&(*(omNonLinearFilter_PF*)filter)._particle_x[i],0)* ( (1.0 + dq_0_i_w ) / f);
		double dq_0_i_y = om_vector_getValue(&(*(omNonLinearFilter_PF*)filter)._particle_x[i],1)* ( (1.0 + dq_0_i_w ) / f);
		double dq_0_i_z = om_vector_getValue(&(*(omNonLinearFilter_PF*)filter)._particle_x[i],2)* ( (1.0 + dq_0_i_w ) / f);

		om_quat_create(&dq_0_i,dq_0_i_w,dq_0_i_x,dq_0_i_y,dq_0_i_z);

		// compute q_0_i
		om_operator_quat_mul(&dq_0_i,&(*(omNonLinearFilter_PF*)filter)._q_est,&(*(omNonLinearFilter_PF*)filter)._particle_q[i]);

		// init weights at time 0
		(*(omNonLinearFilter_PF*)filter)._particle_w[i] = w_0;

		om_vector_free(&dp_0_i);

	}

	//init resample to false
	(*(omNonLinearFilter_PF*)filter)._resample = 0;

	//init ned_frame
	init_ned_frame();

}


void om_pf_free(void *filter){

	om_vector_free(&(*(omNonLinearFilter_PF*)filter)._x_k);
	om_matrix_free(&(*(omNonLinearFilter_PF*)filter)._L);
	om_matrix_free(&(*(omNonLinearFilter_PF*)filter)._R);


	for(int i=0;i<(*(omNonLinearFilter_PF*)filter)._n;i++){
		om_vector_free(&(*(omNonLinearFilter_PF*)filter)._particle_wn[i]);
		om_vector_free(&(*(omNonLinearFilter_PF*)filter)._particle_x[i]);
	}

	free((*(omNonLinearFilter_PF*)filter)._particle_q);
	free((*(omNonLinearFilter_PF*)filter)._particle_x);
	free((*(omNonLinearFilter_PF*)filter)._particle_w);
	free((*(omNonLinearFilter_PF*)filter)._particle_wn);

	(*(omNonLinearFilter_PF*)filter)._particle_q = 0;
	(*(omNonLinearFilter_PF*)filter)._particle_x = 0;
	(*(omNonLinearFilter_PF*)filter)._particle_w = 0;
	(*(omNonLinearFilter_PF*)filter)._particle_wn = 0;


}


void om_pf_process(struct omSensorFusionManager *manager,void *filter){

	(*(omNonLinearFilter_PF*)filter)._seed = (*(omNonLinearFilter_PF*)filter)._seed > 100000000.0 ? 0.0 : (*(omNonLinearFilter_PF*)filter)._seed;

	om_pf_prediction(manager,(omNonLinearFilter_PF*)filter);
	om_pf_update(manager,(omNonLinearFilter_PF*)filter);
}

void om_pf_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_PF *filter){

	//get measurement values
	omVector z_acc_tmp;
	omVector z_mag_tmp;
	omVector z_k;

	om_vector_create(&z_k,6);
	om_vector_create(&z_acc_tmp,3);
	om_vector_create(&z_mag_tmp,3);
	om_vector_clone(&manager->imu_data.data_accelerometer,&z_acc_tmp);
	om_vector_clone(&manager->imu_data.data_magnetometer,&z_mag_tmp);

	om_vector_normalize(&z_acc_tmp);
	om_vector_normalize(&z_mag_tmp);

	for(int l=0;l<3;l++){
		om_vector_setValue(&z_k,l,om_vector_getValue(&z_acc_tmp,l));
		om_vector_setValue(&z_k,l+3,om_vector_getValue(&z_mag_tmp,l));
	}

	om_vector_free(&z_acc_tmp);
	om_vector_free(&z_mag_tmp);

	// propagate quaternion
	omVector b;
	omVector omega;

	om_vector_create(&omega,3);
	om_vector_create(&b,3,om_vector_getValue(&filter->_x_k,3),
 						  om_vector_getValue(&filter->_x_k,4),
					      om_vector_getValue(&filter->_x_k,5));

	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b,&omega);
	om_kinematics_quaternion(&filter->_q_est,&omega,&filter->_q_pred);

	om_vector_free(&b);
	om_vector_free(&omega);

	//set sum_w to 0.0
	filter->_sum_w = 0.0;

	//compute
	omQuaternion q_pred_inv;
	om_quat_inverse(&filter->_q_pred,&q_pred_inv);

	// compute N particle of angular velocity and propagate _q_k_i to _q_(k+1)_i
	for(int i=0;i<filter->_n;++i){

		//om_vector_display(&(*(omNonLinearFilter_PF*)filter)._particle_x[i]);

		// generate normal random number
		omVector noise_dp;
		omVector noise_b;

		om_vector_create(&noise_dp,3,om_vector_getValue(&filter->_particle_wn[i],0),
									 om_vector_getValue(&filter->_particle_wn[i],1),
									 om_vector_getValue(&filter->_particle_wn[i],2));

		om_vector_create(&noise_b,3,om_vector_getValue(&filter->_particle_wn[i],3),
									om_vector_getValue(&filter->_particle_wn[i],4),
									om_vector_getValue(&filter->_particle_wn[i],5));
		// compute angular velocity
		omVector b_i;
		omVector angular_velocity_i;

		om_vector_create(&angular_velocity_i,3);
		om_vector_create(&b_i,3,om_vector_getValue(&filter->_particle_x[i],3),
								om_vector_getValue(&filter->_particle_x[i],4),
								om_vector_getValue(&filter->_particle_x[i],5));

		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);
		om_operator_vector_add(&angular_velocity_i,&noise_dp,&angular_velocity_i);

		// propagate _q_k_i to _q_(k+1)_i
		om_kinematics_quaternion(&filter->_particle_q[i],&angular_velocity_i,&filter->_particle_q[i]);

		// propagate bias_k_i to bias_(k+1)_i
		om_operator_vector_add(&b_i,&noise_b,&b_i);

		// compute d_q_k_i
		omQuaternion dq_i;
		om_operator_quat_mul(&filter->_particle_q[i],&q_pred_inv,&dq_i);

		// propagate d_p_k_i to d_p_(k+1)_i
		omVector dp_i;
		double tmp = (filter->_f * (dq_i._qw < 0.0 ? -1.0 : 1.0) )/(1.0 + fabs(dq_i._qw));
		om_vector_create(&dp_i,3,dq_i._qx*tmp,dq_i._qy*tmp,dq_i._qz*tmp);

		// update _x_k_i
		for(int l=0;l<3;l++){
			om_vector_setValue(&filter->_particle_x[i],l,om_vector_getValue(&dp_i,l));
			om_vector_setValue(&filter->_particle_x[i],l+3,om_vector_getValue(&b_i,l));
		}
		// compute z_k_i
		omVector z_acc;
		omVector z_mag;
		omVector z_k_i;

		om_vector_create(&z_acc,3);
		om_vector_create(&z_mag,3);
		om_vector_create(&z_k_i,6);

		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_gravity,&z_acc);
		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_geographic_north,&z_mag);


		for(int l=0;l<3;l++){
			om_vector_setValue(&z_k_i,l,om_vector_getValue(&z_acc,l));
			om_vector_setValue(&z_k_i,l+3,om_vector_getValue(&z_mag,l));
		}
		// likehood function
		omVector s_k_i;
		om_vector_create(&s_k_i,6);
		om_operator_vector_sub(&z_k,&z_k_i,&s_k_i);

		//om_vector_display(&s_k_i);

		double tmp_acc = 0.0;
		double tmp_mag = 0.0;

		for(int l=0;l<3;l++){
			tmp_acc += om_vector_getValue(&s_k_i,l)*om_vector_getValue(&s_k_i,l);
			tmp_mag += om_vector_getValue(&s_k_i,l+3)*om_vector_getValue(&s_k_i,l+3);
		}

		tmp_acc *= (1.0/manager->imu_params.variance_accelerometer);
		tmp_mag *= (1.0/manager->imu_params.variance_magnetometer);

		double L_k_i = exp( (tmp_acc+tmp_mag)*(-0.5) );

		// update weights
		filter->_particle_w[i] *= L_k_i;

		//printf("L_k_i = %.*f\n",30,L_k_i);
		//printf("tmp_mag = %.*f\n",30,tmp_mag);


		// accumulator to normalize weights
		filter->_sum_w += filter->_particle_w[i];

		//free space
		om_vector_free(&b_i);
		om_vector_free(&angular_velocity_i);
		om_vector_free(&noise_b);
		om_vector_free(&noise_dp);
		om_vector_free(&dp_i);

		om_vector_free(&z_acc);
		om_vector_free(&z_mag);
		om_vector_free(&z_k_i);
		om_vector_free(&s_k_i);

	}
	om_vector_free(&z_k);


}

void om_pf_update(struct omSensorFusionManager *manager,omNonLinearFilter_PF *filter){

	// compute mean x_k^p
	omVector x_k_sigma;

	om_vector_setValues(&filter->_x_k,6,0.0,0.0,0.0,0.0,0.0,0.0);
	om_vector_create(&x_k_sigma,6,0.0,0.0,0.0,0.0,0.0,0.0);

	filter->_sum_n_eff = 0.0;

	//upadate particles
	for(int i=0;i<filter->_n;++i){

		//update weights
		filter->_particle_w[i] /= filter->_sum_w;

		//compute N_eff for resampling step
		filter->_sum_n_eff += (filter->_particle_w[i]*filter->_particle_w[i]);

		//update state vector
		omVector tmp_x;
		om_vector_create(&tmp_x,6);
		om_operator_vector_scal_mul(&filter->_particle_x[i],filter->_particle_w[i],&tmp_x);

		om_operator_vector_add(&filter->_x_k,&tmp_x,&filter->_x_k);
		om_operator_vector_add(&x_k_sigma,&filter->_particle_x[i],&x_k_sigma);

		om_vector_free(&tmp_x);
	}


	om_operator_vector_scal_div(&x_k_sigma,(double)(filter->_n),&x_k_sigma);

	//om_vector_display(&filter->_x_k);
	//om_vector_display(&x_k_sigma);

	// compute dq
	omVector dp;
	omQuaternion dq;

	om_vector_create(&dp,3,om_vector_getValue(&filter->_x_k,0),
						   om_vector_getValue(&filter->_x_k,1),
						   om_vector_getValue(&filter->_x_k,2));

	double dp_norm = om_vector_norm(&dp);
	double dp_norm_square = dp_norm*dp_norm;

	double dq_w = ( (filter->_f*filter->_f) - dp_norm_square)/ ( (filter->_f*filter->_f) + dp_norm_square);
	double dq_x = om_vector_getValue(&dp,0)* ( (1.0 + dq_w ) / filter->_f);
	double dq_y = om_vector_getValue(&dp,1)* ( (1.0 + dq_w ) / filter->_f);
	double dq_z = om_vector_getValue(&dp,2)* ( (1.0 + dq_w ) / filter->_f);

	om_quat_create(&dq,dq_w,dq_x,dq_y,dq_z);

	// update q_est
	om_operator_quat_mul(&dq,&filter->_q_pred,&filter->_q_est);


	//set output
	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}

	// resampling step
	 om_pf_resampling(filter);


	// if resample step has been performed
	 if(filter->_resample == 1){

		 //compute convariance matrix S
		 omMatrix S;
		 om_matrix_create(&S,6,6);

		 for(int i=0;i<filter->_n;++i){

			 omVector dx;
			 om_vector_create(&dx,6);
			 om_operator_vector_sub(&filter->_particle_x[i],&x_k_sigma,&dx);

			 omMatrix tmp_S_i_t;
			 omMatrix tmp_S_i;
			 om_matrix_create(&tmp_S_i,6,1);
			 om_matrix_create(&tmp_S_i_t,1,6);
			 om_matrix_setColumn(&tmp_S_i,0,&dx);
			 om_matrix_transpose(&tmp_S_i,&tmp_S_i_t);

			 omMatrix S_i;
			 om_matrix_create(&S_i,6,6);
			 om_operator_matrix_mul(&tmp_S_i,&tmp_S_i_t,&S_i);

			 om_operator_matrix_add(&S,&S_i,&S);


			 om_matrix_free(&S_i);
			 om_matrix_free(&tmp_S_i);
			 om_matrix_free(&tmp_S_i_t);
			 om_vector_free(&dx);
		 }

		 om_operator_matrix_scal_div(&S,filter->_n-1,&S);
		 om_operator_matrix_scal_mul(&S,(filter->_h*filter->_h),&S);

		 // compute cholesky decomposition of S
		 omMatrix S_cho;
		 om_matrix_create(&S_cho,6,6);
		 om_matrix_choleskyDecomposition(&S,&S_cho);

		 /*
		 if(om_matrix_containsNaN(&S_cho)){
			 for(int i = 0;i<6;i++)
				 for(int j = 0;j<6;j++)
					 om_matrix_setValue(&S_cho,i,j,0.0);

			 om_matrix_squareRoot(&S,&S_cho,1);
		 }*/

		 //om_matrix_display(&S);
		 //printf("\n");
		 //om_matrix_display(&S_cho);

		 // perturb particles with noise
		 for(int i=0;i<filter->_n;++i){

			 // generate gaussian noise
			 omVector noise;
			 om_vector_create(&noise,6);
			 om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&S,filter->_seed,&noise);
			 filter->_seed += 6.0;

			 // perturb x^k_i
			 om_operator_vector_add(&filter->_particle_x[i],&noise,&filter->_particle_x[i]);

			 // "adding" pertupation to _particles_q_k
			 // compute dq
			 omVector dp_i;
			 omQuaternion dq_i;

			 //om_vector_display(&filter->_particle_x[i]);

			 om_vector_create(&dp_i,3,om_vector_getValue(&filter->_particle_x[i],0),
								   	  om_vector_getValue(&filter->_particle_x[i],1),
									  om_vector_getValue(&filter->_particle_x[i],2));

			 double dp_i_norm = om_vector_norm(&dp_i);
			 double dp_i_norm_square = dp_i_norm*dp_i_norm;

			 double dq_i_w = ( (filter->_f*filter->_f) - dp_i_norm_square)/ ( (filter->_f*filter->_f) + dp_i_norm_square);
			 double dq_i_x = om_vector_getValue(&dp_i,0)* ( (1.0 + dq_i_w ) / filter->_f);
			 double dq_i_y = om_vector_getValue(&dp_i,1)* ( (1.0 + dq_i_w ) / filter->_f);
			 double dq_i_z = om_vector_getValue(&dp_i,2)* ( (1.0 + dq_i_w ) / filter->_f);

			 om_quat_create(&dq_i,dq_i_w,dq_i_x,dq_i_y,dq_i_z);
			 om_operator_quat_mul(&dq_i,&filter->_q_pred,&filter->_particle_q[i]);

			 om_vector_free(&dp_i);
			 om_vector_free(&noise);

		 }

		 om_matrix_free(&S);
		 om_matrix_free(&S_cho);
	 }


	 om_vector_free(&x_k_sigma);
	 om_vector_free(&dp);
}

void om_pf_resampling(omNonLinearFilter_PF *filter){

	filter->_resample = 0;

	double N_eff = floor(1.0/filter->_sum_n_eff);


	if(N_eff < filter->_threeshold){


		int n = (int)N_eff;
		filter->_resample = 1;

		om_pf_quicksort(filter,0,filter->_n-1);

		omVector* particle_chosen = (omVector*)malloc(n*sizeof(omVector));

		for (int l = 0; l < n; ++l){
			int ind = (filter->_n - l - 1);
			om_vector_create(&particle_chosen[l],6);
			om_vector_clone(&filter->_particle_x[ind],&particle_chosen[l]);
		}



		for (int k = 0; k < (filter->_n - n); ++k) {
			int index = k%n;

			for (int j = 0; j < 6; ++j)
				om_vector_setValue(&filter->_particle_x[k],j,om_vector_getValue(&particle_chosen[index],j));

			filter->_particle_w[k] = 1.0/(double)(filter->_n);
		}

		for (int l = 0; l < n; ++l)
			om_vector_free(&particle_chosen[l]);

		free(particle_chosen);
		particle_chosen=0;

	}

}

void om_pf_swap(omNonLinearFilter_PF *filter,int i,int j){

	//swap weight
	double tmp_w = filter->_particle_w[i];
	filter->_particle_w[i] = filter->_particle_w[j];
	filter->_particle_w[j] = tmp_w;

	//swap state vector
	for (int l = 0; l < 6; ++l){
		double tmp_x = om_vector_getValue(&filter->_particle_x[i],l);
		om_vector_setValue(&filter->_particle_x[i],l,om_vector_getValue(&filter->_particle_x[j],l));
		om_vector_setValue(&filter->_particle_x[j],l,tmp_x);
	}

}


void om_pf_quicksort(omNonLinearFilter_PF *filter,int left, int right){

	int min = (left+right)/2;

	int i = left;
	int j = right;
    double pivot = filter->_particle_w[min];

    while(left<j || i<right)
    {
        while(filter->_particle_w[i]<pivot)
        i++;
        while(filter->_particle_w[j]>pivot)
        j--;

        if(i<=j){
        	om_pf_swap(filter,i,j);
            i++;
            j--;
        }
        else{
            if(left<j)
            	om_pf_quicksort(filter,left, j);
            if(i<right)
            	om_pf_quicksort(filter,i,right);
            return;
        }
    }

}



///////////////////////////////////////////////////////
/////           NonLinearFilter GDOF              /////
///////////////////////////////////////////////////////



/* initialization of all component used by the nonlinear filter CGO */
void om_gdof_initialization(struct omSensorFusionManager *manager,void *filter){

	// constant initialization
	(*(omNonLinearFilter_GDOF*)filter)._eta = 0.003;
	(*(omNonLinearFilter_GDOF*)filter)._beta = 0.005;


	// attitude initialization
	om_quat_create(&(*(omNonLinearFilter_GDOF*)filter)._q_est,1.0,0.0,0.0,0.0);
	switch(manager->type){

	// convertion according to user choices
	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_GDOF*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_GDOF*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_GDOF*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_GDOF*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}


	om_vector_create(&(*(omNonLinearFilter_GDOF*)filter)._bias_est,3);
	om_vector_create(&(*(omNonLinearFilter_GDOF*)filter)._e_b,3);
	om_vector_create(&(*(omNonLinearFilter_GDOF*)filter)._e_m,3);
	om_vector_create(&(*(omNonLinearFilter_GDOF*)filter)._f_a,3);
	om_vector_create(&(*(omNonLinearFilter_GDOF*)filter)._f_b,3);

	om_matrix_create(&(*(omNonLinearFilter_GDOF*)filter)._J_a,3,4);
	om_matrix_create(&(*(omNonLinearFilter_GDOF*)filter)._J_b,3,4);

	om_vector_clone(&manager->imu_params.bias_gyroscope,&(*(omNonLinearFilter_GDOF*)filter)._bias_est);

	// initialization of North East Down frame vector
	init_ned_frame();

}

/* process function of the nonlinear filter CGO */
void om_gdof_process(struct omSensorFusionManager *manager,void *filter){

	om_gdof_prediction(manager,(omNonLinearFilter_GDOF*)filter);
	om_gdof_update(manager,(omNonLinearFilter_GDOF*)filter);

}



void om_gdof_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter){


	omVector acc;
	omVector mag;

	om_vector_create(&acc,3);
	om_vector_create(&mag,3);

	om_vector_clone(&manager->imu_data.data_accelerometer,&acc);
	om_vector_clone(&manager->imu_data.data_magnetometer,&mag);

	om_vector_normalize(&acc);
	om_vector_normalize(&mag);


	omQuaternion q_inv;
	om_quat_inverse(&filter->_q_est,&q_inv);

	om_rotate_vector_quaternion(&q_inv,&manager->imu_data.data_magnetometer,&filter->_e_m);

	double b_x = sqrt( (om_vector_getValue(&filter->_e_m,0)*om_vector_getValue(&filter->_e_m,0) ) + (om_vector_getValue(&filter->_e_m,1)*om_vector_getValue(&filter->_e_m,1) ));
	double b_z = om_vector_getValue(&filter->_e_m,2);

	om_vector_setValues(&filter->_e_b,3,b_x,0.0,b_z);

	double q1 = filter->_q_est._qw;
	double q2 = filter->_q_est._qx;
	double q3 = filter->_q_est._qy;
	double q4 = filter->_q_est._qz;

	double f_a_x = (2.0*( (q2*q4) - (q1*q3) ))       - om_vector_getValue(&acc,0);
	double f_a_y = (2.0*( (q1*q2) + (q4*q3) ))       - om_vector_getValue(&acc,1);
	double f_a_z = (2.0*( 0.5 - (q2*q2) - (q3*q3) )) - om_vector_getValue(&acc,2);

	double f_b_x = (2.0*b_x*(0.5 - (q3*q3) -(q4*q4))) + (2.0*b_z*( (q2*q4) - (q1*q3) )) - om_vector_getValue(&mag,0);
	double f_b_y = (2.0*b_x*( (q2*q3) - (q1*q4) )) + (2.0*b_z*( (q1*q2) + (q4*q3) ))    - om_vector_getValue(&mag,1);
	double f_b_z = (2.0*b_x*( (q1*q3) + (q2*q4))) + (2.0*b_z*(0.5 - (q3*q3) -(q2*q2)))  - om_vector_getValue(&mag,2);

	om_vector_setValues(&filter->_f_a,3,f_a_x,f_a_y,f_a_z);
	om_vector_setValues(&filter->_f_b,3,f_b_x,f_b_y,f_b_z);

	om_matrix_setValue(&filter->_J_a,0,0,-2.0*q3);
	om_matrix_setValue(&filter->_J_a,0,1, 2.0*q4);
	om_matrix_setValue(&filter->_J_a,0,2,-2.0*q1);
	om_matrix_setValue(&filter->_J_a,0,3, 2.0*q2);

	om_matrix_setValue(&filter->_J_a,1,0, 2.0*q2);
	om_matrix_setValue(&filter->_J_a,1,1, 2.0*q1);
	om_matrix_setValue(&filter->_J_a,1,2, 2.0*q4);
	om_matrix_setValue(&filter->_J_a,1,3, 2.0*q3);

	om_matrix_setValue(&filter->_J_a,2,0, 0.0);
	om_matrix_setValue(&filter->_J_a,2,1, -4.0*q2);
	om_matrix_setValue(&filter->_J_a,2,2, -4.0*q3);
	om_matrix_setValue(&filter->_J_a,2,3, 0.0);


	om_matrix_setValue(&filter->_J_b,0,0, -2.0*b_z*q3);
	om_matrix_setValue(&filter->_J_b,0,1,  2.0*b_z*q4);
	om_matrix_setValue(&filter->_J_b,0,2, (-4.0*b_x*q3) - (2.0*b_z*q1));
	om_matrix_setValue(&filter->_J_b,0,3, (-4.0*b_x*q4) + (2.0*b_z*q2));

	om_matrix_setValue(&filter->_J_b,1,0, (-2.0*b_x*q4) + (2.0*b_z*q2));
	om_matrix_setValue(&filter->_J_b,1,1, (2.0*b_x*q3) + (2.0*b_z*q1));
	om_matrix_setValue(&filter->_J_b,1,2, (2.0*b_x*q2) + (2.0*b_z*q4));
	om_matrix_setValue(&filter->_J_b,1,3, (-2.0*b_x*q1) + (2.0*b_z*q3));

	om_matrix_setValue(&filter->_J_b,2,0, 2.0*b_x*q3);
	om_matrix_setValue(&filter->_J_b,2,1, (2.0*b_x*q4) - (4.0*b_z*q2));
	om_matrix_setValue(&filter->_J_b,2,2, (2.0*b_x*q1) - (4.0*b_z*q3));
	om_matrix_setValue(&filter->_J_b,2,3, 2.0*b_x*q2);


}


void om_gdof_update(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter){

	// variables
	omVector gain_a;
	omVector gain_b;
	omVector gain;
	omQuaternion dq;
	omMatrix J_a_t;
	omMatrix J_b_t;
	omQuaternion q_bias;
	omVector angular_velocity;

	// allocation
	om_vector_create(&gain_a,4);
	om_vector_create(&gain_b,4);
	om_vector_create(&gain,4);
	om_matrix_create(&J_a_t,4,3);
	om_matrix_create(&J_b_t,4,3);
	om_vector_create(&angular_velocity,3);


	// compute gain

	om_matrix_transpose(&filter->_J_a,&J_a_t);
	om_matrix_transpose(&filter->_J_b,&J_b_t);

	om_operator_matrix_vector_mul(&J_a_t,&filter->_f_a,&gain_a);
	om_operator_matrix_vector_mul(&J_b_t,&filter->_f_b,&gain_b);
	om_operator_vector_add(&gain_a,&gain_b,&gain);

	om_quat_create(&dq,om_vector_getValue(&gain,0),om_vector_getValue(&gain,1),om_vector_getValue(&gain,2),om_vector_getValue(&gain,3));

	//compute gyro bias estimation

	om_operator_quat_mul(&filter->_q_est,&dq,&q_bias);
	om_operator_quat_scal_mul(&q_bias,2.0,&q_bias);
	om_quat_imaginary(&q_bias,&filter->_bias_est);
	om_operator_vector_scal_mul(&filter->_bias_est,filter->_eta,&filter->_bias_est);

	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&filter->_bias_est,&angular_velocity);

	// normalize gain
	om_quat_normalize(&dq);
	om_operator_quat_scal_mul(&dq,filter->_beta,&dq);

	// correction
	om_kinematics_quaternion(&filter->_q_est,&angular_velocity,&filter->_q_est);
	om_operator_quat_sub(&filter->_q_est,&dq,&filter->_q_est);
	om_quat_normalize(&filter->_q_est);


	// set output
	double qw = filter->_q_est._qw;
	double qx = filter->_q_est._qx;
	double qy = filter->_q_est._qy;
	double qz = filter->_q_est._qz;

	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,qw,qx,qy,qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&filter->_q_est,1.0,0.0,0.0,0.0);
		break;

	}

	// free memory
	om_vector_free(&gain_a);
	om_vector_free(&gain_b);
	om_vector_free(&gain);
	om_vector_free(&angular_velocity);
	om_matrix_free(&J_a_t);
	om_matrix_free(&J_b_t);

}



/* release all component used for the nonlinear filter GDOF */
void om_gdof_free(void *filter){

	//free memory
	om_vector_free(&(*(omNonLinearFilter_GDOF*)filter)._bias_est);
	om_vector_free(&(*(omNonLinearFilter_GDOF*)filter)._e_b);
	om_vector_free(&(*(omNonLinearFilter_GDOF*)filter)._e_m);
	om_vector_free(&(*(omNonLinearFilter_GDOF*)filter)._f_b);
	om_vector_free(&(*(omNonLinearFilter_GDOF*)filter)._f_a);
	om_matrix_free(&(*(omNonLinearFilter_GDOF*)filter)._J_a);
	om_matrix_free(&(*(omNonLinearFilter_GDOF*)filter)._J_b);

}
///////////////////////////////////////////////////////
/////           NonLinearFilter CFA               /////
///////////////////////////////////////////////////////



/* initialization of all component used by the nonlinear filter CGO */
void om_cfa_initialization(struct omSensorFusionManager *manager,void *filter){

	// constant initialization
	(*(omNonLinearFilter_CFA*)filter)._lambda = 0.00001;
	(*(omNonLinearFilter_CFA*)filter)._beta = 0.5;

	om_quat_create(&(*(omNonLinearFilter_CFA*)filter)._q_est,1.0,0.0,0.0,0.0);
	// attitude initialization
	switch(manager->type){

	// convertion according to user choices
	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_CFA*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_CFA*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_CFA*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_CFA*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	om_vector_create(&(*(omNonLinearFilter_CFA*)filter)._v_acc_pred,3);
	om_vector_create(&(*(omNonLinearFilter_CFA*)filter)._v_mag_pred,3);

	// initialization of North East Down frame vector
	init_ned_frame();

}

/* process function of the nonlinear filter CGO */
void om_cfa_process(struct omSensorFusionManager *manager,void *filter){

	om_cfa_prediction(manager,(omNonLinearFilter_CFA*)filter);
	om_cfa_update(manager,(omNonLinearFilter_CFA*)filter);

}



void om_cfa_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CFA *filter){

	if(  om_vector_norm(&manager->imu_data.data_accelerometer) - G < 0.03 )
		filter->_beta = 0.05;
	else
		filter->_beta = 0.5;

	om_kinematics_quaternion(&filter->_q_est,&manager->imu_data.data_gyroscope,&filter->_q_pred);

	om_rotate_vector_quaternion(&filter->_q_pred,&ned_gravity,&filter->_v_acc_pred);
	om_rotate_vector_quaternion(&filter->_q_pred,&ned_geographic_north,&filter->_v_mag_pred);

	om_operator_vector_scal_mul(&filter->_v_acc_pred,G,&filter->_v_acc_pred);


}


void om_cfa_update(struct omSensorFusionManager *manager,omNonLinearFilter_CFA *filter){

	omMatrix X;
	omMatrix K;
	omMatrix K_tmp;
	omMatrix K_inv;
	omMatrix I;
	omMatrix X_t;
	omMatrix S_mag;
	omMatrix S_acc;
	omVector s_acc;
	omVector s_mag;
	omVector s;
	omVector gain;
	omQuaternion dq;

	om_vector_create(&s_acc,3);
	om_vector_create(&s_mag,3);
	om_vector_create(&s,6);
	om_vector_create(&gain,3);

	om_matrix_create(&X,3,6);
	om_matrix_create(&K,3,6);
	om_matrix_createIdentity(&I,3);
	om_matrix_create(&X_t,6,3);
	om_matrix_create(&S_acc,3,3);
	om_matrix_create(&S_mag,3,3);
	om_matrix_create(&K_tmp,3,3);
	om_matrix_create(&K_inv,3,3);

	om_matrix_skewSymmetricMatrix(&filter->_v_acc_pred,&S_acc);
	om_matrix_skewSymmetricMatrix(&filter->_v_mag_pred,&S_mag);

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
				om_matrix_setValue(&X,i,j,om_matrix_getValue(&S_acc,i,j)*(-2.0));
				om_matrix_setValue(&X,i,j+3,om_matrix_getValue(&S_mag,i,j)*(-2.0));
		}


	// compute K = ( (X*X.transpose()) + (identity(3)*_lambda)).inverse()*X;
	om_matrix_transpose(&X,&X_t);
	om_operator_matrix_mul(&X,&X_t,&K_tmp);
	om_operator_matrix_scal_mul(&I,filter->_lambda,&I);
	om_operator_matrix_add(&K_tmp,&I,&K_tmp);
	om_matrix_inverse(&K_tmp,&K_inv);
	om_operator_matrix_mul(&K_inv,&X,&K);

	//compute gain
	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&filter->_v_acc_pred,&s_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&filter->_v_mag_pred,&s_mag);

	for(int i=0;i<3;i++){
		om_vector_setValue(&s,i,om_vector_getValue(&s_acc,i));
		om_vector_setValue(&s,i+3,om_vector_getValue(&s_mag,i));
	}

	om_operator_matrix_vector_mul(&K,&s,&gain);
	om_quat_create(&dq,1.0,om_vector_getValue(&gain,0),om_vector_getValue(&gain,1),om_vector_getValue(&gain,2));

	om_operator_quat_scal_mul(&dq,filter->_beta,&dq);
	om_operator_quat_mul(&filter->_q_pred,&dq,&filter->_q_est);
	om_quat_normalize(&filter->_q_est);

	// set output
	double qw = filter->_q_est._qw;
	double qx = filter->_q_est._qx;
	double qy = filter->_q_est._qy;
	double qz = filter->_q_est._qz;

	switch(manager->type){

	case Quarternion:
		om_quat_create(&manager->output.quaternion,qw,qx,qy,qz);
		break;

	case Matrix:
		om_convert_quaternion2matrix(&filter->_q_est,&manager->output.matrix);
		break;

	case EulerAngle:
		om_convert_quaternion2euler(&filter->_q_est,&manager->output.euler);
		break;
	case AxisAngle:
		om_convert_quaternion2axisAngle(&filter->_q_est,&manager->output.axis_angle);
		break;
	default:
		om_quat_create(&filter->_q_est,1.0,0.0,0.0,0.0);
		break;

	}


	// free memory
	om_matrix_free(&X);
	om_matrix_free(&K);
	om_matrix_free(&K_inv);
	om_matrix_free(&K_tmp);
	om_matrix_free(&X_t);
	om_matrix_free(&I);
	om_matrix_free(&S_acc);
	om_matrix_free(&S_mag);
	om_vector_free(&s_mag);
	om_vector_free(&s_acc);
	om_vector_free(&s);
	om_vector_free(&gain);

}



/* release all component used for the nonlinear filter GDOF */
void om_cfa_free(void *filter){

	//free memory
	om_vector_free(&(*(omNonLinearFilter_CFA*)filter)._v_acc_pred);
	om_vector_free(&(*(omNonLinearFilter_CFA*)filter)._v_mag_pred);

}






