/*
 * om.c
 *
 *  Created on: 25 May, 2016
 *      Author: Thomas BRAUD, Nizar OUARTI
 */

#include "om.h"

// global variables
omVector ned_gravity;
omVector ned_magnetic_field;
omVector ned_gravity_normalized;
omVector ned_magnetic_field_normalized;



// This function is done to set the initial orientation where "pseudonorth" is chosen as the initial heading orientation. It is done to suit to many Ground Truth dataset
//@ author Nizar Ouarti
void init_frame (double *data){
        
    omVector acc;
    omVector magn;
    omVector magnXacc;
    omVector pseudonorth;
    om_vector_create(&acc,3);
    om_vector_create(&magn,3);
    om_vector_create(&magnXacc,3);
    om_vector_create(&pseudonorth,3);

    //First sample t0
    om_vector_setValues(&acc,3,data[0], data[1], data[2]);
    om_vector_setValues(&magn,3,data[3], data[4], data[5]);
    om_vector_setValues(&magnXacc,3,0.0, 0.0, 0.0);
    om_vector_setValues(&pseudonorth,3,0.0, 0.0, 0.0);

    //normalization
    om_vector_normalize(&magn);
    om_vector_normalize(&acc);


    //compute the initial "pseudonorth"
    om_vector_crossProduct(&magn,&acc,&magnXacc);
    om_vector_normalize(&magnXacc);
    om_vector_crossProduct(&acc,&magnXacc,&pseudonorth);
    om_vector_normalize(&pseudonorth);


    /*
    printf("acc = ");om_vector_display(&acc);
    printf("mag = ");om_vector_display(&magn);
    printf("magnXacc = ");om_vector_display(&magnXacc);
    printf("pseudonorth = ");om_vector_display(&pseudonorth);
	*/

    //*/
    ned_magnetic_field_normalized._values[0]=(double)pseudonorth._values[0];
    ned_magnetic_field_normalized._values[1]=(double)pseudonorth._values[1];
    ned_magnetic_field_normalized._values[2]=(double)pseudonorth._values[2];
    /*/
    ned_gravity_normalized._values[0]=(double)pseudonorth._values[0];
    ned_gravity_normalized._values[1]=(double)pseudonorth._values[1];
    ned_gravity_normalized._values[2]=(double)pseudonorth._values[2];
    om_operator_vector_scal_mul(&ned_gravity_normalized, G, &ned_gravity);
    //*/

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
	(*(omNonLinearFilter_CGO*)filter)._k_acc = 2.0;
	(*(omNonLinearFilter_CGO*)filter)._k_mag = 2.0;

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
	om_rotate_vector_quaternion(&filter->_q_pred,&ned_gravity_normalized,&v_acc_pred);

	// compute v_mag = R(q_(k+1)^-)m^a
	om_rotate_vector_quaternion(&filter->_q_pred,&ned_magnetic_field_normalized,&v_mag_pred);

	// removing
	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	// normalization
	om_vector_normalize(&z_acc);
	om_vector_normalize(&z_mag);

	// compute cross product
	om_vector_crossProduct(&z_acc,&v_acc_pred,&cross_acc);
	om_vector_crossProduct(&z_mag,&v_mag_pred,&cross_mag);


	om_operator_vector_scal_mul(&cross_acc,(filter->_k_acc)/2.0,&cross_acc);
	om_operator_vector_scal_mul(&cross_mag,(filter->_k_mag)/2.0,&cross_mag);

	// compute gain omega
	om_operator_vector_add(&cross_acc,&cross_mag,&omega);

	// compute bias_est =  (_omega*(_k_I)*_delta_t*(-1.0));

	om_operator_vector_scal_mul(&omega,filter->_k_I*DELTA_T*(-1.0),&filter->_bias_est);
	//om_operator_vector_add(&filter->_bias_est,&filter->_bias_pred,&filter->_bias_est);

	// compute angular_velocity = y_gyro - b_est + omega*k_P
	omVector angular_velocity;
	om_vector_create(&angular_velocity,3);
	om_operator_vector_scal_mul(&omega,filter->_k_P,&omega_q);
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&filter->_bias_est,&angular_velocity);
	om_operator_vector_add(&angular_velocity,&omega_q,&angular_velocity);

	// compute q_(k+1)^+ = Omega(Y_gyro - b_pred + omega*k_P)q_(k)^+
	om_kinematics_quaternion(&filter->_q_est,&angular_velocity,&filter->_q_est);
	om_quat_normalize(&filter->_q_est);

	// set output
	double signQw = filter->_q_est._qw < 0.0? -1.0 : 1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
	om_vector_free(&angular_velocity);



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
	//_variance_u = 0.000031623;
	//_variance_v = 0.0031623;

	om_quat_create(&(*(omNonLinearFilter_USQUE*)filter)._q_pred,1.0,0.0,0.0,0.0);

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


	//prediction step
	om_usque_prediction(manager,(omNonLinearFilter_USQUE*)filter);

	//correction step
	om_usque_update(manager,(omNonLinearFilter_USQUE*)filter);

}


void om_usque_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter){

	int bool_start = 0;

	if(om_vector_getValue(&filter->_x_k, 0) == 0.0 )
		bool_start = 1;


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

	om_matrix_squareRoot(&T,&S,40);


	// generate sigma point at time k
	filter->sigma_points = (omVector*)malloc(13 * sizeof(omVector));
	om_vector_create(&filter->sigma_points[0],6);
	om_vector_clone(&filter->_x_k,&filter->sigma_points[0]);

	for(int i=1;i<=6;++i){

		omVector S_row;
		om_vector_create(&S_row,6);

		om_vector_create(&filter->sigma_points[i],6);
		om_vector_create(&filter->sigma_points[i+6],6);

		om_matrix_getRow(&S,(i-1),&S_row);
		om_operator_vector_add(&filter->_x_k,&S_row,&filter->sigma_points[i]);
		om_operator_vector_sub(&filter->_x_k,&S_row,&filter->sigma_points[i+6]);

		om_vector_free(&S_row);

	}


	// generate sigma quaternion at time k
	filter->sigma_quaternion = (omQuaternion*)malloc(13 * sizeof(omQuaternion));
	for(unsigned int i=0;i<13;++i){

		omVector dp_i;
		om_vector_create(&dp_i,3,om_vector_getValue(&filter->sigma_points[i],0),om_vector_getValue(&filter->sigma_points[i],1),om_vector_getValue(&filter->sigma_points[i],2));

		if(i==0){

			double qw = filter->_q_est._qw;
			double qx = filter->_q_est._qx;
			double qy = filter->_q_est._qy;
			double qz = filter->_q_est._qz;
			om_quat_create(&filter->sigma_quaternion[i],qw,qx,qy,qz);

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
			om_operator_quat_mul(&dq_i,&filter->_q_est,&filter->sigma_quaternion[i]);
			om_vector_free(&dqv_i);

		}

		om_vector_free(&dp_i);


	}


	// propagation of sigma quaternion and sigma points
	omQuaternion s_q_0_inv;
	for(int i=0;i<13;++i){

		omVector b_i;
		omVector angular_velocity_i;

		om_vector_create(&b_i,3,om_vector_getValue(&filter->sigma_points[i],3),om_vector_getValue(&filter->sigma_points[i],4),om_vector_getValue(&filter->sigma_points[i],5));
		om_vector_create(&angular_velocity_i,3);

		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);

		// propagation of sigma quaternion
		omQuaternion q_pred_i;
		om_kinematics_quaternion(&filter->sigma_quaternion[i],&angular_velocity_i,&q_pred_i);

		if(i==0){

			om_vector_setValues(&filter->sigma_points[0],6,0.0,0.0,0.0,om_vector_getValue(&b_i,0),om_vector_getValue(&b_i,1),om_vector_getValue(&b_i,2));
			om_quat_inverse(&q_pred_i,&s_q_0_inv);

		}else{

			omQuaternion dq_i;

			om_operator_quat_mul(&q_pred_i,&s_q_0_inv,&dq_i);

			double tmp = (f/ (a + dq_i._qw));
			double dp_i_x = dq_i._qx * tmp;
			double dp_i_y = dq_i._qy * tmp;
			double dp_i_z = dq_i._qz * tmp;

			om_vector_setValues(&filter->sigma_points[i],6,dp_i_x,dp_i_y,dp_i_z,
					om_vector_getValue(&b_i,0),om_vector_getValue(&b_i,1),om_vector_getValue(&b_i,2));

		}


		om_quat_create(&filter->sigma_quaternion[i],q_pred_i._qw,q_pred_i._qx,q_pred_i._qy,q_pred_i._qz);

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
		om_operator_vector_add(&sum_x,&filter->sigma_points[i],&sum_x);
	}

	om_operator_vector_scal_mul(&filter->sigma_points[0],lambda/(6.0+lambda),&sum_a);
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
		om_operator_vector_sub(&filter->sigma_points[i],&filter->_x_k_pred,&var_x_i);
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
	om_operator_vector_sub(&filter->sigma_points[0],&filter->_x_k_pred,&var_x_0);

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
	om_vector_free(&b_pred);
	om_vector_free(&omega);

}

void om_usque_update(struct omSensorFusionManager *manager,omNonLinearFilter_USQUE *filter){

	int bool_start = 0;

	if(om_vector_getValue(&filter->_x_k, 0) == 0.0 )
		bool_start = 1;

	//measurement vector;
	omVector z_acc;
	omVector z_mag;
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);

	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_mag);
	//om_vector_normalize(&z_acc);

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

		om_rotate_vector_quaternion(&filter->sigma_quaternion[i],&ned_gravity,&tmp_res_a);
		om_rotate_vector_quaternion(&filter->sigma_quaternion[i],&ned_magnetic_field_normalized,&tmp_res_m);
		//om_rotate_vector_quaternion(&filter->sigma_quaternion[i],&ned_magnetic_field,&tmp_res_m);


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
	om_operator_vector_sub(&filter->sigma_points[0],&filter->_x_k_pred,&var_x_0);
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
		om_operator_vector_sub(&filter->sigma_points[i],&filter->_x_k_pred,&var_x_i);
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
	om_quat_normalize(&filter->_q_est);

	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
	om_vector_free(&sum_z);
	om_vector_free(&sum_a);
	om_vector_free(&sum_b);
	om_vector_free(&z_acc);
	om_vector_free(&z_mag);

	om_matrix_free(&P_zz_0);
	om_matrix_free(&P_xz_0);
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
		om_vector_free(&filter->sigma_points[i]);
	}

	free(filter->sigma_points);
	filter->sigma_points = 0;
	free(filter->sigma_quaternion);
	filter->sigma_quaternion = 0;

	free(sigma_z);
	sigma_z = 0;

}



///////////////////////////////////////////////////////
/////           NonLinearFilter CSP               /////
///////////////////////////////////////////////////////

/* initialization of all component used by the nonlinear filter USQUE */
void om_csp_initialization(struct omSensorFusionManager *manager,void *filter){

	// constant
	(*(omNonLinearFilter_CSP*)filter)._seed = 0;
	(*(omNonLinearFilter_CSP*)filter)._a = 1.0;
	(*(omNonLinearFilter_CSP*)filter)._f = 4.0;
	(*(omNonLinearFilter_CSP*)filter)._lambda = 0.5;
	(*(omNonLinearFilter_CSP*)filter)._k_I = 0.3;
	(*(omNonLinearFilter_CSP*)filter)._k_P = 1.0;
	(*(omNonLinearFilter_CSP*)filter)._k_acc = 1.0;
	(*(omNonLinearFilter_CSP*)filter)._k_mag = 1.0;

	om_quat_create(&(*(omNonLinearFilter_CSP*)filter)._q_pred,1.0,0.0,0.0,0.0);

	switch(manager->type){

	case Quarternion:

		om_quat_create(&(*(omNonLinearFilter_CSP*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_CSP*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_CSP*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_CSP*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}


	double var_u = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double var_v = manager->imu_params.variance_gyroscope;
	//var_u = 0.000031623;
	//var_v = 0.0031623;




	double delta_t = DELTA_T;
	om_matrix_create(&(*(omNonLinearFilter_CSP*)filter)._Q,6,6);
	for(int i=0;i<3;++i){

		om_matrix_setValue(&(*(omNonLinearFilter_CSP*)filter)._Q,i,i,(var_v*delta_t + 0.33333*var_u*(pow(delta_t,3.0))));
		om_matrix_setValue(&(*(omNonLinearFilter_CSP*)filter)._Q,i+3,i,-(0.5*var_u*delta_t*delta_t));
		om_matrix_setValue(&(*(omNonLinearFilter_CSP*)filter)._Q,i,i+3,-(0.5*var_u*delta_t*delta_t));
		om_matrix_setValue(&(*(omNonLinearFilter_CSP*)filter)._Q,i+3,i+3,(var_u*delta_t));

	}

	om_matrix_create(&(*(omNonLinearFilter_CSP*)filter)._R,6,6);

	for(int i=0;i<3;++i){
		om_matrix_setValue(&(*(omNonLinearFilter_CSP*)filter)._R,i,i,manager->imu_params.variance_accelerometer);
		om_matrix_setValue(&(*(omNonLinearFilter_CSP*)filter)._R,i+3,i+3,manager->imu_params.variance_magnetometer);
	}


	double b_x = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double b_y = om_vector_getValue(&manager->imu_params.bias_gyroscope,1);
	double b_z = om_vector_getValue(&manager->imu_params.bias_gyroscope,2);

	om_matrix_create(&(*(omNonLinearFilter_CSP*)filter)._P_k,6,6);
	om_matrix_create(&(*(omNonLinearFilter_CSP*)filter)._P_k_pred,6,6);
	om_vector_create(&(*(omNonLinearFilter_CSP*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);
	om_vector_create(&(*(omNonLinearFilter_CSP*)filter)._x_k_pred,6);




}


void om_csp_free(void *filter){


	om_vector_free(&(*(omNonLinearFilter_CSP*)filter)._x_k);
	om_vector_free(&(*(omNonLinearFilter_CSP*)filter)._x_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_CSP*)filter)._P_k_pred);
	om_matrix_free(&(*(omNonLinearFilter_CSP*)filter)._P_k);
	om_matrix_free(&(*(omNonLinearFilter_CSP*)filter)._Q);
	om_matrix_free(&(*(omNonLinearFilter_CSP*)filter)._R);

}



void om_csp_process(struct omSensorFusionManager *manager,void *filter){

	om_csp_prediction(manager,(omNonLinearFilter_CSP*)filter);
	om_csp_update(manager,(omNonLinearFilter_CSP*)filter);
}




void om_csp_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CSP *filter){


	double a = filter->_a;
	double f = filter->_f;
	double lambda = filter->_lambda;

	omMatrix S;
	omMatrix T;

	om_matrix_create(&S,6,6);
	om_matrix_create(&T,6,6);

	// compute S = sqrt( (P + Q)*(n + lambda )
	om_operator_matrix_add(&filter->_P_k,&filter->_Q,&T);
	//om_matrix_clone(&filter->_P_k,&T);
	om_operator_matrix_scal_mul(&T,lambda + 6.0,&T);
	//om_matrix_choleskyDecomposition(&T,&S);

	om_matrix_squareRoot(&T,&S,40);

	// generate sigma point at time k
	filter->sigma_points = (omVector*)malloc(13 * sizeof(omVector));
	om_vector_create(&filter->sigma_points[0],6);
	om_vector_clone(&filter->_x_k,&filter->sigma_points[0]);

	for(int i=1;i<=6;++i){

		omVector S_row;
		om_vector_create(&S_row,6);

		om_vector_create(&filter->sigma_points[i],6);
		om_vector_create(&filter->sigma_points[i+6],6);

		om_matrix_getRow(&S,(i-1),&S_row);
		om_operator_vector_add(&filter->_x_k,&S_row,&filter->sigma_points[i]);
		om_operator_vector_sub(&filter->_x_k,&S_row,&filter->sigma_points[i+6]);

		om_vector_free(&S_row);

	}


	// generate sigma quaternion at time k
	filter->sigma_quaternion = (omQuaternion*)malloc(13 * sizeof(omQuaternion));
	for(unsigned int i=0;i<13;++i){

		omVector dp_i;
		om_vector_create(&dp_i,3,om_vector_getValue(&filter->sigma_points[i],0),om_vector_getValue(&filter->sigma_points[i],1),om_vector_getValue(&filter->sigma_points[i],2));

		if(i==0){

			double qw = filter->_q_est._qw;
			double qx = filter->_q_est._qx;
			double qy = filter->_q_est._qy;
			double qz = filter->_q_est._qz;
			om_quat_create(&filter->sigma_quaternion[i],qw,qx,qy,qz);

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
			om_operator_quat_mul(&dq_i,&filter->_q_est,&filter->sigma_quaternion[i]);
			om_vector_free(&dqv_i);

		}

		om_vector_free(&dp_i);


	}


	// variable
	omVector z_acc;
	omVector z_mag;
	omVector v_acc_pred_i;
	omVector v_mag_pred_i;
	omVector cross_acc_i;
	omVector cross_mag_i;
	omVector omega_i;
	omVector omega_q_i;

	// allocation
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);
	om_vector_create(&v_acc_pred_i,3);
	om_vector_create(&v_mag_pred_i,3);
	om_vector_create(&omega_i,3);
	om_vector_create(&omega_q_i,3);
	om_vector_create(&cross_acc_i,3);
	om_vector_create(&cross_mag_i,3);

	// removing
	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	// normalization
	om_vector_normalize(&z_acc);
	om_vector_normalize(&z_mag);


	// propagation of sigma quaternion and sigma points
	omQuaternion s_q_0_inv;
	for(int i=0;i<13;++i){

		omVector b_i;
		omVector b_i_tmp;
		omVector angular_velocity_i;

		om_vector_create(&b_i,3,om_vector_getValue(&filter->sigma_points[i],3),om_vector_getValue(&filter->sigma_points[i],4),om_vector_getValue(&filter->sigma_points[i],5));
		om_vector_create(&b_i_tmp,3);
		om_vector_create(&angular_velocity_i,3);
		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);

		// propagation of sigma quaternion
		omQuaternion q_pred_i;
		om_kinematics_quaternion(&filter->sigma_quaternion[i],&angular_velocity_i,&q_pred_i);

		/////
		// constraint computation
		/////

		// compute v_acc = R(q_(k+1)^-)g^a
		om_rotate_vector_quaternion(&q_pred_i,&ned_gravity_normalized,&v_acc_pred_i);

		// compute v_mag = R(q_(k+1)^-)m^a
		om_rotate_vector_quaternion(&q_pred_i,&ned_magnetic_field_normalized,&v_mag_pred_i);

		om_vector_normalize(&v_acc_pred_i);
		om_vector_normalize(&v_mag_pred_i);


		// compute cross product
		om_vector_crossProduct(&z_acc,&v_acc_pred_i,&cross_acc_i);
		om_vector_crossProduct(&z_mag,&v_mag_pred_i,&cross_mag_i);


		om_operator_vector_scal_mul(&cross_acc_i,(filter->_k_acc)/2.0,&cross_acc_i);
		om_operator_vector_scal_mul(&cross_mag_i,(filter->_k_mag)/2.0,&cross_mag_i);

		// compute gain omega
		om_operator_vector_add(&cross_acc_i,&cross_mag_i,&omega_i);

		// compute bias_est =  (_omega*(_k_I)*_delta_t*(-1.0));
		//printf("b_i       = ");om_vector_display(&b_i);
		om_operator_vector_scal_mul(&omega_i,filter->_k_I*DELTA_T*(-1.0),&b_i);
		//double toto = filter->_k_I*DELTA_T*(-1.0);
		//om_operator_vector_scal_mul(&omega_i,toto,&b_i_tmp);
		//om_operator_vector_add(&b_i,&b_i_tmp,&b_i);


		//printf("\nfilter->_k_I = %f\n",filter->_k_I);
		//printf("DELTA_T = %f\n",DELTA_T);
		//printf("toto = %f\n",toto);
		//printf("omega_i = ");om_vector_display(&omega_i);
		//printf("after b_i = ");om_vector_display(&b_i_tmp);


		// compute angular_velocity = y_gyro - b_est + omega*k_P
		om_operator_vector_scal_mul(&omega_i,filter->_k_P,&omega_q_i);
		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);
		om_operator_vector_add(&angular_velocity_i,&omega_q_i,&angular_velocity_i);

		// compute q_(k+1)^+ = Omega(Y_gyro - b_pred + omega*k_P)q_(k)^+
		om_kinematics_quaternion(&filter->sigma_quaternion[i],&angular_velocity_i,&q_pred_i);

		if(i==0){

			om_vector_setValues(&filter->sigma_points[0],6,0.0,0.0,0.0,om_vector_getValue(&b_i,0),om_vector_getValue(&b_i,1),om_vector_getValue(&b_i,2));
			om_quat_inverse(&q_pred_i,&s_q_0_inv);

		}else{


			omQuaternion dq_i;

			om_operator_quat_mul(&q_pred_i,&s_q_0_inv,&dq_i);

			double tmp = (f/ (a + dq_i._qw));
			double dp_i_x = dq_i._qx * tmp;
			double dp_i_y = dq_i._qy * tmp;
			double dp_i_z = dq_i._qz * tmp;

			om_vector_setValues(&filter->sigma_points[i],6,dp_i_x,dp_i_y,dp_i_z,om_vector_getValue(&b_i,0),om_vector_getValue(&b_i,1),om_vector_getValue(&b_i,2));

		}

		om_quat_create(&filter->sigma_quaternion[i],q_pred_i._qw,q_pred_i._qx,q_pred_i._qy,q_pred_i._qz);

		om_vector_free(&b_i);
		om_vector_free(&b_i_tmp);
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
		om_operator_vector_add(&sum_x,&filter->sigma_points[i],&sum_x);
	}

	om_operator_vector_scal_mul(&filter->sigma_points[0],lambda/(6.0+lambda),&sum_a);
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
		om_operator_vector_sub(&filter->sigma_points[i],&filter->_x_k_pred,&var_x_i);
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
	om_operator_vector_sub(&filter->sigma_points[0],&filter->_x_k_pred,&var_x_0);

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

	/*/
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_pred,&omega);
	om_kinematics_quaternion(&filter->_q_est,&omega,&filter->_q_pred);
	/*/
	om_kinematics_quaternion(&filter->_q_est,&manager->imu_data.data_gyroscope,&filter->_q_pred);
	//*/


	//free memory
	om_vector_free(&b_pred);
	om_vector_free(&omega);
	om_vector_free(&z_acc);
	om_vector_free(&z_mag);
	om_vector_free(&cross_acc_i);
	om_vector_free(&cross_mag_i);
	om_vector_free(&omega_i);
	om_vector_free(&omega_q_i);
	om_vector_free(&v_mag_pred_i);
	om_vector_free(&v_acc_pred_i);

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




}

void om_csp_update(struct omSensorFusionManager *manager,omNonLinearFilter_CSP *filter){



	//measurement vector;
	omVector z_acc;
	omVector z_mag;
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);

	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_mag);
	//om_vector_normalize(&z_acc);

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

		om_rotate_vector_quaternion(&filter->sigma_quaternion[i],&ned_gravity,&tmp_res_a);
		om_rotate_vector_quaternion(&filter->sigma_quaternion[i],&ned_magnetic_field_normalized,&tmp_res_m);

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
	om_operator_vector_sub(&filter->sigma_points[0],&filter->_x_k_pred,&var_x_0);
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
		om_operator_vector_sub(&filter->sigma_points[i],&filter->_x_k_pred,&var_x_i);
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
	om_quat_normalize(&filter->_q_est);

	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}


	//free memory
	om_vector_free(&z_acc);
	om_vector_free(&z_mag);
	om_vector_free(&z_k);
	om_vector_free(&dp_est);
	om_vector_free(&mean_z);
	om_vector_free(&s_z_k);
	om_vector_free(&corr_x);
	om_vector_free(&var_x_0);
	om_vector_free(&var_z_0);

	om_vector_free(&sum_a);
	om_vector_free(&sum_b);
	om_vector_free(&sum_z);


	om_matrix_free(&P_zz_0);
	om_matrix_free(&P_xz_0);
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
		om_vector_free(&filter->sigma_points[i]);
	}

	free(filter->sigma_points);
	filter->sigma_points = 0;
	free(filter->sigma_quaternion);
	filter->sigma_quaternion = 0;


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

	om_rotate_vector_quaternion(&q_var,&ned_magnetic_field_normalized,&tmp_res_m);
	om_rotate_vector_quaternion(&q_var,&ned_gravity_normalized,&tmp_res_a);
	om_operator_vector_scal_mul(&tmp_res_a,9.81,&tmp_res_a);


	for(int i=0;i<3;i++){
		om_vector_setValue(h_x,i,om_vector_getValue(&tmp_res_a,i)); //+ om_vector_getValue(&filter->_v_k,i));
		om_vector_setValue(h_x,i+3,om_vector_getValue(&tmp_res_m,i)); // + om_vector_getValue(&filter->_v_k,i+3) );
	}

	om_vector_free(&dp);
	om_vector_free(&tmp_res_a);
	om_vector_free(&tmp_res_m);


}

void om_mekf_f_jacobian(struct omSensorFusionManager *manager,omNonLinearFilter_MEKF *filter){

	double h = DELTA_T;

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


	double h = DELTA_T;

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


	//om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&filter->_Q_cho,filter->_seed,&filter->_w_k);
	//filter->_seed += 6.0;
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

	//measurement vector;
	omVector z_acc;
	omVector z_mag;
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);

	om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_mag);
	//om_vector_normalize(&z_acc);

	omVector z;
	om_vector_create(&z,6,om_vector_getValue(&z_acc,0),om_vector_getValue(&z_acc,1),om_vector_getValue(&z_acc,2),
							om_vector_getValue(&z_mag,0),om_vector_getValue(&z_mag,1),om_vector_getValue(&z_mag,2));


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
	om_quat_normalize(&filter->_q_est);

	// reset error dp to null vector
	for(int i=0;i<3;++i)
		om_vector_setValue(&filter->_x_k,i,0.0);



	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
/////           NonLinearFilter MEKF              /////
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

	(*(omNonLinearFilter_REQUEST*)filter)._n = 2;
	(*(omNonLinearFilter_REQUEST*)filter)._index = 0;
	(*(omNonLinearFilter_REQUEST*)filter)._bool_start = 0;

	// init weights
	(*(omNonLinearFilter_REQUEST*)filter)._a = (double*)malloc((*(omNonLinearFilter_REQUEST*)filter)._n*sizeof(double));

	// init measuerment vector & reference vector
	(*(omNonLinearFilter_REQUEST*)filter)._b = (omVector*)malloc((*(omNonLinearFilter_REQUEST*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_REQUEST*)filter)._r = (omVector*)malloc((*(omNonLinearFilter_REQUEST*)filter)._n*sizeof(omVector));


	for(int i=0;i<(*(omNonLinearFilter_REQUEST*)filter)._n;i+=2){
		(*(omNonLinearFilter_REQUEST*)filter)._a[i] = 1.0 / manager->imu_params.variance_accelerometer;
		(*(omNonLinearFilter_REQUEST*)filter)._a[i+1] = 1.0 / manager->imu_params.variance_magnetometer;

		om_vector_create(&(*(omNonLinearFilter_REQUEST*)filter)._r[i],3);
		om_vector_create(&(*(omNonLinearFilter_REQUEST*)filter)._r[i+1],3);

		om_vector_clone(&ned_gravity_normalized, &(*(omNonLinearFilter_REQUEST*)filter)._r[i]);
		om_vector_clone(&ned_magnetic_field_normalized, &(*(omNonLinearFilter_REQUEST*)filter)._r[i+1]);

	}

	//normalize weight
	double sum = 0.0;
	for(int i=0;i<(*(omNonLinearFilter_REQUEST*)filter)._n;i++){
		//printf("a[i] = %f\n",(*(omNonLinearFilter_QUEST*)filter)._a[i]);
		sum +=  (*(omNonLinearFilter_REQUEST*)filter)._a[i];
	}


	for(int i=0;i<(*(omNonLinearFilter_REQUEST*)filter)._n;i++){

		(*(omNonLinearFilter_REQUEST*)filter)._a[i] /= sum;
		om_vector_create(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],3);


		//printf("after a[i] = %f\n",(*(omNonLinearFilter_QUEST*)filter)._a[i]);
	}
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._K_k,4,4);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._K_k_pred,4,4);

}

void om_request_free(void *filter){


	for(int i=0;i<(*(omNonLinearFilter_REQUEST*)filter)._n;i++){
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



}



void om_request_process(struct omSensorFusionManager *manager,void *filter){

	om_request_preprocess(manager,(omNonLinearFilter_REQUEST*)filter);

	if((*(omNonLinearFilter_REQUEST*)filter)._bool_start){
		om_request_prediction(manager,(omNonLinearFilter_REQUEST*)filter);
		om_request_update(manager,(omNonLinearFilter_REQUEST*)filter);
	}
}

void om_request_preprocess(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){

	/* update measurement vector */
	omVector z_acc;
	omVector z_mag;

	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);

	om_vector_clone(&manager->imu_data.data_accelerometer, &z_acc);
	om_vector_clone(&manager->imu_data.data_magnetometer, &z_mag);
	//om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	//om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_acc);
	om_vector_normalize(&z_mag);

	om_vector_setValues(&filter->_b[filter->_index*2],3,om_vector_getValue(&z_acc,0),om_vector_getValue(&z_acc,1),om_vector_getValue(&z_acc,2));
	om_vector_setValues(&filter->_b[filter->_index*2 + 1],3,om_vector_getValue(&z_mag,0),om_vector_getValue(&z_mag,1),om_vector_getValue(&z_mag,2));

	filter->_index = (filter->_index + 1) % (filter->_n/2);
	filter->_bool_start = ( filter->_bool_start == 1 || (filter->_index == (filter->_n/2 - 1)) )? 1 : 0;

	om_vector_free(&z_acc);
	om_vector_free(&z_mag);

	//printf("index = %d\n",filter->_index);

}

void om_request_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){


	omMatrix B;
	omMatrix B_t;
	omVector z;
	omMatrix S;
	omMatrix S_tmp;
	omMatrix d_K;
	omMatrix I;

	om_matrix_create(&B,3,3);
	om_matrix_create(&B_t,3,3);
	om_matrix_create(&S,3,3);
	om_matrix_create(&S_tmp,3,3);
	om_matrix_create(&d_K,4,4);
	om_matrix_createIdentity(&I,3);
	om_vector_create(&z,3,0.0,0.0,0.0);


	for(int i=0;i<filter->_n;i++){

		/*/
		printf("b[i] = ");om_vector_display(&filter->_b[i]);
		printf("r[i] = ");om_vector_display(&filter->_r[i]);
		//printf("a[i] = %f\n",(*(omNonLinearFilter_QUEST*)filter)._a[i]);
		//*/

		omVector cross;
		omMatrix m_bi;
		omMatrix m_ri_t;
		omMatrix B_tmp;

		om_vector_create(&cross,3);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&B_tmp,3,3);

		// compute _d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		om_vector_crossProduct(&filter->_b[i],&filter->_r[i],&cross);
		om_operator_vector_scal_mul(&cross,filter->_a[i],&cross);
		om_operator_vector_add(&z,&cross,&z);


		//compute B = B + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&filter->_b[i],j));
			om_matrix_setValue(&m_ri_t,0,j,om_vector_getValue(&filter->_r[i],j));
		}

		om_operator_matrix_mul(&m_bi,&m_ri_t,&B_tmp);
		om_operator_matrix_scal_mul(&B_tmp,filter->_a[i],&B_tmp);
		om_operator_matrix_add(&B,&B_tmp,&B);

		//free memory
		om_vector_free(&cross);
		om_matrix_free(&m_bi);
		om_matrix_free(&m_ri_t);
		om_matrix_free(&B_tmp);

	}

	//printf("\nMatrix B\n");
	//om_matrix_display(&B);

	// get transpose B
	om_matrix_transpose(&B,&B_t);

	// compute S = B + B_t
	om_operator_matrix_add(&B,&B_t,&S);

	// compute S - I*tr(B)
	double trB = om_matrix_trace(&B);
	om_operator_matrix_scal_mul(&I,trB,&I);
	om_operator_matrix_sub(&S,&I,&S_tmp);


	// compute matrix K
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
			om_matrix_setValue(&d_K,i,j,om_matrix_getValue(&S_tmp,i,j));

		om_matrix_setValue(&d_K,i,3,om_vector_getValue(&z,i));
		om_matrix_setValue(&d_K,3,i,om_vector_getValue(&z,i));
	}
	om_matrix_setValue(&d_K,3,3,trB);

	// compute state transition matrix Phy
	omMatrix Phi;
	omMatrix Phi_t;
	omMatrix tmp_K;
	omVector angular_velocity;

	om_matrix_create(&tmp_K,4,4);
	om_matrix_create(&Phi_t,4,4);
	om_matrix_create(&Phi,4,4);

	om_vector_create(&angular_velocity,3);
	//om_operator_vector_sub(&manager->imu_data.data_gyroscope,&manager->imu_params.bias_gyroscope,&angular_velocity);
	om_vector_clone(&manager->imu_data.data_gyroscope,&angular_velocity);


	om_operator_omega_kinematics(&angular_velocity,&Phi);
	om_matrix_transpose(&Phi,&Phi_t);

	// propagate K_k
	om_operator_matrix_mul(&Phi,&filter->_K_k,&tmp_K);
	om_operator_matrix_mul(&tmp_K,&Phi_t,&filter->_K_k_pred);


	// propagate m_k
	double rho_k_opt = 0.1;
	double d_m_kp1 = 0.0;

	for(int i=0;i<2;i++){
		d_m_kp1 += filter->_a[filter->_index*2 + i];
	}

	double m_kp1 = (filter->_m_k) + (d_m_kp1);

	// some variable
	double tmp_a = (rho_k_opt)*(filter->_m_k/m_kp1);
	double tmp_b = (1.0/m_kp1);

	filter->_m_k = m_kp1;

	// update K
	om_operator_matrix_scal_mul(&filter->_K_k_pred,tmp_a,&filter->_K_k_pred);
	om_operator_matrix_scal_mul(&d_K,tmp_b,&d_K);
	om_operator_matrix_add(&filter->_K_k_pred,&d_K,&filter->_K_k);


	//printf("\nMatrix K\n");
	//om_matrix_display(&filter->_K_k);
	//printf("\nd_m_kp1 = %.*f\n",10,d_m_kp1);

	// free memory
	om_matrix_free(&Phi);
	om_matrix_free(&Phi_t);
	om_matrix_free(&tmp_K);
	om_vector_free(&angular_velocity);



}

void om_request_update(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){

	// get eigen values of K
	omVector* eigen_vector;
	double* eigen_values;
	om_matrix_getEingenValues(&filter->_K_k,&eigen_vector,&eigen_values,25);





	// calcul of lambda_max
	double lambda=0.0;
	for(int i=0;i<4;i++){

		//printf("l_i = %.*f v_i = ",10,eigen_values[i]);om_vector_display(&eigen_vector[i]);

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


	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}



	// free memory
	free(eigen_values);
	free(eigen_vector);
	eigen_values = 0;
	eigen_vector = 0;



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
	(*(omNonLinearFilter_PF*)filter)._n = 2000;
	(*(omNonLinearFilter_PF*)filter)._threeshold = ( (double)(*(omNonLinearFilter_PF*)filter)._n)*(1.0/3.0) ;

	// init state
	double b_x = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double b_y = om_vector_getValue(&manager->imu_params.bias_gyroscope,1);
	double b_z = om_vector_getValue(&manager->imu_params.bias_gyroscope,2);

	om_vector_create(&(*(omNonLinearFilter_PF*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);

	//compute Cholesky Decomposition of matrix Q in order to generate noise

	//double var_u = 0.000031623;
	//double var_v = 0.0031623;

	double var_u = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double var_v = manager->imu_params.variance_gyroscope;

	om_matrix_create(&(*(omNonLinearFilter_PF*)filter)._cov_L,6,6);
	double rate = 0.1*DELTA_T;
	for(int i=0;i<3;i++){
		double tmp = ((var_v*var_v)/rate) + ((var_u*var_u*rate)/12.0);

		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._cov_L,i,i,sqrt(tmp));
		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._cov_L,i+3,i+3, var_u*sqrt(rate));
		om_matrix_setValue(&(*(omNonLinearFilter_PF*)filter)._cov_L,i+3,i,(-0.5)*var_u*sqrt(rate));
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
		om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&(*(omNonLinearFilter_PF*)filter)._cov_L,(*(omNonLinearFilter_PF*)filter)._seed,&(*(omNonLinearFilter_PF*)filter)._particle_wn[i]);
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


}


void om_pf_free(void *filter){

	om_vector_free(&(*(omNonLinearFilter_PF*)filter)._x_k);
	om_matrix_free(&(*(omNonLinearFilter_PF*)filter)._cov_L);
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
		//om_kinematics_quaternion(&filter->_q_est,&angular_velocity_i,&filter->_particle_q[i]);

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

		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_gravity_normalized,&z_acc);
		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_magnetic_field_normalized,&z_mag);


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
			tmp_acc += om_vector_getValue(&s_k_i,l)*om_vector_getValue(&s_k_i,l)*(1.0/manager->imu_params.variance_accelerometer);
			tmp_mag += om_vector_getValue(&s_k_i,l+3)*om_vector_getValue(&s_k_i,l+3)*(1.0/manager->imu_params.variance_magnetometer);
		}

		//tmp_acc *= (1.0/manager->imu_params.variance_accelerometer);
		//tmp_mag *= (1.0/manager->imu_params.variance_magnetometer);

		double L_k_i = exp( (tmp_acc+tmp_mag)*(-0.5) );

		// update weights
		filter->_particle_w[i] *= L_k_i;

		//printf("L_k_i = %.*f\n",30,L_k_i);
		//if(i == 0){
		//	printf("tmp_acc = %.*f\n",30,tmp_acc);
		//	printf("tmp_mag = %.*f\n",30,tmp_mag);
		//}

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

		/*/
		om_vector_free(&v_acc_pred_i);
		om_vector_free(&v_mag_pred_i);
		om_vector_free(&cross_acc_i);
		om_vector_free(&cross_mag_i);
		om_vector_free(&omega_i);
		om_vector_free(&omega_q_i);
		om_vector_free(&b_i_tmp);
		//*/



	}
	om_vector_free(&z_k);
	om_vector_free(&z_acc_tmp);
	om_vector_free(&z_mag_tmp);


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
		//if(isnan(filter->_particle_w[i]))
		//	filter->_particle_w[i] = 1.0/(double)(filter->_n);

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
	om_quat_normalize(&filter->_q_est);


	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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

		//*/
		//om_enkpf_quicksort(filter,0,filter->_n-1);

		double* c = (double*)malloc(filter->_n*sizeof(double));
		c[0] = 0;

		for (int i = 1; i < (filter->_n); ++i) {
			c[i] = c[i-1] + filter->_particle_w[i];
		}

		double u1 = om_random_uniformDistribution(filter->_seed)*(1.0/filter->_n);
		filter->_seed+=1.0;

		int index = 0;
		for (int j = 1; j <= (filter->_n); ++j) {

			double uj = u1 + (j - 1.0)/(filter->_n);
			while(index < (filter->_n) && uj > c[index])
				index++;

			if(index < (filter->_n) ){
				//om_enkpf_swap(filter, index, j);
				om_vector_clone(&filter->_particle_x[index],&filter->_particle_x[j-1]);
			}
			filter->_particle_w[j-1] = (1.0/filter->_n);
		}
		/*/

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
		//*/
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



/* initialization of all component used by the nonlinear filter GDOF */
void om_gdof_initialization(struct omSensorFusionManager *manager,void *filter){

	// constant initialization
	(*(omNonLinearFilter_GDOF*)filter)._eta = 0.003;
	(*(omNonLinearFilter_GDOF*)filter)._beta = 0.08;

	//(*(omNonLinearFilter_GDOF*)filter)._eta = 0.0001;
	//(*(omNonLinearFilter_GDOF*)filter)._beta = 0.002;



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


}

/* process function of the nonlinear filter GDOF */
void om_gdof_process(struct omSensorFusionManager *manager,void *filter){

	om_gdof_prediction(manager,(omNonLinearFilter_GDOF*)filter);
	om_gdof_update(manager,(omNonLinearFilter_GDOF*)filter);

}





void om_gdof_f_b_function(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter,omVector *x,omVector *f_x){

	omVector mag;
	omVector v_mag;
	omQuaternion q;

	om_vector_create(&mag,3);
	om_vector_create(&v_mag,3);
	om_quat_create(&q,om_vector_getValue(x,0),om_vector_getValue(x,1),om_vector_getValue(x,2),om_vector_getValue(x,3));


	om_vector_clone(&manager->imu_data.data_magnetometer,&mag);
	om_rotate_vector_quaternion(&q,&filter->_e_b,&v_mag);

	om_vector_normalize(&v_mag);
	om_vector_normalize(&mag);

	double f_b_x = om_vector_getValue(&v_mag,0) - om_vector_getValue(&mag,0);
	double f_b_y = om_vector_getValue(&v_mag,1) - om_vector_getValue(&mag,1);
	double f_b_z = om_vector_getValue(&v_mag,2) - om_vector_getValue(&mag,2);

	om_vector_setValues(f_x,f_b_x,f_b_y,f_b_z);

	om_vector_free(&mag);
	om_vector_free(&v_mag);

}

void om_gdof_J_b_jacobian(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter){

	double h = 0.001;

	for(int i=0;i<3;i++)
		for(int j=0;j<4;j++)
			om_matrix_setValue(&filter->_J_b,i,j,0.0);


	omVector x;
	om_vector_create(&x,4,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);


	for (int j = 0; j < 4; j++){

		omVector xhp;
		omVector xhm;
		omVector f_xhm;
		omVector f_xhp;
		omVector value;

		om_vector_create(&value,3);
		om_vector_create(&xhp,4);
		om_vector_create(&xhm,4);
		om_vector_create(&f_xhm,3);
		om_vector_create(&f_xhp,3);


		om_vector_clone(&x,&xhp);
		om_vector_clone(&x,&xhm);

		om_vector_setValue(&xhp,j,om_vector_getValue(&xhp,j)+h);
		om_vector_setValue(&xhm,j,om_vector_getValue(&xhm,j)-h);

		om_gdof_f_b_function(manager,filter,&xhp,&f_xhp);
		om_gdof_f_b_function(manager,filter,&xhm,&f_xhm);

		om_operator_vector_sub(&f_xhp,&f_xhm,&value);
		om_operator_vector_scal_div(&value,2.0*h,&value);

		for (int i = 0; i < 3; i++)
			om_matrix_setValue(&filter->_J_b,i,j,om_vector_getValue(&value,i));

		om_vector_free(&xhp);
		om_vector_free(&xhm);
		om_vector_free(&f_xhp);
		om_vector_free(&f_xhm);
		om_vector_free(&value);


	}

}


void om_gdof_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_GDOF *filter){


	omVector acc;
	omVector mag;

	omVector v_acc;
	omVector v_mag;

	om_vector_create(&acc,3);
	om_vector_create(&mag,3);
	om_vector_create(&v_acc,3);
	om_vector_create(&v_mag,3);


	om_vector_clone(&manager->imu_data.data_accelerometer,&acc);
	om_vector_clone(&manager->imu_data.data_magnetometer,&mag);

	om_vector_normalize(&acc);
	om_vector_normalize(&mag);


	//*/
	omQuaternion q_inv;
	om_quat_inverse(&filter->_q_est,&q_inv);
	om_rotate_vector_quaternion(&q_inv,&mag,&filter->_e_m);

	double b_x = sqrt( (om_vector_getValue(&filter->_e_m,0)*om_vector_getValue(&filter->_e_m,0) ) + (om_vector_getValue(&filter->_e_m,1)*om_vector_getValue(&filter->_e_m,1) ) );
	double b_z = om_vector_getValue(&filter->_e_m,2);
	/*/
	double b_x = sqrt( (om_vector_getValue(&ned_magnetic_field_normalized,0)*om_vector_getValue(&ned_magnetic_field_normalized,0) ) + (om_vector_getValue(&ned_magnetic_field_normalized,1)*om_vector_getValue(&ned_magnetic_field_normalized,1) ));
	double b_z = om_vector_getValue(&ned_magnetic_field_normalized,2);
	//*/

	//om_vector_clone(&ned_magnetic_field_normalized,&filter->_e_b);
	om_vector_setValues(&filter->_e_b,3,b_x,0.0,b_z);

	double q1 = filter->_q_est._qw;
	double q2 = filter->_q_est._qx;
	double q3 = filter->_q_est._qy;
	double q4 = filter->_q_est._qz;

	//*/
	double f_a_x = (2.0*( (q2*q4) - (q1*q3) ))       - om_vector_getValue(&acc,0);
	double f_a_y = (2.0*( (q1*q2) + (q4*q3) ))       - om_vector_getValue(&acc,1);
	double f_a_z = (2.0*( 0.5 - (q2*q2) - (q3*q3) )) - om_vector_getValue(&acc,2);

	double f_b_x = (2.0*b_x*(0.5 - (q3*q3)  -(q4*q4)) ) + (2.0*b_z*( (q2*q4) - (q1*q3) )) - om_vector_getValue(&mag,0);
	double f_b_y = (2.0*b_x*( (q2*q3) - (q1*q4) )) + (2.0*b_z*( (q1*q2) + (q4*q3) ))      - om_vector_getValue(&mag,1);
	double f_b_z = (2.0*b_x*( (q1*q3) + (q2*q4))) + (2.0*b_z*(0.5 - (q3*q3) -(q2*q2)))    - om_vector_getValue(&mag,2);

	/*/
	om_rotate_vector_quaternion(&filter->_q_est,&ned_gravity_normalized,&v_acc);
	om_rotate_vector_quaternion(&filter->_q_est,&filter->_e_b,&v_mag);


	om_vector_normalize(&v_acc);
	om_vector_normalize(&v_mag);

	double f_a_x = om_vector_getValue(&v_acc,0) - om_vector_getValue(&acc,0);
	double f_a_y = om_vector_getValue(&v_acc,1) - om_vector_getValue(&acc,1);
	double f_a_z = om_vector_getValue(&v_acc,2) - om_vector_getValue(&acc,2);

	double f_b_x = om_vector_getValue(&v_mag,0) - om_vector_getValue(&mag,0);
	double f_b_y = om_vector_getValue(&v_mag,1) - om_vector_getValue(&mag,1);
	double f_b_z = om_vector_getValue(&v_mag,2) - om_vector_getValue(&mag,2);

	//*/

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

	//*/
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
	/*/
	//om_vector_setValues(&filter->_e_b,3,0.0,b_x,b_z);
	om_gdof_J_b_jacobian(manager,filter);
	//*/

	om_vector_free(&acc);
	om_vector_free(&mag);
	om_vector_free(&v_acc);
	om_vector_free(&v_mag);


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
	omVector delta_bias;

	// allocation
	om_vector_create(&gain_a,4);
	om_vector_create(&gain_b,4);
	om_vector_create(&gain,4);
	om_vector_create(&delta_bias,3);
	om_matrix_create(&J_a_t,4,3);
	om_matrix_create(&J_b_t,4,3);
	om_vector_create(&angular_velocity,3);


	// compute gain

	om_matrix_transpose(&filter->_J_a,&J_a_t);
	om_matrix_transpose(&filter->_J_b,&J_b_t);

	om_operator_matrix_vector_mul(&J_a_t,&filter->_f_a,&gain_a);
	om_operator_matrix_vector_mul(&J_b_t,&filter->_f_b,&gain_b);

	//om_operator_vector_scal_mul(&gain_b,0.0,&gain_b);

	om_operator_vector_add(&gain_a,&gain_b,&gain);
	om_quat_create(&dq,om_vector_getValue(&gain,0),om_vector_getValue(&gain,1),om_vector_getValue(&gain,2),om_vector_getValue(&gain,3));

	//compute gyro bias estimation


	om_operator_quat_mul(&filter->_q_est,&dq,&q_bias);

	/*/
	om_operator_quat_scal_mul(&q_bias,2.0,&q_bias);
	om_quat_imaginary(&q_bias,&filter->_bias_est);
	/*/
	om_operator_quat_scal_mul(&q_bias,2.0*DELTA_T,&q_bias);
	om_quat_imaginary(&q_bias,&delta_bias);
	om_operator_vector_add(&filter->_bias_est,&delta_bias,&filter->_bias_est);
	//*/


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
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
	om_vector_free(&delta_bias);
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



/* initialization of all component used by the nonlinear filter CFA */
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


}

/* process function of the nonlinear filter CFA */
void om_cfa_process(struct omSensorFusionManager *manager,void *filter){

	om_cfa_prediction(manager,(omNonLinearFilter_CFA*)filter);
	om_cfa_update(manager,(omNonLinearFilter_CFA*)filter);

}



void om_cfa_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_CFA *filter){

	//*/
	//printf("COUCOU CONNARD!!!! %f  \n",om_vector_norm(&manager->imu_data.data_accelerometer) - G);
	if(  om_vector_norm(&manager->imu_data.data_accelerometer) - G < 0.03 )
		filter->_beta = 5.0;
	else
		filter->_beta = 0.5;
	//*/

	if(om_vector_norm(&manager->imu_data.data_gyroscope) > 0.0)
		om_kinematics_quaternion(&filter->_q_est,&manager->imu_data.data_gyroscope,&filter->_q_pred);
	else
		om_quat_create(&filter->_q_pred,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);

	//om_vector_normalize(&manager->imu_data.data_accelerometer);
	om_vector_normalize(&manager->imu_data.data_magnetometer);

	om_rotate_vector_quaternion(&filter->_q_pred,&ned_gravity,&filter->_v_acc_pred);
	//om_rotate_vector_quaternion(&filter->_q_pred,&ned_gravity_normalized,&filter->_v_acc_pred);
	om_rotate_vector_quaternion(&filter->_q_pred,&ned_magnetic_field_normalized,&filter->_v_mag_pred);
	//om_rotate_vector_quaternion(&filter->_q_pred,&ned_magnetic_field,&filter->_v_mag_pred);

	//om_operator_vector_scal_mul(&filter->_v_acc_pred,G,&filter->_v_acc_pred);



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
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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



/* release all component used for the nonlinear filter CFA */
void om_cfa_free(void *filter){

	//free memory
	om_vector_free(&(*(omNonLinearFilter_CFA*)filter)._v_acc_pred);
	om_vector_free(&(*(omNonLinearFilter_CFA*)filter)._v_mag_pred);

}


///////////////////////////////////////////////////////
/////           NonLinearFilter QUEST             /////
///////////////////////////////////////////////////////


/* initialization of all component used by the nonlinear filter CFA */
void om_quest_initialization(struct omSensorFusionManager *manager,void *filter){


	// attitude initialization
	switch(manager->type){

	// convertion according to user choices
	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_QUEST*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_QUEST*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_QUEST*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_QUEST*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	(*(omNonLinearFilter_QUEST*)filter)._n = 2;
	(*(omNonLinearFilter_QUEST*)filter)._index = 0;
	(*(omNonLinearFilter_QUEST*)filter)._bool_start = 0;

	// init weights measuerment vector & reference vector
	(*(omNonLinearFilter_QUEST*)filter)._a = (double*)malloc((*(omNonLinearFilter_QUEST*)filter)._n*sizeof(double));
	(*(omNonLinearFilter_QUEST*)filter)._b = (omVector*)malloc((*(omNonLinearFilter_QUEST*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_QUEST*)filter)._r = (omVector*)malloc((*(omNonLinearFilter_QUEST*)filter)._n*sizeof(omVector));


	for(int i=0;i<(*(omNonLinearFilter_QUEST*)filter)._n;i+=2){
		(*(omNonLinearFilter_QUEST*)filter)._a[i] = 1.0 / manager->imu_params.variance_accelerometer;
		(*(omNonLinearFilter_QUEST*)filter)._a[i+1] = 1.0 / manager->imu_params.variance_magnetometer;

		om_vector_create(&(*(omNonLinearFilter_QUEST*)filter)._r[i],3);
		om_vector_create(&(*(omNonLinearFilter_QUEST*)filter)._r[i+1],3);

		om_vector_clone(&ned_gravity_normalized, &(*(omNonLinearFilter_QUEST*)filter)._r[i]);
		om_vector_clone(&ned_magnetic_field_normalized, &(*(omNonLinearFilter_QUEST*)filter)._r[i+1]);

	}

	//normalize weight
	double sum = 0.0;
	for(int i=0;i<(*(omNonLinearFilter_QUEST*)filter)._n;i++){
		//printf("a[i] = %f\n",(*(omNonLinearFilter_QUEST*)filter)._a[i]);
		sum +=  (*(omNonLinearFilter_QUEST*)filter)._a[i];
	}


	for(int i=0;i<(*(omNonLinearFilter_QUEST*)filter)._n;i++){

		(*(omNonLinearFilter_QUEST*)filter)._a[i] /= sum;

		//(*(omNonLinearFilter_QUEST*)filter)._a[i] = 1.0/((double)(*(omNonLinearFilter_QUEST*)filter)._n);
		om_vector_create(&(*(omNonLinearFilter_QUEST*)filter)._b[i],3);


		//printf("after a[i] = %f\n",(*(omNonLinearFilter_QUEST*)filter)._a[i]);
	}

	//printf("\ntoto = %f\n",toto);
	//printf("sum = %f\n",sum);
}

/* process function of the nonlinear filter CFA */
void om_quest_process(struct omSensorFusionManager *manager,void *filter){

	om_quest_preprocess(manager,(omNonLinearFilter_QUEST*)filter);

	if((*(omNonLinearFilter_QUEST*)filter)._bool_start){
		om_quest_prediction(manager,(omNonLinearFilter_QUEST*)filter);
		om_quest_update(manager,(omNonLinearFilter_QUEST*)filter);
	}

}



void om_quest_preprocess(struct omSensorFusionManager *manager,omNonLinearFilter_QUEST *filter){

	/* update measurement vector */
	omVector z_acc;
	omVector z_mag;

	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);

	om_vector_clone(&manager->imu_data.data_accelerometer, &z_acc);
	om_vector_clone(&manager->imu_data.data_magnetometer, &z_mag);
	//om_operator_vector_sub(&manager->imu_data.data_accelerometer,&manager->imu_params.bias_accelerometer,&z_acc);
	//om_operator_vector_sub(&manager->imu_data.data_magnetometer,&manager->imu_params.bias_magnetometer,&z_mag);

	om_vector_normalize(&z_acc);
	om_vector_normalize(&z_mag);

	om_vector_setValues(&filter->_b[filter->_index*2],3,om_vector_getValue(&z_acc,0),om_vector_getValue(&z_acc,1),om_vector_getValue(&z_acc,2));
	om_vector_setValues(&filter->_b[filter->_index*2 + 1],3,om_vector_getValue(&z_mag,0),om_vector_getValue(&z_mag,1),om_vector_getValue(&z_mag,2));

	filter->_index = (filter->_index + 1) % (filter->_n/2);
	filter->_bool_start = ( filter->_bool_start == 1 || (filter->_index == (filter->_n/2 - 1)) )? 1 : 0;

	om_vector_free(&z_acc);
	om_vector_free(&z_mag);

	//printf("index = %d\n",filter->_index);

}

void om_quest_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_QUEST *filter){

	omMatrix B;
	omMatrix B_t;
	omVector z;
	omMatrix S;
	omMatrix S_tmp;
	omMatrix K;
	omMatrix I;

	om_matrix_create(&B,3,3);
	om_matrix_create(&B_t,3,3);
	om_matrix_create(&S,3,3);
	om_matrix_create(&S_tmp,3,3);
	om_matrix_create(&K,4,4);
	om_matrix_createIdentity(&I,3);
	om_vector_create(&z,3,0.0,0.0,0.0);

	for(int i=0;i<filter->_n;i++){

		/*/
		printf("b[i] = ");om_vector_display(&filter->_b[i]);
		printf("r[i] = ");om_vector_display(&filter->_r[i]);
		//printf("a[i] = %f\n",(*(omNonLinearFilter_QUEST*)filter)._a[i]);
		//*/

		omVector cross;
		omMatrix m_bi;
		omMatrix m_ri_t;
		omMatrix B_tmp;

		om_vector_create(&cross,3);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&B_tmp,3,3);

		// compute _d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		om_vector_crossProduct(&filter->_b[i],&filter->_r[i],&cross);
		om_operator_vector_scal_mul(&cross,filter->_a[i],&cross);
		om_operator_vector_add(&z,&cross,&z);


		//compute B = B + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&filter->_b[i],j));
			om_matrix_setValue(&m_ri_t,0,j,om_vector_getValue(&filter->_r[i],j));
		}

		om_operator_matrix_mul(&m_bi,&m_ri_t,&B_tmp);
		om_operator_matrix_scal_mul(&B_tmp,filter->_a[i],&B_tmp);
		om_operator_matrix_add(&B,&B_tmp,&B);

		//free memory
		om_vector_free(&cross);
		om_matrix_free(&m_bi);
		om_matrix_free(&m_ri_t);
		om_matrix_free(&B_tmp);

	}

	//printf("\nMatrix B\n");
	//om_matrix_display(&B);

	// get transpose B
	om_matrix_transpose(&B,&B_t);

	// compute S = B + B_t
	om_operator_matrix_add(&B,&B_t,&S);

	// compute S - I*tr(B)
	double trB = om_matrix_trace(&B);
	om_operator_matrix_scal_mul(&I,trB,&I);
	om_operator_matrix_sub(&S,&I,&S_tmp);


	// compute matrix K
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
			om_matrix_setValue(&K,i,j,om_matrix_getValue(&S_tmp,i,j));

		om_matrix_setValue(&K,i,3,om_vector_getValue(&z,i));
		om_matrix_setValue(&K,3,i,om_vector_getValue(&z,i));
	}
	om_matrix_setValue(&K,3,3,trB);


	// get eigen values of K
	omVector* eigen_vector;
	double* eigen_values;
	om_matrix_getEingenValues(&K,&eigen_vector,&eigen_values,5);

	//*/
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
	om_quat_normalize(&filter->_q_est);
	/*/

	// calcul of lambda_max
	double lambda=0.0;
	for(int i=0;i<4;i++){

		//get the max values
		if(lambda < eigen_values[i]){
			lambda = eigen_values[i];
		}

		// free memory
		om_vector_free(&eigen_vector[i]);

	}


	omVector x;
	omMatrix adj_S;
	omMatrix square_S;
	omMatrix tot_S;

	om_vector_create(&x,3);
	om_matrix_create(&adj_S,3,3);
	om_matrix_create(&square_S,3,3);
	om_matrix_create(&tot_S,3,3);

	om_matrix_adjugate(&S,&adj_S);

	double alpha = (lambda*lambda) - (trB*trB) + om_matrix_trace(&adj_S);
	double gamma = (alpha*( lambda + trB )) - om_matrix_determinant(&S);

	for(int i=0;i<3;i++)
		om_matrix_setValue(&I,i,i,alpha);

	om_operator_matrix_scal_mul(&S,(lambda - trB),&S_tmp);
	om_operator_matrix_mul(&S,&S,&square_S);
	om_operator_matrix_add(&I,&S_tmp,&tot_S);
	om_operator_matrix_add(&tot_S,&square_S,&tot_S);

	om_operator_matrix_vector_mul(&tot_S,&z,&x);

	double qw = gamma;
	double qx = om_vector_getValue(&x,0);
	double qy = om_vector_getValue(&x,1);
	double qz = om_vector_getValue(&x,2);


	// update q_est
	if( !isnan(qw) && !isnan(qx)&& !isnan(qy)&& !isnan(qz) )
		om_quat_create(&filter->_q_est,qw,qx,qy,qz);


	//om_matrix_free(&adj_S);
	//om_matrix_free(&square_S);
	//om_matrix_free(&tot_S);
	//om_vector_free(&x);


	//*/

	// free memory
	om_matrix_free(&B);
	om_matrix_free(&B_t);
	om_matrix_free(&S);
	om_matrix_free(&S_tmp);
	om_matrix_free(&K);
	om_vector_free(&z);



	free(eigen_values);
	free(eigen_vector);
	eigen_values = 0;
	eigen_vector = 0;


}


void om_quest_update(struct omSensorFusionManager *manager,omNonLinearFilter_QUEST *filter){


	om_quat_normalize(&filter->_q_est);

	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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


}



/* release all component used for the nonlinear filter CFA */
void om_quest_free(void *filter){


	for(int i=0;i<(*(omNonLinearFilter_QUEST*)filter)._n;i++){
		om_vector_free(&(*(omNonLinearFilter_QUEST*)filter)._r[i]);
		om_vector_free(&(*(omNonLinearFilter_QUEST*)filter)._b[i]);
	}

	free((*(omNonLinearFilter_QUEST*)filter)._r);
	free((*(omNonLinearFilter_QUEST*)filter)._b);
	free((*(omNonLinearFilter_QUEST*)filter)._a);

	(*(omNonLinearFilter_QUEST*)filter)._r = 0;
	(*(omNonLinearFilter_QUEST*)filter)._b = 0;
	(*(omNonLinearFilter_QUEST*)filter)._a = 0;


}



///////////////////////////////////////////////////////
/////           NonLinearFilter EnKPF             /////
///////////////////////////////////////////////////////



void om_enkpf_initialization(struct omSensorFusionManager *manager,void *filter){

	// init q_0
	switch(manager->type){

	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_EnKPF*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
		break;

	case Matrix:
		om_convert_matrix2quaternion(&manager->output.matrix,&(*(omNonLinearFilter_EnKPF*)filter)._q_est);
		break;

	case EulerAngle:
		om_convert_euler2quaternion(&manager->output.euler,&(*(omNonLinearFilter_EnKPF*)filter)._q_est);
		break;

	default:
		om_quat_create(&(*(omNonLinearFilter_EnKPF*)filter)._q_est,1.0,0.0,0.0,0.0);
		break;

	}

	(*(omNonLinearFilter_EnKPF*)filter)._seed = 0.0;
	(*(omNonLinearFilter_EnKPF*)filter)._h = 0.1;
	(*(omNonLinearFilter_EnKPF*)filter)._f = 4.0;
	(*(omNonLinearFilter_EnKPF*)filter)._lambda = 5.0;
	(*(omNonLinearFilter_EnKPF*)filter)._n = 50.0;
	(*(omNonLinearFilter_EnKPF*)filter)._k_I = 0.3;
	(*(omNonLinearFilter_EnKPF*)filter)._k_P = 1.0;
	(*(omNonLinearFilter_EnKPF*)filter)._k_acc = 2.0;
	(*(omNonLinearFilter_EnKPF*)filter)._k_mag = 2.0;
	(*(omNonLinearFilter_EnKPF*)filter)._threeshold = ( (double)(*(omNonLinearFilter_EnKPF*)filter)._n)/2.0 ;

	// init state
	double b_x = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double b_y = om_vector_getValue(&manager->imu_params.bias_gyroscope,1);
	double b_z = om_vector_getValue(&manager->imu_params.bias_gyroscope,2);

	om_vector_create(&(*(omNonLinearFilter_EnKPF*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);

	//compute Cholesky Decomposition of matrix Q in order to generate noise
	double var_u = om_vector_getValue(&manager->imu_params.bias_gyroscope,0);
	double var_v = manager->imu_params.variance_gyroscope;

	//var_u = 0.000031623;
	//var_v = 0.0031623;

	//*/
	om_matrix_create(&(*(omNonLinearFilter_EnKPF*)filter)._cov_L,6,6);
	double rate = DELTA_T;//*0.0001;
	for(int i=0;i<3;i++){
		double tmp = ((var_v*var_v)/rate) + ((var_u*var_u*rate)/12.0);

		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._cov_L,i,i,sqrt(tmp));
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._cov_L,i+3,i+3, var_u*sqrt(rate));
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._cov_L,i,i+3,(-0.5)*var_u*sqrt(rate));
	}
	/*/


	double delta_t = DELTA_T;//*0.0001;
	om_matrix_create(&(*(omNonLinearFilter_EnKPF*)filter)._L,6,6);
	for(int i=0;i<3;++i){

		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._L,i,i,(var_v*delta_t + 0.33333*var_u*(pow(delta_t,3.0))));
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._L,i+3,i,-(0.5*var_u*delta_t*delta_t));
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._L,i,i+3,-(0.5*var_u*delta_t*delta_t));
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._L,i+3,i+3,(var_u*delta_t));
	}
	//*/

	// compute covariance matrix R
	om_matrix_create(&(*(omNonLinearFilter_EnKPF*)filter)._R,6,6);
	for(int i=0;i<3;++i){
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._R,i,i,manager->imu_params.variance_accelerometer);
		om_matrix_setValue(&(*(omNonLinearFilter_EnKPF*)filter)._R,i+3,i+3,manager->imu_params.variance_magnetometer);
	}

	(*(omNonLinearFilter_EnKPF*)filter)._particle_q = (omQuaternion*)malloc((*(omNonLinearFilter_EnKPF*)filter)._n*sizeof(omQuaternion));
	(*(omNonLinearFilter_EnKPF*)filter)._particle_x = (omVector*)malloc((*(omNonLinearFilter_EnKPF*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_EnKPF*)filter)._particle_wn = (omVector*)malloc((*(omNonLinearFilter_EnKPF*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_EnKPF*)filter)._particle_z = (omVector*)malloc((*(omNonLinearFilter_EnKPF*)filter)._n*sizeof(omVector));
	(*(omNonLinearFilter_EnKPF*)filter)._particle_w = (double*)malloc((*(omNonLinearFilter_EnKPF*)filter)._n*sizeof(double));

	double w_0 = 1.0/ ((double)(*(omNonLinearFilter_EnKPF*)filter)._n);
	double f = (*(omNonLinearFilter_EnKPF*)filter)._f;

	// generation of N particles a time t=0
	for(int i=0;i<(*(omNonLinearFilter_EnKPF*)filter)._n;i++){

		om_vector_create(&(*(omNonLinearFilter_EnKPF*)filter)._particle_z[i],6,0.0,0.0,0.0,0.0,0.0,0.0);

		// generate normal random number
		om_vector_create(&(*(omNonLinearFilter_EnKPF*)filter)._particle_wn[i],6,0.0,0.0,0.0,0.0,0.0,0.0);
		om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&(*(omNonLinearFilter_EnKPF*)filter)._cov_L,(*(omNonLinearFilter_EnKPF*)filter)._seed,&(*(omNonLinearFilter_EnKPF*)filter)._particle_wn[i]);
		(*(omNonLinearFilter_EnKPF*)filter)._seed += 6.0;

		// perturbation of x_0
		om_vector_create(&(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i],6,0.0,0.0,0.0,0.0,0.0,0.0);
		om_operator_vector_add(&(*(omNonLinearFilter_EnKPF*)filter)._x_k,&(*(omNonLinearFilter_EnKPF*)filter)._particle_wn[i], &(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i]);

		// get d_p_0_i
		omVector dp_0_i;
		om_vector_create(&dp_0_i,3);
		for(int l=0;l<3;l++)
			om_vector_setValue(&dp_0_i,l,om_vector_getValue(&(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i],l));

		// compute d_q_0_i
		omQuaternion dq_0_i;

		double dp_0_i_norm = om_vector_norm(&dp_0_i);
		double dp_0_i_norm_square = dp_0_i_norm*dp_0_i_norm;

		double dq_0_i_w = ( (f*f) - dp_0_i_norm_square)/ ( (f*f) + dp_0_i_norm_square);
		double dq_0_i_x = om_vector_getValue(&(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i],0)* ( (1.0 + dq_0_i_w ) / f);
		double dq_0_i_y = om_vector_getValue(&(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i],1)* ( (1.0 + dq_0_i_w ) / f);
		double dq_0_i_z = om_vector_getValue(&(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i],2)* ( (1.0 + dq_0_i_w ) / f);

		om_quat_create(&dq_0_i,dq_0_i_w,dq_0_i_x,dq_0_i_y,dq_0_i_z);

		// compute q_0_i
		om_operator_quat_mul(&dq_0_i,&(*(omNonLinearFilter_EnKPF*)filter)._q_est,&(*(omNonLinearFilter_EnKPF*)filter)._particle_q[i]);
		om_quat_normalize(&(*(omNonLinearFilter_EnKPF*)filter)._particle_q[i]);

		// init weights at time 0
		(*(omNonLinearFilter_EnKPF*)filter)._particle_w[i] = w_0;

		om_vector_free(&dp_0_i);

	}

	//init resample to false
	(*(omNonLinearFilter_EnKPF*)filter)._resample = 0;


}


void om_enkpf_free(void *filter){

	om_vector_free(&(*(omNonLinearFilter_EnKPF*)filter)._x_k);
	om_matrix_free(&(*(omNonLinearFilter_EnKPF*)filter)._cov_L);
	om_matrix_free(&(*(omNonLinearFilter_EnKPF*)filter)._R);


	for(int i=0;i<(*(omNonLinearFilter_EnKPF*)filter)._n;i++){
		om_vector_free(&(*(omNonLinearFilter_EnKPF*)filter)._particle_wn[i]);
		om_vector_free(&(*(omNonLinearFilter_EnKPF*)filter)._particle_x[i]);
		om_vector_free(&(*(omNonLinearFilter_EnKPF*)filter)._particle_z[i]);
	}

	free((*(omNonLinearFilter_EnKPF*)filter)._particle_q);
	free((*(omNonLinearFilter_EnKPF*)filter)._particle_x);
	free((*(omNonLinearFilter_EnKPF*)filter)._particle_w);
	free((*(omNonLinearFilter_EnKPF*)filter)._particle_z);
	free((*(omNonLinearFilter_EnKPF*)filter)._particle_wn);

	(*(omNonLinearFilter_EnKPF*)filter)._particle_q = 0;
	(*(omNonLinearFilter_EnKPF*)filter)._particle_x = 0;
	(*(omNonLinearFilter_EnKPF*)filter)._particle_w = 0;
	(*(omNonLinearFilter_EnKPF*)filter)._particle_z = 0;
	(*(omNonLinearFilter_EnKPF*)filter)._particle_wn = 0;


}


void om_enkpf_process(struct omSensorFusionManager *manager,void *filter){

	(*(omNonLinearFilter_EnKPF*)filter)._seed = (*(omNonLinearFilter_EnKPF*)filter)._seed > 100000000.0 ? 0.0 : (*(omNonLinearFilter_EnKPF*)filter)._seed;

	om_enkpf_prediction(manager,(omNonLinearFilter_EnKPF*)filter);
	om_enkpf_update(manager,(omNonLinearFilter_EnKPF*)filter);
}

void om_enkpf_prediction(struct omSensorFusionManager *manager,omNonLinearFilter_EnKPF *filter){


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


	// propagate quaternion
	omVector b;
	omVector omega;

	om_vector_create(&omega,3);
	om_vector_create(&b,3,om_vector_getValue(&filter->_x_k,3),
 						  om_vector_getValue(&filter->_x_k,4),
					      om_vector_getValue(&filter->_x_k,5));

	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b,&omega);
	om_kinematics_quaternion(&filter->_q_est,&omega,&filter->_q_pred);


	//set sum_w to 0.0
	filter->_sum_w = 0.0;

	//compute
	omQuaternion q_pred_inv;
	om_quat_inverse(&filter->_q_pred,&q_pred_inv);

	omVector x_mean;
	omVector z_mean;

	om_vector_create(&x_mean,6,0.0,0.0,0.0,0.0,0.0,0.0);
	om_vector_create(&z_mean,6,0.0,0.0,0.0,0.0,0.0,0.0);

	// compute N particle of angular velocity and propagate _q_k_i to _q_(k+1)_i
	for(int i=0;i<filter->_n;++i){

		// generate normal random number
		omVector noise_dp;
		omVector noise_b;

		//om_random_generateWhiteNoiseFromCovarianceMatrix(0.0,&filter->_L,filter->_seed,&filter->_particle_wn[i]);
		//filter->_seed += 6.0;

		om_vector_create(&noise_dp,3,om_vector_getValue(&filter->_particle_wn[i],0),
									 om_vector_getValue(&filter->_particle_wn[i],1),
									 om_vector_getValue(&filter->_particle_wn[i],2));

		om_vector_create(&noise_b,3,om_vector_getValue(&filter->_particle_wn[i],3),
									om_vector_getValue(&filter->_particle_wn[i],4),
									om_vector_getValue(&filter->_particle_wn[i],5));
		// compute angular velocity
		omVector b_i;
		omVector angular_velocity_i;
		omQuaternion q_pred_i;
		om_vector_create(&angular_velocity_i,3);
		om_vector_create(&b_i,3,om_vector_getValue(&filter->_particle_x[i],3),
								om_vector_getValue(&filter->_particle_x[i],4),
								om_vector_getValue(&filter->_particle_x[i],5));

		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);
		om_operator_vector_add(&angular_velocity_i,&noise_dp,&angular_velocity_i);

		// propagate _q_k_i to _q_(k+1)_i
		//om_kinematics_quaternion(&filter->_particle_q[i],&angular_velocity_i,&filter->_particle_q[i]);
		om_kinematics_quaternion(&filter->_particle_q[i],&angular_velocity_i,&q_pred_i);
		//om_kinematics_quaternion(&filter->_q_est,&angular_velocity_i,&q_pred_i);

		//*/
		////
		// constraint computation
		/////

		omVector cross_acc_i;
		omVector cross_mag_i;
		omVector omega_i;
		omVector omega_q_i;
		omVector b_i_tmp;

		omVector v_acc_pred_i;
		omVector v_mag_pred_i;

		om_vector_create(&v_acc_pred_i,3);
		om_vector_create(&v_mag_pred_i,3);
		om_vector_create(&cross_acc_i,3);
		om_vector_create(&cross_mag_i,3);
		om_vector_create(&omega_i,3);
		om_vector_create(&omega_q_i,3);
		om_vector_create(&b_i_tmp,3);


		// compute v_acc = R(q_(k+1)^-)g^a
		om_rotate_vector_quaternion(&q_pred_i,&ned_gravity_normalized,&v_acc_pred_i);

		// compute v_mag = R(q_(k+1)^-)m^a
		om_rotate_vector_quaternion(&q_pred_i,&ned_magnetic_field_normalized,&v_mag_pred_i);

		om_vector_normalize(&v_acc_pred_i);
		om_vector_normalize(&v_mag_pred_i);


		// compute cross product
		om_vector_crossProduct(&z_acc_tmp,&v_acc_pred_i,&cross_acc_i);
		om_vector_crossProduct(&z_mag_tmp,&v_mag_pred_i,&cross_mag_i);


		om_operator_vector_scal_mul(&cross_acc_i,filter->_k_acc/2.0,&cross_acc_i);
		om_operator_vector_scal_mul(&cross_mag_i,filter->_k_mag/2.0,&cross_mag_i);

		// compute gain omega
		om_operator_vector_add(&cross_acc_i,&cross_mag_i,&omega_i);


		// compute bias_est =  (_omega*(_k_I)*_delta_t*(-1.0));
		om_operator_vector_scal_mul(&omega_i,filter->_k_I*DELTA_T*(-1.0),&b_i_tmp);
		//om_operator_vector_scal_mul(&omega_i,filter->_k_I*DELTA_T*(-1.0),&b_i);
		om_operator_vector_add(&b_i,&b_i_tmp,&b_i);



		// compute angular_velocity = y_gyro - b_est + omega*k_P
		om_operator_vector_scal_mul(&omega_i,filter->_k_P,&omega_q_i);
		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity_i);
		om_operator_vector_add(&angular_velocity_i,&omega_q_i,&angular_velocity_i);

		// compute q_(k+1)^+ = Omega(Y_gyro - b_pred + omega*k_P)q_(k)^+
		om_kinematics_quaternion(&filter->_particle_q[i],&angular_velocity_i,&filter->_particle_q[i]);
		//om_kinematics_quaternion(&filter->_q_est,&angular_velocity_i,&filter->_particle_q[i]);
		om_quat_normalize(&filter->_particle_q[i]);
		//*/

		// propagate bias_k_i to bias_(k+1)_i
		om_operator_vector_add(&b_i,&noise_b,&b_i);
		//printf("\nnoise_b = ");om_vector_display(&noise_b);
		//printf("b_i = ");om_vector_display(&b_i);

		// compute d_q_k_i
		omQuaternion dq_i;
		om_operator_quat_mul(&filter->_particle_q[i],&q_pred_inv,&dq_i);
		//printf("filter->_particle_q[i] = ");om_quat_display(&filter->_particle_q[i]);
		//printf("dq_i = ");om_quat_display(&dq_i);


		// propagate d_p_k_i to d_p_(k+1)_i
		omVector dp_i;
		double tmp = (filter->_f * (dq_i._qw < 0.0 ? -1.0 : 1.0) )/(1.0 + fabs(dq_i._qw));
		om_vector_create(&dp_i,3,dq_i._qx * tmp,dq_i._qy * tmp,dq_i._qz * tmp);

		// update _x_k_i
		for(int l=0;l<3;l++){
			om_vector_setValue(&filter->_particle_x[i],l,om_vector_getValue(&dp_i,l));
			om_vector_setValue(&filter->_particle_x[i],l+3,om_vector_getValue(&b_i,l));
		}
		//printf("particle_w[i] = %f   particle_x[i] = ",filter->_particle_w[i]);om_vector_display(&filter->_particle_x[i]);

		// compute z_k_i
		omVector z_acc;
		omVector z_mag;
		om_vector_create(&z_acc,3);
		om_vector_create(&z_mag,3);

		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_gravity_normalized,&z_acc);
		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_magnetic_field_normalized,&z_mag);

		om_vector_normalize(&z_acc);
		om_vector_normalize(&z_mag);

		for(int l=0;l<3;l++){
			om_vector_setValue(&filter->_particle_z[i],l,om_vector_getValue(&z_acc,l));
			om_vector_setValue(&filter->_particle_z[i],l+3,om_vector_getValue(&z_mag,l));
		}

		// compute mean
		om_operator_vector_add(&x_mean,&filter->_particle_x[i],&x_mean);
		om_operator_vector_add(&z_mean,&filter->_particle_z[i],&z_mean);

		//free memory
		om_vector_free(&b_i);
		om_vector_free(&angular_velocity_i);
		om_vector_free(&noise_b);
		om_vector_free(&noise_dp);
		om_vector_free(&dp_i);
		om_vector_free(&z_acc);
		om_vector_free(&z_mag);

		//*/
		om_vector_free(&v_acc_pred_i);
		om_vector_free(&v_mag_pred_i);
		om_vector_free(&cross_acc_i);
		om_vector_free(&cross_mag_i);
		om_vector_free(&omega_i);
		om_vector_free(&omega_q_i);
		om_vector_free(&b_i_tmp);
		//*/



	}


	om_operator_vector_scal_div(&x_mean,filter->_n,&x_mean);
	om_operator_vector_scal_div(&z_mean,filter->_n,&z_mean);

	/*/
	printf("x_mean = ");
	om_vector_display(&x_mean);
	printf("z_mean = ");
	om_vector_display(&z_mean);
	//*/

	// covariance and cross covariance
	omMatrix Pxz;
	omMatrix Pzz;

	om_matrix_create(&Pxz,6,6);
	om_matrix_create(&Pzz,6,6);

	for(int i=0;i<filter->_n;++i){

		omMatrix m_x_i;
		omMatrix m_z_i;
		omMatrix m_z_i_t;
		omMatrix P_xz_i;
		omMatrix P_zz_i;

		om_matrix_create(&P_xz_i,6,6);
		om_matrix_create(&P_zz_i,6,6);

		om_matrix_create(&m_x_i,6,1);
		om_matrix_create(&m_z_i,6,1);
		om_matrix_create(&m_z_i_t,1,6);

		for(int j=0;j<6;++j){
			om_matrix_setValue(&m_x_i,j,0,om_vector_getValue(&filter->_particle_x[i],j) - om_vector_getValue(&x_mean,j) );
			om_matrix_setValue(&m_z_i,j,0,om_vector_getValue(&filter->_particle_z[i],j) - om_vector_getValue(&z_mean,j));
			om_matrix_setValue(&m_z_i_t,0,j,om_vector_getValue(&filter->_particle_z[i],j) - om_vector_getValue(&z_mean,j));

		}

		om_operator_matrix_mul(&m_x_i,&m_z_i_t,&P_xz_i);
		om_operator_matrix_mul(&m_z_i,&m_z_i_t,&P_zz_i);

		om_operator_matrix_add(&Pxz,&P_xz_i,&Pxz);
		om_operator_matrix_add(&Pzz,&P_zz_i,&Pzz);

		om_matrix_free(&m_x_i);
		om_matrix_free(&m_z_i_t);
		om_matrix_free(&m_z_i);
		om_matrix_free(&P_xz_i);
		om_matrix_free(&P_zz_i);

	}

	om_operator_matrix_scal_div(&Pzz,(double)(filter->_n - 1) ,&Pzz);
	om_operator_matrix_scal_div(&Pxz,(double)(filter->_n - 1) ,&Pxz);


	for(int j=0;j<6;++j){
		double tmp = om_matrix_getValue(&Pzz,j,j) + (om_matrix_getValue(&filter->_R,j,j)/filter->_lambda);
		om_matrix_setValue(&Pzz,j,j,tmp);
	}

	// kalman gain
	omMatrix K;
	omMatrix P_zz_inv;

	om_matrix_create(&K,6,6);
	om_matrix_create(&P_zz_inv,6,6);

	om_matrix_inverse(&Pzz,&P_zz_inv);
	om_operator_matrix_mul(&Pxz,&P_zz_inv,&K);



	for(int i=0;i<filter->_n;++i){

		// correction x
		omVector s_k_i;
		omVector x_corr;
		om_vector_create(&s_k_i,6);
		om_vector_create(&x_corr,6);

		om_operator_vector_sub(&z_k,&filter->_particle_z[i],&s_k_i);
		om_operator_matrix_vector_mul(&K,&s_k_i,&x_corr);
		om_operator_vector_add(&filter->_particle_x[i],&x_corr,&filter->_particle_x[i]);

		// compute dq
		omVector dp_i;
		omQuaternion dq_i;

		om_vector_create(&dp_i,3,om_vector_getValue(&filter->_particle_x[i],0),
							     om_vector_getValue(&filter->_particle_x[i],1),
							     om_vector_getValue(&filter->_particle_x[i],2));

		double dp_norm = om_vector_norm(&dp_i);
		double dp_norm_square = dp_norm*dp_norm;

		double dq_w = ( (filter->_f*filter->_f) - dp_norm_square)/ ( (filter->_f*filter->_f) + dp_norm_square);
		double dq_x = om_vector_getValue(&dp_i,0)* ( (1.0 + dq_w ) / filter->_f);
		double dq_y = om_vector_getValue(&dp_i,1)* ( (1.0 + dq_w ) / filter->_f);
		double dq_z = om_vector_getValue(&dp_i,2)* ( (1.0 + dq_w ) / filter->_f);

		om_quat_create(&dq_i,dq_w,dq_x,dq_y,dq_z);

		// update q_est
		om_operator_quat_mul(&dq_i,&filter->_q_pred,&filter->_particle_q[i]);

		// compute z_k_i
		omVector z_acc;
		omVector z_mag;
		omVector z_k_i;

		om_vector_create(&z_acc,3);
		om_vector_create(&z_mag,3);
		om_vector_create(&z_k_i,6);

		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_gravity_normalized,&z_acc);
		om_rotate_vector_quaternion(&filter->_particle_q[i],&ned_magnetic_field_normalized,&z_mag);

		om_vector_normalize(&z_mag);
		om_vector_normalize(&z_acc);

		for(int l=0;l<3;l++){
			om_vector_setValue(&z_k_i,l,om_vector_getValue(&z_acc,l));
			om_vector_setValue(&z_k_i,l+3,om_vector_getValue(&z_mag,l));
		}
		// likehood function
		om_operator_vector_sub(&z_k,&z_k_i,&s_k_i);

		double tmp_acc = 0.0;
		double tmp_mag = 0.0;

		for(int l=0;l<3;l++){
			tmp_acc += om_vector_getValue(&s_k_i,l)*om_vector_getValue(&s_k_i,l)/om_matrix_getValue(&filter->_R, l, l);
			tmp_mag += om_vector_getValue(&s_k_i,l+3)*om_vector_getValue(&s_k_i,l+3)/om_matrix_getValue(&filter->_R, l+3, l+3);
		}

		//tmp_acc *= (1.0/manager->imu_params.variance_accelerometer);
		//tmp_mag *= (1.0/manager->imu_params.variance_magnetometer);

		double L_k_i = exp( (tmp_acc+tmp_mag)*(-0.5) );

		// update weights
		filter->_particle_w[i] *= L_k_i;

		// accumulator to normalize weights
		filter->_sum_w += filter->_particle_w[i];



		//free memory
		om_vector_free(&dp_i);
		om_vector_free(&x_corr);
		om_vector_free(&z_acc);
		om_vector_free(&z_mag);
		om_vector_free(&z_k_i);
		om_vector_free(&s_k_i);

	}


	//free memory
	om_matrix_free(&Pxz);
	om_matrix_free(&Pzz);
	om_matrix_free(&K);
	om_matrix_free(&P_zz_inv);

	om_vector_free(&z_acc_tmp);
	om_vector_free(&z_mag_tmp);
	om_vector_free(&b);
	om_vector_free(&omega);
	om_vector_free(&z_k);
	om_vector_free(&x_mean);
	om_vector_free(&z_mean);



}

void om_enkpf_update(struct omSensorFusionManager *manager,omNonLinearFilter_EnKPF *filter){

	// compute mean x_k^p
	omVector x_k_sigma;

	om_vector_setValues(&filter->_x_k,6,0.0,0.0,0.0,0.0,0.0,0.0);
	om_vector_create(&x_k_sigma,6,0.0,0.0,0.0,0.0,0.0,0.0);

	filter->_sum_n_eff = 0.0;


	//upadate particles
	for(int i=0;i<filter->_n;++i){

		//update weights
		if(filter->_sum_w > 0.0)
			filter->_particle_w[i] /= filter->_sum_w;
		else{
			filter->_particle_w[i] = 1.0/((double)filter->_n);
			//printf("WRONG\n");
			//exit(1);
		}
		//printf("particle_w[i] = %f   particle_x[i] = ",filter->_particle_w[i]);om_vector_display(&filter->_particle_x[i]);

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

	//printf("\n");

	om_operator_vector_scal_div(&x_k_sigma,(double)(filter->_n),&x_k_sigma);

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


	// set output
	double signQw = filter->_q_est._qw < 0.0 ? -1.0 :1.0;
	double qw = filter->_q_est._qw * signQw;
	double qx = filter->_q_est._qx * signQw;
	double qy = filter->_q_est._qy * signQw;
	double qz = filter->_q_est._qz * signQw;

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
		om_quat_create(&manager->output.quaternion,filter->_q_est._qw,filter->_q_est._qx,filter->_q_est._qy,filter->_q_est._qz);
		break;

	}

	// resampling step
	om_enkpf_resampling(filter);

	// if resample step has been performed
	 if(filter->_resample == 1){

		 //printf("resampling \n");

		 //compute convariance matrix S
		 omMatrix S;
		 omMatrix L;
		 om_matrix_create(&L,6,6);
		 om_matrix_create(&S,6,6);

		 for(int i=0;i<filter->_n;++i){
			 //printf("after particle_w[i] = %f   particle_x[i] = ",filter->_particle_w[i]);om_vector_display(&filter->_particle_x[i]);
			 omVector dx;
			 om_vector_create(&dx,6);
			 om_operator_vector_sub(&filter->_particle_x[i],&x_k_sigma,&dx);

			 omMatrix tmp_S_i_t;
			 omMatrix tmp_S_i;
			 om_matrix_create(&tmp_S_i,6,1);
			 om_matrix_create(&tmp_S_i_t,1,6);

			 for(int l=0;l<6;l++)
				 om_matrix_setValue(&tmp_S_i,l,0,om_vector_getValue(&dx,l));

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

		 //printf("\n");

		 om_operator_matrix_scal_div(&S,filter->_n - 1,&S);
		 om_operator_matrix_scal_mul(&S,(filter->_h*filter->_h),&S);

		 //om_matrix_choleskyDecomposition(&S, &L);
		 om_matrix_squareRoot(&S, &L, 50);

		 //printf("Matrix S\n");om_matrix_display(&S);
		 //printf("Matrix L\n");om_matrix_display(&L);
		// exit(1);

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
			 om_quat_normalize(&filter->_particle_q[i]);

			 om_vector_free(&dp_i);
			 om_vector_free(&noise);
			 //printf("pertu particle_w[i] = %f   _particle_q[i] = ",filter->_particle_w[i]);om_quat_display(&filter->_particle_q[i]);

		 }


		 om_matrix_free(&S);

	 }

	 om_vector_free(&x_k_sigma);
	 om_vector_free(&dp);
}

void om_enkpf_resampling(omNonLinearFilter_EnKPF *filter){

	filter->_resample = 0;

	double N_eff = floor(1.0/filter->_sum_n_eff);
	//printf("N_eff = %f\n",N_eff);

	if(N_eff < filter->_threeshold){


		int n = (int)N_eff;
		filter->_resample = 1;

		//*/
		//om_enkpf_quicksort(filter,0,filter->_n-1);

		double* c = (double*)malloc(filter->_n*sizeof(double));
		c[0] = 0;

		for (int i = 1; i < (filter->_n); ++i) {
			c[i] = c[i-1] + filter->_particle_w[i];
		}

		double u1 = om_random_uniformDistribution(filter->_seed)*(1.0/filter->_n);
		filter->_seed+=1.0;

		int index = 0;
		for (int j = 1; j <= (filter->_n); ++j) {

			double uj = u1 + (j - 1.0)/(filter->_n);
			while(index < (filter->_n) && uj > c[index])
				index++;

			if(index < (filter->_n) ){
				//om_enkpf_swap(filter, index, j);
				om_vector_clone(&filter->_particle_x[index],&filter->_particle_x[j-1]);
			}
			filter->_particle_w[j-1] = (1.0/filter->_n);
		}

		free(c);

		/*/
		om_enkpf_quicksort(filter,0,filter->_n-1);

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
		//*/


	}

}

void om_enkpf_swap(omNonLinearFilter_EnKPF *filter,int i,int j){

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


void om_enkpf_quicksort(omNonLinearFilter_EnKPF *filter,int left, int right){

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
        	om_enkpf_swap(filter,i,j);
            i++;
            j--;
        }
        else{
            if(left<j)
            	om_enkpf_quicksort(filter,left, j);
            if(i<right)
            	om_enkpf_quicksort(filter,i,right);
            return;
        }
    }

}


