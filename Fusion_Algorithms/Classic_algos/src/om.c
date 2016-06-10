/*
 * om.c
 *
 *  Created on: 25 May, 2016
 *      Author: thomas
 */

#include "om.h"

omVector ned_magnetic_field;
omVector ned_gravity;
omVector ned_geographic_north;

void init_ned_frame(){

	om_vector_create(&ned_geographic_north,3,1.0,0.0,0.0);
	om_vector_create(&ned_magnetic_field,3,0.9687632386,0.003889250977,-0.2479569747);
	om_vector_create(&ned_gravity,3,0.0,0.0,1.0);

}


///////////////////////////////////////////////////////
/////           NonLinearFilter CGO               /////
///////////////////////////////////////////////////////


void om_cgo_initialization(struct omSensorFusionManager *manager,void *filter){

	(*(omNonLinearFilter_CGO*)filter)._seed = 0;
	(*(omNonLinearFilter_CGO*)filter)._k_I = 0.3;
	(*(omNonLinearFilter_CGO*)filter)._k_P = 1.0;
	(*(omNonLinearFilter_CGO*)filter)._k_acc = 1.0;
	(*(omNonLinearFilter_CGO*)filter)._k_mag = 1.0;

	om_quat_create(&(*(omNonLinearFilter_CGO*)filter)._q_est,1.0,0.0,0.0,0.0);
	switch(manager->type){

	case Quarternion:
		om_quat_create(&(*(omNonLinearFilter_CGO*)filter)._q_est,manager->output.quaternion._qw,manager->output.quaternion._qx,manager->output.quaternion._qy,manager->output.quaternion._qz);
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

	om_vector_clone(&manager->imu_params.bias_gyroscope,&(*(omNonLinearFilter_CGO*)filter)._bias_est);
	om_vector_clone(&manager->imu_params.bias_gyroscope,&(*(omNonLinearFilter_CGO*)filter)._bias_pred);

	init_ned_frame();

}

void om_cgo_process(struct omSensorFusionManager *manager,void *filter){

	//////
	// Prediction step
	/////



	// compute b_(k+1)^- = b_(k)^+
	for(int i=0;i<3;++i)
		om_vector_setValue(&(*(omNonLinearFilter_CGO*)filter)._bias_pred,i,(*(omNonLinearFilter_CGO*)filter)._bias_est._values[i]);

	// compute angular_velocity = y_gyro - b_pred
	omVector angular_velocity;
	om_vector_create(&angular_velocity,3);
	om_operator_vector_sub(&manager->imu_data.data_gyroscope,&(*(omNonLinearFilter_CGO*)filter)._bias_pred,&angular_velocity);

	// compute q_(k+1)^- = Omega(Y_gyro - b_pred)q_(k)^+
	om_kinematics_quaternion(&(*(omNonLinearFilter_CGO*)filter)._q_est,&angular_velocity,&(*(omNonLinearFilter_CGO*)filter)._q_pred);

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
	omVector omega_b;
	omVector q_tmp;

	omMatrix Xi;
	omMatrix S_acc;
	omMatrix S_mag;
	omMatrix S;

	omQuaternion dq;

	// allocation
	om_vector_create(&z_acc,3);
	om_vector_create(&z_mag,3);
	om_vector_create(&omega_b,3);
	om_vector_create(&omega_q,3);
	om_vector_create(&omega_b,3);
	om_vector_create(&q_tmp,4);
	om_matrix_create(&S,3,3);

	// compute v_acc = R(q_(k+1)^-)g^a
	om_rotate_vector_quaternion(&(*(omNonLinearFilter_CGO*)filter)._q_pred,&ned_gravity,&v_acc_pred);

	// compute v_mag = R(q_(k+1)^-)m^a
	om_rotate_vector_quaternion(&(*(omNonLinearFilter_CGO*)filter)._q_pred,&ned_geographic_north,&v_mag_pred);

	// toto
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

	// compute S_acc and S_mag
	om_matrix_skewSymetricMatrix(&cross_acc,&S_acc);
	om_matrix_skewSymetricMatrix(&cross_mag,&S_mag);

	om_operator_matrix_const_mul(&S_acc,(-1.0) * (*(omNonLinearFilter_CGO*)filter)._k_acc /2.0,&S_acc);
	om_operator_matrix_const_mul(&S_mag,(-1.0) * (*(omNonLinearFilter_CGO*)filter)._k_mag /2.0,&S_mag);

	// compute omega = - vex (S_acc+S_mag)
	om_operator_matrix_add(&S_acc,&S_mag,&S);
	om_operator_vex(&S,&omega);
	om_operator_vector_const_mul(&omega,(-1.0),&omega);

	// compute omega_q = omega * k_P * 0.5 * DELTA_T
	om_operator_vector_const_mul(&omega,(*(omNonLinearFilter_CGO*)filter)._k_P * (0.5*DELTA_T),&omega_q);

	// compute omega_b = omega * k_I * DELTA_T
	om_operator_vector_const_mul(&omega,(*(omNonLinearFilter_CGO*)filter)._k_I * DELTA_T,&omega_b);

	// compute Xi(q_(k+1)^-)
	om_operator_xi(&(*(omNonLinearFilter_CGO*)filter)._q_pred , &Xi);

	// compute dq = Xi(q_(k+1)^-)*omega_q
	om_operator_matrix_vector_mul(&Xi,&omega_q,&q_tmp);
	om_quat_create(&dq,q_tmp._values[3],q_tmp._values[0],q_tmp._values[1],q_tmp._values[2]);


	//compute q_(k+1)^+ = q_(k+1)^- + dq
	om_operator_quat_add(&(*(omNonLinearFilter_CGO*)filter)._q_pred, &dq , &(*(omNonLinearFilter_CGO*)filter)._q_est);

	//compute b_(k+1)^+ = b_(k+1)^- - omega_b
	om_operator_vector_sub(&(*(omNonLinearFilter_CGO*)filter)._bias_pred,&omega_b,&(*(omNonLinearFilter_CGO*)filter)._bias_est);

	// set output
	double qw = (*(omNonLinearFilter_CGO*)filter)._q_est._qw;
	double qx = (*(omNonLinearFilter_CGO*)filter)._q_est._qx;
	double qy = (*(omNonLinearFilter_CGO*)filter)._q_est._qy;
	double qz = (*(omNonLinearFilter_CGO*)filter)._q_est._qz;


	om_quat_create(&manager->output.quaternion,qw,qx,qy,qz);



	//dispose object
	om_vector_dispose(&z_acc);
	om_vector_dispose(&z_mag);
	om_vector_dispose(&cross_acc);
	om_vector_dispose(&cross_mag);
	om_vector_dispose(&omega);
	om_vector_dispose(&omega_q);
	om_vector_dispose(&omega_b);
	om_vector_dispose(&q_tmp);
	om_vector_dispose(&v_mag_pred);
	om_vector_dispose(&v_acc_pred);

	om_matrix_dispose(&S_acc);
	om_matrix_dispose(&S_mag);
	om_matrix_dispose(&S);
	om_matrix_dispose(&Xi);





}


///////////////////////////////////////////////////////
/////           NonLinearFilter USQUE             /////
///////////////////////////////////////////////////////


void om_usque_initialization(struct omSensorFusionManager *manager,void *filter){

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
	om_operator_matrix_const_mul(&T,lambda + 6.0,&T);
	om_matrix_choleskyDecomposition(&T,&S);

	// generate sigma point at time k
	omVector* sigma_points = (omVector*)malloc(13 * sizeof(omVector));
	om_vector_clone(&filter->_x_k,&sigma_points[0]);

	for(int i=1;i<=6;++i){

		omVector S_row;
		om_vector_create(&sigma_points[i],6);
		om_vector_create(&sigma_points[i+6],6);


		om_matrix_getRow(&S,(i-1),&S_row);
		om_operator_vector_add(&filter->_x_k,&S_row,&sigma_points[i]);
		om_operator_vector_sub(&filter->_x_k,&S_row,&sigma_points[i+6]);

		om_vector_dispose(&S_row);

	}


	om_matrix_dispose(&S);
	om_matrix_dispose(&T);

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
			om_operator_vector_const_mul(&dp_i,( (1.0/f) * (a + dqw_i )),&dqv_i);
			omQuaternion dq_i;
			om_quat_create(&dq_i,dqw_i,om_vector_getValue(&dqv_i,0),om_vector_getValue(&dqv_i,1),om_vector_getValue(&dqv_i,2) );

			// compute sigma quaternion i
			om_operator_quat_mul(&dq_i,&filter->_q_est,&sigma_quaternion[i]);
			om_vector_dispose(&dqv_i);

		}

		om_vector_dispose(&dp_i);


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

		om_vector_dispose(&b_i);
		om_vector_dispose(&angular_velocity_i);
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

	om_operator_vector_const_mul(&sigma_points[0],lambda/(6.0+lambda),&sum_a);
	om_operator_vector_const_mul(&sum_x,1.0/(2.0*(6.0+lambda)),&sum_b);
	om_operator_vector_add(&sum_a,&sum_b,&filter->_x_k_pred);

	om_vector_dispose(&sum_a);
	om_vector_dispose(&sum_b);
	om_vector_dispose(&sum_x);


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

		//compute  X_k(i) - x_k^-
		om_operator_vector_sub(&sigma_points[i],&filter->_x_k_pred,&var_x_i);
		for(int i=0;i<6;i++)
			om_matrix_setValue(&tmp_i,i,0,om_vector_getValue(&var_x_i,i));

		//compute P_k^- +=  (x_k^- - X_k(i))(x_k^- - X_k(i))^T
		om_matrix_transpose(&tmp_i,&tmp_i_t);
		om_operator_matrix_mul(&tmp_i,&tmp_i_t,&P_i);
		om_operator_matrix_add(&sum_P,&P_i,&sum_P);

		//free memory
		om_vector_dispose(&var_x_i);
		om_matrix_dispose(&tmp_i);
		om_matrix_dispose(&tmp_i_t);
		om_matrix_dispose(&P_i);

	}


	omVector var_x_0;
	omMatrix tmp_0;
	omMatrix tmp_0_t;
	omMatrix P_0;

	om_vector_create(&var_x_0,6);
	om_matrix_create(&tmp_0,6,1);
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
	om_operator_matrix_const_mul(&P_0,lambda/(6.0+lambda),&sum_P_a);
	om_operator_matrix_const_mul(&sum_P,1.0/(2.0*(6.0+lambda)),&sum_P_b);
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
	om_vector_dispose(&var_x_0);
	om_matrix_dispose(&tmp_0);
	om_matrix_dispose(&tmp_0_t);
	om_matrix_dispose(&P_0);
	om_matrix_dispose(&sum_P);
	om_matrix_dispose(&sum_P_a);
	om_matrix_dispose(&sum_P_b);

	for(int i=0;i<13;++i){
		om_vector_dispose(&sigma_points[i]);
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

		om_rotate_vector_quaternion(&sigma_quaternion[i],&ned_gravity,&tmp_res_a);
		om_rotate_vector_quaternion(&sigma_quaternion[i],&ned_geographic_north,&tmp_res_m);

		om_operator_vector_const_mul(&tmp_res_a,G,&tmp_res_a);

		//om_vector_normalize(&tmp_res_m);

		om_vector_create(&sigma_z[i],6,om_vector_getValue(&tmp_res_a,0),om_vector_getValue(&tmp_res_a,1),om_vector_getValue(&tmp_res_a,2),
									   om_vector_getValue(&tmp_res_m,0),om_vector_getValue(&tmp_res_m,1),om_vector_getValue(&tmp_res_m,2));

        om_vector_dispose(&tmp_res_m);
        om_vector_dispose(&tmp_res_a);



		//printf("sigma_z[i] = ");
		//om_vector_display(&sigma_z[i]);

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

	om_operator_vector_const_mul(&sigma_z[0],lambda/(6.0+lambda),&sum_a);
	om_operator_vector_const_mul(&sum_z,1.0/(2.0*(6.0+lambda)),&sum_b);
	om_operator_vector_add(&sum_a,&sum_b,&mean_z);

	om_vector_dispose(&sum_a);
	om_vector_dispose(&sum_b);
	om_vector_dispose(&sum_z);

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
		om_vector_dispose(&var_x_i);
		om_vector_dispose(&var_z_i);

		om_matrix_dispose(&tmp_z_i);
		om_matrix_dispose(&tmp_x_i);
		om_matrix_dispose(&tmp_z_i_t);
		om_matrix_dispose(&P_xz_i);
		om_matrix_dispose(&P_zz_i);



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
	om_operator_matrix_const_mul(&P_xz_0,lambda/(6.0+lambda),&sum_P_xz_a);
	om_operator_matrix_const_mul(&sum_P_xz,1.0/(2.0*(6.0+lambda)),&sum_P_xz_b);
	om_operator_matrix_const_mul(&P_zz_0,lambda/(6.0+lambda),&sum_P_zz_a);
	om_operator_matrix_const_mul(&sum_P_zz,1.0/(2.0*(6.0+lambda)),&sum_P_zz_b);

	om_operator_matrix_add(&sum_P_xz_a,&sum_P_xz_b,&P_xz);
	om_operator_matrix_add(&sum_P_zz_a,&sum_P_zz_b,&P_zz);
	om_operator_matrix_add(&filter->_R,&P_zz,&P_zz);


	//Kalman gain
	omMatrix K;
	omMatrix P_zz_inv;

	om_matrix_create(&K,6,6);
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
	om_operator_quat_mul(&dq,&filter->_q_pred,&manager->output.quaternion);

	//free memory
	om_vector_dispose(&z_k);
	om_vector_dispose(&dp_est);
	om_vector_dispose(&mean_z);
	om_vector_dispose(&s_z_k);
	om_vector_dispose(&corr_x);
	om_vector_dispose(&var_x_0);
	om_vector_dispose(&var_z_0);


	om_matrix_dispose(&tmp_x_0);
	om_matrix_dispose(&tmp_z_0);
	om_matrix_dispose(&tmp_z_0_t);
	om_matrix_dispose(&K);
	om_matrix_dispose(&K_t);
	om_matrix_dispose(&S);
	om_matrix_dispose(&S_tmp);
	om_matrix_dispose(&sum_P_xz);
	om_matrix_dispose(&sum_P_zz);
	om_matrix_dispose(&P_zz);
	om_matrix_dispose(&P_zz_inv);
	om_matrix_dispose(&P_xz);
	om_matrix_dispose(&sum_P_xz_a);
	om_matrix_dispose(&sum_P_xz_b);
	om_matrix_dispose(&sum_P_zz_a);
	om_matrix_dispose(&sum_P_zz_b);

	for(int i=0;i<13;++i){
		om_vector_dispose(&sigma_z[i]);
	}

	free(sigma_z);
	sigma_z = 0;

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
	(*(omNonLinearFilter_REQUEST*)filter)._a[0]=0.34;
	(*(omNonLinearFilter_REQUEST*)filter)._a[1]=0.33;
	(*(omNonLinearFilter_REQUEST*)filter)._a[2]=0.33;


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

		(*(omNonLinearFilter_REQUEST*)filter)._d_m_k += (*(omNonLinearFilter_REQUEST*)filter)._a[i];
		(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k += om_vector_dotProduct(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],&(*(omNonLinearFilter_REQUEST*)filter)._r[i])*(*(omNonLinearFilter_REQUEST*)filter)._a[i];

		// compute _d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		om_vector_crossProduct(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],&(*(omNonLinearFilter_REQUEST*)filter)._r[i],&cross);
		om_operator_vector_const_mul(&cross,(*(omNonLinearFilter_REQUEST*)filter)._a[i],&cross);
		om_operator_vector_add(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,&cross,&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k);


		//compute _d_B_k = _d_B_k + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&M,3,3);

		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&(*(omNonLinearFilter_REQUEST*)filter)._b[i],j));
			om_matrix_setValue(&m_ri_t,0,j,om_vector_getValue(&(*(omNonLinearFilter_REQUEST*)filter)._r[i],j));
		}

		om_operator_matrix_mul(&m_bi,&m_ri_t,&M);
		om_operator_matrix_const_mul(&M,(*(omNonLinearFilter_REQUEST*)filter)._a[i],&M);
		om_operator_matrix_add(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,&M,&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k);

		//free memory
		om_vector_dispose(&cross);
		om_matrix_dispose(&m_bi);
		om_matrix_dispose(&m_ri_t);
		om_matrix_dispose(&M);


	}

	(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k /= (*(omNonLinearFilter_REQUEST*)filter)._d_m_k;
	om_operator_matrix_const_div(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,(*(omNonLinearFilter_REQUEST*)filter)._d_m_k,&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k);
	om_operator_vector_const_div(&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k,(*(omNonLinearFilter_REQUEST*)filter)._d_m_k,&(*(omNonLinearFilter_REQUEST*)filter)._d_z_k);


	// init S_0
	omMatrix d_B_k_t;

	om_matrix_transpose(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,&d_B_k_t);
	om_matrix_create(&(*(omNonLinearFilter_REQUEST*)filter)._d_S_k,3,3);

	//compute _d_S_k =_d_B_k + _d_B_k.transpose();
	om_operator_matrix_add(&(*(omNonLinearFilter_REQUEST*)filter)._d_B_k,&d_B_k_t,&(*(omNonLinearFilter_REQUEST*)filter)._d_S_k);

	om_matrix_dispose(&d_B_k_t);

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
	om_operator_matrix_const_mul(&I,(*(omNonLinearFilter_REQUEST*)filter)._d_sigma_k,&I);
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


	om_matrix_dispose(&T);
	om_matrix_dispose(&I);

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


	om_vector_dispose(&z_acc);
	om_vector_dispose(&z_mag);
	om_vector_dispose(&y_axis);

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

		filter->_d_m_k += filter->_a[i];
		filter->_d_sigma_k += om_vector_dotProduct(&filter->_b[i],&filter->_r[i])*filter->_a[i];

		// compute _d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		om_vector_crossProduct(&filter->_b[i],&filter->_r[i],&cross);
		om_operator_vector_const_mul(&cross,filter->_a[i],&cross);
		om_operator_vector_add(&filter->_d_z_k,&cross,&filter->_d_z_k);


		//compute _d_B_k = _d_B_k + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri_t,1,3);
		om_matrix_create(&M,3,3);

		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&filter->_b[i],j));
			om_matrix_setValue(&m_ri_t,0,j,om_vector_getValue(&filter->_r[i],j));
		}

		om_operator_matrix_mul(&m_bi,&m_ri_t,&M);
		om_operator_matrix_const_mul(&M,filter->_a[i],&M);
		om_operator_matrix_add(&filter->_d_B_k,&M,&filter->_d_B_k);

		//free memory
		om_vector_dispose(&cross);
		om_matrix_dispose(&m_bi);
		om_matrix_dispose(&m_ri_t);
		om_matrix_dispose(&M);


	}



	filter->_d_sigma_k /= filter->_d_m_k;
	om_operator_matrix_const_div(&filter->_d_B_k,filter->_d_m_k,&filter->_d_B_k);
	om_operator_vector_const_div(&filter->_d_z_k,filter->_d_m_k,&filter->_d_z_k);


	//compute _d_S_k =_d_B_k + _d_B_k.transpose();
	omMatrix d_B_k_t;
	om_matrix_transpose(&filter->_d_B_k,&d_B_k_t);
	om_operator_matrix_add(&filter->_d_B_k,&d_B_k_t,&filter->_d_S_k);
	om_matrix_dispose(&d_B_k_t);

	// compute d_K
	omMatrix I;
	omMatrix T;

	om_matrix_create(&T,3,3);

	for(int i=0;i<4;++i)
		for(int j=0;j<4;++j)
			om_matrix_setValue(&filter->_d_K_k,i,j,0.0);

	om_matrix_createIdentity(&I,3);
	om_operator_matrix_const_mul(&I,filter->_d_sigma_k,&I);
	om_operator_matrix_sub(&filter->_d_S_k,&I,&T);


	//set d_K_k values
	for(int i=0;i<3;i++){

		om_matrix_setValue(&filter->_d_K_k,i,3,om_vector_getValue(&filter->_d_z_k,i));
		om_matrix_setValue(&filter->_d_K_k,3,i,om_vector_getValue(&filter->_d_z_k,i));

		for(int j=0;j<3;j++)
			om_matrix_setValue(&filter->_d_K_k,i,j,om_matrix_getValue(&T,i,j));

	}

	om_matrix_setValue(&filter->_d_K_k,3,3,filter->_d_sigma_k);

	om_matrix_dispose(&T);
	om_matrix_dispose(&I);

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
	om_matrix_dispose(&Phi);
	om_matrix_dispose(&Phi_t);
	om_matrix_dispose(&tmp_K);
	om_matrix_dispose(&tmp_P);
	om_vector_dispose(&angular_velocity);



}

void om_request_update(struct omSensorFusionManager *manager,omNonLinearFilter_REQUEST *filter){

	// compute optimal gain
	double trace_P = om_matrix_trace(&filter->_P_k_pred);
	double trace_R = om_matrix_trace(&filter->_R);
	double a = (filter->_m_k*filter->_m_k)*trace_P;
	double rho_k_opt = ( a )/( a + ((filter->_d_m_k*filter->_d_m_k)*trace_R)  );
	rho_k_opt=0.001;


	// propagate m_k
	double m_kp1 = ((1.0 - rho_k_opt)*filter->_m_k) + (rho_k_opt*filter->_d_m_k);

	// some variable
	double tmp_a = (1.0 - rho_k_opt)*(filter->_m_k/m_kp1);
	double tmp_b = (rho_k_opt)*(filter->_d_m_k/m_kp1);

	// update K
	om_operator_matrix_const_mul(&filter->_K_k_pred,tmp_a,&filter->_K_k_pred);
	om_operator_matrix_const_mul(&filter->_d_K_k,tmp_b,&filter->_d_K_k);
	om_operator_matrix_add(&filter->_K_k_pred,&filter->_d_K_k,&filter->_K_k);

	// update P
	om_operator_matrix_const_mul(&filter->_P_k_pred,tmp_a*tmp_a,&filter->_P_k_pred);
	om_operator_matrix_const_mul(&filter->_R,tmp_b*tmp_b,&filter->_R);
	om_operator_matrix_add(&filter->_P_k_pred,&filter->_R,&filter->_P_k);

	/*
	printf("\n matrix _K_k\n");
	om_matrix_display(&filter->_K_k);

	printf("\n matrix _P_k\n");
	om_matrix_display(&filter->_P_k);
	*/


	// get eigen values of K
	omVector* eigen_vector;
	double* eigen_values;
	om_matrix_getEingenValues(&filter->_K_k,&eigen_vector,&eigen_values,50);

	// calcul of lambda_max
	double lambda=0.0;
	for(int i=0;i<4;i++){


		if(lambda < eigen_values[i]){

			lambda = eigen_values[i];

			double qw = om_vector_getValue(&eigen_vector[i],3);
			double qx = om_vector_getValue(&eigen_vector[i],0);
			double qy = om_vector_getValue(&eigen_vector[i],1);
			double qz = om_vector_getValue(&eigen_vector[i],2);

			om_quat_create(&filter->_q_est,qw,qx,qy,qz);
			om_quat_create(&manager->output.quaternion,qw,qx,qy,qz);

		}

		om_vector_dispose(&eigen_vector[i]);

	}

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


		om_matrix_create(&m_bi,3,1);
		om_matrix_create(&m_ri,3,1);

		//vectorToMatrix
		for(int j=0;j<3;j++){
			om_matrix_setValue(&m_bi,j,0,om_vector_getValue(&filter->_b[i],j));
			om_matrix_setValue(&m_ri,j,0,om_vector_getValue(&filter->_r[i],j));
		}

		// compute S_ri = skewSymetricMatrix(_r[i]);
		om_matrix_skewSymetricMatrix(&filter->_r[i],&S_ri);

		// compute I = identity(3)*(3.0 -(dotProduct(_r[i],_b[i])*dotProduct(_r[i],_b[i])));
		om_matrix_createIdentity(&I,3);
		om_operator_matrix_const_mul(&I,3.0 - (dot*dot),&I);

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
		om_operator_matrix_const_mul(&R11_i_a,dot,&R11_i_a);
		om_operator_matrix_add(&R11_i_a,&I,&R11_i_a);


		// compute (S_ri*(m_bi*m_bi.transpose())*S_ri.transpose())  )
		om_operator_matrix_mul(&S_ri,&tmp_bi_bi,&tmp);
		om_operator_matrix_mul(&tmp,&S_ri_t,&R11_i_b);

		//compute R11 = R11 + ( I + ( ( (m_bi*m_ri.transpose()) +(m_ri*m_bi.transpose())   )*dotProduct(_b[i],_r[i]) )  + (S_ri*(m_bi*m_bi.transpose())*S_ri.transpose())  );
		om_operator_matrix_add(&R11_i_a,&R11_i_b,&tmp);
		om_operator_matrix_add(&R11,&tmp,&R11);



		//free memory
		om_matrix_dispose(&m_bi);
		om_matrix_dispose(&m_ri);
		om_matrix_dispose(&S_ri);
		om_matrix_dispose(&S_ri_t);
		om_matrix_dispose(&I);
		om_matrix_dispose(&tmp_bi_bi);
		om_matrix_dispose(&tmp_bi_ri);
		om_matrix_dispose(&tmp_ri_bi);
		om_matrix_dispose(&m_bi_t);
		om_matrix_dispose(&m_ri_t);
		om_matrix_dispose(&R11_i_a);
		om_matrix_dispose(&R11_i_b);
		om_matrix_dispose(&tmp);


	}

	//compute R11 = R11*(_mu_k/n_k);
	om_operator_matrix_const_mul(&R11,(filter->_mu_k/n_k),&R11);


	//set R values
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			om_matrix_setValue(&filter->_R,i,j,om_matrix_getValue(&R11,i,j));

	om_matrix_setValue(&filter->_R,3,3,R22);

	om_matrix_dispose(&R11);



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
	om_operator_matrix_const_mul(&I,tmp,&I);

	om_operator_matrix_sub(&tmp_d_B_k_2,&tmp_d_B_k_square,&T);
	om_operator_matrix_sub(&T,&tmp_d_B_k_t_square,&T);
	om_operator_matrix_const_mul(&T,2.0,&T);

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
	om_matrix_create(&D,3,3);
	om_matrix_createIdentity(&I2,3);
	om_operator_matrix_const_mul(&I2,filter->_d_sigma_k,&I2);

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
	om_operator_vector_const_mul(&Q12,eta_k*(-1.0),&Q12);

	//set Q values
	for(int i=0;i<3;i++){

		om_matrix_setValue(&filter->_Q,i,3,om_vector_getValue(&Q12,i));
		om_matrix_setValue(&filter->_Q,3,i,om_vector_getValue(&Q12,i));

		for(int j=0;j<3;j++)
			om_matrix_setValue(&filter->_Q,i,j,om_matrix_getValue(&Q11,i,j));

	}

	om_matrix_setValue(&filter->_Q,3,3,Q22);
	om_operator_matrix_const_mul(&filter->_Q,DELTA_T*DELTA_T,&filter->_Q);


	//free memory
	om_matrix_dispose(&d_B_k_t);
	om_matrix_dispose(&tmp_d_B_k);
	om_matrix_dispose(&tmp_d_B_k_2);
	om_matrix_dispose(&tmp_d_B_k_square);
	om_matrix_dispose(&tmp_d_B_k_t_square);
	om_matrix_dispose(&I);
	om_matrix_dispose(&T);
	om_matrix_dispose(&Q11);
	om_matrix_dispose(&M);
	om_matrix_dispose(&M_t);
	om_matrix_dispose(&D);
	om_matrix_dispose(&I2);
	om_matrix_dispose(&A);

	om_vector_dispose(&y);
	om_vector_dispose(&x);
	om_vector_dispose(&Q12);

}






