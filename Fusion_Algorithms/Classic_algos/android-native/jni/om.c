/*
 * om.c
 *
 *  Created on: 25 May, 2016
 *      Author: thomas
 */

#include "om.h"



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
	om_quat_create(&(*(omNonLinearFilter_CGO*)filter)._q_pred,1.0,0.0,0.0,0.0);

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
	om_rotate_vector_quaternion(&(*(omNonLinearFilter_CGO*)filter)._q_pred,&ned_magnetic_field,&v_mag_pred);

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
	om_quat_create(&manager->output.quaternion,(*(omNonLinearFilter_CGO*)filter)._q_est._qw,(*(omNonLinearFilter_CGO*)filter)._q_est._qx,(*(omNonLinearFilter_CGO*)filter)._q_est._qy,(*(omNonLinearFilter_CGO*)filter)._q_est._qz);

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
	om_vector_create(&(*(omNonLinearFilter_USQUE*)filter)._x_k,6,0.0,0.0,0.0,b_x,b_y,b_z);

	init_ned_frame();

}


void om_usque_process(struct omSensorFusionManager *manager,void *filter){



}


void om_usque_prediction(struct omSensorFusionManager *manager,void *filter){

	omMatrix S;
	omMatrix T;

	// compute S = sqrt( (P + Q)*(n + lambda )
	om_operator_matrix_add(&(*(omNonLinearFilter_USQUE*)filter)._P_k,&(*(omNonLinearFilter_USQUE*)filter)._Q,&T);
	om_operator_matrix_const_mul(&T,(*(omNonLinearFilter_USQUE*)filter)._lambda + 6.0,&T);
	om_matrix_squareRoot(&T,&S);


	// generate sigma point at time k
	omVector* sigma_points = (omVector*)malloc(13 * sizeof(omVector));

	for(int i=0;i<=6;++i){

		if(i==0){

			om_vector_clone(&(*(omNonLinearFilter_USQUE*)filter)._x_k,&sigma_points[0]);

		}else{

			omVector S_row;

			om_matrix_getRow(&S,(i-1),&S_row);
			om_operator_vector_add(&(*(omNonLinearFilter_USQUE*)filter)._x_k,&S_row,&sigma_points[i-1]);
			om_operator_vector_add(&(*(omNonLinearFilter_USQUE*)filter)._x_k,&S_row,&sigma_points[6 + i-1]);
		}
	}

	// generate sigma quaternion at time k
	omQuaternion* sigma_quaternion = (omQuaternion*)malloc(13 * sizeof(omQuaternion));

	for(unsigned int i=0;i<13;++i){

		if(i==0){

			double qw = (*(omNonLinearFilter_USQUE*)filter)._q_est._qw;
			double qx = (*(omNonLinearFilter_USQUE*)filter)._q_est._qx;
			double qy = (*(omNonLinearFilter_USQUE*)filter)._q_est._qy;
			double qz = (*(omNonLinearFilter_USQUE*)filter)._q_est._qz;

			om_quat_create(&sigma_quaternion[0],qw,qx,qy,qz);

		}else{

			// compute error quaternion i
			omVector dp_i;
			om_vector_create(&dp_i,3,om_vector_getValue(&sigma_points[i],0),om_vector_getValue(&sigma_points[i],1),om_vector_getValue(&sigma_points[i],2));
			double dp_i_norm = om_vector_norm(&dp_i);
			double dp_i_norm_sq = dp_i_norm*dp_i_norm;
			double a = (*(omNonLinearFilter_USQUE*)filter)._a;
			double f = (*(omNonLinearFilter_USQUE*)filter)._f;
			double dqw_i = ( ( a*(-1.0)*dp_i_norm_sq ) + (f*sqrt( f*f + (1.0 - a*a)*dp_i_norm_sq  ))) / ( (f*f) + dp_i_norm_sq );

			omVector dqv_i;
			om_operator_vector_const_mul(&dp_i,( (1.0/f) * (a + dqw_i )),&dqv_i);

			omQuaternion dq_i;
			om_quat_create(&dq_i,dqw_i,om_vector_getValue(&dqv_i,0),om_vector_getValue(&dqv_i,1),om_vector_getValue(&dqv_i,2) );

			// compute sigma quaternion i
			om_operator_quat_mul(&dq_i,&(*(omNonLinearFilter_USQUE*)filter)._q_est,&sigma_quaternion[i]);


		}
	}

	// propagation of sigma quaternion and sigma points
	for(unsigned int i=0;i<13;++i){


		omVector b_i;
		om_vector_create(&b_i,3,om_vector_getValue(&sigma_points[i],3),om_vector_getValue(&sigma_points[i],4),om_vector_getValue(&sigma_points[i],5));

		//
		omVector angular_velocity;
		om_operator_vector_sub(&manager->imu_data.data_gyroscope,&b_i,&angular_velocity);

		// propagation of sigma quaternion
		omQuaternion q_pred_i;
		om_kinematics_quaternion(&sigma_quaternion[i],&angular_velocity,&q_pred_i);
		/*
		if(i==0){

			_sigma_quaternion_k[0] = q_pred_i;
			_sigma_points_x_kp1.push_back(Vector(6,0.0,0.0,0.0,b_i.getValue(0),b_i.getValue(1),b_i.getValue(2)));

		}else{

			_sigma_quaternion_k[i] = q_pred_i;

			Quaternion dq_i = q_pred_i*_sigma_quaternion_k[0].inverse();

			Vector dqv_i (3,dq_i.getQx(),dq_i.getQy(),dq_i.getQz());
			Vector dp_i = dqv_i*(_f/(_a+dq_i.getQw()));

			_sigma_points_x_kp1.push_back(Vector(6,dp_i.getValue(0),dp_i.getValue(1),dp_i.getValue(2),b_i.getValue(0),b_i.getValue(1),b_i.getValue(2)));

		}*/

	}





}

void om_usque_update(struct omSensorFusionManager *manager,void *filter){

}








