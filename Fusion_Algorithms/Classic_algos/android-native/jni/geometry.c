/*
 * geometry.c
 *
 *  Created on: 17 May, 2016
 *      Author: thomas
 */


#include "geometry.h"

///////////////////////////////////////////////////////
/////           Kinematics - Rotation             /////
///////////////////////////////////////////////////////


void om_rotate_vector_quaternion(struct omQuaternion *q,struct omVector *in,struct omVector *out){

	omQuaternion q_inv;
	omQuaternion q_vec;

	om_quat_inverse(q,&q_inv);
	om_quat_create(&q_vec,0.0,in->_values[0],in->_values[1],in->_values[2]);

	om_operator_quat_mul(&q_inv,&q_vec,&q_vec);
	om_operator_quat_mul(&q_vec,q,&q_vec);

	om_vector_create(out,3,q_vec._qx,q_vec._qy,q_vec._qz);

}

void om_rotate_vector_matrix(struct omMatrix *M,struct omVector *in,struct omVector *out){

	om_operator_matrix_vector_mul(M,in,out);

}

void om_rotate_vector_axisAngle(struct omAxisAngle *aa,struct omVector *in,struct omVector *out){

	omQuaternion q;
	om_convert_axisAngle2quaternion(aa,&q);
	om_rotate_vector_quaternion(&q,in,out);

}


void om_kinematics_quaternion(struct omQuaternion *q_t,struct omVector *angular_velocity,struct omQuaternion *q_tp1){

	omMatrix Omega;
	omMatrix Omega_tmp;
	omMatrix S_phi;
	omMatrix I;
	omVector phi;
	omVector q_t_vec;
	omVector q_tp1_vec;

	om_matrix_create(&Omega_tmp,3,3);
	om_matrix_create(&Omega,4,4);
	om_matrix_createIdentity(&I,3);
	om_vector_create(&q_t_vec,4,q_t->_qx,q_t->_qy,q_t->_qz,q_t->_qw);
	om_vector_create(&q_tp1_vec,4);

	double norm = om_vector_norm(angular_velocity);
	double cos_tmp = cosf(0.5*norm*DELTA_T);
	double sin_tmp = sinf(0.5*norm*DELTA_T);

	om_vector_create(&phi,3);
	om_vector_clone(angular_velocity,&phi);
	om_operator_vector_const_mul(&phi,( sin_tmp/norm ),&phi);

	om_matrix_skewSymetricMatrix(&phi,&S_phi);
	om_operator_matrix_const_mul(&I,cos_tmp,&I);

	om_operator_matrix_sub(&I,&S_phi,&Omega_tmp);

	for(int i=0;i<3;++i){
		om_matrix_setValue(&Omega,i,3,phi._values[i]);
		om_matrix_setValue(&Omega,3,i,-phi._values[i]);

		for(int j=0;j<3;++j)
			om_matrix_setValue(&Omega,i,j,om_matrix_getValue(&Omega_tmp,i,j));
	}

	om_matrix_setValue(&Omega,3,3,cos_tmp);

	om_operator_matrix_vector_mul(&Omega,&q_t_vec,&q_tp1_vec);
	om_quat_create(q_tp1,q_tp1_vec._values[3],q_tp1_vec._values[0],q_tp1_vec._values[1],q_tp1_vec._values[2]);

	om_matrix_dispose(&Omega);
	om_matrix_dispose(&Omega_tmp);
	om_matrix_dispose(&I);
	om_matrix_dispose(&S_phi);
	om_vector_dispose(&phi);
	om_vector_dispose(&q_t_vec);
	om_vector_dispose(&q_tp1_vec);

}


///////////////////////////////////////////////////////
/////                Conversion                   /////
///////////////////////////////////////////////////////

void om_convert_quaternion2matrix(struct omQuaternion *in,struct omMatrix* out){

	omMatrix I;
	omMatrix S;
	omMatrix S_tmp;
	omVector in_img;

	om_matrix_create(out,3,3);
	om_matrix_createIdentity(&I,3);
	om_quat_imaginary(in,&in_img);
	om_matrix_skewSymetricMatrix(&in_img,&S);

	om_operator_matrix_mul(&S,&S,&S_tmp);
	om_operator_matrix_const_mul(&S_tmp,2.0,&S_tmp);

	om_operator_matrix_const_mul(&S,2.0*in->_qw,&S);

	om_operator_matrix_add(&I,&S,out);
	om_operator_matrix_add(out,&S_tmp,out);

	om_matrix_dispose(&I);
	om_matrix_dispose(&S);
	om_matrix_dispose(&S_tmp);
	om_vector_dispose(&in_img);

}

void om_convert_quaternion2euler(struct omQuaternion *in,struct omEulerAngle* out){

	double sqx = in->_qx*in->_qx;
	double sqy = in->_qy*in->_qy;
	double sqz = in->_qz*in->_qz;

	double test = 2.0* ( (in->_qx*in->_qz) - (in->_qw*in->_qy)   );
	double phi,theta,psy;

	/* north pole singularity detected */
	if ( fabs(test + 1.0) < EPSILON ) {

		psy =  atan2(  (in->_qy*in->_qx - in->_qz*in->_qw) , (in->_qx*in->_qz + in->_qw*in->_qy)  );
		theta = PI/2.0;
		phi = 0;

	/* south pole singularity detected */
	} else if (fabs(test - 1.0) < EPSILON ) {

		psy =  - atan2(  (in->_qy*in->_qx - in->_qz*in->_qw) , (in->_qx*in->_qz + in->_qw*in->_qy)  );
		theta = -PI/2.0;
		phi = 0;

	}else{
		phi =  atan2( (in->_qx*in->_qw + in->_qy*in->_qz) , ( 0.5 - (sqx+sqy )));
		theta = asin(- 2.0* ( in->_qx*in->_qz - in->_qw*in->_qy   ) );
		psy =  atan2( (in->_qz*in->_qw + in->_qy*in->_qx) , ( 0.5 - (sqz+sqy )));

	}

	out->_pitch = phi*RAD_TO_DEG;
	out->_roll = theta * RAD_TO_DEG;
	out->_yaw = psy * RAD_TO_DEG;



}

void om_convert_quaternion2axisAngle(struct omQuaternion *in,struct omAxisAngle *out){


	if (  (fabs(in->_qw - 1.0) < EPSILON)  ){

		out->_angle = 0.0;
		om_vector_create(&out->_axis,3,1.0,0.0,0.0);


	}else{

		om_quat_imaginary(in,&out->_axis);

		out->_angle = 2.0*atan2( om_vector_norm(&out->_axis),in->_qw ) *RAD_TO_DEG;
		om_vector_normalize(&out->_axis);


	}

}

void om_convert_matrix2quaternion(struct omMatrix *in,struct omQuaternion *out){


	double tr = om_matrix_trace(in);
	double qw,qx,qy,qz;

	if (tr > 0) {

	  double S = sqrt(tr + 1.0) * 2.0;
	  qw = 0.25 * S;
	  qx = (om_matrix_getValue(in,2,1) - om_matrix_getValue(in,1,2)) / S;
	  qy = (om_matrix_getValue(in,0,2) - om_matrix_getValue(in,2,0)) / S;
	  qz = (om_matrix_getValue(in,1,0) - om_matrix_getValue(in,0,1)) / S;

	} else if ( (om_matrix_getValue(in,0,0) > om_matrix_getValue(in,1,1)) && (om_matrix_getValue(in,0,0) > om_matrix_getValue(in,2,2)) ) {

	  float S = sqrt( 1.0 + om_matrix_getValue(in,0,0) - om_matrix_getValue(in,1,1) - om_matrix_getValue(in,2,2) ) * 2.0; // S=4*qx

	  qw = (om_matrix_getValue(in,2,1) - om_matrix_getValue(in,1,2)) / S;
	  qx = 0.25 * S;
	  qy = (om_matrix_getValue(in,0,1) + om_matrix_getValue(in,1,0)) / S;
	  qz = (om_matrix_getValue(in,0,2) + om_matrix_getValue(in,2,0)) / S;

	} else if ( ( om_matrix_getValue(in,1,1) > om_matrix_getValue(in,2,2)) ) {

	float S = sqrt(1.0 + om_matrix_getValue(in,1,1) - om_matrix_getValue(in,0,0) - om_matrix_getValue(in,2,2)) * 2.0; // S=4*qy
	  qw = (om_matrix_getValue(in,0,2) - om_matrix_getValue(in,2,0)) / S;
	  qx = (om_matrix_getValue(in,0,1) + om_matrix_getValue(in,1,0)) / S;
	  qy = 0.25 * S;
	  qz = (om_matrix_getValue(in,1,2)  + om_matrix_getValue(in,2,1)) / S;


	} else {

	  float S = sqrt(1.0 + om_matrix_getValue(in,2,2) - om_matrix_getValue(in,0,0) - om_matrix_getValue(in,1,1)) * 2.0; // S=4*qz

	  qw = (om_matrix_getValue(in,1,0) - om_matrix_getValue(in,0,1)) / S;
	  qx = (om_matrix_getValue(in,0,2) + om_matrix_getValue(in,2,0)) / S;
	  qy = (om_matrix_getValue(in,1,2) + om_matrix_getValue(in,2,1)) / S;
	  qz = 0.25 * S;

	}

	om_quat_create(out,qw,qx,qy,qz);

}





void om_convert_axisAngle2quaternion(struct omAxisAngle *in,struct omQuaternion *out){

	// angle in rad
	double rad_angle = (in->_angle/2.0) * DEG_TO_RAD;

	double cos_a = cos(rad_angle);
	double sin_a = sin(rad_angle);

	om_quat_create(out,cos_a,in->_axis._values[0]*sin_a,in->_axis._values[1]*sin_a,in->_axis._values[2]*sin_a);

}

void om_convert_axisAngle2matrix(struct omAxisAngle *in,struct omMatrix *out){

	omMatrix I;
	omMatrix S;
	omMatrix S_tmp;

	om_matrix_createIdentity(&I,3);
	om_matrix_skewSymetricMatrix(&in->_axis,&S);


	double rad_angle = (in->_angle) * DEG_TO_RAD;
	double cos_a = cos(rad_angle);
	double sin_a = sin(rad_angle);

	// I + S*sin_a + S*S*( 1.0 - cos_a )
	om_operator_matrix_mul(&S,&S,&S_tmp);
	om_operator_matrix_const_mul(&S_tmp,(1.0 - cos_a),&S_tmp);
	om_operator_matrix_const_mul(&S,sin_a,&S);
	om_operator_matrix_add(&S,&S_tmp,out);
	om_operator_matrix_add(out,&I,out);

	om_matrix_dispose(&I);
	om_matrix_dispose(&S);
	om_matrix_dispose(&S_tmp);

}


void om_convert_euler2quaternion(struct omEulerAngle *in,struct omQuaternion *out){

	//angle in radian
	double phi = in->_pitch*DEG_TO_RAD / 2.0;
	double theta = in->_roll*DEG_TO_RAD / 2.0;
	double psy = in->_yaw*DEG_TO_RAD / 2.0;

	double qw = (cos(phi)*cos(theta)*cos(psy)) + (sin(phi)*sin(theta)*sin(psy));
	double qx = (sin(phi)*cos(theta)*cos(psy)) - (cos(phi)*sin(theta)*sin(psy));
	double qy = (cos(phi)*sin(theta)*cos(psy)) + (sin(phi)*cos(theta)*sin(psy));
	double qz= (cos(phi)*cos(theta)*sin(psy)) - (sin(phi)*sin(theta)*cos(psy));

	om_quat_create(out,qw,qx,qy,qz);

}




///////////////////////////////////////////////////////
/////                Operators                    /////
///////////////////////////////////////////////////////


void om_operator_vex(struct omMatrix *in,struct omVector *out){


	double x = om_matrix_getValue(in,2,1);
	double y = om_matrix_getValue(in,0,2);
	double z = om_matrix_getValue(in,1,0);

	om_vector_create(out,3,x,y,z);

}

void om_operator_omega(struct omVector *in,struct omMatrix *out){

	om_matrix_create(out,4,4);

	om_matrix_setValue(out,0,0,0.0);
	om_matrix_setValue(out,0,1,in->_values[2]);
	om_matrix_setValue(out,0,2,in->_values[1]*(-1.0));
	om_matrix_setValue(out,0,3,in->_values[0]);

	om_matrix_setValue(out,1,0,in->_values[2]*(-1.0));
	om_matrix_setValue(out,1,1,0.0);
	om_matrix_setValue(out,1,2,in->_values[0]);
	om_matrix_setValue(out,1,3,in->_values[1]);

	om_matrix_setValue(out,2,0,in->_values[1]);
	om_matrix_setValue(out,2,1,in->_values[0]*(-1.0));
	om_matrix_setValue(out,2,2,0.0);
	om_matrix_setValue(out,2,3,in->_values[2]);

	om_matrix_setValue(out,3,0,in->_values[0]*(-1.0));
	om_matrix_setValue(out,3,1,in->_values[1]*(-1.0));
	om_matrix_setValue(out,3,2,in->_values[2]*(-1.0));
	om_matrix_setValue(out,3,3,0.0);

}

void om_operator_gamma(struct omVector *in,struct omMatrix *out){

	om_matrix_create(out,4,4);

	om_matrix_setValue(out,0,0,0.0);
	om_matrix_setValue(out,0,1,in->_values[2]*(-1.0));
	om_matrix_setValue(out,0,2,in->_values[1]);
	om_matrix_setValue(out,0,3,in->_values[0]);

	om_matrix_setValue(out,1,0,in->_values[2]);
	om_matrix_setValue(out,1,1,0.0);
	om_matrix_setValue(out,1,2,in->_values[0]*(-1.0));
	om_matrix_setValue(out,1,3,in->_values[1]);

	om_matrix_setValue(out,2,0,in->_values[1]*(-1.0));
	om_matrix_setValue(out,2,1,in->_values[0]);
	om_matrix_setValue(out,2,2,0.0);
	om_matrix_setValue(out,2,3,in->_values[2]);

	om_matrix_setValue(out,3,0,in->_values[0]*(-1.0));
	om_matrix_setValue(out,3,1,in->_values[1]*(-1.0));
	om_matrix_setValue(out,3,2,in->_values[2]*(-1.0));
	om_matrix_setValue(out,3,3,0.0);

}

void om_operator_xi(struct omQuaternion *in,struct omMatrix *out){

	om_matrix_create(out,4,3);


	om_matrix_setValue(out,0,0,in->_qw);
	om_matrix_setValue(out,0,1,-in->_qz);
	om_matrix_setValue(out,0,2,in->_qy);

	om_matrix_setValue(out,1,0,in->_qz);
	om_matrix_setValue(out,1,1,in->_qw);
	om_matrix_setValue(out,1,2,-in->_qx);

	om_matrix_setValue(out,2,0,-in->_qy);
	om_matrix_setValue(out,2,1,in->_qx);
	om_matrix_setValue(out,2,2,in->_qw);

	om_matrix_setValue(out,3,0,-in->_qx);
	om_matrix_setValue(out,3,1,-in->_qy);
	om_matrix_setValue(out,3,2,-in->_qz);


}

void om_operator_psy(struct omQuaternion *in,struct omMatrix *out){

	om_matrix_create(out,4,3);

	om_matrix_setValue(out,0,0,in->_qw);
	om_matrix_setValue(out,0,1,in->_qz);
	om_matrix_setValue(out,0,2,in->_qy*(-1.0));

	om_matrix_setValue(out,1,0,in->_qz*(-1.0));
	om_matrix_setValue(out,1,1,in->_qw);
	om_matrix_setValue(out,1,2,in->_qx);

	om_matrix_setValue(out,2,0,in->_qy);
	om_matrix_setValue(out,2,1,in->_qx*(-1.0));
	om_matrix_setValue(out,2,2,in->_qw);

	om_matrix_setValue(out,3,0,in->_qx*(-1.0));
	om_matrix_setValue(out,3,1,in->_qy*(-1.0));
	om_matrix_setValue(out,3,2,in->_qz*(-1.0));

}










