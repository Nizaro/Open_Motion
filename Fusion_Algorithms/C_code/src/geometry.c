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


/* proceed a rotation of a 3D vector define in a reference frame A into a reference frame B with a quaternion q representing the orentation of B relative to A */
void om_rotate_vector_quaternion(struct omQuaternion *q,struct omVector *in,struct omVector *out){

	// variables
	omQuaternion q_inv;
	omQuaternion q_vec;
	omQuaternion q_tmp;

	// get inverse of q
	om_quat_inverse(q,&q_inv);

	// transform the vector into quaternion
	om_quat_create(&q_vec,0.0,in->_values[0],in->_values[1],in->_values[2]);

	// compute v_b = imaginary( q^{-1} * q_va * q)
	om_operator_quat_mul(&q_inv,&q_vec,&q_tmp);
	om_operator_quat_mul(&q_tmp,q,&q_vec);

	// set output
	om_vector_setValues(out,3,q_vec._qx,q_vec._qy,q_vec._qz);

}

/* proceed a rotation of a 3D vector define in a reference frame A into a reference frame B with a rotation matrix M representing the orentation of B relative to A */
void om_rotate_vector_matrix(struct omMatrix *M,struct omVector *in,struct omVector *out){

	// compute v_b = M * v_a
	om_operator_matrix_vector_mul(M,in,out);

}

/* proceed a rotation of a 3D vector define in a reference frame A into a reference frame B with an axis n and an angle a representing the orentation of B relative to A */
void om_rotate_vector_axisAngle(struct omAxisAngle *aa,struct omVector *in,struct omVector *out){

	// convert the axis angle into a quaternion
	omQuaternion q;
	om_convert_axisAngle2quaternion(aa,&q);

	// proceed rotation
	om_rotate_vector_quaternion(&q,in,out);

}

/* proceed the kinematics equation on a quaternion with an angular velocity */
void om_kinematics_quaternion(struct omQuaternion *q_t,struct omVector *angular_velocity,struct omQuaternion *q_tp1){

	// variable
	omMatrix Omega;
	omVector q_t_vec;
	omVector q_tp1_vec;

	// allocation of variables
	om_vector_create(&q_t_vec,4,q_t->_qx,q_t->_qy,q_t->_qz,q_t->_qw);
	om_vector_create(&q_tp1_vec,4);


	// compute transition matrix Omega
	om_matrix_create(&Omega,4,4);
	om_operator_omega_kinematics(angular_velocity,&Omega);

	// compute q_tp1 = Omega*q_t
	om_operator_matrix_vector_mul(&Omega,&q_t_vec,&q_tp1_vec);

	// set output
	om_quat_create(q_tp1,q_tp1_vec._values[3],q_tp1_vec._values[0],q_tp1_vec._values[1],q_tp1_vec._values[2]);

	// free memory
	om_matrix_free(&Omega);
	om_vector_free(&q_t_vec);
	om_vector_free(&q_tp1_vec);

}


///////////////////////////////////////////////////////
/////                Conversion                   /////
///////////////////////////////////////////////////////

/* convert a quaternion into a rotation matrix */
void om_convert_quaternion2matrix(struct omQuaternion *in,struct omMatrix* out){

	// variables
	omMatrix I;
	omMatrix S;
	omMatrix S_tmp;
	omVector in_img;

	// allocation
	om_matrix_create(&S,3,3);
	om_matrix_create(&S_tmp,3,3);
	om_matrix_createIdentity(&I,3);
	om_vector_create(&in_img,3);

	// get imaginary of q
	om_quat_imaginary(in,&in_img);

	// compute skew symmetric matrix of imaginary of q
	om_matrix_skewSymmetricMatrix(&in_img,&S);

	// compute 	R = I + S*(2.0*q_w) + S*S*2.0;
	om_operator_matrix_mul(&S,&S,&S_tmp);
	om_operator_matrix_scal_mul(&S_tmp,2.0,&S_tmp);
	om_operator_matrix_scal_mul(&S,2.0*in->_qw,&S);
	om_operator_matrix_add(&I,&S,out);
	om_operator_matrix_add(out,&S_tmp,out);

	// release memory
	om_matrix_free(&I);
	om_matrix_free(&S);
	om_matrix_free(&S_tmp);
	om_vector_free(&in_img);

}

/* convert a quaternion into euler angle representation */
void om_convert_quaternion2euler(struct omQuaternion *in,struct omEulerAngle* out){

	double sqx = in->_qx*in->_qx;
	double sqy = in->_qy*in->_qy;
	double sqz = in->_qz*in->_qz;
	double phi,theta,psy;

	double test = 2.0* ( (in->_qx*in->_qz) - (in->_qw*in->_qy)   );

	// if north pole singularity detected
	if ( fabs(test + 1.0) < EPSILON ) {

		psy =  atan2(  (in->_qy*in->_qx - in->_qz*in->_qw) , (in->_qx*in->_qz + in->_qw*in->_qy)  );
		theta = PI/2.0;
		phi = 0;

	// if south pole singularity detected
	} else if (fabs(test - 1.0) < EPSILON ) {

		psy =  - atan2(  (in->_qy*in->_qx - in->_qz*in->_qw) , (in->_qx*in->_qz + in->_qw*in->_qy)  );
		theta = -PI/2.0;
		phi = 0;

	}else{
		phi =  atan2( (in->_qx*in->_qw + in->_qy*in->_qz) , ( 0.5 - (sqx+sqy )));
		theta = asin(- 2.0* ( in->_qx*in->_qz - in->_qw*in->_qy   ) );
		psy =  atan2( (in->_qz*in->_qw + in->_qy*in->_qx) , ( 0.5 - (sqz+sqy )));

	}

	// set output in radian
	out->_pitch = phi;
	out->_roll = theta ;
	out->_yaw = psy;



}

/* convert a quaternion into axis angle representation  */
void om_convert_quaternion2axisAngle(struct omQuaternion *in,struct omAxisAngle *out){

	// if singularity detected
	if (  (fabs(in->_qw - 1.0) < EPSILON)  ){

		// set angle to 0
		out->_angle = 0.0;

		// the axis can't be defined in this situation.
		// then a default vector is set.
		om_vector_create(&out->_axis,3,1.0,0.0,0.0);

	}else{

		// axis n is the normalized imaginary part of the quaternion
		om_vector_create(&out->_axis,3);
		om_quat_imaginary(in,&out->_axis);

		// compute angle theta in radian such that
		// theta = 2*arctan2( ||n|| , q_w )
		out->_angle = 2.0*atan2( om_vector_norm(&out->_axis),in->_qw );

		// normalization
		om_vector_normalize(&out->_axis);


	}

}

/* convert a rotation matrix into quaternion representation  */
void om_convert_matrix2quaternion(struct omMatrix *in,struct omQuaternion *out){

	// variables
	double qw,qx,qy,qz;
	double tr;

	//get trace of M
	tr = om_matrix_trace(in);

	if (tr > 0) {

		double S = sqrt(tr + 1.0) * 2.0;
		qw = 0.25 * S;
		qx = (om_matrix_getValue(in,2,1) - om_matrix_getValue(in,1,2)) / S;
		qy = (om_matrix_getValue(in,0,2) - om_matrix_getValue(in,2,0)) / S;
		qz = (om_matrix_getValue(in,1,0) - om_matrix_getValue(in,0,1)) / S;

	}
	else if ( (om_matrix_getValue(in,0,0) > om_matrix_getValue(in,1,1)) && (om_matrix_getValue(in,0,0) > om_matrix_getValue(in,2,2)) ) {

	  float S = sqrt( 1.0 + om_matrix_getValue(in,0,0) - om_matrix_getValue(in,1,1) - om_matrix_getValue(in,2,2) ) * 2.0; // S=4*qx

	  qw = (om_matrix_getValue(in,2,1) - om_matrix_getValue(in,1,2)) / S;
	  qx = 0.25 * S;
	  qy = (om_matrix_getValue(in,0,1) + om_matrix_getValue(in,1,0)) / S;
	  qz = (om_matrix_getValue(in,0,2) + om_matrix_getValue(in,2,0)) / S;

	}
	else if ( ( om_matrix_getValue(in,1,1) > om_matrix_getValue(in,2,2)) ) {

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

	// set output
	om_quat_create(out,qw,qx,qy,qz);

}

/* convert an axis angle into quaternion representation  */
void om_convert_axisAngle2quaternion(struct omAxisAngle *in,struct omQuaternion *out){

	// variables
	double cos_a;
	double sin_a;

	cos_a = cos(in->_angle/2.0);
	sin_a = sin(in->_angle/2.0);

	// compute output such that
	// q_w = cos(angle/2)
	// q_x = n_x*sin(angle/2)
	// q_y = n_y*sin(angle/2)
	// q_z = n_z*sin(angle/2)
	om_quat_create(out,cos_a,in->_axis._values[0]*sin_a,in->_axis._values[1]*sin_a,in->_axis._values[2]*sin_a);

}

/* convert an axis angle into rotation matrix representation   */
void om_convert_axisAngle2matrix(struct omAxisAngle *in,struct omMatrix *out){

	// variables
	omMatrix I;
	omMatrix S;
	omMatrix S_tmp;
	double cos_a = cos(in->_angle);
	double sin_a = sin(in->_angle);

	// allocation
	om_matrix_createIdentity(&I,3);
	om_matrix_create(&S,3,3);
	om_matrix_create(&S_tmp,3,3);
	om_matrix_skewSymmetricMatrix(&in->_axis,&S);

	// compute R = I + S*sin(angle) + S*S*( 1.0 - cos(angle) )
	om_operator_matrix_mul(&S,&S,&S_tmp);
	om_operator_matrix_scal_mul(&S_tmp,(1.0 - cos_a),&S_tmp);
	om_operator_matrix_scal_mul(&S,sin_a,&S);
	om_operator_matrix_add(&S,&S_tmp,out);
	om_operator_matrix_add(out,&I,out);

	// free memory
	om_matrix_free(&I);
	om_matrix_free(&S);
	om_matrix_free(&S_tmp);

}

/* convert Euler angles into rotation matrix representation   */
void om_convert_euler2quaternion(struct omEulerAngle *in,struct omQuaternion *out){

	// angle in radian
	double phi = in->_pitch / 2.0;
	double theta = in->_roll / 2.0;
	double psy = in->_yaw / 2.0;

	double qw = (cos(phi)*cos(theta)*cos(psy)) + (sin(phi)*sin(theta)*sin(psy));
	double qx = (sin(phi)*cos(theta)*cos(psy)) - (cos(phi)*sin(theta)*sin(psy));
	double qy = (cos(phi)*sin(theta)*cos(psy)) + (sin(phi)*cos(theta)*sin(psy));
	double qz= (cos(phi)*cos(theta)*sin(psy)) - (sin(phi)*sin(theta)*cos(psy));

	// set output
	om_quat_create(out,qw,qx,qy,qz);

}




///////////////////////////////////////////////////////
/////                Operators                    /////
///////////////////////////////////////////////////////

/* operator vex */
void om_operator_vex(struct omMatrix *in,struct omVector *out){

	// compute vex( [ a x ] ) = a
	// where [ a x ] is skew symmetric matrix of a
	double x = om_matrix_getValue(in,2,1);
	double y = om_matrix_getValue(in,0,2);
	double z = om_matrix_getValue(in,1,0);

	// set output
	om_vector_create(out,3,x,y,z);

}

/* operator omega */
void om_operator_omega(struct omVector *in,struct omMatrix *out){

	// set values such as:
	//
	// v = ( v_x , v_y , v_z )^T
	//
	//         /  0.0   v_z -v_y  v_x \
	// Omega = | -v_z   0.0  v_x  v_y |
	//         |  v_y  -v_x  0.0  v_z |
	//         \ -v_x  -v_y -v_z  0.0 /
	//
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

/* operator gamma */
void om_operator_gamma(struct omVector *in,struct omMatrix *out){

	// set values such as:
	//
	// v = ( v_x , v_y , v_z )^T
	//
	//         /  0.0  -v_z  v_y  v_x \
	// Gamma = |  v_z   0.0 -v_x  v_y |
	//         | -v_y   v_x  0.0  v_z |
	//         \ -v_x  -v_y -v_z  0.0 /
	//
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

/* operator xi */
void om_operator_xi(struct omQuaternion *in,struct omMatrix *out){

	// set values such as:
	//
	// q = ( q_w , q_x , q_y , q_z )^T
	//
	//         /  q_w  -q_z  q_y \
	//   Xi  = |  q_z   q_w -q_x |
	//         | -q_y   q_x  q_w |
	//         \ -q_x  -q_y -q_z /
	//
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

/* operator psy */
void om_operator_psy(struct omQuaternion *in,struct omMatrix *out){

	// set values such as:
	//
	// q = ( q_w , q_x , q_y , q_z )^T
	//
	//         /  q_w   q_z -q_y \
	//   Xi  = | -q_z   q_w  q_x |
	//         |  q_y  -q_x  q_w |
	//         \ -q_x  -q_y -q_z /
	//
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

/* compute transition matrix for kinematics quaternion equation */
void om_operator_omega_kinematics(struct omVector *angular_velocity,struct omMatrix *out){

	// variables
	omMatrix Omega_tmp;
	omMatrix S_phi;
	omMatrix I;
	omVector phi;

	// allocation
	om_matrix_create(&Omega_tmp,3,3);
	om_matrix_create(&S_phi,3,3);
	om_matrix_createIdentity(&I,3);
	om_vector_create(&phi,3);


	// set values such as:
	//
	// w = ( w_x , w_y , w_z )^T
	// phy = cos(0.5 ||w|| Delta_T)*(w/||w||)
	//
	// [ w x ] is the skew symmetric matrix of w
	//
	//         / cos(0.5 ||w|| Delta_T)*I_3 - [ w x ]              phy             \
	//   M  =  \         - phy^T                          cos(0.5 ||w|| Delta_T)   /
	//
	double norm = om_vector_norm(angular_velocity);
	double cos_tmp = cosf(0.5*norm*DELTA_T);
	double sin_tmp = sinf(0.5*norm*DELTA_T);

	om_vector_clone(angular_velocity,&phi);
	om_operator_vector_scal_mul(&phi,( sin_tmp/norm ),&phi);
	om_matrix_skewSymmetricMatrix(&phi,&S_phi);
	om_operator_matrix_scal_mul(&I,cos_tmp,&I);
	om_operator_matrix_sub(&I,&S_phi,&Omega_tmp);

	for(int i=0;i<3;++i){
		om_matrix_setValue(out,i,3,phi._values[i]);
		om_matrix_setValue(out,3,i,-phi._values[i]);

		for(int j=0;j<3;++j)
			om_matrix_setValue(out,i,j,om_matrix_getValue(&Omega_tmp,i,j));
	}
	om_matrix_setValue(out,3,3,cos_tmp);

	// free memory
	om_matrix_free(&Omega_tmp);
	om_matrix_free(&I);
	om_matrix_free(&S_phi);
	om_vector_free(&phi);

}








