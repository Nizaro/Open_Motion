/**
 * \file geometry.h
 * \author Thomas BRAUD, Nizar OUARTI
 * \date 10 june 2016
 * \brief File containing linear algebra methods
 *
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "random.h"

///////////////////////////////////////////////////////
/////                Structure                   /////
///////////////////////////////////////////////////////

typedef struct omAxisAngle{

	double _angle;
	omVector _axis;

}omAxisAngle;


typedef struct omEulerAngle{

	double _pitch;
	double _roll;
	double _yaw;

}omEulerAngle;


///////////////////////////////////////////////////////
/////           Kinematics - Rotation             /////
///////////////////////////////////////////////////////


void om_rotate_vector_quaternion(struct omQuaternion *q,struct omVector *in,struct omVector *out);
void om_rotate_vector_matrix(struct omMatrix *M,struct omVector *in,struct omVector *out);
void om_rotate_vector_axisAngle(struct omAxisAngle *aa,struct omVector *in,struct omVector *out);

void om_kinematics_quaternion(struct omQuaternion *q_t,struct omVector *angular_velocity,struct omQuaternion *q_tp1);


///////////////////////////////////////////////////////
/////                Conversion                   /////
///////////////////////////////////////////////////////

void om_convert_quaternion2matrix(struct omQuaternion *in,struct omMatrix* out);
void om_convert_quaternion2euler(struct omQuaternion *in,struct omEulerAngle* out);
void om_convert_quaternion2axisAngle(struct omQuaternion *in,struct omAxisAngle *out);

void om_convert_matrix2quaternion(struct omMatrix *in,struct omQuaternion *out);

void om_convert_axisAngle2quaternion(struct omAxisAngle *in,struct omQuaternion *out);
void om_convert_axisAngle2matrix(struct omAxisAngle *in,struct omMatrix *out);

void om_convert_euler2quaternion(struct omEulerAngle *in,struct omQuaternion *out);

///////////////////////////////////////////////////////
/////                Operators                    /////
///////////////////////////////////////////////////////

void om_operator_omega_kinematics(struct omVector *angular_velocity,struct omMatrix *out);
void om_operator_vex(struct omMatrix *in,struct omVector *out);
void om_operator_omega(struct omVector *in,struct omMatrix *out);
void om_operator_gamma(struct omVector *in,struct omMatrix *out);
void om_operator_xi(struct omQuaternion *in,struct omMatrix *out);
void om_operator_psy(struct omQuaternion *in,struct omMatrix *out);

#endif /* GEOMETRY_H_ */
