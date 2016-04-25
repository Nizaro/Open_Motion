

/**
 * \file main.c
 * \brief Contains all geometric functions
 * \author Thomas BRAUD
 * \version 0.1
 * \date 11/01/2016
 *
 *	This file contains all functions simulating euclidian geometry. Conversion between attitude representation, and rotation kinematics.
 *
 */



#ifndef GEOMETRY_H_
#define GEOMETRY_H_


#include "random.h"

// namespace isf
namespace om {

// pointer to function with a vector in parameter
typedef Vector (*numerical_function)(const Vector& x);


/**
 * \struct NumericalFunction
 * \brief structure representing a numerical function
 *
 */
struct NumericalFunction{

	// dimension of the output vector
	int m;

	// function pointer
	numerical_function function;

};

/**
 * \struct AxisAngle
 * \brief Axis-Angle representation
 *
 * Axis-Angle representation, dispose an angle and a 3D unit vector representing the axis of rotation.
 */
struct AxisAngle{

	//angle in radian
	double angle;

	// 3D unit vector
	Vector axis;

};





///////////////////////////////////////////////////////
/////               Matrix tools                  /////
///////////////////////////////////////////////////////


/**
 * \fn Matrix identity(int n)
 * \brief create the identity matrix of size n
 * \author Thomas BRAUD
 *
 * \param n size of the square matrix
 * \return identity matrix of size n
 */
Matrix identity(int n);

/**
 * \fn Matrix skewSymetricMatrix(const Vector& x)
 * \brief create the skew symetric matrix associate to the 3D vector x
 * \author Thomas BRAUD
 *
 * \param x the 3D vector
 * \return skew symetric matrix associate to the 3D vector x
 */
Matrix skewSymetricMatrix(const Vector& x);

/**
 * \fn Matrix jacobianMatrix(Vector x, NumericalFunction f )
 * \brief compute the jacobian matrix of f(x)
 * \author Thomas BRAUD
 *
 * \param x the 3D vector
 * \return the jacobian matrix
 */
Matrix jacobianMatrix(Vector x, NumericalFunction f );


/**
 * \fn Matrix rotationMatrixZYX(double phi,double theta,double psy)
 * \brief compute the rotation matrix from the euler angle phi, theta, psy
 * \author Thomas BRAUD
 *
 * \param phi pitch
 * \param theta roll
 * \param psy yaw
 * \return the rotation matrix
 */
Matrix rotationMatrixZYX(double phi,double theta,double psy);

/**
 * \fn Matrix rotationMatrixXYZ(double phi,double theta,double psy)
 * \brief compute the rotation matrix from the euler angle phi, theta, psy
 * \author Thomas BRAUD
 *
 * \param phi pitch
 * \param theta roll
 * \param psy yaw
 * \return the rotation matrix
 */
Matrix rotationMatrixXYZ(double phi,double theta,double psy);

/**
 * \fn Matrix matExponential(const Matrix& A)
 * \brief compute the exponential matrix from the matrix A
 * \author Thomas BRAUD
 *
 * \param A the matrix, must be square
 * \return the rotation matrix
 */
Matrix matExponential(const Matrix& A);

/**
 * \fn Matrix matSquareRoot1x1(Matrix A)
 * \brief compute the square root matrix from the matrix A
 * \author Thomas BRAUD
 *
 * \param A the matrix of size 1x1
 * \return the square root matrix
 */
Matrix matSquareRoot1x1(Matrix A);

/**
 * \fn Matrix matSquareRoot2x2(Matrix A)
 * \brief compute the square root matrix from the matrix A
 * \author Thomas BRAUD
 *
 * \param A the matrix of size 2x2, must be semi-definite positif
 * \return the square root matrix
 */
Matrix matSquareRoot2x2( Matrix A);

/**
 * \fn Matrix matSquareRoot2x2(Matrix A)
 * \brief compute the square root matrix from the matrix A
 * \author Thomas BRAUD
 *
 * \param A the matrix of size nxn with n>2, must be semi-definite positif
 * \return the square root matrix
 */
Matrix matSquareRoot(Matrix A);


/**
 * \fn Matrix operatorOmega(const Vector& v)
 * \brief
 * \author Thomas BRAUD
 *
 * \param v a 3D vector
 * \return a matrix
 */
Matrix operatorOmega(const Vector& v);

/**
 * \fn Matrix operatorGamma(const Vector& v)
 * \brief
 * \author Thomas BRAUD
 *
 * \param v a 3D vector
 * \return a matrix
 */
Matrix operatorGamma(const Vector& v);

/**
 * \fn Matrix operatorXi(Quaternion q)
 * \brief
 * \author Thomas BRAUD
 *
 * \param q a quaternion
 * \return a matrix
 */
Matrix operatorXi(Quaternion q);

/**
 * \fn Matrix operatorPsy(Quaternion q)
 * \brief
 * \author Thomas BRAUD
 *
 * \param q a quaternion
 * \return a matrix
 */
Matrix operatorPsy(Quaternion q);



Matrix transitionOmega(const Vector& gyro);


/**
 * \fn Matrix operatorPsy(Quaternion q)
 * \brief compute the eigen values and their associated eigen vector of a matrix
 * \author Thomas BRAUD
 *
 * \param A the matrix
 * \return an array of real representing the eigen values and their associated eigen vector
 */
void getEigenValuesVector(const Matrix& A,vector<double>& eigenValues,vector<Vector>& eigenVector,int N);



///////////////////////////////////////////////////////
/////               Vector tools                  /////
///////////////////////////////////////////////////////


/**
 * \fn Matrix operatorPsy(Quaternion q)
 * \brief compute the eigen values and their associated eigen vector of a matrix
 * \author Thomas BRAUD
 *
 * \param A the matrix
 * \return an array of real representing the eigen values and their associated eigen vector
 */

Vector operatorVex(const Matrix& S);


/**
 * \fn Matrix operatorPsy(Quaternion q)
 * \brief compute the eigen values and their associated eigen vector of a matrix
 * \author Thomas BRAUD
 *
 * \param A the matrix
 * \return an array of real representing the eigen values and their associated eigen vector
 */

Vector operatorVex(const Vector& a,const Vector& b);


/**
 * \fn Matrix operatorPsy(Quaternion q)
 * \brief compute the eigen values and their associated eigen vector of a matrix
 * \author Thomas BRAUD
 *
 * \param A the matrix
 * \return an array of real representing the eigen values and their associated eigen vector
 */

Matrix vectorToMatrix(const Vector& v);

Vector projection(Vector v1, Vector v2 );

Matrix outerProduct(const Vector& v1,const Vector& v2 );

Vector crossProduct(const Vector& v1,const Vector& v2 );

double dotProduct(const Vector& v1,const Vector& v2 );

double euclidianDistance(const Vector& v1,const Vector& v2 );


Vector rotationVector(const Vector& v,Quaternion q);

///////////////////////////////////////////////////////
/////               Angles tools                  /////
///////////////////////////////////////////////////////


double angularSum(double alpha,double beta,double coeff_alpha,double coeff_beta);

double angularDistance(const Vector& v1,const Vector& v2 );

Vector angularVelocity(Quaternion q_t,Quaternion q_tm1);


///////////////////////////////////////////////////////
/////             Conversion tools                /////
///////////////////////////////////////////////////////

AxisAngle quatToAxisAngle(Quaternion q);

Vector rotationMatrixToEuler (const Matrix& m);

Matrix quatToRotationMatrix (Quaternion q);

AxisAngle eulerToAxisAngle( double phi, double theta, double psy);

Vector quatToEuler(Quaternion q);

Matrix axisAngleToRotationMatrix (double angle,const Vector& axis);

Matrix axisAngleToRotationMatrix (AxisAngle axis_angle);

Quaternion rotationMatrixToQuat(Matrix m);


///////////////////////////////////////////////////////
/////               Quaternion tools              /////
///////////////////////////////////////////////////////

Quaternion quatExponential(Quaternion q);

Quaternion quatLogarithm(Quaternion q);

Quaternion quatPower(Quaternion q,double p);

double dotProduct(Quaternion q1,Quaternion q2);

Quaternion slerp(Quaternion q1,Quaternion q2,double t);

Quaternion nlerp(Quaternion q1,Quaternion q2,double t);

///////////////////////////////////////////////////////
/////               Others tools                  /////
///////////////////////////////////////////////////////


double calculErrorOrientation(Quaternion q_real,Quaternion q_est);




} /* namespace isf */

#endif /* GEOMETRY_H_ */
