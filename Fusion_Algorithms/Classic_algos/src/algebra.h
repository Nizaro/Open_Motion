/**
 * @file algebra.h
 * @author Thomas BRAUD, Nizar OUARTI
 * @date 10 june 2016
 * @brief File containing linear algebra methods
 *
 * Here typically goes a more extensive explanation of what the header
 * defines. Doxygens tags are words preceeded by either a backslash @\
 */

#ifndef ALGEBRA_H_
#define ALGEBRA_H_


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <float.h>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif


#if defined(_OPENCL)
#include <CL/cl.h>
#endif


///////////////////////////////////////////////////////
/////            Global constants                 /////
///////////////////////////////////////////////////////

// PI Variable
#define PI 3.14159265359

// Time variation
#define DELTA_T 0.001

// Gravitational constant.
#define G 9.81

// Convert radian to degree
#define RAD_TO_DEG 180.0/PI

// Convert degree to radian
#define DEG_TO_RAD PI/180

#define EPSILON 0.000000001




///////////////////////////////////////////////////////
/////             Vector class                    /////
///////////////////////////////////////////////////////

/**
 * @brief Use brief, otherwise the index won't have a brief explanation.
 *
 * A vector of size n containing
 */
typedef struct omVector{

	int _length;
	double* _values;

}omVector;

/**
 * update the length of the vector (normally it will be unecessary
 * @param length the new length of the vector
 */
void om_vector_create(struct omVector *vector,int size,...);
void om_vector_clone(struct omVector *in,struct omVector *out);
void om_vector_setValue(struct omVector *vector,int index,double value);
void om_vector_setValues(struct omVector *vector,int size,...);
double om_vector_getValue(struct omVector *vector,int index);
double om_vector_rms(struct omVector *vector);
double om_vector_norm(struct omVector *vector);
void om_vector_normalize(struct omVector *vector);
void om_vector_display(struct omVector *vector);
void om_vector_dispose(struct omVector *vector);


///////////////////////////////////////////////////////
/////             Matrix class                    /////
///////////////////////////////////////////////////////

typedef struct omMatrix{

	int _columns;
	int _rows;
	double** _values;

}omMatrix;

void om_matrix_create(struct omMatrix *matrix,int rows,int columns);

double om_matrix_getValue(struct omMatrix *matrix,int i,int j);
void om_matrix_getColumn(struct omMatrix *matrix,int column,struct omVector *out);
void om_matrix_getRow(struct omMatrix *matrix,int row,struct omVector *out);

void om_matrix_setValue(struct omMatrix *matrix,int i,int j,double value);
void om_matrix_setColumn(struct omMatrix *matrix,int column,struct omVector *in);
void om_matrix_setRow(struct omMatrix *matrix,int row,struct omVector *in);

double om_matrix_norm(struct omMatrix *matrix);
double om_matrix_determinant(struct omMatrix *matrix);
double om_matrix_trace(struct omMatrix *matrix);

void om_matrix_exponantial(struct omMatrix *matrix,struct omMatrix *m_exp,int N);
void om_matrix_squareRoot(struct omMatrix *matrix,struct omMatrix *m_sqrt);

int om_matrix_isSquare(struct omMatrix *matrix);
int om_matrix_containsNaN(struct omMatrix *matrix);
int om_matrix_isNull(struct omMatrix *matrix);

void om_matrix_submatrix(struct omMatrix *matrix,struct omMatrix *submatrix,int row, int column,int n, int m);
void om_matrix_comatrix(struct omMatrix *matrix,struct omMatrix *comatrix);
void om_matrix_transpose(struct omMatrix *matrix,struct omMatrix *transpose);
void om_matrix_inverse(struct omMatrix *matrix,struct omMatrix *inverse);
void om_matrix_adjugate(struct omMatrix *matrix,struct omMatrix *adjugate);
void om_matrix_factorizationQR(struct omMatrix *matrix,struct omMatrix *Q,struct omMatrix *R);
void om_matrix_factorizationLU(struct omMatrix *matrix,struct omMatrix *L,struct omMatrix *U);
void om_matrix_choleskyDecomposition(struct omMatrix *matrix,struct omMatrix *L);
void om_matrix_schurDecomposition(struct omMatrix *matrix,struct omMatrix *T,struct omMatrix *U);

void om_matrix_display(struct omMatrix *matrix);
void om_matrix_dispose(struct omMatrix *matrix);
void om_matrix_clone(struct omMatrix *in,struct omMatrix *out);

void om_matrix_getEingenValues(struct omMatrix *matrix,struct omVector **eigen_vectors,double **eigen_values,int N);
void om_matrix_skewSymetricMatrix(struct omVector *in,struct omMatrix *out);
void om_matrix_createIdentity(struct omMatrix *I,int n);

///////////////////////////////////////////////////////
/////             Quaternion class                /////
///////////////////////////////////////////////////////


typedef struct omQuaternion {

	double _qx;
	double _qy;
	double _qz;
	double _qw;

}omQuaternion;


void om_quat_create(struct omQuaternion *quat,double qw,double qx,double qy,double qz);
void om_quat_conjugate(struct omQuaternion *quat,struct omQuaternion *conjugate);
void om_quat_inverse(struct omQuaternion *quat,struct omQuaternion *inverse);
void om_quat_imaginary(struct omQuaternion *quat,struct omVector *imaginary);
void om_quat_normalize(struct omQuaternion *quat);
double om_quat_norm(struct omQuaternion *quat);
void om_quat_display(struct omQuaternion *quat);


///////////////////////////////////////////////////////
/////              Operators                      /////
///////////////////////////////////////////////////////

void om_operator_vector_add(struct omVector *a,struct omVector *b,struct omVector *out);
void om_operator_vector_sub(struct omVector *a,struct omVector *b,struct omVector *out);
void om_operator_vector_const_mul(struct omVector *a,double b,struct omVector *out);
void om_operator_vector_const_div(struct omVector *a,double b,struct omVector *out);

void om_operator_matrix_add(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out);
void om_operator_matrix_sub(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out);
void om_operator_matrix_mul(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out);
void om_operator_matrix_const_mul(struct omMatrix *a,double b,struct omMatrix *out);
void om_operator_matrix_const_div(struct omMatrix *a,double b,struct omMatrix *out);
void om_operator_matrix_vector_mul(struct omMatrix *a,struct omVector *b,struct omVector *out);

void om_operator_quat_add(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out);
void om_operator_quat_sub(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out);
void om_operator_quat_mul(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out);
void om_operator_quat_const_mul(struct omQuaternion *a,double b,struct omQuaternion *out);
void om_operator_quat_const_div(struct omQuaternion *a,double b,struct omQuaternion *out);




///////////////////////////////////////////////////////
/////                Divers                       /////
///////////////////////////////////////////////////////

void om_vector_crossProduct(struct omVector *a,struct omVector *b,struct omVector *cross);
double om_vector_dotProduct(struct omVector *a,struct omVector *b);
double om_quat_dotProduct(struct omQuaternion *a,struct omQuaternion *b);

void om_solvingLinearSystem(struct omMatrix *A,struct omVector *b,struct omVector *x);
void om_solvingLinearSystemLU(struct omMatrix *L,struct omMatrix *U,struct omVector *b,struct omVector *x);

#endif /* ALGEBRA_H_ */
