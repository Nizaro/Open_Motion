/**
 * \file algebra.h
 * \author Thomas BRAUD, Nizar OUARTI
 * \date 10 june 2016
 * \brief File containing linear algebra methods
 *
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

#define CENTRAL_RANGE 0.7

// PI Variable
#define PI 3.14159265359

// Gravitational constant.
#define G 9.80665

// Convert radian to degree
#define RAD_TO_DEG 180.0/PI

// Convert degree to radian
#define DEG_TO_RAD PI/180.0

#define EPSILON 0.000000001


// Time variation
extern double DELTA_T;

///////////////////////////////////////////////////////
/////             Vector class                    /////
///////////////////////////////////////////////////////

/**
 * \struct omVector
 * \brief A vector of size n containing real values
 *
 */
typedef struct omVector{

	int _length;
	double* _values;

}omVector;

/**
 * \fn void om_vector_create(struct omVector *vector,int size,...)
 * \brief create a new vector with a fixed length
 * \param vector the new length of the vector
 * \param size the new length of the vector
 *
 */
void om_vector_create(struct omVector *vector,int size,...);

/**
 * \fn void om_vector_clone(struct omVector *in,struct omVector *out)
 * \brief copy all values of a vector into a new one
 * \param in the vector to copy
 * \param out the new vector
 *
 */
void om_vector_clone(struct omVector *in,struct omVector *out);

/**
 * \fn void om_vector_setValue(struct omVector *vector,int index,double value)
 * \brief change a value of a vector
 * \param vector the vector
 * \param index the index of the value to change
 * \param value the new value
 *
 */
void om_vector_setValue(struct omVector *vector,int index,double value);

/**
 * \fn om_vector_setValue(struct omVector *vector,int index,double value)
 * \brief change all values of a vector
 * \param vector the vector
 * \param size the new length of the vector
 *
 */
void om_vector_setValues(struct omVector *vector,int size,...);

/**
 * \fn double om_vector_getValue(struct omVector *vector,int index)
 * \brief get a value of a vector
 * \param vector the vector
 * \param index the index of the value to get
 *
 */
double om_vector_getValue(struct omVector *vector,int index);


/**
 * \fn double om_vector_getValue(struct omVector *vector,int index)
 * \brief get a value of a vector
 * \param vector the vector
 * \param index the index of the value to get
 *
 */
int om_vector_getLength(struct omVector *vector);

/**
 * \fn double om_vector_rms(struct omVector *vector)
 * \brief get the root mean square of the vector
 * \param vector the vector
 *
 */
double om_vector_rms(struct omVector *vector);

/**
 * \fn double om_vector_norm(struct omVector *vector)
 * \brief get the euclidian norm of the vector
 * \param vector the vector
 *
 */
double om_vector_norm(struct omVector *vector);

/**
 * \fn void om_vector_normalize(struct omVector *vector)
 * \brief normalize a vector
 * \param vector the vector
 *
 */
void om_vector_normalize(struct omVector *vector);

/**
 * \fn void om_vector_normalize(struct omVector *vector)
 * \brief display in a console the values of the vector
 * \param vector the vector
 *
 */
void om_vector_display(struct omVector *vector);

/**
 * \fn double om_vector_mean(struct omVector *vector)
 * \brief display in a console the values of the vector
 * \param vector the vector
 *
 */
double om_vector_mean(struct omVector *vector);

/**
 * \fn double om_vector_median(struct omVector *vector)
 * \brief display in a console the values of the vector
 * \param vector the vector
 *
 */
double om_vector_median(struct omVector *vector);



/**
 * \fn void om_vector_free(struct omVector *vector)
 * \brief free the memory allocated by the vector
 * \param vector the vector
 *
 */
double om_vector_interpolation(omVector *x,omVector *y,double xq);


/**
 * \fn void om_vector_free(struct omVector *vector)
 * \brief free the memory allocated by the vector
 * \param vector the vector
 *
 */
void om_vector_free(struct omVector *vector);

/**
 * \fn void om_vector_crossProduct(struct omVector *vector)
 * \brief free the memory allocated by the vector
 * \param vector the vector
 *
 */
void om_vector_crossProduct(struct omVector *a,struct omVector *b,struct omVector *cross);

/**
 * \fn void om_vector_dotProduct(struct omVector *vector)
 * \brief free the memory allocated by the vector
 * \param vector the vector
 *
 */
double om_vector_dotProduct(struct omVector *a,struct omVector *b);

///////////////////////////////////////////////////////
/////             Matrix class                    /////
///////////////////////////////////////////////////////

/**
 * \struct omMatrix
 * \brief A Matrix of size n x m containing real values
 *
 */
typedef struct omMatrix{

	int _columns;
	int _rows;
	double** _values;

}omMatrix;


/**
 * \fn void om_matrix_create(struct omMatrix *matrix,int rows,int columns)
 * \brief create a matrix with n rows and m columns
 * \param matrix the matrix
 * \param rows the number of rows
 * \param columns the number of columns
 *
 */
void om_matrix_create(struct omMatrix *matrix,int rows,int columns);

/**
 * \fn double om_matrix_getValue(struct omMatrix *matrix,int i,int j)
 * \brief get the (i x j)th value of a matrix
 * \param matrix the matrix
 * \param i the row index of the value to get
 * \param j the column index of the value to get
 *
 */
double om_matrix_getValue(struct omMatrix *matrix,int i,int j);

/**
 * \fn void om_matrix_getColumn(struct omMatrix *matrix,int column,struct omVector *out)
 * \brief get the jth column of a matrix
 * \param matrix the matrix
 * \param column the column index to get
 * \param out the column
 *
 */
void om_matrix_getColumn(struct omMatrix *matrix,int column,struct omVector *out);

/**
 * \fn void om_matrix_getColumn(struct omMatrix *matrix,int column,struct omVector *out)
 * \brief get the ith row of a matrix
 * \param matrix the matrix
 * \param row the row index to get
 * \param out the row
 *
 */
void om_matrix_getRow(struct omMatrix *matrix,int row,struct omVector *out);

/**
 * \fn void om_matrix_setValue(struct omMatrix *matrix,int i,int j,double value)
 * \brief set the (i x j)th value of a matrix
 * \param matrix the matrix
 * \param i the row index of the value to change
 * \param j the column index of the value to change
 * \param value the new value
 *
 */
void om_matrix_setValue(struct omMatrix *matrix,int i,int j,double value);


/**
 *
 */
void om_matrix_convolution2D(omMatrix *A,omMatrix *B,omMatrix *C);


/**
 *
 */
void om_matrix_convolution2D_valid(omMatrix *A,omMatrix *B,omMatrix *C);


/**
 * \fn void om_matrix_setColumn(struct omMatrix *matrix,int column,struct omVector *in)
 * \brief set the jth column of a matrix
 * \param matrix the matrix
 * \param j the column index of the value to change
 * \param in the new values of the column
 *
 */
void om_matrix_setColumn(struct omMatrix *matrix,int column,struct omVector *in);

/**
 * \fn void om_matrix_setRow(struct omMatrix *matrix,int row,struct omVector *in)
 * \brief set the ith row of a matrix
 * \param matrix the matrix
 * \param i the column index of the value to change
 * \param in the new values of the column
 *
 */
void om_matrix_setRow(struct omMatrix *matrix,int row,struct omVector *in);


/**
 * \fn double om_matrix_norm(struct omMatrix *matrix)
 * \brief compute the norm of a square matrix
 * \param matrix the matrix
 * \return the norm of the matrix
 */
double om_matrix_norm(struct omMatrix *matrix);


/**
 * \fn double om_matrix_determinant(struct omMatrix *matrix)
 * \brief compute the determinant of a square matrix
 * \param matrix the matrix
 * \return the determinant of the matrix
 */
double om_matrix_determinant(struct omMatrix *matrix);

/**
 * \fn double om_matrix_trace(struct omMatrix *matrix)
 * \brief compute the trace (sum of diagonal values) of a square matrix
 * \param matrix the matrix
 * \return the norm of the matrix
 */
double om_matrix_trace(struct omMatrix *matrix);

/**
 * \fn void om_matrix_exponantial(struct omMatrix *matrix,struct omMatrix *m_exp,int N)
 * \brief compute the exponantial of a square matrix
 * \param matrix the matrix
 * \param m_exp the exponential matrix
 * \param N number of iteration
 */
void om_matrix_exponential(struct omMatrix *matrix,struct omMatrix *m_exp,int N);

/**
 * \fn void om_matrix_squareRoot(struct omMatrix *matrix,struct omMatrix *m_sqrt)
 * \brief compute the square root of a square matrix
 * \param matrix the matrix
 * \param m_sqrt the square root matrix
 */
void om_matrix_squareRoot(struct omMatrix *matrix,struct omMatrix *m_sqrt,int N);

/**
 * \fn int om_matrix_isSquare(struct omMatrix *matrix)
 * \brief determine is a matrix is square or not. The function will test if the number of rows is equals to the number od
 * \param matrix the matrix
 * \return 1 if the matrix is square
 */
int om_matrix_isSquare(struct omMatrix *matrix);

/**
 * \fn int om_matrix_isSquare(struct omMatrix *matrix)
 * \brief determine is a matrix is square or not. The function will test if the number of rows is equals to the number od
 * \param matrix the matrix
 * \return 1 if the matrix is square
 */
int om_matrix_isSymmetric(struct omMatrix *matrix);



/**
 * \fn int om_matrix_containsNaN(struct omMatrix *matrix)
 * \brief determine is a matrix is square or not. The function will test if the number of rows is equals to the number od
 * \param matrix the matrix
 * \return 1 if the matrix is square
 */
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
void om_matrix_schurDecomposition(struct omMatrix *matrix,struct omMatrix *T,struct omMatrix *U,int N);

void om_matrix_display(struct omMatrix *matrix);
void om_matrix_free(struct omMatrix *matrix);
void om_matrix_clone(struct omMatrix *in,struct omMatrix *out);

void om_matrix_getEingenValues(struct omMatrix *matrix,struct omVector **eigen_vectors,double **eigen_values,int N);
void om_matrix_skewSymmetricMatrix(struct omVector *in,struct omMatrix *out);
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
void om_quat_clone(struct omQuaternion *in,struct omQuaternion *out);

///////////////////////////////////////////////////////
/////              Operators                      /////
///////////////////////////////////////////////////////

void om_operator_vector_add(struct omVector *a,struct omVector *b,struct omVector *out);
void om_operator_vector_sub(struct omVector *a,struct omVector *b,struct omVector *out);
void om_operator_vector_scal_mul(struct omVector *a,double b,struct omVector *out);
void om_operator_vector_scal_div(struct omVector *a,double b,struct omVector *out);
void om_operator_vector_outer_product (struct omVector* a, struct omVector* b, struct omMatrix* out);

void om_operator_matrix_add(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out);
void om_operator_matrix_sub(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out);
void om_operator_matrix_mul(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out);
void om_operator_matrix_scal_mul(struct omMatrix *a,double b,struct omMatrix *out);
void om_operator_matrix_scal_div(struct omMatrix *a,double b,struct omMatrix *out);
void om_operator_matrix_vector_mul(struct omMatrix *a,struct omVector *b,struct omVector *out);

void om_operator_quat_add(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out);
void om_operator_quat_sub(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out);
void om_operator_quat_mul(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out);
void om_operator_quat_scal_mul(struct omQuaternion *a,double b,struct omQuaternion *out);
void om_operator_quat_scal_div(struct omQuaternion *a,double b,struct omQuaternion *out);


///////////////////////////////////////////////////////
/////                Divers                       /////
///////////////////////////////////////////////////////

double om_maths_erfinv( double y);


void om_convert_vector2matrix(struct omVector* a, struct omMatrix* out);

double om_quat_dotProduct(struct omQuaternion *a,struct omQuaternion *b);

void om_solvingLinearSystem(struct omMatrix *A,struct omVector *b,struct omVector *x);
void om_solvingLinearSystemLU(struct omMatrix *L,struct omMatrix *U,struct omVector *b,struct omVector *x);

void om_least_square_method (struct omMatrix *pX,struct omVector *pY,struct omVector *pBeta);
double simpsonadapt(double (*fnct)(double), double a,double b, double mid, double epsilon, double maxh,double minh, double fa, double fb, double fmid, double *bada, double *badb,int *success);
void cdiv(double xr, double xi, double yr, double yi,double* cdivr,double* cdivi);

#endif /* ALGEBRA_H_ */
