#include "algebra.h"

double DELTA_T = 0;

///////////////////////////////////////////////////////
/////             Vector class                    /////
///////////////////////////////////////////////////////

/* initialization of a vector */
void om_vector_create(struct omVector *vector,int length,...){

	// set length
	vector->_length = length;

	// allocation of the array containing vector's values
	vector->_values = (double*)malloc(length*sizeof(double));

	// read the arguments
	va_list args;
	va_start ( args, length );

	//for each arguments
	for ( int i = 0; i < vector->_length ; i++ ){

		// read argument's value
		double value =  va_arg ( args, double );

		//set value in the vector
		if(value)
			vector->_values[i] = value;
		else
			vector->_values[i] = 0.0;
	}

	va_end(args);
}

/* set all values of a vector */
void om_vector_setValues(struct omVector *vector,int length,...){

	// read the arguments
	va_list args;
	va_start ( args, length );

	//for each arguments
	for ( int i = 0; i < vector->_length ; i++ ){

		// read argument's value
		double value =  va_arg ( args, double );

		//set value in the vector
		if(value)
			vector->_values[i] = value;
		else
			vector->_values[i] = 0.0;
	}

	va_end(args);

}

/* set the ith value of a vector */
void om_vector_setValue(struct omVector *vector,int index,double value){

	//set value in the vector
	vector->_values[index]=value;
}


/* get the ith value of a vector */
double om_vector_getValue(struct omVector *vector,int index){

	//return the value
	return vector->_values[index];
}

/* get a copy of the vector */
void om_vector_clone(struct omVector *in,struct omVector *out){


	// set the values of the copy
	for ( int i = 0; i < in->_length ; i++ )
		om_vector_setValue(out,i,in->_values[i]);
}



/* compute the root mean square of the vector */
double om_vector_rms(struct omVector *vector){

	// set rms value to zero
	double rms = 0.0;

	// compute the sum of all vector's square values
	for(int i = 0; i < vector->_length;++i)
		rms += vector->_values[i]*vector->_values[i];

	// divide the sum by the length of the vector
	rms /= (double)(vector->_length);

	//return the root mean square
	return sqrt(rms);

}


/* compute the mean value of the vector */
double om_vector_mean(struct omVector *vector){

	// set mean value to zero
	double mean = 0.0;

	// get sum of all vector values
	for(int i = 0; i < vector->_length;++i)
		mean += vector->_values[i];

	// return the sum divide by the number of elements
	return mean / (double)(vector->_length);
}

/* get the number of elements */
int om_vector_getLength(struct omVector *vector){
	return vector->_length;
}


/* compute the median value of the vector */
double om_vector_median(struct omVector *vector){

	// set mean value to zero
	double median = 0.0;

	// if there is an odd number of elements
	if(vector->_length % 2 == 1){
		double tmp = ((double)(vector->_length) + 1.0)/2.0;
		int upper = (int)(ceil(tmp)) - 1;
		int lower = (int)(floor(tmp)) - 1;

		median = (vector->_values[upper] + vector->_values[lower])/2.0;

	}
	// even number of elements
	else{
		int index = ((vector->_length + 1) / 2) - 1;
		median = vector->_values[index];
	}

	//return median
	return median;
}


/* interpolation */
double om_vector_interpolation(omVector *x,omVector *y,double xq){

	double res;

	for(int i=0;i<x->_length-1;++i){

		double x0 = om_vector_getValue(x, i);
		double x1 = om_vector_getValue(x, i+1);

		if(xq >= x0 && xq <= x1){

			double y0 = om_vector_getValue(y, i);
			double y1 = om_vector_getValue(y, i+1);

			res = y0 + ( ((y1 - y0)/(x1 - x0))*(xq - x0)  );

			break;
		}

	}


	return res;

}




/* compute the norm of a vector */
double om_vector_norm(struct omVector *vector){

	// set rms value to zero
	double norm = 0.0;

	// compute the sum of all vector's square values
	for(int i = 0; i < vector->_length;++i)
		norm += vector->_values[i]*vector->_values[i];

	//return the root square of the sum
	return sqrt(norm);

}

/* normalize a vector */
void om_vector_normalize(struct omVector *vector){

	// get the norm of the vector
	double norm = om_vector_norm(vector);

	// divide all vector's values by the norm
	for(int i = 0; i < vector->_length;++i)
		om_vector_setValue(vector,i,vector->_values[i]/norm);

}

/* display a vector in a terminal */
void om_vector_display(struct omVector *vector){

	int numberofdecimals = 10;

	printf("[");
	for (int i = 0; i < vector->_length; i++)
		printf(" %.*f ", numberofdecimals, vector->_values[i]);
	printf("]\n");

}

/* release memory allocated by the vector */
void om_vector_free(struct omVector *vector){

	free(vector->_values);
	vector->_values = 0;

}


///////////////////////////////////////////////////////
/////             Matrix class                    /////
///////////////////////////////////////////////////////


/* allocation of a new matrix object */
void om_matrix_create(struct omMatrix *matrix,int rows,int columns){

	// set number of rows
	matrix->_rows = rows;

	// set number of rows
	matrix->_columns = columns;

	// allocation of the array
	matrix->_values = (double**)malloc( matrix->_rows* sizeof(double*));

	//set all values to zero
   for (int i=0; i<rows; i++){
	   matrix->_values[i] = (double*)malloc( matrix->_columns* sizeof(double));
	   for (int j=0; j<columns; j++)
		   om_matrix_setValue(matrix,i,j,0.0);

   }
}

/* release memory allocated by the matrix */
void om_matrix_free(struct omMatrix *matrix){


	for (int i=0; i<matrix->_rows; i++){
		free(matrix->_values[i]);
		matrix->_values[i] = 0;
	}

	free(matrix->_values);
	matrix->_values = 0;

}

/* set the (i x j)th value of a matrix */
void om_matrix_setValue(struct omMatrix *matrix,int i,int j,double value){

	//set the value
	matrix->_values[i][j] = value;
}

/* get the (i x j)th value of a matrix */
double om_matrix_getValue(struct omMatrix *matrix,int i,int j){

	//return the value
	return matrix->_values[i][j];
}


/**
 *
 */
void om_matrix_convolution2D(omMatrix *A,omMatrix *B,omMatrix *C){

	int nA = A->_rows;
	int mA = A->_columns;

	int nB = B->_rows;
	int mB = B->_columns;

	int nC = nA + nB - 1;
	int mC = mA + mB - 1;

	//printf("nc = %d mc = %d\n",nC,mC);
	om_matrix_create(C,nC,mC);


	for(int i=0;i<nC;i++){
		for(int j=0;j<mC;j++){

			double value = 0.0;

			for(int p=0;p<nA;p++){
				for(int q=0;q<mA;q++){

					if(i - p >= 0 && j - q >= 0 && i - p < nB && j - q < mB)
						value += om_matrix_getValue(A, p, q)*om_matrix_getValue(B, i - p , j - q );

					//printf("idx_1 = %d idx_2 = %d\n",i - p,j - q);
				}
			}
			om_matrix_setValue(C, i, j, value);
		}
	}

}


/**
 *
 */
void om_matrix_convolution2D_valid(omMatrix *A,omMatrix *B,omMatrix *C){

	int nA = A->_rows;
	int mA = A->_columns;

	int nB = B->_rows;
	int mB = B->_columns;

	int nC = nA + nB - 1;
	int mC = mA + mB - 1;

	om_matrix_create(C,nA - nB + 1,mA - mB + 1);

	for(int i=0;i<nC;i++){
		for(int j=0;j<mC;j++){

			if( i >= nB-1 && i<nC - nB + 1 && j >= mB-1 && j < mC - mB + 1 ){
				double value = 0.0;

				for(int p=0;p<nA;p++){
					for(int q=0;q<mA;q++){

						if(i - p  >= 0 && j - q >= 0 && i - p < nB && j - q < mB)
							value += om_matrix_getValue(A, p, q)*om_matrix_getValue(B, i - p , j - q );
					}
				}

				om_matrix_setValue(C, i - nB + 1, j- mB + 1, value);
			}
		}
	}

}




/* compute the norm of a matrix */
double om_matrix_norm(struct omMatrix *matrix){

	// set max to zero
	double sum=0.0;

	// compute the sum of all values
	for(int j=0;j<matrix->_columns;++j){

		for(int i=0;i<matrix->_rows;++i)
			sum += om_matrix_getValue(matrix,i,j)*om_matrix_getValue(matrix,i,j);

	}

	//return the norm
	return sqrt(sum);


}


int om_matrix_isSymmetric(struct omMatrix *matrix){

	int bool;

	if(om_matrix_isSquare(matrix) == 1){
		omMatrix transpose;
		omMatrix diff;
		om_matrix_create(&transpose,matrix->_columns,matrix->_rows);
		om_matrix_create(&diff,matrix->_columns,matrix->_rows);

		om_matrix_transpose(matrix,&transpose);
		om_operator_matrix_sub(matrix,&transpose,&diff);

		if(om_matrix_norm(&diff) < EPSILON)
			bool=1;
		else
			bool=0;

		om_matrix_free(&transpose);
		om_matrix_free(&diff);

	}else{
		bool = 0;
	}


	return bool;
}


/* compute the determinant of a matrix */
double om_matrix_determinant(struct omMatrix *matrix){

	   double determinant = -1.0;

	   // if the matrix is not square, return -1.0
	   if( om_matrix_isSquare(matrix) == 1 ){

		   // variables
		   omMatrix L;
		   omMatrix U;
		   double detL = 1.0;
		   double detU = 1.0;

		   om_matrix_create(&L,matrix->_rows,matrix->_columns);
		   om_matrix_create(&U,matrix->_rows,matrix->_columns);

		   // get LU factorization of the matrix
		   om_matrix_factorizationLU(matrix,&L,&U);

		   // compute the determinant of L and U.
		   //As they are both triangular, the determinant is the product of diagonal values.
		   for(int i=0;i<matrix->_rows;++i){
			   detL *= om_matrix_getValue(&L,i,i);
			   detU *= om_matrix_getValue(&U,i,i);
		   }

		   // multiply the determinant of L and U to get the determinant of the matrix
		   // Indeed det(AB) = det(A)det(B)
		   determinant = detL*detU;

		   om_matrix_free(&L);
		   om_matrix_free(&U);

	   }

	   // return the determinant
	   return  determinant;

}

/* get the trace of a matrix */
double om_matrix_trace(struct omMatrix *matrix){

	//set trace to zero
	double trace = 0.0;

	//if the matrix is non square, return null
	if( om_matrix_isSquare(matrix) == 1 ){

		// compute the sum of diagonal values
		for (int i = 0; i < matrix->_rows; i++)
			trace += om_matrix_getValue(matrix,i,i);

	}

	// return the trace
	return trace;
}

/* get all eigen values and associated eigen vector of a matrix */
void om_matrix_getEingenValues(struct omMatrix *matrix,struct omVector **eigen_vectors,double **eigen_values,int N){


	// if the matrix  symmetric, apply the Jacobi cyclic method
	if(om_matrix_isSymmetric(matrix)== 1){

		//*/
        // Initialize the eigenvalues to the identity matrix.
		int row, i, j, k, m;
		double *pAk, *pAm, *p_r, *p_e;
		double threshold_norm;
		double threshold;
		double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
		double sin_2phi, cos_2phi, cot_2phi;
		double dum1;
		double dum2;
		double dum3;
		double r;
		double max;
		int n = matrix->_columns;
		double* eigenvectors = (double*)malloc(n*n*sizeof(double));
		double* A = (double*)malloc(n*n*sizeof(double));

		for(int i2=0;i2<n;i2++)
			for(int j2=0;j2<n;j2++)
				A[i2*n + j2]=om_matrix_getValue(matrix, i2, j2);


		for (p_e = eigenvectors, i = 0; i < n; i++)
			for (j = 0; j < n; p_e++, j++)
				if (i == j) *p_e = 1.0; else *p_e = 0.0;

		// Calculate the threshold and threshold_norm.

		for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++)
			for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);


		threshold = sqrt(threshold + threshold);
		threshold_norm = threshold * DBL_EPSILON;
		max = threshold + 1.0;

		while (threshold > threshold_norm) {

			threshold /= 10.0;

			if (max < threshold) continue;
			max = 0.0;

			for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
				for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
					if ( fabs(*(pAk + m)) < threshold ) continue;

					// Calculate the sin and cos of the rotation angle which
					// annihilates A[k][m].

					cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
					dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
					if (cot_2phi < 0.0) dum1 = -dum1;

					tan_phi = -cot_2phi + dum1;
					tan2_phi = tan_phi * tan_phi;
					sin2_phi = tan2_phi / (1.0 + tan2_phi);
					cos2_phi = 1.0 - sin2_phi;
					sin_phi = sqrt(sin2_phi);

					if (tan_phi < 0.0) sin_phi = - sin_phi;

					cos_phi = sqrt(cos2_phi);
					sin_2phi = 2.0 * sin_phi * cos_phi;
					cos_2phi = cos2_phi - sin2_phi;

				   // Rotate columns k and m for both the matrix A
				   //     and the matrix of eigenvectors.

					p_r = A;
					dum1 = *(pAk + k);
					dum2 = *(pAm + m);
					dum3 = *(pAk + m);
					*(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
					*(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
					*(pAk + m) = 0.0;
					*(pAm + k) = 0.0;
					for (i = 0; i < n; p_r += n, i++) {
						if ( (i == k) || (i == m) ) continue;
						if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
						if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
						dum3 = dum1 * cos_phi + dum2 * sin_phi;
						if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
						dum3 = - dum1 * sin_phi + dum2 * cos_phi;
						if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
					}
					for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
						dum1 = *(p_e + k);
						dum2 = *(p_e + m);
						*(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
						*(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
					}
				}
				for (i = 0; i < n; i++)
					if ( i == k ) continue;
					else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
			}
		}

		(*eigen_values) = (double*)malloc(matrix->_columns*sizeof(double));
		(*eigen_vectors) = (omVector*)malloc(matrix->_columns*sizeof(omVector));

		for (pAk = A, k = 0; k < n; pAk += n, k++) (*eigen_values)[k] = *(pAk + k);


		for(int l=0;l<n;++l){

			//eigen vector are the column values of matrix P
			om_vector_create(&(*eigen_vectors)[l],n);
		}


		for(int l=0;l<n;++l)
			for(int k=0;k<n;++k)
				//om_vector_setValue(&(*eigen_vectors)[l],k,  k== 0? eigenvectors[k*n + l]*(-1.0) : eigenvectors[k*n + l] );
				om_vector_setValue(&(*eigen_vectors)[l],k,   eigenvectors[k*n + l] );


		free(A);
		free(eigenvectors);
		A = 0;
		eigenvectors = 0;
		//free(pAk);
		//free(pAm);
		//free(p_r);
		//free(p_e);



		/*/
		// variables
		omMatrix A;
		omMatrix Q;
		omMatrix R;
		omMatrix P;

		int m = matrix->_rows;
		int n = matrix->_columns;

		//allocation
		om_matrix_create(&Q,m,n);
		om_matrix_create(&R,m,n);
		om_matrix_create(&A,m,n);
		om_matrix_createIdentity(&P,m);

		// initialization of D
		om_matrix_clone(matrix,&A);

		// for each iteration
		for(int index=0;index<N;index++){

			//some variables
			omMatrix P_tmp;
			om_matrix_create(&P_tmp,matrix->_rows,matrix->_columns);

			// QR factorization of A
			om_matrix_factorizationQR(&A,&Q,&R);


			// compute A = R*Q
			om_operator_matrix_mul(&R,&Q,&A);

			// compute P_tp1 = P_t*Q
			om_operator_matrix_mul(&P,&Q,&P_tmp);
			om_matrix_clone(&P_tmp,&P);



			//free memory
			om_matrix_free(&P_tmp);
		}

		// allocation of eigen values and eigen vector
		(*eigen_values) = (double*)malloc(matrix->_columns*sizeof(double));
		(*eigen_vectors) = (omVector*)malloc(matrix->_columns*sizeof(omVector));

		for(int l=0;l<matrix->_columns;++l){

			//eigen values are the diagonal values of matrix D
			(*eigen_values)[l] = om_matrix_getValue(&A,l,l);

			//eigen vector are the column values of matrix P
			om_vector_create(&(*eigen_vectors)[l],n);
			om_matrix_getColumn(&P,l,&(*eigen_vectors)[l]);

		}

		// free memory
		om_matrix_free(&A);
		om_matrix_free(&Q);
		om_matrix_free(&R);
		om_matrix_free(&P);
		//*/


	}else{


		int n = matrix->_columns;
		omMatrix H;
		omMatrix V;
		omVector o;
		omVector d;
		omVector e;

		om_vector_create(&o,n);
		om_vector_create(&d,n);
		om_vector_create(&e,n);
		om_matrix_create(&H,n,n);
		om_matrix_create(&V,n,n);
		om_matrix_clone(matrix,&H);


		for(int index=0;index<n;index++){
			om_vector_setValue(&d,index,0.0);
			om_vector_setValue(&e,index,0.0);
		}

		//Nonsymmetric reduction to Hessenberg form
		//  This is derived from the Algol procedures orthes and ortran,
		//  by Martin and Wilkinson, Handbook for Auto. Comp.,
		//  Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutines in EISPACK.

		int low = 0;
		int high = n-1;

		for (int m = low+1; m <= high-1; m++) {

			// Scale column.
			double scale = 0.0;
			for (int i = m; i <= high; i++) {
				 //scale = scale + Math.abs(H[i][m-1]);
				scale = scale + fabs(om_matrix_getValue(&H,i,m-1));

			}
			if (scale != 0.0) {

				// Compute Householder transformation.
				double h = 0.0;
				for (int i = high; i >= m; i--) {

					//ort[i] = H[i][m-1]/scale;
					om_vector_setValue(&o,i,om_matrix_getValue(&H,i,m-1)/scale);

					//h += ort[i] * ort[i];
					h += om_vector_getValue(&o,i)*om_vector_getValue(&o,i);

				}

				double g = sqrt(h);
				if (om_vector_getValue(&o,m) > 0) {
					g *= -1.0;
				}

				//h = h - ort[m] * g;
				h = h - (om_vector_getValue(&o,m) * g);

				//ort[m] = ort[m] - g;
				om_vector_setValue(&o,m,om_vector_getValue(&o,m) - g);

				// Apply Householder similarity transformation
				// H = (I-u*u'/h)*H*(I-u*u')/h)

				for (int j = m; j < n; j++) {
					double f = 0.0;

					for (int i = high; i >= m; i--) {
						//f += ort[i]*H[i][j];
						f += om_vector_getValue(&o,i)*om_matrix_getValue(&H,i,j);
					}

					f /= h;

					for (int i = m; i <= high; i++) {
						//H[i][j] -= f*ort[i];
						om_matrix_setValue(&H,i,j,om_matrix_getValue(&H,i,j) - (f*om_vector_getValue(&o,i)));

					}
				}

				for (int i = 0; i <= high; i++) {
					double f = 0.0;
					for (int j = high; j >= m; j--) {
						//f += ort[j]*H[i][j];
						f += om_vector_getValue(&o,j)*om_matrix_getValue(&H,i,j);
					}
					f /= h;
					for (int j = m; j <= high; j++) {
						//H[i][j] -= f*ort[j];
						om_matrix_setValue(&H,i,j,om_matrix_getValue(&H,i,j) - (f*om_vector_getValue(&o,j)));

					}
				}
				//ort[m] = scale*ort[m];
				om_vector_setValue(&o,m,scale*om_vector_getValue(&o,m));

				//H[m][m-1] = scale*g;
				om_matrix_setValue(&H,m,m-1,scale*g);
			}
		}

		// Accumulate transformations (Algol's ortran).
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				//V[i][j] = (i == j ? 1.0 : 0.0);
				om_matrix_setValue(&V,i,j,(i == j ? 1.0 : 0.0));
			}
		}

		for (int m = high-1; m >= low+1; m--) {

			//if (H[m][m-1] != 0.0) {
			if (om_matrix_getValue(&H,m,m-1) != 0.0) {

				for (int i = m+1; i <= high; i++) {
					//ort[i] = H[i][m-1];
					om_vector_setValue(&o,i,om_matrix_getValue(&H,i,m-1));
				}

				for (int j = m; j <= high; j++) {
					double g = 0.0;
					for (int i = m; i <= high; i++) {
						//g += ort[i] * V[i][j];
						g += om_vector_getValue(&o,i)*om_matrix_getValue(&V,i,j);
					}

					// Double division avoids possible underflow
					//g = (g / ort[m]) / H[m][m-1];
					g = (g / om_vector_getValue(&o,m))/om_matrix_getValue(&H,m,m-1);
					for (int i = m; i <= high; i++) {
						//V[i][j] += g * ort[i];
						om_matrix_setValue(&V,i,j,om_matrix_getValue(&V,i,j) + (g * om_vector_getValue(&o,i) ));
					}

				}
			}
		}





		//  This is derived from the Algol procedure hqr2,
		//  by Martin and Wilkinson, Handbook for Auto. Comp.,
		//  Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		// Initialize

		int nn = n;
		n = nn-1;
		low = 0;
		high = nn - 1;

		double eps = pow(2.0,-52.0);
		double exshift = 0.0;
		double p=0,q=0,r=0,s=0,z=0,t,w,x,y;


		// Store roots isolated by balanc and compute matrix norm

		double norm = 0.0;
		for (int i = 0; i < nn; i++) {
			if ( (i < low) | (i > high) ) {
				//d[i] = H[i][i];
				om_vector_setValue(&d,i,om_matrix_getValue(&H,i,i));

				//e[i] = 0.0;
				om_vector_setValue(&e,i,0.0);
			}

			//for (int j = max(i-1, 0) ; j < nn; j++) {
			for (int j = (i-1 > 0 ? i-1 : 0) ; j < nn; j++) {
				//norm = norm + fabs(H[i][j]);
				norm = norm + fabs(om_matrix_getValue(&H,i,j));
			}
		}



		// Outer loop over eigenvalue index

		int iter = 0;
		while (n >= low) {

			// Look for single small sub-diagonal element
			int l = n;
			while (l > low) {
				//s = fabs(H[l-1][l-1]) + fabs(H[l][l]);
				s = fabs(om_matrix_getValue(&H,l-1,l-1)) + fabs(om_matrix_getValue(&H,l,l));

				if (s == 0.0) {
					s = norm;
				}
				//if (fabs(H[l][l-1]) < eps * s) {
				if (fabs(om_matrix_getValue(&H,l,l-1)) < eps * s) {
					break;
				}
				l--;
			}



			// Check for convergence
			// One root found

			if (l == n) {



				//H[n][n] = H[n][n] + exshift;
				om_matrix_setValue(&H,n,n,om_matrix_getValue(&H,n,n) + exshift);

				//d[n] = H[n][n];
				om_vector_setValue(&d,n,om_matrix_getValue(&H,n,n));

				//e[n] = 0.0;
				om_vector_setValue(&e,n,0.0);

				n--;
				iter = 0;

			// Two roots found
			} else if (l == n-1) {

				//w = H[n][n-1] * H[n-1][n];
				w = om_matrix_getValue(&H,n,n-1)*om_matrix_getValue(&H,n-1,n);

				//p = (H[n-1][n-1] - H[n][n]) / 2.0;
				p = (om_matrix_getValue(&H,n-1,n-1) - om_matrix_getValue(&H,n,n)) / 2.0;

				q = (p * p) + w;

				z = sqrt(fabs(q));

				//H[n][n] = H[n][n] + exshift;
				om_matrix_setValue(&H,n,n,om_matrix_getValue(&H,n,n) + exshift);

				//H[n-1][n-1] = H[n-1][n-1] + exshift;
				om_matrix_setValue(&H,n-1,n-1,om_matrix_getValue(&H,n-1,n-1) + exshift);

				//x = H[n][n];
				x = om_matrix_getValue(&H,n,n);

				// Real pair
				if (q >= 0) {
					if (p >= 0) {
						z = p + z;
					} else {
						z = p - z;
					}

					//d[n-1] = x + z;
					om_vector_setValue(&d,n-1,x + z);

					//d[n] = d[n-1];
					om_vector_setValue(&d,n,x + z);

					if (z != 0.0) {
						//d[n] = x - w / z;
						om_vector_setValue(&d,n,x - w / z);
					}

					//e[n-1] = 0.0;
					om_vector_setValue(&e,n-1,0.0);

					//e[n] = 0.0;
					om_vector_setValue(&e,n,0.0);

					//x = H[n][n-1];
					x = om_matrix_getValue(&H,n,n-1);

					s = fabs(x) + fabs(z);
					q = z / s;

					r = sqrt(p * p + q * q);
					p = p / r;
					q = q / r;

					// Row modification

					for (int j = n-1; j < nn; j++) {
						//z = H[n-1][j];
						z = om_matrix_getValue(&H,n-1,j);

						//H[n-1][j] = q * z + p * H[n][j];
						om_matrix_setValue(&H,n-1,j, (q * z) + (p * om_matrix_getValue(&H,n,j)));

						//H[n][j] = q * H[n][j] - p * z;
						om_matrix_setValue(&H,n,j,(q * om_matrix_getValue(&H,n,j)) - (p * z));
					}

					// Column modification

					for (int i = 0; i <= n; i++) {
						//z = H[i][n-1];
						z = om_matrix_getValue(&H,i,n-1);

						//H[i][n-1] = q * z + p * H[i][n];
						om_matrix_setValue(&H,i,n-1, (q * z) + (p * om_matrix_getValue(&H,i,n)));

						//H[i][n] = q * H[i][n] - p * z;
						om_matrix_setValue(&H,i,n, (q * om_matrix_getValue(&H,i,n)) - (p * z));
					}

					// Accumulate transformations

					for (int i = low; i <= high; i++) {

						//z = V[i][n-1];
						z = om_matrix_getValue(&V,i,n-1);

						//V[i][n-1] = q * z + p * V[i][n];
						om_matrix_setValue(&V,i,n-1, (q * z) + (p * om_matrix_getValue(&V,i,n)));

						//V[i][n] = q * V[i][n] - p * z;
						om_matrix_setValue(&V,i,n, (q * om_matrix_getValue(&V,i,n)) - (p * z));
					}

				// Complex pair
				} else {

					//d[n-1] = x + p;
					om_vector_setValue(&d,n-1, x + p);

					//d[n] = x + p;
					om_vector_setValue(&d,n, x + p);

					//e[n-1] = z;
					om_vector_setValue(&e,n-1,z);

					//e[n] = -z;
					om_vector_setValue(&e,n,-z);
				}
				n = n - 2;
				iter = 0;




			// No convergence yet
			} else {
				// Form shift

				//x = H[n][n];
				x = om_matrix_getValue(&H,n,n);
				y = 0.0;
				w = 0.0;

				if (l < n) {
					//y = H[n-1][n-1];
					y = om_matrix_getValue(&H,n-1,n-1);

					//w = H[n][n-1] * H[n-1][n];
					w = om_matrix_getValue(&H,n,n-1) * om_matrix_getValue(&H,n-1,n);
				}


				// Wilkinson's original ad hoc shift

				if (iter == 10) {
					exshift += x;
					for (int i = low; i <= n; i++) {
						//H[i][i] -= x;
						om_matrix_setValue(&H,i,i, om_matrix_getValue(&H,i,i) - x);

					}
					//s = fabs(H[n][n-1]) + fabs(H[n-1][n-2]);
					s = fabs(om_matrix_getValue(&H,n,n-1)) + fabs(om_matrix_getValue(&H,n-1,n-2));

					x = y = 0.75 * s;
					w = -0.4375 * s * s;
				}


				// MATLAB's new ad hoc shift
				if (iter == 30) {
					s = (y - x) / 2.0;
					s = s * s + w;
					if (s > 0) {
						s = sqrt(s);
						if (y < x) {
							s = -s;
						}
						s = x - w / ((y - x) / 2.0 + s);
						for (int i = low; i <= n; i++) {
							//H[i][i] -= s;
							om_matrix_setValue(&H,i,i, om_matrix_getValue(&H,i,i) - s);
						}
						exshift += s;
						x = y = w = 0.964;
					}
				}

				iter = iter + 1;   // (Could check iteration count here.)

				// Look for two consecutive small sub-diagonal elements

				int m = n-2;
				while (m >= l) {
					//z = H[m][m];
					z = om_matrix_getValue(&H,m,m);

					r = x - z;
					s = y - z;
					//p = (r * s - w) / H[m+1][m] + H[m][m+1];
					p = (r * s - w) / om_matrix_getValue(&H,m+1,m) + om_matrix_getValue(&H,m,m+1);

					//q = H[m+1][m+1] - z - r - s;
					q = om_matrix_getValue(&H,m+1,m+1) - z - r - s;

					//r = H[m+2][m+1];
					r = om_matrix_getValue(&H,m+2,m+1);

					s = fabs(p) + fabs(q) + fabs(r);

					p = p / s;
					q = q / s;
					r = r / s;

					if (m == l) {
						break;
					}

					double tmp = fabs( om_matrix_getValue(&H,m,m-1) ) * (fabs(q) + fabs(r));
					double tmp2 = (fabs(p) * ( fabs( om_matrix_getValue(&H,m-1,m-1) ) + fabs(z) + fabs(om_matrix_getValue(&H,m+1,m+1)) ) );
					if (  tmp < eps * tmp2 ) {
						break;
					}
					m--;
				}

				for (int i = m+2; i <= n; i++) {
					//H[i][i-2] = 0.0;
					om_matrix_setValue(&H,i,i-2,0.0);

					if (i > m+2) {
						//H[i][i-3] = 0.0;
						om_matrix_setValue(&H,i,i-3,0.0);
					}


				}

					// Double QR step involving rows l:n and columns m:n

					for (int k = m; k <= n-1; k++) {
						int notlast = (k != n-1) ? 1 : 0 ;

						if (k != m) {
							//p = H[k][k-1];
							p = om_matrix_getValue(&H,k,k-1);

		                	//q = H[k+1][k-1];
		                	q = om_matrix_getValue(&H,k+1,k-1);

		                	//r = (notlast ? H[k+2][k-1] : 0.0);
		                	r = (notlast == 1 ? om_matrix_getValue(&H,k+2,k-1) : 0.0);

		                	x = fabs(p) + fabs(q) + fabs(r);
		                	if (x != 0.0) {
		                		p = p / x;
		                		q = q / x;
		                		r = r / x;
		                	}
						}
		                if (x == 0.0) {
		                   break;
		                }

		                s = sqrt(p * p + q * q + r * r);

		                if (p < 0) {
		                   s = -s;
		                }

		                if (s != 0) {
		                   if (k != m) {
		                      //H[k][k-1] = -s * x;
		                	   om_matrix_setValue(&H,k,k-1,-s * x);

		                   } else if (l != m) {
		                      //H[k][k-1] = -H[k][k-1];
		                	   om_matrix_setValue(&H,k,k-1, om_matrix_getValue(&H,k,k-1)*(-1.0) );
		                   }
		                   p = p + s;
		                   x = p / s;
		                   y = q / s;
		                   z = r / s;
		                   q = q / p;
		                   r = r / p;

		                   // Row modification

		                   for (int j = k; j < nn; j++) {
		                      //p = H[k][j] + q * H[k+1][j];
		                      p = om_matrix_getValue(&H,k,j) + (q * om_matrix_getValue(&H,k+1,j));

		                      if (notlast == 1) {
		                         //p = p + r * H[k+2][j];
		                    	  p = p + (r * om_matrix_getValue(&H,k+2,j));

		                         //H[k+2][j] = H[k+2][j] - p * z;
		                         om_matrix_setValue(&H,k+2,j,om_matrix_getValue(&H,k+2,j) - (p*z));
		                      }

		                      //H[k][j] = H[k][j] - p * x;
		                      om_matrix_setValue(&H,k,j,om_matrix_getValue(&H,k,j) - (p*x));

		                      //H[k+1][j] = H[k+1][j] - p * y;
		                      om_matrix_setValue(&H,k+1,j,om_matrix_getValue(&H,k+1,j) - (p*y));
		                   }

		                   // Column modification

		                   //for (int i = 0; i <= min(n,k+3); i++) {
		                   for (int i = 0; i <= (n < k+3 ? n : k+3); i++) {

		                	   //p = x * H[i][k] + y * H[i][k+1];
		                      p = (x * om_matrix_getValue(&H,i,k)) + (y * om_matrix_getValue(&H,i,k+1));

		                      if (notlast == 1) {
		                         //p = p + z * H[i][k+2];
		                    	  p = p + (z * om_matrix_getValue(&H,i,k+2));

		                    	  //H[i][k+2] = H[i][k+2] - p * r;
		                    	  om_matrix_setValue(&H,i,k+2,om_matrix_getValue(&H,i,k+2) - (p*r));
		                      }

		                      //H[i][k] = H[i][k] - p;
		                      om_matrix_setValue(&H,i,k,om_matrix_getValue(&H,i,k) - p);

		                      //H[i][k+1] = H[i][k+1] - p * q;
		                      om_matrix_setValue(&H,i,k+1,om_matrix_getValue(&H,i,k+1) - (p*q));
		                   }

		                   // Accumulate transformations

		                   for (int i = low; i <= high; i++) {
		                      //p = x * V[i][k] + y * V[i][k+1];
		                	   p = (x * om_matrix_getValue(&V,i,k)) + (y * om_matrix_getValue(&V,i,k+1));

		                	   if (notlast == 1) {
		                         //p = p + z * V[i][k+2];
		                		   p = p + (z * om_matrix_getValue(&V,i,k+2));

		                         //V[i][k+2] = V[i][k+2] - p * r;
		                		   om_matrix_setValue(&V,i,k+2,om_matrix_getValue(&V,i,k+2) - (p*r));

		                      }

		                      //V[i][k] = V[i][k] - p;
		                	   om_matrix_setValue(&V,i,k,om_matrix_getValue(&V,i,k) - p);

		                      //V[i][k+1] = V[i][k+1] - p * q;
		                      om_matrix_setValue(&V,i,k+1,om_matrix_getValue(&V,i,k+1) - (p*q));
		                   }
		                }  // (s != 0)
		             }  // k loop
		          }  // check convergence
		       }  // while (n >= low)

		       // Backsubstitute to find vectors of upper triangular form

		       if (norm == 0.0) {
		          return;
		       }

		       for (n = nn-1; n >= 0; n--) {
		          //p = d[n];
		    	  p = om_vector_getValue(&d,n);
		    	  //q = e[n];
		    	  q = om_vector_getValue(&e,n);

		    	  // Real vector

		          if (q == 0.0) {
		             int l = n;

		             //H[n][n] = 1.0;
		             om_matrix_setValue(&H,n,n,1.0);

		             for (int i = n-1; i >= 0; i--) {
		                //w = H[i][i] - p;
		            	w = om_matrix_getValue(&H,i,i) - p;
		                r = 0.0;

		                for (int j = l; j <= n; j++) {
		                   //r = r + H[i][j] * H[j][n];
		                	r = r + (om_matrix_getValue(&H,i,j) * om_matrix_getValue(&H,j,n));

		                }
		                //if (e[i] < 0.0) {
		                if (om_vector_getValue(&e,i) < 0.0) {
		                   z = w;
		                   s = r;
		                } else {
		                   l = i;
		                   //if (e[i] == 0.0) {
		                   if (om_vector_getValue(&e,i) == 0.0) {

		                      if (w != 0.0) {
		                         //H[i][n] = -r / w;
		                    	 om_matrix_setValue(&H,i,n,-r / w);
		                      } else {
		                         //H[i][n] = -r / (eps * norm);
		                    	  om_matrix_setValue(&H,i,n,-r / (eps*norm));
		                      }

		                   // Solve real equations

		                   } else {
		                      //x = H[i][i+1];
		                	  x = om_matrix_getValue(&H,i,i+1);

		                      //y = H[i+1][i];
		                	  y = om_matrix_getValue(&H,i+1,i);

		                      //q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
		                	  q = (om_vector_getValue(&d,i) - p) * (om_vector_getValue(&d,i) - p) + (om_vector_getValue(&e,i) * om_vector_getValue(&e,i));

		                      t = (x * s - z * r) / q;

		                      //H[i][n] = t;
		                      om_matrix_setValue(&H,i,n,t);

		                      if (fabs(x) > fabs(z)) {
		                         //H[i+1][n] = (-r - w * t) / x;
		                    	  om_matrix_setValue(&H,i+1,n,(-r - w * t) / x);
		                      } else {
		                         //H[i+1][n] = (-s - y * t) / z;
		                    	  om_matrix_setValue(&H,i+1,n,(-s - y * t) / z);
		                      }
		                   }

		                   // Overflow control

		                   //t = fabs(H[i][n]);
		                   t = fabs(om_matrix_getValue(&H,i,n));

		                   if ((eps * t) * t > 1.0) {
		                      for (int j = i; j <= n; j++) {
		                    	  //H[j][n] = H[j][n] / t;
		                    	  om_matrix_setValue(&H,j,n,om_matrix_getValue(&H,j,n)/t);
		                      }
		                   }
		                }
		             }

		          // Complex vector

		          } else if (q < 0.0) {
		             int l = n-1;

		             // Last vector component imaginary so matrix is triangular

		             //if (fabs(H[n][n-1]) > fabs(H[n-1][n])) {
		             if (fabs(om_matrix_getValue(&H,n,n-1)) > fabs(om_matrix_getValue(&H,n-1,n))) {

		            	 //H[n-1][n-1] = q / H[n][n-1];
		            	 om_matrix_setValue(&H,n-1,n-1,q / om_matrix_getValue(&H,n,n-1));

		            	 //H[n-1][n] = -(H[n][n] - p) / H[n][n-1];
		            	 om_matrix_setValue(&H,n-1,n, (om_matrix_getValue(&H,n,n) - p) / om_matrix_getValue(&H,n,n-1));
		             } else {
		            	 double cdivr,cdivi;

		                //cdiv(0.0,-H[n-1][n]  ,H[n-1][n-1]-p,q,&cdivr,&cdivi);
		            	 cdiv(0.0,-om_matrix_getValue(&H,n-1,n) , om_matrix_getValue(&H,n-1,n-1) - p ,q,&cdivr,&cdivi);

		            	 //H[n-1][n-1] = cdivr;
		            	 om_matrix_setValue(&H,n-1,n-1,cdivr);

		            	 //H[n-1][n] = cdivi;
		            	 om_matrix_setValue(&H,n-1,n,cdivi);
		             }

		             //H[n][n-1] = 0.0;
		             om_matrix_setValue(&H,n,n-1,0.0);

		             //H[n][n] = 1.0;
		             om_matrix_setValue(&H,n,n,1.0);

		             for (int i = n-2; i >= 0; i--) {
		                double ra,sa,vr,vi;
		                ra = 0.0;
		                sa = 0.0;
		                for (int j = l; j <= n; j++) {
		                   //ra = ra + H[i][j] * H[j][n-1];
		                	ra = ra + (om_matrix_getValue(&H,i,j) * om_matrix_getValue(&H,j,n-1));

		                	//sa = sa + H[i][j] * H[j][n];
		                	sa = sa + (om_matrix_getValue(&H,i,j) * om_matrix_getValue(&H,j,n));
		                }

		                //w = H[i][i] - p;
		                w = om_matrix_getValue(&H,i,i) - p;

		                //if (e[i] < 0.0) {
		                if (om_vector_getValue(&e,i) < 0.0) {
		                   z = w;
		                   r = ra;
		                   s = sa;
		                } else {
		                   l = i;

		                   //if (e[i] == 0) {
		                   if (om_vector_getValue(&e,i) == 0.0) {
		                	   double cdivr,cdivi;
		                	  cdiv(-ra,-sa,w,q,&cdivr,&cdivi);

		                	  //H[i][n-1] = cdivr;
		                	  om_matrix_setValue(&H,i,n-1,cdivr);

		                	  //H[i][n] = cdivi;
		                	  om_matrix_setValue(&H,i,n,cdivi);
		                   } else {

		                      // Solve complex equations

		                      //x = H[i][i+1];
		                      x =  om_matrix_getValue(&H,i,i+1);

		                	  //y = H[i+1][i];
		                      y = om_matrix_getValue(&H,i+1,i);

		                	  //vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
		                      vr = (om_vector_getValue(&d,i) - p) * (om_vector_getValue(&d,i) - p) + (om_vector_getValue(&e,i) * om_vector_getValue(&e,i)) - (q * q);

		                	  //vi = (d[i] - p) * 2.0 * q;
		                      vi = (om_vector_getValue(&d,i) - p) * 2.0 * q;


		                	  if ( (vr == 0.0) & (vi == 0.0)) {
		                         vr = eps * norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z));
		                      }

		                	  double cdivr,cdivi;
		                	  cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi,&cdivr,&cdivi);

		                	  //H[i][n-1] = cdivr;
		                	  om_matrix_setValue(&H,i,n-1,cdivr);

		                	  //H[i][n] = cdivi;
		                	  om_matrix_setValue(&H,i,n,cdivi);

		                      if (fabs(x) > (fabs(z) + fabs(q))) {

		                    	 //H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x;
		                    	 double tmp =  (-ra - w *om_matrix_getValue(&H,i,n-1) + q * om_matrix_getValue(&H,i,n)) / x;
		                    	 om_matrix_setValue(&H,i+1,n-1,tmp);

		                         //H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x;
		                    	 tmp =  (-sa - w *om_matrix_getValue(&H,i,n) + q * om_matrix_getValue(&H,i,n-1)) / x;
		                    	 om_matrix_setValue(&H,i+1,n,tmp);

		                      } else {

		                         //cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q,&cdivr,&cdivi);
		                    	 cdiv(-r - y * om_matrix_getValue(&H,i,n-1) ,-s-y * om_matrix_getValue(&H,i,n),z,q,&cdivr,&cdivi);

		                    	 //H[i+1][n-1] = cdivr;
			                	 om_matrix_setValue(&H,i+1,n-1,cdivr);

			                	 //H[i+1][n] = cdivi;
			                	 om_matrix_setValue(&H,i+1,n,cdivi);
		                      }
		                   }

		                   // Overflow control

		                   //t = Math.max(fabs(H[i][n-1]),fabs(H[i][n]));
		                   t = fabs(om_matrix_getValue(&H,i,n-1)) > fabs(om_matrix_getValue(&H,i,n)) ? fabs(om_matrix_getValue(&H,i,n-1)) : fabs(om_matrix_getValue(&H,i,n));


		                   if ((eps * t) * t > 1) {
		                      for (int j = i; j <= n; j++) {
		                         //H[j][n-1] = H[j][n-1] / t;
		                    	 om_matrix_setValue(&H,j,n-1,om_matrix_getValue(&H,j,n-1)/t);

		                    	 //H[j][n] = H[j][n] / t;
		                    	 om_matrix_setValue(&H,j,n,om_matrix_getValue(&H,j,n)/t);
		                      }
		                   }
		                }
		             }
		          }
		       }

		       // Vectors of isolated roots

		       for (int i = 0; i < nn; i++) {
		          if ( (i < low) | (i > high) ) {
		             for (int j = i; j < nn; j++) {
		                //V[i][j] = H[i][j];
		            	om_matrix_setValue(&V,i,j,om_matrix_getValue(&H,i,j));
		             }
		          }
		       }

		       // Back transformation to get eigenvectors of original matrix

		       for (int j = nn-1; j >= low; j--) {
		          for (int i = low; i <= high; i++) {
		             z = 0.0;
		             //for (int k = low; k <= min(j,high); k++) {
		             for (int k = low; k <= (j < high ? j : high); k++) {
		                //z = z + V[i][k] * H[k][j];
		            	 z = z + (om_matrix_getValue(&V,i,k) * om_matrix_getValue(&H,k,j));
		             }
		             //V[i][j] = z;
		             om_matrix_setValue(&V,i,j,z);
		          }
		       }
		

		// allocation of eigen values and eigen vector
		(*eigen_values) = (double*)malloc(matrix->_columns*sizeof(double));
		(*eigen_vectors) = (omVector*)malloc(matrix->_columns*sizeof(omVector));

		for(int l=0;l<matrix->_columns;++l){

			//eigen values are the diagonal values of matrix D
			(*eigen_values)[l] = om_vector_getValue(&d,l);

			//eigen vector are the column values of matrix P
			om_vector_create(&(*eigen_vectors)[l],matrix->_columns);
			om_matrix_getColumn(&V,l,&(*eigen_vectors)[l]);

		}

		om_vector_free(&o);
		om_vector_free(&d);
		om_vector_free(&e);

		om_matrix_free(&H);
		om_matrix_free(&V);



	}

	// test
	/*/
	for(int l=0;l<matrix->_columns;l++){
		omMatrix lL;

		om_matrix_createIdentity(&lL,3);
		om_operator_matrix_scal_mul(&lL,(*eigen_values)[l],&lL);
		om_operator_matrix_sub(matrix,&lL,&lL);

		omVector u_test_a;
		omVector u_test_b;

		om_vector_create(&u_test_a,3);
		om_vector_create(&u_test_b,3);

		om_operator_matrix_vector_mul(matrix,&(*eigen_vectors)[l],&u_test_a);
		om_operator_vector_scal_mul(&(*eigen_vectors)[l],(*eigen_values)[l],&u_test_b);

		printf("\n\nl = %d values test a = ",l) ;
		om_vector_display(&u_test_a);

		printf("l = %d values test b = ",l) ;
		om_vector_display(&u_test_b);

		printf("l = %d det(N - lL) = %f\n",l,om_matrix_determinant(&lL));
		printf("l = %d lambda = %f\n",l,(*eigen_values)[l]);

		om_matrix_free(&lL);
		om_vector_free(&u_test_b);
		om_vector_free(&u_test_a);
	}
	//*/





}

/* get the skew symmetric matrix of a 3d vector */
void om_matrix_skewSymmetricMatrix(struct omVector *in,struct omMatrix *out){


	// set values such as:
	//
	// v = ( v_x , v_y , v_z )^T
	//
	//     /  0.0  -v_z  v_y \
	// S = |  v_z   0.0 -v_x |
	//     \ -v_y   v_x  0.0 /
	//
	om_matrix_setValue(out,0,0,0.0);
	om_matrix_setValue(out,1,1,0.0);
	om_matrix_setValue(out,2,2,0.0);
	om_matrix_setValue(out,0,1,-in->_values[2]);
	om_matrix_setValue(out,0,2,in->_values[1]);
	om_matrix_setValue(out,1,0,in->_values[2]);
	om_matrix_setValue(out,1,2,-in->_values[0]);
	om_matrix_setValue(out,2,0,-in->_values[1]);
	om_matrix_setValue(out,2,1,in->_values[0]);

}

/* create an identity matrix */
void om_matrix_createIdentity(struct omMatrix *I,int n){

	//allocation
	om_matrix_create(I,n,n);

	//set all diagonal values to 1
	for(int i=0;i<n;i++)
		om_matrix_setValue(I,i,i,1.0);

}


/* compute the exponential of a matrix */
void om_matrix_exponential(struct omMatrix *matrix,struct omMatrix *exp,int N){

	// if the matrix is square
	if(matrix->_columns == matrix->_rows){

		// some variables
		omMatrix acc;
		double acc_n = 1.0;

		// the idea is to use Taylor's expansion for exponential such as:
		// exp(A) = I + A + (A^2)*(1/facto(2)) + (A^3)*(1/facto(3)) + ... +(A^N)*(1/facto(N))
		//
		// initialization by computing the first two iteration
		// exp(A) = I + A
		om_matrix_create(&acc,matrix->_rows,matrix->_columns);
		om_matrix_clone(matrix,&acc);
		om_matrix_createIdentity(exp,matrix->_rows);
		om_operator_matrix_add(exp,matrix,exp);

		// compute the rest of Taylor's expansion
		for (int i = 2 ; i<N; ++i){

			acc_n = acc_n*(double)(i);
			omMatrix tmp;

			om_matrix_create(&tmp,matrix->_rows,matrix->_columns);
			om_operator_matrix_mul(&acc,matrix,&tmp);
			om_matrix_clone(&tmp,&acc);

			om_operator_matrix_scal_div(&tmp,acc_n,&tmp);
			om_operator_matrix_add(exp,&tmp,exp);

			om_matrix_free(&tmp);
		}

		// free memory
		om_matrix_free(&acc);

	}


}

/* compute the square root of a matrix */
void om_matrix_squareRoot(struct omMatrix *matrix,struct omMatrix *m_sqrt,int N){

	// if the matrix is square
	if(matrix->_rows == matrix->_columns){

		// if the matrix is (1 x 1)
		if(matrix->_rows == 1){

			om_matrix_setValue(m_sqrt,0,0,sqrt(  om_matrix_getValue(matrix,0,0)));
		}
		// if the matrix is (2 x 2)
		else if(matrix->_rows == 2){

			omMatrix I;
			om_matrix_createIdentity(&I,2);

			double det = sqrt(om_matrix_determinant(matrix));
			double tmp = sqrt(om_matrix_trace(matrix) + (2.0*det));

			om_operator_matrix_scal_mul(&I,det,&I);
			om_operator_matrix_add(matrix,&I,m_sqrt);

			om_operator_matrix_scal_div(m_sqrt,tmp,m_sqrt);

			om_matrix_free(&I);

		}
		// otherwise, for all matrices (n x n)
		else{

			// variables
			omMatrix D;
			omMatrix Q;
			omMatrix R;
			omMatrix P;
			omMatrix P_inv;
			omMatrix S_tmp;
			omMatrix squareD;
			omMatrix P_tmp;

			// allocation
			om_matrix_create(&D,matrix->_rows,matrix->_rows);
			om_matrix_create(&P,matrix->_rows,matrix->_rows);
			om_matrix_create(&Q,matrix->_rows,matrix->_rows);
			om_matrix_create(&R,matrix->_rows,matrix->_rows);
			om_matrix_create(&P_tmp,matrix->_rows,matrix->_rows);
			om_matrix_create(&P_inv,matrix->_rows,matrix->_rows);
			om_matrix_create(&squareD,D._rows,D._columns);
			om_matrix_create(&S_tmp,matrix->_rows,matrix->_rows);

			// initialization
			// D = A
			om_matrix_clone(matrix,&D);

			// (Q,R) <- factoQR(D)
			om_matrix_factorizationQR(&D,&Q,&R);

			// P = Q
			om_matrix_clone(&Q,&P);

			//for each iteration
			for(int i=0;i<N;i++){

				// D = R*Q
				om_operator_matrix_mul(&R,&Q,&D);

				// (Q,R) <- factoQR(D)
				om_matrix_factorizationQR(&D,&Q,&R);

				// P = P*Q
				om_operator_matrix_mul(&P,&Q,&P_tmp);
				om_matrix_clone(&P_tmp,&P);

			}

			// as D is a diagonal matrix
			// sqrt(D) = ( sqrt(d_{ii}) );
			for(int i=0;i< squareD._rows;i++)
				om_matrix_setValue(&squareD, i, i, sqrt(om_matrix_getValue(&D, i, i)));

			// get inverse of P
			om_matrix_inverse(&P,&P_inv);

			// sqrt(A) = P*sqrt(D)*P^{-1}
			om_operator_matrix_mul(&P,&squareD,&S_tmp);
			om_operator_matrix_mul(&S_tmp,&P_inv,m_sqrt);




			// free memory
			om_matrix_free(&P);
			om_matrix_free(&P_inv);
			om_matrix_free(&S_tmp);
			om_matrix_free(&D);
			om_matrix_free(&R);
			om_matrix_free(&Q);
			om_matrix_free(&squareD);
			om_matrix_free(&P_tmp);

		}

	}

}

/* return 1 if the matrix is square, 0 otherwise */
int om_matrix_isSquare(struct omMatrix *matrix){

	// return true if the number of rows is equals to the number of columns
	if(matrix->_rows == matrix->_columns)
		return 1;
	else
		return 0;

}

/* return 1 if the matrix contains NaN values, 0 otherwise */
int om_matrix_containsNaN(struct omMatrix *matrix){


	//variables
	int i = 0;
	int j = 0;
	int bool = 0;

	// test if at least one value is NaN
	while(bool != 1 && i < matrix->_rows){
		while(bool != 1 && j < matrix->_columns){

			if(isnan(matrix->_values[i][j]))
				bool = 1;

			j++;
		}
		j=0;
		i++;
	}

	return bool;
}

/* return 1 if the matrix is null, 0 otherwise */
int om_matrix_isNull(struct omMatrix *matrix){

	//variables
	int i = 0;
	int j = 0;
	int bool = 1;

	// test if at least one value is different than 0.0
	while(bool != 0 && i < matrix->_rows){
		while(bool != 0 && j < matrix->_columns){
			if(fabs(matrix->_values[i][j]) > 0.0)
				bool = 0;

			j++;
		}
		i++;
	}

	return bool;

}

/* get the comatrix of a matrix */
void om_matrix_comatrix(struct omMatrix *matrix,struct omMatrix *comatrix){

	// if the matrix is square
	if( om_matrix_isSquare(matrix) == 1 ){

		int NMAX=matrix->_rows;

		 omMatrix b;

		 om_matrix_create(&b,NMAX-1,NMAX-1);

		 for (int q = 0; q < NMAX; q++) {
			 for (int p = 0; p < NMAX; p++) {

				 int m = 0;
				 int n = 0;

				 // Sub matrix
				 for (int i = 0; i < NMAX; i++) {
					 for (int j = 0; j < NMAX; j++) {


						 if (i != q && j != p) {
							 om_matrix_setValue(&b,m,n,om_matrix_getValue(matrix,i,j) );

							 if (n < (NMAX - 2))
								 n++;
							 else {
								 n = 0;
								 m++;
							 }
						 }
					 }
				 }

				 om_matrix_setValue(comatrix,q,p,pow(-1.0,(double)( q + p ) ) * om_matrix_determinant(&b));

			 }
		 }

		 // free submatrix
		 om_matrix_free(&b);

	}


}

/* get the transpose of a matrix */
void om_matrix_transpose(struct omMatrix *matrix,struct omMatrix *transpose){

	// ( a_{ij} )^T = ( a_{ji} )
	for (int i = 0; i < matrix->_rows; i++)
		for (int j = 0; j < matrix->_columns; j++)
			om_matrix_setValue(transpose,j,i,om_matrix_getValue(matrix,i,j));

}

/* get the ith column of a matrix (return a vector) */
void om_matrix_getColumn(struct omMatrix *matrix,int column,struct omVector *out){

	for(int i=0;i<matrix->_rows;i++)
		om_vector_setValue(out,i,om_matrix_getValue(matrix,i,column));

}

/* get the ith row of a matrix (return a vector) */
void om_matrix_getRow(struct omMatrix *matrix,int row,struct omVector *out){


	for(int i=0;i<matrix->_columns;i++)
		om_vector_setValue(out,i,om_matrix_getValue(matrix,row,i));

}

/* set the ith column of a matrix (return a vector) */
void om_matrix_setColumn(struct omMatrix *matrix,int column,struct omVector *in){

	for(int i=0;i<matrix->_rows;i++)
		om_matrix_setValue(matrix,i,column,in->_values[i]);

}

/* set the ith row of a matrix (return a vector) */
void om_matrix_setRow(struct omMatrix *matrix,int row,struct omVector *in){

	for(int i=0;i<matrix->_columns;i++)
		om_matrix_setValue(matrix,row,i,in->_values[i]);

}




/* get the inverse of a matrix if possible */
void om_matrix_inverse(struct omMatrix *matrix,struct omMatrix *inverse){

	// if the matrix is square
	if( om_matrix_isSquare(matrix) == 1 ){

		//variables
		omMatrix L;
		omMatrix U;
		omMatrix I;

		// allocation
		om_matrix_create(&L,matrix->_rows,matrix->_columns);
		om_matrix_create(&U,matrix->_rows,matrix->_columns);
		om_matrix_createIdentity(&I,matrix->_rows);

		// (L,U) <- factoLU(A)
		om_matrix_factorizationLU(matrix,&L,&U);

		for(int i=0;i<matrix->_rows;++i){

			omVector x_i;
			omVector I_i;

			om_vector_create(&x_i,matrix->_rows);
			om_vector_create(&I_i,matrix->_rows);

			om_matrix_getColumn(&I,i,&I_i);
			om_solvingLinearSystemLU(&L,&U,&I_i,&x_i);
			om_matrix_setColumn(inverse,i,&x_i);

			om_vector_free(&x_i);
			om_vector_free(&I_i);
		}


		om_matrix_free(&L);
		om_matrix_free(&U);
		om_matrix_free(&I);
	}



}

/* get the adjugate of a matrix */
void om_matrix_adjugate(struct omMatrix *matrix,struct omMatrix *adjugate){

	omMatrix comatrix;

	// get the comatrix
	om_matrix_create(&comatrix,matrix->_rows,matrix->_columns);
	om_matrix_comatrix(matrix,&comatrix);

	// adj(A) = comatrix(A)^T
	om_matrix_transpose(&comatrix,adjugate);

	// free memory
	om_matrix_free(&comatrix);
}

/* get the LU decomposition of a matrix */
void om_matrix_factorizationLU(struct omMatrix *matrix,struct omMatrix *L,struct omMatrix *U){

	// if the matrix is square
	if( om_matrix_isSquare(matrix) == 1 ){

		//variables
		 int n = matrix->_rows;
		 omMatrix A;

		 // allocation
		 om_matrix_create(&A,n,n);
		 om_matrix_clone(matrix,&A);

		 //*/
		 for(int k=0;k<n;k++){
			 om_matrix_setValue(L,k,k,1.0);

			 for(int i=k+1;i<n;i++){
				 om_matrix_setValue(L,i,k,om_matrix_getValue(&A,i,k)/om_matrix_getValue(&A,k,k));

				 for(int j=k+1;j<n;j++)
					 om_matrix_setValue(&A,i,j,om_matrix_getValue(&A,i,j) - (om_matrix_getValue(L,i,k)*om_matrix_getValue(&A,k,j)) );
			 }

			 for(int j=k;j<n;j++)
				 om_matrix_setValue(U,k,j,om_matrix_getValue(&A,k,j));
		 }
		 /*/
		for(int j=0; j<n; j++)
		{
			for(int i=0; i<n; i++)
			{
			if(i <= j)
			{
				//U[i][j]=A[i][j];
				om_matrix_setValue(U,i,j,om_matrix_getValue(&A,i,j));

				for(int k=0; k<i-1; k++){
					//U[i][j] -= L[i][k]*U[k][j];
					double tmp = om_matrix_getValue(L,i,k)*om_matrix_getValue(U,k,j);
					om_matrix_setValue(U,i,j,om_matrix_getValue(U,i,j) - tmp);

				}

				if(i==j){
					//L[i][j]=1;
					om_matrix_setValue(L,i,j,1.0);
				}

				else{
					//L[i][j]=0;
					om_matrix_setValue(L,i,j,0.0);
				}


			}
			else
			{

				//L[i][j]=A[i][j];
				om_matrix_setValue(L,i,j,om_matrix_getValue(&A,i,j));

				for(int k=0; k <= j-1; k++){

					//L[i][j]-= L[i][k]*U[k][j];
					double tmp = om_matrix_getValue(L,i,k)*om_matrix_getValue(U,k,j);
					om_matrix_setValue(L,i,j,om_matrix_getValue(L,i,j) - tmp);

				}

				if(om_matrix_getValue(U,j,j) == 0.0){
					printf("WTF\n");
				}

				//L[i][j]/=U[j][j];
				om_matrix_setValue(L,i,j,om_matrix_getValue(L,i,j)/om_matrix_getValue(U,j,j));

				//U[i][j]=0;
				om_matrix_setValue(U,i,j,0.0);
			}
		}
	}

		 //*/

		 // free memory
		 om_matrix_free(&A);

	}

}


/* get the LU decomposition of a matrix */
void om_matrix_factorizationQR(struct omMatrix *matrix,struct omMatrix *Q,struct omMatrix *R){

	int m=matrix->_rows;
	int n=matrix->_columns;
	omMatrix A;

	om_matrix_create(&A,m,n);
	om_matrix_clone(matrix,&A);

	for(int k=0;k<n;++k){

		double s = 0.0;

		for(int j=0;j<m;++j)
			s += om_matrix_getValue(&A,j,k)*om_matrix_getValue(&A,j,k);

		om_matrix_setValue(R,k,k,sqrt(s));

		for(int j=0;j<m;++j)
			om_matrix_setValue(Q,j,k, om_matrix_getValue(&A,j,k)/om_matrix_getValue(R,k,k));

		for(int i=k+1;i<n;i++){
			s = 0.0;

			for(int j=0;j<m;++j)
				s += om_matrix_getValue(&A,j,i)*om_matrix_getValue(Q,j,k);

			om_matrix_setValue(R,k,i,s);

			for(int j=0;j<m;++j){
				double value = om_matrix_getValue(&A,j,i) - (om_matrix_getValue(R,k,i) *om_matrix_getValue(Q,j,k) );
				om_matrix_setValue(&A,j,i, value);
			}


		}

	}

	om_matrix_free(&A);

}


/* Cholesky decomposition of a matrix */
void om_matrix_choleskyDecomposition(struct omMatrix *matrix,struct omMatrix *L){

	if( om_matrix_isSquare(matrix) == 1 ){

		int n = matrix->_rows;


		double L00 = sqrt(om_matrix_getValue(matrix,0,0));
		om_matrix_setValue(L,0,0,L00);

		for(int j=1;j<n;++j)
			om_matrix_setValue(L,j,0,om_matrix_getValue(matrix,j,0)/L00);

		for(int i=1;i<n;++i){

			double sum1 = 0.0;

			for(int k=0;k<i;++k)
				sum1 += (om_matrix_getValue(L,i,k)*om_matrix_getValue(L,i,k));

			om_matrix_setValue(L,i,i,sqrt(om_matrix_getValue(matrix,i,i) - sum1));

			for(int j=i+1;j<n;++j){

				double sum2 = 0.0;

				for(int k=0;k<i;++k)
					sum2 += om_matrix_getValue(L,i,k)*om_matrix_getValue(L,j,k);

				om_matrix_setValue(L,j,i,( om_matrix_getValue(matrix,i,j) - sum2)/om_matrix_getValue(L,i,i));

			}
		}
	}
}


/* Schur decomposition of a matrix */
void om_matrix_schurDecomposition(struct omMatrix *matrix,struct omMatrix *T,struct omMatrix *U,int N){


	// first iteration T = A
	om_matrix_clone(matrix,T);

	// first iteration U = I
	for(int i=0;i<matrix->_rows;++i)
		om_matrix_setValue(U,i,i,1.0);

	// for each iteration
	for(int k=0;k < N ;k++){

		//printf("index %d \n",k);

		// variable
		omMatrix Q;
		omMatrix R;
		omMatrix U_tmp;


		// allocation
		om_matrix_create(&U_tmp,matrix->_rows,matrix->_columns);
		om_matrix_create(&Q,matrix->_rows,matrix->_columns);
		om_matrix_create(&R,matrix->_rows,matrix->_columns);

		// QR factorization of T
		om_matrix_factorizationQR(T,&Q,&R);

		// T = R*Q
		om_operator_matrix_mul(&R,&Q,T);

		// U = U*Q
		om_operator_matrix_mul(U,&Q,&U_tmp);
		om_matrix_clone(&U_tmp,U);

		// free memory
		om_matrix_free(&U_tmp);
		om_matrix_free(&R);
		om_matrix_free(&Q);
	}

}


/* create a clone of a matrix */
void om_matrix_clone(struct omMatrix *in,struct omMatrix *out){

	// copy all values in the new matrix
	for(int i=0;i<in->_rows;++i)
		for(int j=0;j<in->_columns;++j)
			om_matrix_setValue(out,i,j,om_matrix_getValue(in,i,j));

}



/* display the matrix values in a terminal */
void om_matrix_display(struct omMatrix *matrix){

	int numberofdecimals = 10;

	for (int i = 0; i < matrix->_rows; i++) {

		if (i == 0)
			printf("/");
		else if (i == matrix->_rows -1)
			printf("\\");
		else
			printf("|");

		for (int j = 0; j < matrix->_columns; j++)
	        	printf(" %.*f ", numberofdecimals, om_matrix_getValue(matrix,i,j));

		if (i == 0)
			printf("\\\n");
		else if (i == matrix->_rows -1)
			printf("/\n");
		else
			printf("|\n");

	}

}



///////////////////////////////////////////////////////
/////             Quaternion class                /////
///////////////////////////////////////////////////////

/* create a new quaternion */
void om_quat_create(struct omQuaternion *quat,double qw,double qx,double qy,double qz){

	quat->_qw = qw;
	quat->_qx = qx;
	quat->_qy = qy;
	quat->_qz = qz;

}

/* get the conjugate of a quaternion */
void om_quat_conjugate(struct omQuaternion *quat,struct omQuaternion *conjugate){

	// q = [ qw , qx , qy , qz ]
	// conj(q) = [ qw , -qx , -qy , -qz ]
	om_quat_create(conjugate,quat->_qw,quat->_qx * (-1.0),quat->_qy * (-1.0),quat->_qz * (-1.0));

}

/* get the inverse of a quaternion */
void om_quat_inverse(struct omQuaternion *quat,struct omQuaternion *inverse){

	// q = [ qw , qx , qy , qz ]
	// inv(q) = conj(q) / ||conj(q)||
	om_quat_conjugate(quat,inverse);
	om_quat_normalize(inverse);

}

/* get the imaginary part of a quaternion */
void om_quat_imaginary(struct omQuaternion *quat,struct omVector *imaginary){

	// q = [ qw , qx , qy , qz ]
	// img(q_ = [ qx , qy , qz ]
	om_vector_setValue(imaginary,0,quat->_qx);
	om_vector_setValue(imaginary,1,quat->_qy);
	om_vector_setValue(imaginary,2,quat->_qz);


}

/* normalize a quaternion */
void om_quat_normalize(struct omQuaternion *quat){

	// get the norm
	double norm = om_quat_norm(quat);

	// divide all values by the norm
	quat->_qw /= norm;
	quat->_qx /= norm;
	quat->_qy /= norm;
	quat->_qz /= norm;


}

/* get the norm of a quaternion */
double om_quat_norm(struct omQuaternion *quat){

	double norm = (quat->_qw*quat->_qw) + (quat->_qx*quat->_qx) + (quat->_qy*quat->_qy) + (quat->_qz*quat->_qz);

	return sqrt(norm);
}


/* display a quaternion into the terminal */
void om_quat_display(struct omQuaternion *quat){

	int numberofdecimals = 10;

	printf("Quaternion : [ %.*f ",numberofdecimals,quat->_qw);
	printf(", %.*f ",numberofdecimals,quat->_qx);
	printf(", %.*f ",numberofdecimals,quat->_qy);
	printf(", %.*f ]\n",numberofdecimals,quat->_qz);

}


/* display a quaternion into the terminal */
void om_quat_clone(struct omQuaternion *in,struct omQuaternion *out){
	om_quat_create(out,in->_qw,in->_qx,in->_qy,in->_qz);
}

///////////////////////////////////////////////////////
/////              Operators                      /////
///////////////////////////////////////////////////////

/* allows to add 2 vectors */
void om_operator_vector_add(struct omVector *a,struct omVector *b,struct omVector *out){

	// if the 2 vectors have the same length
	if (a->_length == b->_length){

		// for all values
		// c_i = a_i +  b_i
		for(int i=0;i<out->_length;i++)
			om_vector_setValue(out,i,a->_values[i] + b->_values[i]);

	}

}

/* allows to subtract 2 vectors */
void om_operator_vector_sub(struct omVector *a,struct omVector *b,struct omVector *out){

	// if the 2 vectors have the same length
	if (a->_length == b->_length){

		// for all values
		// c_i = a_i -  b_i
		for(int i=0;i<out->_length;i++)
			om_vector_setValue(out,i,a->_values[i] - b->_values[i]);

	}


}

/* allows to multiply a vector with a scalar */
void om_operator_vector_scal_mul(struct omVector *a,double b,struct omVector *out){

	// for all values
	// c_i = a_i * b
	for(int i=0;i<out->_length;i++)
		om_vector_setValue(out,i,a->_values[i] * b);

}

/* allows to divide a vector with a scalar */
void om_operator_vector_scal_div(struct omVector *a,double b,struct omVector *out){

	// for all values
	// c_i = a_i / b
	for(int i=0;i<out->_length;i++)
		om_vector_setValue(out,i,a->_values[i] / b);

}

/* compute the outer product of two vector */
void om_operator_vector_outer_product (struct omVector* a, struct omVector* b, struct omMatrix* out)
{
    int i,j;
    int m = a->_length;
    int n = b->_length;

	if (a->_length == b->_length){

	    for(i = 0; i < m; i++) {
	        for (j = 0; j < n; j++) {
	            om_matrix_setValue(out,i,j,om_vector_getValue(a,i) * om_vector_getValue(b,j));
	        }
	    }

    }
}

/* allows to add 2 matrices */
void om_operator_matrix_add(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out){

	// if the 2 matrices have the same number of rows and columns
	if( (a->_rows == b->_rows) && (a->_columns == b->_columns) ){

		// for all values
		// c_ij = a_ij + b_ij
		for(int i=0;i<out->_rows;i++)
			for(int j=0;j<out->_columns;j++)
				om_matrix_setValue(out,i,j,om_matrix_getValue(a,i,j) + om_matrix_getValue(b,i,j));

	}

}

/* allows to subtract 2 matrices */
void om_operator_matrix_sub(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out){

	// if the 2 matrices have the same number of rows and columns
	if( (a->_rows == b->_rows) && (a->_columns == b->_columns) ){

		// for all values
		// c_ij = a_ij - b_ij
		for(int i=0;i<out->_rows;i++)
			for(int j=0;j<out->_columns;j++)
				om_matrix_setValue(out,i,j,om_matrix_getValue(a,i,j) - om_matrix_getValue(b,i,j));

	}


}

/* allows to multiply 2 matrices */
void om_operator_matrix_mul(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out){

	// if the 2 matrices have the same number of rows and columns
	if( a->_columns == b->_rows){

		for(int i=0;i<out->_rows;++i)
			for(int j=0;j<out->_columns;++j){

				double sum=0.0;

				// for all values
				// c_ij = sum (a_ik - b_kj)
				for(int k=0;k<a->_columns;++k)
					 sum += om_matrix_getValue(a,i,k)*om_matrix_getValue(b,k,j);

				om_matrix_setValue(out,i,j,sum);

			}
	}

}


/* allows to multiply a matrix with a vector */
void om_operator_matrix_vector_mul(struct omMatrix *a,struct omVector *b,struct omVector *out){

	if (b->_length == a->_columns){

		for(int i=0;i<out->_length;++i){

			double sum=0.0;

			for(int k=0;k<a->_columns;++k)
				sum += (b->_values[k] * om_matrix_getValue(a,i,k));

			om_vector_setValue(out,i,sum);
		}

	}


}


/* allows to multiply a matrix with a scalar */
void om_operator_matrix_scal_mul(struct omMatrix *a,double b,struct omMatrix *out){

	// for all values
	// c_ij = a_ij * b
	for(int i=0;i<out->_rows;i++)
		for(int j=0;j<out->_columns;j++)
			om_matrix_setValue(out,i,j,om_matrix_getValue(a,i,j)*b);

}

/* allows to divide a matrix with a scalar */
void om_operator_matrix_scal_div(struct omMatrix *a,double b,struct omMatrix *out){

	// for all values
	// c_ij = a_ij / b
	for(int i=0;i<out->_rows;i++)
		for(int j=0;j<out->_columns;j++)
			om_matrix_setValue(out,i,j, om_matrix_getValue(a,i,j)/b);


}

/* allows to add 2 quaternions  */
void om_operator_quat_add(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out){

	om_quat_create(out,a->_qw + b->_qw,a->_qx + b->_qx,a->_qy + b->_qy,a->_qz + b->_qz);


}

/* allows to subtract 2 quaternions  */
void om_operator_quat_sub(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out){

	om_quat_create(out,a->_qw - b->_qw,a->_qx - b->_qx,a->_qy - b->_qy,a->_qz - b->_qz);

}

/* allows to multiply 2 quaternions  */
void om_operator_quat_mul(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out){

	double qw = ( a->_qw*b->_qw - a->_qx*b->_qx - a->_qy*b->_qy - a->_qz*b->_qz);
	double qx = ( a->_qw*b->_qx + a->_qx*b->_qw + a->_qy*b->_qz - a->_qz*b->_qy);
	double qy = ( a->_qw*b->_qy - a->_qx*b->_qz + a->_qy*b->_qw + a->_qz*b->_qx);
	double qz = ( a->_qw*b->_qz + a->_qx*b->_qy - a->_qy*b->_qx + a->_qz*b->_qw);

	om_quat_create(out,qw,qx,qy,qz);

}

/* allows to multiply a quaternion with a scalar */
void om_operator_quat_scal_mul(struct omQuaternion *a,double b,struct omQuaternion *out){

	om_quat_create(out,a->_qw * b,a->_qx * b,a->_qy * b,a->_qz * b);

}

/* allows to divide a quaternion with a scalar */
void om_operator_quat_scal_div(struct omQuaternion *a,double b,struct omQuaternion *out){

	om_quat_create(out,a->_qw / b,a->_qx / b,a->_qy / b,a->_qz / b);
}



///////////////////////////////////////////////////////
/////                Divers                       /////
///////////////////////////////////////////////////////

void om_convert_vector2matrix(struct omVector* a, struct omMatrix* out)
{
    if(out->_columns==1){
          om_matrix_setColumn(out,0,a);
    }
    else  printf("The size of the matrix does not match a vector\n");
}


/* get the cross product of 2 vectors */
void om_vector_crossProduct(struct omVector *a,struct omVector *b,struct omVector *cross){

	// if both vectors have 3 dimensions
	if(a->_length == 3 && b->_length == a->_length){

		double cross_x = (a->_values[1] * b->_values[2]) - (a->_values[2] * b->_values[1]);
		double cross_y = (a->_values[2] * b->_values[0]) - (a->_values[0] * b->_values[2]);
		double cross_z = (a->_values[0] * b->_values[1]) - (a->_values[1] * b->_values[0]);

		om_vector_setValue(cross,0,cross_x);
		om_vector_setValue(cross,1,cross_y);
		om_vector_setValue(cross,2,cross_z);

	}

}

/* get the dot product of 2 vectors */
double om_vector_dotProduct(struct omVector *a,struct omVector *b){

	double dot = -1.0;

	if( b->_length == a->_length ){

		dot = 0.0;

		for(int i =0;i<a->_length;++i)
			dot += b->_values[i]*a->_values[i];

	}

	return dot;
}


/* get the dot product of 2 quaternions */
double om_quat_dotProduct(struct omQuaternion *a,struct omQuaternion *b){

	double dot = 0.0;

	dot += a->_qw*b->_qw;
	dot += a->_qx*b->_qx;
	dot += a->_qy*b->_qy;
	dot += a->_qz*b->_qz;

	return dot;

}

/* allows to solve the linear A*x = b */
void om_solvingLinearSystem(struct omMatrix *A,struct omVector *b,struct omVector *x){

	// variables
	omMatrix L;
	omMatrix U;

	// allocation
	om_matrix_create(&L,A->_rows,A->_columns);
	om_matrix_create(&U,A->_rows,A->_columns);

	// get the LU decomposition of A
	om_matrix_factorizationLU(A,&L,&U);

	// solve the system with the LU decomposition
	om_solvingLinearSystemLU(&L,&U,b,x);

	// free memory
	om_matrix_free(&L);
	om_matrix_free(&U);

}

/* allows to solve the linear A*x = b  with the LU decomposition */
void om_solvingLinearSystemLU(struct omMatrix *L,struct omMatrix *U,struct omVector *b,struct omVector *x){

	// variables
	int n = b->_length;
	omVector d;

	// allocation
	om_vector_create(&d,n);

	// solving Ld = b
	om_vector_setValue(&d,0,b->_values[0]);

	for(int i=1;i<n;++i){
		double acc = 0.0;

		for(int j=0;j<i;++j)
			acc += om_matrix_getValue(L,i,j)*d._values[j];

		om_vector_setValue(&d,i,b->_values[i]-acc);

	}

	// solving Ux = d
	om_vector_setValue(x,n-1,d._values[n-1]/ om_matrix_getValue(U,(n-1),(n-1)) );

	for(int i=n-2;i>=0;--i){

		double acc = 0.0;

		for(int j=i+1;j<n;++j)
			acc += om_matrix_getValue(U,i,j)*x->_values[j];

		om_vector_setValue(x,i, (d._values[i]-acc)/om_matrix_getValue(U,i,i));

	}

	// free memory
	om_vector_free(&d);

}


/*
    This function is an implementation of the least square method.

    Parameters : - struct omMatrix *pX : data matrix
                        - struct omVector *pY : response variable vector (Y=f(X))
                        - struct omVector *pBeta : output vector of the estimation of the best coefficients

    Algorithm : Beta = (X_Transposed . X)^(-1) . X_Transposed . Y

*/

void om_least_square_method (struct omMatrix *pX,struct omVector *pY,struct omVector *pBeta)
{

     /*Declaration of pointers*/

    omMatrix* pX_transposed;
    omMatrix* pX_transposed_Product_X;
    omMatrix* p_Inverse_of_X_transposed_product_X;
    omMatrix* p_Inverse_of_X_transposed_product_X_Multiplied_by_X_transposed;

    /*Check the size of the vectors and matrices for the product*/

    if (pBeta->_length != pY->_length){
        printf("Error : Beta length and Y length are not the same size\n");
        return;
    }

    if (pY->_length != pX->_rows){
        printf("Error : Y length and X number of rows are not the same size\n");
        return;
    }

    /*Dynamic creation of matrices*/

    om_matrix_create(pX_transposed,pX->_columns,pX->_rows); // transposed of the X matrix
    om_matrix_create(pX_transposed_Product_X,pX->_columns,pX->_columns); //  product of the X matrix and its transposed
    om_matrix_create(p_Inverse_of_X_transposed_product_X,pX->_columns,pX->_columns); // Inverse of the product of the X matrix and its transposed
    om_matrix_create(p_Inverse_of_X_transposed_product_X_Multiplied_by_X_transposed, pX->_columns, pX->_rows); // Product of the previous matrix and the transposed of X

   /*Computation of intermediate matrices*/

    om_matrix_transpose(pX, pX_transposed);
    om_operator_matrix_mul(pX_transposed, pX, pX_transposed_Product_X);
    om_matrix_inverse(pX_transposed_Product_X, p_Inverse_of_X_transposed_product_X);
    om_operator_matrix_mul(p_Inverse_of_X_transposed_product_X, pX_transposed, p_Inverse_of_X_transposed_product_X_Multiplied_by_X_transposed);

    /*Computing of the estimation of the coefficients of the Beta matrix*/

    om_operator_matrix_vector_mul(p_Inverse_of_X_transposed_product_X_Multiplied_by_X_transposed, pY, pBeta);

    /*Destruction of dynamically  created  structures*/

    om_matrix_free(pX_transposed);
    om_matrix_free(pX_transposed_Product_X);
    om_matrix_free(p_Inverse_of_X_transposed_product_X);
    om_matrix_free(p_Inverse_of_X_transposed_product_X_Multiplied_by_X_transposed);
}

/*
   This function performs the integration of the function fnct between a and b applying the adaptive Simpson method.

    Parameters : - (*fnct) : pointer to the function to integrate
                        - a, b : lower and upper  limits of the integral
                        - mid : midpoint of the interval integration (has to be computed outside the function because the simpsonadapt call itself recursively
                        - epsilon : maximum allowable error
                        - maxh, minh : maximum and minimum lengths of the subdivision
                        - fa, fb, fmid : fnct values at a,b and mid
                        - *bada, *badb : endpoints of the subinterval on which the calculation failed
                        - success : variable equal to 1 if successful or 0 otherwise

    Source : INTRODUCTION TO NUMERICAL ANALYSIS WITH C PROGRAMS, Attila MATE, Brooklyn College of the City University of New York
*/

double simpsonadapt(double (*fnct)(double), double a, double b, double mid, double epsilon, double maxh, double minh, double fa, double fb, double fmid, double *bada, double *badb, int *success)
{
   double integr1, integr2, integr = 0.0, mid1, mid2, fmid1, fmid2, h, s1, s2;

    h = b - a;

    if ( h >= minh )
    {
        mid1 = (a + mid)/2.0; fmid1 = (*fnct)(mid1);
        mid2 = (mid + b)/2.0; fmid2 = (*fnct)(mid2);
        s1 = (fa+4.0*fmid+fb)*h/6.0;
        s2 = (fa+4.0*fmid1+2.0*fmid+4.0*fmid2+fb)*h/12.0;
        if ( h<=maxh && fabs(s2-s1)<= 15.0*epsilon )
        {
            integr = s2;
            *success = 1;
        }
        else
        {
            integr1 = simpsonadapt(fnct, a, mid, mid1, epsilon/2.0,
            maxh, minh, fa, fmid, fmid1, bada, badb, success);
            if ( *success )
            {
                integr2 = simpsonadapt(fnct, mid, b, mid2, epsilon/2.0,
                maxh, minh, fmid, fb, fmid2, bada, badb, success);

            }
            if ( *success )
            {
                  integr = integr1+integr2;
            }
        }
    }
    else
    {
        *success = 0; *bada = a; *badb = b;
    }
     return integr;
}

// Complex scalar division.
void cdiv(double xr, double xi, double yr, double yi,double* cdivr,double* cdivi) {
      double r,d;
      if (fabs(yr) > fabs(yi)) {
         r = yi/yr;
         d = yr + r*yi;
         (*cdivr) = (xr + r*xi)/d;
         (*cdivi) = (xi - r*xr)/d;
      } else {
         r = yr/yi;
         d = yi + r*yr;
         (*cdivr) = (r*xr + xi)/d;
         (*cdivi) = (r*xi - xr)/d;
      }
   }




double om_maths_erfinv( double y)
{
        double x,z,num,dem; /*working variables */
        /* coefficients in rational expansion */
        double a[4]={ 0.886226899, -1.645349621,  0.914624893, -0.140543331};
        double b[4]={-2.118377725,  1.442710462, -0.329097515,  0.012229801};
        double c[4]={-1.970840454, -1.624906493,  3.429567803,  1.641345311};
        double d[2]={ 3.543889200,  1.637067800};
        if(fabs(y) > 1.0) return (atof("NaN"));  /* This needs IEEE constant*/
        if(fabs(y) == 1.0) return((copysign(1.0,y))*DBL_MAX);
        if( fabs(y) <= CENTRAL_RANGE )
        {
                z = y*y;
                num = (((a[3]*z + a[2])*z + a[1])*z + a[0]);
                dem = ((((b[3]*z + b[2])*z + b[1])*z +b[0])*z + 1.0);
                x = y*num/dem;
        }
        else if( (fabs(y) > CENTRAL_RANGE) && (fabs(y) < 1.0) )
        {
                z = sqrt(-log((1.0-fabs(y))/2.0));
                num = ((c[3]*z + c[2])*z + c[1])*z + c[0];
                dem = (d[1]*z + d[0])*z + 1.0;
                x = (copysign(1.0,y))*num/dem;
        }
        /* Two steps of Newton-Raphson correction */
        x = x - (erf(x) - y)/( (2.0/sqrt(PI))*exp(-x*x));
        x = x - (erf(x) - y)/( (2.0/sqrt(PI))*exp(-x*x));

        return(x);
}

