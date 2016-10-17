#include "algebra.h"


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

/* compute the norm of a matrix */
double om_matrix_norm(struct omMatrix *matrix){

	// set max to zero
	double max = 0.0;

	for(int j=0;j<matrix->_columns;++j){

		// compute the sum of the ith row values
		double sum=0.0;
		for(int i=0;i<matrix->_rows;++i)
			sum += om_matrix_getValue(matrix,i,j);

		// get the max
		max = max < sum ? sum : max;

	}

	//return the norm
	return max;


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

	// variables
	omMatrix D;
	omMatrix Q;
	omMatrix R;
	omMatrix P;
	int m=matrix->_rows;
	int n=matrix->_columns;

	//allocation
	om_matrix_create(&Q,m,n);
	om_matrix_create(&R,m,n);
	om_matrix_create(&D,m,n);
	om_matrix_create(&P,m,n);

	// initialization of D
	om_matrix_clone(matrix,&D);

	// initialization of Q,R and P
	om_matrix_factorizationQR(&D,&Q,&R);
	om_matrix_clone(&Q,&P);

	// for each iteration
	for(int index=0;index<N;index++){

		//some variables
		omMatrix P_tmp;
		om_matrix_create(&P_tmp,matrix->_rows,matrix->_columns);

		// compute D = R*Q
		om_operator_matrix_mul(&R,&Q,&D);

		// QR factorization of D
		om_matrix_factorizationQR(&D,&Q,&R);

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
		(*eigen_values)[l] = om_matrix_getValue(&D,l,l);

		//eigen vector are the column values of matrix P
		om_vector_create(&(*eigen_vectors)[l],4);
		om_matrix_getColumn(&P,l,&(*eigen_vectors)[l]);

	}

	// free memory
	om_matrix_free(&D);
	om_matrix_free(&Q);
	om_matrix_free(&R);
	om_matrix_free(&P);



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

		// allocation
		//om_matrix_create(m_sqrt,matrix->_rows,matrix->_columns);

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
				om_matrix_setValue(&squareD,i,i,sqrt(om_matrix_getValue(&D,i,i)));

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
			if(abs(matrix->_values[i][j]) > 0.0)
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
			om_matrix_setValue(L,j,0,om_matrix_getValue(matrix,0,j)/L00);

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
        if ( h<=maxh && absval(s2-s1)<= 15.0*epsilon )
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
