#include "algebra.h"


///////////////////////////////////////////////////////
/////             Vector class                    /////
///////////////////////////////////////////////////////

/* initialization of a vector */
void om_vector_create(struct omVector *vector,int length,...){

	vector->_length = length;
	vector->_values = (double*)malloc(length*sizeof(double));

	va_list args;

	va_start ( args, length );

	for ( int i = 0; i < vector->_length ; i++ ){

		double value =  va_arg ( args, double );
		if(value)
			vector->_values[i] = value;
		else
			vector->_values[i] = 0.0;
	}

	va_end(args);
}


void om_vector_setValues(struct omVector *vector,int length,...){

	va_list args;

	va_start ( args, length );

	for ( int i = 0; i < vector->_length ; i++ ){

		double value =  va_arg ( args, double );
		if(value)
			vector->_values[i] = value;
		else
			vector->_values[i] = 0.0;
	}

	va_end(args);

}


void om_vector_setValue(struct omVector *vector,int index,double value){

	vector->_values[index]=value;
}

double om_vector_getValue(struct omVector *vector,int index){
	return vector->_values[index];
}


void om_vector_clone(struct omVector *in,struct omVector *out){
	om_vector_create(out,in->_length);

	for ( int i = 0; i < in->_length ; i++ )
		om_vector_setValue(out,i,in->_values[i]);
}




double om_vector_rms(struct omVector *vector){

	double rms = 0.0;

	for(int i = 0; i < vector->_length;++i)
		rms += vector->_values[i]*vector->_values[i];

	rms /= (double)(vector->_length);

	return sqrt(rms);

}

double om_vector_norm(struct omVector *vector){

	double norm = 0.0;

	for(int i = 0; i < vector->_length;++i)
		norm += vector->_values[i]*vector->_values[i];

	return sqrt(norm);

}

void om_vector_normalize(struct omVector *vector){

	double norm = om_vector_norm(vector);

	for(int i = 0; i < vector->_length;++i)
		om_vector_setValue(vector,i,vector->_values[i]/norm);

}

void om_vector_display(struct omVector *vector){

	int numberofdecimals = 10;

	printf("[");
	for (int i = 0; i < vector->_length; i++)
		printf(" %.*f ", numberofdecimals, vector->_values[i]);
	printf("]\n");
}

void om_vector_dispose(struct omVector *vector){

	if(vector->_values != 0){
		free(vector->_values);
		vector->_values = 0;
	}

}


///////////////////////////////////////////////////////
/////             Matrix class                    /////
///////////////////////////////////////////////////////

void om_matrix_create(struct omMatrix *matrix,int rows,int columns){

	matrix->_rows = rows;
	matrix->_columns = columns;
	matrix->_values = (double**)malloc( matrix->_rows* sizeof(double*));


   for (int i=0; i<rows; i++){
	   matrix->_values[i] = (double*)malloc( matrix->_columns* sizeof(double));
	   for (int j=0; j<columns; j++)
		   om_matrix_setValue(matrix,i,j,0.0);

   }
}


void om_matrix_dispose(struct omMatrix *matrix){

	if(matrix->_values != 0){

	    for (int i=0; i<matrix->_rows; i++){
	    	free(matrix->_values[i]);
	    	matrix->_values[i] = 0;
	    }

	    free(matrix->_values);
	    matrix->_values = 0;
	}

}



void om_matrix_setValue(struct omMatrix *matrix,int i,int j,double value){
	matrix->_values[i][j] = value;
}

double om_matrix_getValue(struct omMatrix *matrix,int i,int j){
	return matrix->_values[i][j];
}

double om_matrix_norm(struct omMatrix *matrix){

	double max = 0.0;

	for(int j=0;j<matrix->_columns;++j){

		double sum=0.0;
		for(int i=0;i<matrix->_rows;++i)
			sum += om_matrix_getValue(matrix,i,j);

		max = max < sum ? sum : max;

	}


	return max;


}

double om_matrix_determinant(struct omMatrix *matrix){

	   double determinant = -1.0;

	   if( om_matrix_isSquare(matrix) == TRUE ){

			omMatrix L;
			omMatrix U;

			om_matrix_factorizationLU(matrix,&L,&U);

			double detL = 1.0;
			double detU = 1.0;

		   for(int i=0;i<matrix->_rows;++i){
			   detL *= om_matrix_getValue(&L,i,i);
			   detU *= om_matrix_getValue(&U,i,i);
		   }

		   determinant = detL*detU;

	   }

	   return  determinant;

}

double om_matrix_trace(struct omMatrix *matrix){

	double trace = 0.0;

	if( om_matrix_isSquare(matrix) == TRUE ){

		for (int i = 0; i < matrix->_rows; i++)
			trace += om_matrix_getValue(matrix,i,i);

	}

	return trace;
}


void om_matrix_getEingenValues(struct omMatrix *matrix,struct omVector **eigen_vectors,double **eigen_values,int N){


	omMatrix D;
	omMatrix Q;
	omMatrix R;
	omMatrix P;

	om_matrix_clone(matrix,&D);
	om_matrix_factorizationQR(&D,&Q,&R);

	om_matrix_clone(&Q,&P);

	omMatrix P_tmp;

	for(int i=0;i<N;i++){

		om_operator_matrix_mul(&R,&Q,&D);
		om_matrix_factorizationQR(&D,&Q,&R);
		om_operator_matrix_mul(&P,&Q,&P_tmp);
		om_matrix_clone(&P_tmp,&P);

	}

	(*eigen_values) = (double*)malloc(matrix->_columns*sizeof(double));
	(*eigen_vectors) = (omVector*)malloc(matrix->_columns*sizeof(omVector));

	for(int l=0;l<matrix->_columns;++l){

		(*eigen_values)[l] = om_matrix_getValue(&D,l,l);
		om_matrix_getColumn(&P,l,&(*eigen_vectors)[l]);

	}

	om_matrix_dispose(&D);
	om_matrix_dispose(&Q);
	om_matrix_dispose(&R);
	om_matrix_dispose(&P);
	om_matrix_dispose(&P_tmp);


}

void om_matrix_skewSymetricMatrix(struct omVector *in,struct omMatrix *out){

	om_matrix_create(out,3,3);

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

void om_matrix_createIdentity(struct omMatrix *I,int n){

	om_matrix_create(I,n,n);

	for(int i=0;i<n;i++)
		om_matrix_setValue(I,i,i,1.0);

}



void om_matrix_exponantial(struct omMatrix *matrix,struct omMatrix *exp,int N){


	if(matrix->_columns == matrix->_rows){

		omMatrix acc;
		double acc_n = 1.0;

		om_matrix_clone(matrix,&acc);

		om_matrix_createIdentity(exp,matrix->_rows);
		om_operator_matrix_add(exp,matrix,exp);

		for (int i = 2 ; i<N; ++i){

			acc_n = acc_n*(double)(i);
			omMatrix tmp;

			om_operator_matrix_mul(&acc,matrix,&tmp);
			om_matrix_clone(&tmp,&acc);
			om_operator_matrix_const_div(&tmp,acc_n,&tmp);
			om_operator_matrix_add(exp,&tmp,exp);

			om_matrix_dispose(&tmp);
		}


		om_matrix_dispose(&acc);

	}


}

void om_matrix_squareRoot(struct omMatrix *matrix,struct omMatrix *m_sqrt){

	if(matrix->_rows == matrix->_columns){

		om_matrix_create(m_sqrt,matrix->_rows,matrix->_columns);

		if(matrix->_rows == 1){

			om_matrix_setValue(m_sqrt,0,0,sqrt(  om_matrix_getValue(matrix,0,0)));

		}else if(matrix->_rows == 2){

			omMatrix I;
			om_matrix_createIdentity(&I,2);

			double det = sqrt(om_matrix_determinant(matrix));
			double tmp = sqrt(om_matrix_trace(matrix) + (2.0*det));

			om_operator_matrix_const_mul(&I,det,&I);
			om_operator_matrix_add(matrix,&I,m_sqrt);

			om_operator_matrix_const_div(m_sqrt,tmp,m_sqrt);

			om_matrix_dispose(&I);

		}else{

			int N=20;

			omMatrix D;
			omMatrix Q;
			omMatrix R;
			omMatrix P;
			omMatrix P_inv;
			omMatrix S_tmp;
			omMatrix squareD;
			omMatrix P_tmp;

			om_matrix_clone(matrix,&D);
			om_matrix_factorizationQR(&D,&Q,&R);
			om_matrix_clone(&Q,&P);
			om_matrix_create(&P_tmp,matrix->_rows,matrix->_rows);

			for(int i=0;i<N;i++){

				om_operator_matrix_mul(&R,&Q,&D);
				om_matrix_factorizationQR(&D,&Q,&R);
				om_operator_matrix_mul(&P,&Q,&P_tmp);

				om_matrix_clone(&P_tmp,&P);

			}

			om_matrix_create(&squareD,D._rows,D._columns);
			om_matrix_create(&S_tmp,matrix->_rows,matrix->_rows);

			for(int i=0;i< squareD._rows;i++)
				om_matrix_setValue(&squareD,i,i,sqrt(om_matrix_getValue(&D,i,i)));

			om_matrix_inverse(&P,&P_inv);

			om_operator_matrix_mul(&P,&squareD,&S_tmp);
			om_operator_matrix_mul(&S_tmp,&P_inv,m_sqrt);

			om_matrix_dispose(&P);
			om_matrix_dispose(&P_inv);
			om_matrix_dispose(&S_tmp);
			om_matrix_dispose(&D);
			om_matrix_dispose(&R);
			om_matrix_dispose(&Q);
			om_matrix_dispose(&squareD);
			om_matrix_dispose(&P_tmp);

		}

	}

}

bool om_matrix_isSquare(struct omMatrix *matrix){

	if(matrix->_rows == matrix->_columns)
		return TRUE;
	else
		return FALSE;

}

bool om_matrix_containsNaN(struct omMatrix *matrix){

	return FALSE;
}

bool om_matrix_isNull(struct omMatrix *matrix){

	return FALSE;

}

void om_matrix_comatrix(struct omMatrix *matrix,struct omMatrix *comatrix){

	if( om_matrix_isSquare(matrix) == TRUE ){
		 int NMAX=matrix->_rows;

		 omMatrix b;
		 om_matrix_create(comatrix,NMAX,NMAX);
		 om_matrix_create(&b,NMAX-1,NMAX-1);

		 for (int q = 0; q < NMAX; q++) {
			 for (int p = 0; p < NMAX; p++) {

				 int m = 0;
				 int n = 0;

				 //sous matrice
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

		 om_matrix_dispose(&b);

	}


}

void om_matrix_transpose(struct omMatrix *matrix,struct omMatrix *transpose){

	om_matrix_create(transpose,matrix->_columns,matrix->_rows);

	for (int i = 0; i < matrix->_rows; i++)
		for (int j = 0; j < matrix->_columns; j++)
			om_matrix_setValue(transpose,j,i,om_matrix_getValue(matrix,i,j));

}


void om_matrix_getColumn(struct omMatrix *matrix,int column,struct omVector *out){

	om_vector_create(out,matrix->_rows);

	for(int i=0;i<matrix->_rows;i++)
		om_vector_setValue(out,i,om_matrix_getValue(matrix,i,column));

}

void om_matrix_getRow(struct omMatrix *matrix,int row,struct omVector *out){

	om_vector_create(out,matrix->_columns);

	for(int i=0;i<matrix->_columns;i++)
		om_vector_setValue(out,i,om_matrix_getValue(matrix,row,i));

}

void om_matrix_setColumn(struct omMatrix *matrix,int column,struct omVector *in){

	for(int i=0;i<matrix->_rows;i++)
		om_matrix_setValue(matrix,i,column,in->_values[i]);

}

void om_matrix_setRow(struct omMatrix *matrix,int row,struct omVector *in){

	for(int i=0;i<matrix->_columns;i++)
		om_matrix_setValue(matrix,row,i,in->_values[i]);

}





void om_matrix_inverse(struct omMatrix *matrix,struct omMatrix *inverse){

	if( om_matrix_isSquare(matrix) == TRUE ){

		omMatrix L;
		omMatrix U;
		omMatrix I;

		om_matrix_create(inverse,matrix->_rows,matrix->_columns);
		om_matrix_createIdentity(&I,matrix->_rows);
		om_matrix_factorizationLU(matrix,&L,&U);

		/*
		printf("\nmatrix L \n");
		om_matrix_display(&L);

		printf("\nmatrix U \n");
		om_matrix_display(&U);
		 */

		for(int i=0;i<matrix->_rows;++i){

			omVector x_i;
			omVector I_i;

			om_matrix_getColumn(&I,i,&I_i);
			om_solvingLinearSystemLU(&L,&U,&I_i,&x_i);
			om_matrix_setColumn(inverse,i,&x_i);
			om_vector_dispose(&x_i);
			om_vector_dispose(&I_i);
		}


		om_matrix_dispose(&L);
		om_matrix_dispose(&U);
		om_matrix_dispose(&I);
	}



}

void om_matrix_adjugate(struct omMatrix *matrix,struct omMatrix *adjugate){

	omMatrix comatrix;
	om_matrix_comatrix(matrix,&comatrix);

	om_matrix_transpose(&comatrix,adjugate);

	om_matrix_dispose(&comatrix);
}

void om_matrix_factorizationLU(struct omMatrix *matrix,struct omMatrix *L,struct omMatrix *U){

	if( om_matrix_isSquare(matrix) == TRUE ){

		 int n = matrix->_rows;

		 omMatrix A;
		 om_matrix_create(L,n,n);
		 om_matrix_create(U,n,n);
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

		 om_matrix_dispose(&A);

	}

}


void om_matrix_factorizationQR(struct omMatrix *matrix,struct omMatrix *Q,struct omMatrix *R){

	int m=matrix->_rows;
	int n=matrix->_columns;

	omMatrix A;

	om_matrix_create(Q,m,n);
	om_matrix_create(R,m,n);
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

			for(int j=0;j<m;++j)
				om_matrix_setValue(&A,j,i, om_matrix_getValue(&A,j,i) - (om_matrix_getValue(R,k,i) *om_matrix_getValue(R,j,k) ) );

		}

	}

	om_matrix_dispose(&A);

}


/* Cholesky decomposition of a matrix */
void om_matrix_choleskyDecomposition(struct omMatrix *matrix,struct omMatrix *L){

	if( om_matrix_isSquare(matrix) == TRUE ){

		int n = matrix->_rows;

		om_matrix_create(L,n,n);

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
void om_matrix_schurDecomposition(struct omMatrix *matrix,struct omMatrix *T,struct omMatrix *U){

	int N=75;


}


/* create a clone of a matrix */
void om_matrix_clone(struct omMatrix *in,struct omMatrix *out){

	om_matrix_create(out,in->_rows,in->_columns);

	for(int i=0;i<in->_rows;++i)
		for(int j=0;j<in->_columns;++j)
			om_matrix_setValue(out,i,j,om_matrix_getValue(in,i,j));

}




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


void om_quat_create(struct omQuaternion *quat,double qw,double qx,double qy,double qz){

	quat->_qw = qw;
	quat->_qx = qx;
	quat->_qy = qy;
	quat->_qz = qz;

}

void om_quat_conjugate(struct omQuaternion *quat,struct omQuaternion *conjugate){
	om_quat_create(conjugate,quat->_qw,quat->_qx * (-1.0),quat->_qy * (-1.0),quat->_qz * (-1.0));
}

void om_quat_inverse(struct omQuaternion *quat,struct omQuaternion *inverse){

	om_quat_conjugate(quat,inverse);
	om_quat_normalize(inverse);

}

void om_quat_imaginary(struct omQuaternion *quat,struct omVector *imaginary){

	om_vector_create(imaginary,3);

	om_vector_setValue(imaginary,0,quat->_qx);
	om_vector_setValue(imaginary,1,quat->_qy);
	om_vector_setValue(imaginary,2,quat->_qz);


}

void om_quat_normalize(struct omQuaternion *quat){

	double norm = om_quat_norm(quat);

	quat->_qw /= norm;
	quat->_qx /= norm;
	quat->_qy /= norm;
	quat->_qz /= norm;


}

double om_quat_norm(struct omQuaternion *quat){

	double norm = (quat->_qw*quat->_qw) + (quat->_qx*quat->_qx) + (quat->_qy*quat->_qy) + (quat->_qz*quat->_qz);
	return sqrt(norm);
}


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

void om_operator_vector_add(struct omVector *a,struct omVector *b,struct omVector *out){

	if (a->_length == b->_length){

		if(out->_values == NULL)
			om_vector_create(out,a->_length);

		for(int i=0;i<out->_length;i++)
			om_vector_setValue(out,i,a->_values[i] + b->_values[i]);

	}

}


void om_operator_vector_sub(struct omVector *a,struct omVector *b,struct omVector *out){

	if (a->_length == b->_length){

		if(out->_values == NULL)
			om_vector_create(out,a->_length);

		for(int i=0;i<out->_length;i++)
			om_vector_setValue(out,i,a->_values[i] - b->_values[i]);

	}


}

void om_operator_vector_const_mul(struct omVector *a,double b,struct omVector *out){

	if(out->_values == NULL)
		om_vector_create(out,a->_length);

	for(int i=0;i<out->_length;i++)
		om_vector_setValue(out,i,a->_values[i] * b);

}

void om_operator_vector_const_div(struct omVector *a,double b,struct omVector *out){

	if(out->_values == NULL)
		om_vector_create(out,a->_length);

	for(int i=0;i<out->_length;i++)
		om_vector_setValue(out,i,a->_values[i] / b);

}


void om_operator_matrix_add(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out){

	if( (a->_rows == b->_rows) && (a->_columns == b->_columns) ){

		if(out->_values == NULL){
			om_matrix_create(out,a->_rows,a->_columns);
		}


		for(int i=0;i<out->_rows;i++)
			for(int j=0;j<out->_columns;j++)
				om_matrix_setValue(out,i,j,om_matrix_getValue(a,i,j) + om_matrix_getValue(b,i,j));

	}

}

void om_operator_matrix_sub(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out){

	if( (a->_rows == b->_rows) && (a->_columns == b->_columns) ){

		if(out->_values == NULL){
			om_matrix_create(out,a->_rows,a->_columns);
		}

		for(int i=0;i<out->_rows;i++)
			for(int j=0;j<out->_columns;j++)
				om_matrix_setValue(out,i,j,om_matrix_getValue(a,i,j) - om_matrix_getValue(b,i,j));

	}


}

void om_operator_matrix_mul(struct omMatrix *a,struct omMatrix *b,struct omMatrix *out){


	if( a->_columns == b->_rows){

		if(out->_values == NULL)
			om_matrix_create(out,a->_rows,b->_columns);

		for(int i=0;i<out->_rows;++i)
			for(int j=0;j<out->_columns;++j){

				double sum=0.0;


				for(int k=0;k<a->_columns;++k)
					 sum += om_matrix_getValue(a,i,k)*om_matrix_getValue(b,k,j);


				om_matrix_setValue(out,i,j,sum);

			}
	}

}

void om_operator_matrix_const_mul(struct omMatrix *a,double b,struct omMatrix *out){

	if(out->_values == NULL)
		om_matrix_create(out,a->_rows,a->_columns);

	for(int i=0;i<out->_rows;i++)
		for(int j=0;j<out->_columns;j++)
			om_matrix_setValue(out,i,j,om_matrix_getValue(a,i,j)*b);

}

void om_operator_matrix_const_div(struct omMatrix *a,double b,struct omMatrix *out){

	if(out->_values == NULL)
		om_matrix_create(out,a->_rows,a->_columns);

	for(int i=0;i<out->_rows;i++)
		for(int j=0;j<out->_columns;j++)
			om_matrix_setValue(out,i,j, om_matrix_getValue(a,i,j)/b);


}

void om_operator_matrix_vector_mul(struct omMatrix *a,struct omVector *b,struct omVector *out){

	if (b->_length == a->_columns){

		if(out->_values == NULL)
			om_vector_create(out,b->_length);

		for(int i=0;i<out->_length;++i){

			double sum=0.0;

			for(int k=0;k<a->_columns;++k)
				sum += (b->_values[k] * om_matrix_getValue(a,i,k));

			om_vector_setValue(out,i,sum);
		}

	}


}



void om_operator_quat_add(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out){

	om_quat_create(out,a->_qw + b->_qw,a->_qx + b->_qx,a->_qy + b->_qy,a->_qz + b->_qz);


}

void om_operator_quat_sub(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out){

	om_quat_create(out,a->_qw - b->_qw,a->_qx - b->_qx,a->_qy - b->_qy,a->_qz - b->_qz);

}

void om_operator_quat_mul(struct omQuaternion *a,struct omQuaternion *b,struct omQuaternion *out){

	double qw = ( a->_qw*b->_qw - a->_qx*b->_qx - a->_qy*b->_qy - a->_qz*b->_qz);
	double qx = ( a->_qw*b->_qx + a->_qx*b->_qw + a->_qy*b->_qz - a->_qz*b->_qy);
	double qy = ( a->_qw*b->_qy - a->_qx*b->_qz + a->_qy*b->_qw + a->_qz*b->_qx);
	double qz = ( a->_qw*b->_qz + a->_qx*b->_qy - a->_qy*b->_qx + a->_qz*b->_qw);

	om_quat_create(out,qw,qx,qy,qz);

}

void om_operator_quat_const_mul(struct omQuaternion *a,double b,struct omQuaternion *out){

	om_quat_create(out,a->_qw * b,a->_qx * b,a->_qy * b,a->_qz * b);

}

void om_operator_quat_const_div(struct omQuaternion *a,double b,struct omQuaternion *out){

	om_quat_create(out,a->_qw / b,a->_qx / b,a->_qy / b,a->_qz / b);
}




///////////////////////////////////////////////////////
/////                Divers                       /////
///////////////////////////////////////////////////////


void om_vector_crossProduct(struct omVector *a,struct omVector *b,struct omVector *cross){

	if(a->_length == 3 && b->_length == a->_length){

		om_vector_create(cross,3);

		double cross_x = (a->_values[1] * b->_values[2]) - (a->_values[2] * b->_values[1]);
		double cross_y = (a->_values[2] * b->_values[0]) - (a->_values[0] * b->_values[2]);
		double cross_z = (a->_values[0] * b->_values[1]) - (a->_values[1] * b->_values[0]);

		om_vector_setValue(cross,0,cross_x);
		om_vector_setValue(cross,1,cross_y);
		om_vector_setValue(cross,2,cross_z);

	}

}

double om_vector_dotProduct(struct omVector *a,struct omVector *b){

	double dot = -1.0;

	if( b->_length == a->_length ){

		dot = 0.0;

		for(int i =0;i<a->_length;++i)
			dot += b->_values[i]*a->_values[i];

	}

	return dot;
}

double om_quat_dotProduct(struct omQuaternion *a,struct omQuaternion *b){

	double dot = 0.0;

	dot += a->_qw*b->_qw;
	dot += a->_qx*b->_qx;
	dot += a->_qy*b->_qy;
	dot += a->_qz*b->_qz;

	return dot;

}


void om_solvingLinearSystem(struct omMatrix *A,struct omVector *b,struct omVector *x){

	omMatrix L;
	omMatrix U;

	om_matrix_factorizationLU(A,&L,&U);
	om_solvingLinearSystemLU(&L,&U,b,x);

	om_matrix_dispose(&L);
	om_matrix_dispose(&U);

}

void om_solvingLinearSystemLU(struct omMatrix *L,struct omMatrix *U,struct omVector *b,struct omVector *x){

	int n = b->_length;
	omVector d;

	om_vector_create(x,n);
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

	om_vector_dispose(&d);

}



