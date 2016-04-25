/*
 * Matrix.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#include "Matrix.h"

namespace om {


///////////////////////////////////////////////////////
/////        Constructors-Destructors             /////
///////////////////////////////////////////////////////

Matrix::Matrix () {

   _rows = -1;
   _columns = -1;
   _values = 0;

}


Matrix::Matrix (const Matrix& m){

   _rows = m._rows;
   _columns = m._columns;

   _values = (double**)malloc(_rows * sizeof(double*));

   for (int i=0; i<_rows; i++)
       _values[i] = (double*)malloc ( _columns * sizeof(double) );


   for (int i=0; i<_rows; i++)
      for (int j=0; j<_columns; j++)
         _values[i][j] = m._values[i][j];

}


Matrix::Matrix ( int rows, int columns){

   _rows = rows;
   _columns = columns;

   _values = (double**)malloc(_rows * sizeof(double*));

   for (int i=0; i<_rows; i++)
       _values[i] = (double*)malloc ( _columns * sizeof(double) );


   for (int i=0; i<_rows; i++)
      for (int j=0; j<_columns; j++)
         _values[i][j] = 0.0;

}




Matrix::~Matrix (){


	if(_values != 0){
	    for (int i=0; i<_rows; i++){
	    	free(_values[i]);
	       _values[i] = 0;
	    }

	    free(_values);
	    _values = 0;
	}


}

///////////////////////////////////////////////////////
/////              Getters-Setters                /////
///////////////////////////////////////////////////////


Vector Matrix::getColumn(int index){

	Vector column_i(static_cast<double>(_rows));


	for(int i=0;i<_rows;++i)
		column_i.setValue(i,_values[i][index]);

	return column_i;
}



Vector Matrix::getRow(int index){

	Vector row_i(static_cast<double>(_columns));


	for(int i=0;i<_columns;++i)
		row_i.setValue(i,_values[index][i]);

	return row_i;
}


void Matrix::setColumn(int index,const Vector& column){


	if(column.getLength() == _rows){

		for(int i=0;i<_rows;i++)
			_values[i][index] = column.getValue(i);

	}else
		cerr << "ERROR size!" << endl;


}


void Matrix::setRow(int index,const Vector& row){


	if(row.getLength() == _columns){

		for(int i=0;i<_columns;i++)
			_values[index][i] = row.getValue(i);

	}else
		cerr << "ERROR size!" << endl;

}


///////////////////////////////////////////////////////
/////              Methods                        /////
///////////////////////////////////////////////////////


void Matrix::dispose(){

	if(_values != 0){
	    for (int i=0; i<_rows; i++){
	    	free(_values[i]);
	       _values[i] = 0;
	    }

	    free(_values);
	    _values = 0;
	}

}


void Matrix::setValues(int row,int column, Matrix submatrix){

	int n = submatrix.getRows();
	int m = submatrix.getColumns();

    for (int i=0; i<n; i++)
      for (int j=0; j<m; j++)
    	  _values[row+i][column+j] = submatrix._values[i][j];



}


Matrix Matrix::transpose (){

   Matrix res(_columns,_rows);

    for (int i=0; i<_rows; i++)
      for (int j=0; j<_columns; j++)
    	  res.setValue(j,i,_values[i][j]);

   return res;

}


bool Matrix::isNull(){

	bool isNull = true;

	for(int i=0;i<_rows && isNull;++i)
		for(int j=0;j<_columns && isNull;++j){
			isNull = abs(_values[i][j]) < EPSILON;
		}

	return isNull;
}

bool Matrix::isSquare() const{
    return _rows == _columns;
}

Matrix Matrix::clone(){

	Matrix clone = Matrix(_rows,_columns);

    for (int i=0; i<_rows; i++)
      for (int j=0; j<_columns; j++)
    	  clone._values[i][j] = _values[i][j];

    return clone;

}


double Matrix::trace() const{

	double res;
	if (isSquare()){
		res = 0.0;
		for(int i=0; i<_rows;++i)
			res += _values[i][i];
	}else{
		res = -1.0;
	}

	return res;

}



double Matrix::norm(){


	double max = 0.0;

	for(int j=0;j<_columns;++j){

		double sum=0.0;
		for(int i=0;i<_rows;++i)
			sum += _values[i][j];

		max = max < sum ? sum : max;

	}


	return max;

}



Vector solvingSystem( Matrix A,const Vector& b){

	pair<Matrix,Matrix> LU = A.factorisationLU();

	return solvingSystem(LU.first,LU.second,b);


}

Vector solvingSystem(const Matrix& L,const Matrix& U,const Vector& b){

	int n = b.getLength();
	Vector x(n);
	Vector d(n);

	// solving Ld = b
	d.setValue(0,b.getValue(0));

	for(int i=1;i<n;++i){
		double acc = 0.0;
		for(int j=0;j<i;++j){
			acc += L.getValue(i,j)*d.getValue(j);
		}
		d.setValue(i, b.getValue(i) - acc);
	}

	// solving Ux = d
	x.setValue(n-1,d.getValue(n-1)/U.getValue(n-1,n-1));
	for(int i=n-2;i>=0;--i){

		double acc = 0.0;

		for(int j=i+1;j<n;++j){
			acc += U.getValue(i,j)*x.getValue(j);
		}

		x.setValue(i, (d.getValue(i) - acc)/U.getValue(i,i) );
	}


	return x;
}



Matrix Matrix::submatrix(int row, int column,int n, int m){


	Matrix S(n,m);

	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			S._values[i][j] = _values[row+i][column+j];

	return S;
}


double Matrix::determinant() {

   double res;
   if (!isSquare()){
      res = -1.0;
      cerr << "Error !! Impossible operation" <<endl;
   }else{

	   pair<Matrix,Matrix> LU = factorisationLU();

	   double detL = 1.0;
	   double detU = 1.0;


	   for(int i=0;i<_rows;++i){
		   detL *= LU.first.getValue(i,i);
		   detU *= LU.second.getValue(i,i);
	   }

	   res = detL*detU;

   }

   return  res;

}


bool Matrix::containsNan(){

	bool res=false;
	int i=0,j=0;

	while(i<_rows && !res){
		j=0;
		while(j<_columns && !res){
			res = isnan(_values[i][j]) || isinf(_values[i][j]);
			j++;
		}
		i++;
	}


	return res;

}



void Matrix::display(){

	for(int x=0;x<_rows;x++)  // loop 3 times for three lines
	 {
		if (x == 0)
			 std::cout<<"/ ";
		else if (x == _rows -1)
			std::cout<<"\\ ";
		else
			std::cout<<"| ";

		for(int y=0;y<_columns;y++)  // loop for the three elements on the line
		{
			std::cout << prd(_values[x][y],9,13) <<" ";   // display the current element out of the array
		}

		if (x == 0)
			std::cout<<"\\ " <<std::endl;
		else if (x == _rows -1)
			std::cout<<"/ " <<std::endl;
		else
			std::cout<<"|" <<std::endl;
	 }

}

/* Cholesky decomposition of the matrix */
Matrix Matrix::choleskyDecomposition(){

	Matrix L(_rows,_columns);


	if( isSquare() ){
		int n = _rows;

		double sum ;
		L.setValue(0,0, sqrt(abs(_values[0][0])));

		//#pragma omp parallel for
		for(int j=1;j<n;++j)
			L.setValue(j,0,(_values[j][0]/L.getValue(0,0)));

		for(int i=1;i<n;++i){

			sum = 0.0;

			//#pragma omp parallel for default(shared) reduction(+:sum)
			for(int k=0;k<i;++k)
				sum += L.getValue(i,k)*L.getValue(i,k);

			L.setValue(i,i,sqrt(abs(_values[i][i] - sum)));


			for(int j=i+1;j<n;++j){

				sum = 0.0;
				//#pragma omp parallel for default(shared) reduction(+:sum)
				for(int k=0;k<i;++k)
					sum += L.getValue(i,k)*L.getValue(j,k);

				L.setValue(j,i,(_values[i][j]-sum)/L.getValue(i,i)  );

			}
		}


	}else{
		cerr << "Error !! Non-Square Matrix" <<endl;
	}
	return L;
}


/* Schur decomposition of the matrix */
pair<Matrix,Matrix> Matrix::schurDecomposition(){

	pair<Matrix,Matrix> TU;

	int N=75;

	Matrix T(*this);
	Matrix U(_rows,_columns);

	for(int i=0;i<_rows;++i)
		U.setValue(i,i,1.0);

	for(int i=0;i<N;++i){
		pair<Matrix,Matrix> QR = T.factorisationQR();
		T = QR.second*QR.first;
		U = U*QR.first;

	}

	TU = std::make_pair(T,U);

	return TU;
}

/* LU factorisation of the matrix */
pair<Matrix,Matrix> Matrix::factorisationLU(){


	 int n = _rows;

	 Matrix L(_rows,_columns);
	 Matrix U(_rows,_columns);
	 Matrix A_clone = clone();

	 pair<Matrix,Matrix> LU;

	 for(int k=0;k<n;k++){
		 L.setValue(k,k,1.0);

		 for(int i=k+1;i<n;i++){
			 L.setValue(i,k,A_clone.getValue(i,k)/A_clone.getValue(k,k) );

			 for(int j=k+1;j<n;j++)
				 A_clone.setValue(i,j, A_clone.getValue(i,j) - L.getValue(i,k)*A_clone.getValue(k,j));

		 }

		 for(int j=k;j<n;j++){
			 U.setValue(k,j,A_clone.getValue(k,j));
		 }
	 }

	LU = std::make_pair(L,U);

	return LU;

}

/* QR decomposition of the matrix */
pair<Matrix,Matrix> Matrix::factorisationQR(){

	pair<Matrix,Matrix> QR;

	Matrix Q(_rows,_columns);
	Matrix R(_rows,_columns);
	Matrix A(clone());

	int m=_rows;
	int n=_columns;

	for(int k=0;k<n;++k){

		double s = 0.0;

		for(int j=0;j<m;++j)
			s += A._values[j][k]* A._values[j][k];

		R._values[k][k] = sqrt(s);

		for(int j=0;j<m;++j)
			Q._values[j][k] = A._values[j][k]/R._values[k][k];

		for(int i=k+1;i<n;i++){
			s = 0.0;

			for(int j=0;j<m;++j)
				s += A._values[j][i]*Q._values[j][k];

			R._values[k][i] = s;

			for(int j=0;j<m;++j)
				A._values[j][i] -= R._values[k][i]*Q._values[j][k];

		}

	}


	QR = std::make_pair(Q,R);

	return QR;

}

/* compute the inverse of the matrix */
Matrix Matrix::inverse(){

	 Matrix inv(_rows,_columns);

	 if(isSquare()){

		pair<Matrix,Matrix> LU = factorisationLU();

		Matrix L(LU.first);
		Matrix U(LU.second);

		Matrix I(_rows,_columns);
		for(int i=0;i<_rows;++i)
			I.setValue(i,i,1.0);

		for(int i=0;i<_rows;i++){

			Vector x_i = solvingSystem(L,U,I.getColumn(i));
			inv.setColumn(i,x_i);

		}

	 }else{

		 cerr << "Error !! Impossible inversion" <<endl;

	 }


	 return inv;

}


Matrix Matrix::adjugate(){
	return comatrix().transpose();
}


/* allows to compute the cofactor matrix */
 Matrix Matrix::comatrix(){

	 int NMAX=_rows;

	 Matrix comatrix(NMAX,NMAX);
	 Matrix b(NMAX-1,NMAX-1);

	 for (int q = 0; q < NMAX; q++) {
		 for (int p = 0; p < NMAX; p++) {

			 int m = 0;
			 int n = 0;

			 //sous matrice
			 for (int i = 0; i < NMAX; i++) {
				 for (int j = 0; j < NMAX; j++) {


					 if (i != q && j != p) {
						 b._values[m][n] = _values[i][j];

						 if (n < (NMAX - 2))
							 n++;
						 else {
							 n = 0;
							 m++;
						 }
					 }
				 }
			 }

			 comatrix._values[q][p] = pow(-1.0,static_cast<double>( q + p ) ) * b.determinant();
		 }
	 }


	 return comatrix;
 }


///////////////////////////////////////////////////////
/////              Operators                      /////
///////////////////////////////////////////////////////


Matrix& Matrix::operator=( const Matrix& m){

    if (this != &m) {

    	if(_rows != -1){

            for (int i=0; i<_rows; i++){
               delete [] _values[i];
               _values[i] = 0;
            }

            delete [] _values;
            _values = 0;

    	}

    	_rows = m._rows;
    	_columns = m._columns;


    	_values = (double**)malloc(_rows * sizeof(double*));

	   for (int i=0; i<_rows; i++)
		   _values[i] = (double*)malloc ( _columns * sizeof(double) );

    	for(int i=0;i<_rows;++i)
    		for(int j=0;j<_columns;++j)
    			_values[i][j] = m._values[i][j];

    }

    return *this;

}


Matrix Matrix::operator+(const Matrix& m){


	Matrix res(_rows,_columns);
	if( _rows == m._rows && _columns == m._columns ){


		for(int i=0;i<_rows;++i)
			for(int j=0;j<_columns;++j){
				res._values[i][j] =  _values[i][j] + m._values[i][j];
			}
	}else{
		 cerr << "Error !! Matrix Impossible addition" <<endl;
	}

	return res;
}

Matrix Matrix::operator-(const Matrix& m){

	Matrix res(_rows,_columns);
	if( _rows == m._rows && _columns == m._columns ){


		for(int i=0;i<_rows;++i)
			for(int j=0;j<_columns;++j)
				res._values[i][j] = _values[i][j] - m._values[i][j];

	}else{
			 cerr << "Error !! Matrix Impossible soustraction" <<endl;
		}

	return res;
}


Matrix Matrix::operator*(const Matrix& m){


	Matrix res(_rows,m._columns);

	if( _columns == m._rows){
		for(int i=0;i<res._rows;++i)
			for(int j=0;j<res._columns;++j){

				double sum=0.0;
				//#pragma omp parallel for default(shared) reduction(+:sum)
				for(int k=0;k<_columns;++k){
					 sum += _values[i][k]*m._values[k][j];
				}

				res._values[i][j] = sum ;
			}

	}else{
			 cerr << "Error !! Matrix Impossible multiplication" <<endl;
		}

	return res;

}


Matrix Matrix::operator*(const double& alpha){


	Matrix res(_rows,_columns);


	for(int i=0;i<res._rows;++i)
		for(int j=0;j<res._columns;++j)
				res._values[i][j] = _values[i][j]*alpha;

	return res;

}

Matrix Matrix::operator/(const double& alpha){

	Matrix res(_rows,_columns);

	for(int i=0;i<res._rows;++i)
		for(int j=0;j<res._columns;++j)
				res._values[i][j] = _values[i][j]/alpha;

	return res;

}






///////////////////////////////////////////////////////
/////       Operator Matrix Vector                /////
///////////////////////////////////////////////////////

Vector operator*(const Vector& v, const Matrix& m){

	Vector res(v.getLength());

	for(int i=0;i<res.getLength();++i)
			for(int k=0;k<m.getRows();++k){
				double tmp = res.getValue(i) + (v.getValue(k) * m.getValue(k,i));
				res.setValue(i, tmp);
			}

	return res;
}


Vector operator*(const Matrix& m,const Vector& v){

	Vector res(m.getRows());



	if (v.getLength() == m.getColumns()){


		for(int i=0;i<res.getLength();++i){
			res.setValue(i, 0.0);

			double sum=0.0;

			for(int k=0;k<m.getColumns();++k){
				sum += (v.getValue(k) * m.getValue(i,k));
			}

			res.setValue(i, sum);
		}

	}else{

		cerr << "ERROR!! Incompatible dimension " << v <<endl;
	}



	return res;

}


} /* namespace isf */
