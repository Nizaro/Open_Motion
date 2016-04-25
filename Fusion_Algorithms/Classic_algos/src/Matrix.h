/*
 * Matrix.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "Vector.h"

namespace om {


class Matrix
{

public:

	//Constructors-Destructor

	/**
	 * default constructor. It's an empty matrix
	 */
   Matrix();

   /**
    * initialize null matrix in M(rows,columns)
    *
    * @param rows number of rows
    * @param columns number of columns
    *
    */
   Matrix(int rows,int columns);

   /**
    * Create a copy of a matrix m
    *
    * @param a matrix m
    */
   Matrix (const Matrix& m);

   ~Matrix();



   //getters
   int getRows() const {return _rows;}
   int getColumns() const {return _columns;}
   double getValue(int row,int column) const { return _values[row][column];}

   Vector getColumn(int index);
   Vector getRow(int index);

   //setters
   void setRows( int rows ) {_rows = rows;}
   void setColumns( int columns ) {_columns = columns;}
   void setValue( int rows,int column, double value ) { _values[rows][column] = value;}

   void setValues(int row,int column, Matrix submatrix);

   void setColumn(int index,const Vector&  column);
   void setRow(int index,const Vector&  row);


   //operators
   Matrix operator+(const Matrix& m);
   Matrix operator-(const Matrix& m);
   Matrix operator*(const Matrix& m);
   Matrix operator*(const double& alpha);
   Matrix operator/(const double& alpha);
   Matrix& operator=(const Matrix& m);


   //Diplayers
   void display() ;

   // methods

   /**
    * If the matrix is inversible (det != 0), it returns the inverse of the matrix
    *@return the inverse of the matrix
    */
   Matrix inverse();

   /**
    * The Schur decomposition of the matrix
    * @return a pair of matrices T  and U such as A=UTU*
    */
   pair<Matrix,Matrix> schurDecomposition();

   /**
    * The Cholesky decomposition of the matrix
    *
    */
   Matrix choleskyDecomposition();

   /**
    * The QR factorisation of the matrix
    *
    */
   pair<Matrix,Matrix> factorisationQR();

   /**
       * The QR factorisation of the matrix
       *
       */
   pair<Matrix,Matrix> factorisationLU();


   /**
    * get the submatrix
    * @param
    */
   Matrix submatrix(int i, int j, int n,int m);

   /**
    * get the norm of the matrix
    */
   double norm();

   /**
    * the transpose of the matrix
    *@return the transpose of the matrix
    */
   Matrix transpose();

   /**
    * the comatrix of the matrix
    * @return the comatrix of the matrix
    */
   Matrix comatrix();


   Matrix adjugate();
   /**
    *
    * @return the determinant of the matrix
    */
   double determinant();

   /**
    * allows to know if the matrix is square or not
    *@return true if it's a square matrix
    */
   bool isSquare() const;

   bool containsNan();

   bool isNull();
   /**
    * create a copy of the matrix
    * @return a copy of the matrix
    */
   Matrix clone();

   /**
    * get the trave value of the matrix
    *@return the sum of the diagonal value if the matrix is square, -1 otherwise
    */
   double trace() const;


   void dispose();

private:

   //attributes

   /**
    * represents the number of rows
    */
   int _rows;

   /**
    * represents the number of column
    */
   int _columns;

   /**
    * the values of the matrix
    */
   double** _values;

};

Vector solvingSystem(const Matrix& L,const Matrix& U,const Vector& b);
Vector solvingSystem(Matrix A,const Vector& b);


///////////////////////////////////////////////////////
/////       Operator Matrix Vector                /////
///////////////////////////////////////////////////////


Vector operator*( const Vector& v, const Matrix& m);

Vector operator*(const Matrix& m,const Vector& v);


} /* namespace isf */

#endif /* MATRIX_H_ */
