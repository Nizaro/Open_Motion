/*
 * Vector.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#ifndef VECTOR_H_
#define VECTOR_H_


#include "annexes.h"


namespace om {

/**
 * Represents a vector of n values
 */
class Vector {

public:

	//Constructors
	Vector();
	Vector(double length);
	Vector(const Vector& v);
	Vector(int length, ...);

	Vector(std::vector<Vector> vec);

	//Destructors
	virtual ~Vector();

	//getters

	/**
	 * get the length of the vector
	 * @return the length of the vector
	 *
	 */
	int getLength() const {return _length;}

	/**
	 * get the i-th value of the vector
	 * @return the i-th value
	 */
	double getValue(int i) const {return _values[i];}

	//setters

	/**
	 * update the length of the vector (normally it will be unecessary
	 * @param length the new length of the vector
	 */
	void setLength(int length) {_length = length;}

	/**
	 * update the i-th value of the vector
	 * @param i position of the wanted value
	 * @param value the new value
	 */
	void setValue(int i,double value) {_values[i] = value;}


	/**
	 * update all values of the vector
	 * @param n length of the vector
	 */
	void setValues(int n,...);


	/**
	 * update all values of the vector
	 * @param v the new value
	 */
	void setValues(const Vector& v);

    //Operators

	/**
	 * allows multiplication between a vector and a constant
	 * @param v the second vector
	 * @return the sum of the two vector
	 */
	Vector operator*(const double& alpha);


	/**
	 * allows multiplication between a vector and a constant
	 * @param v the second vector
	 * @return the sum of the two vector
	 */
	Vector operator/(const double& alpha);

	/**
	 * allows addition between two vector of the same size
	 * @param v the second vector
	 * @return the sum of the two vector
	 */
	Vector operator+(const Vector& v);


	/**
	 * allows subtraction between two vector of the same size
	 * @param v the second vector
	 * @return the difference of the two vector
	 */
	Vector operator-(const Vector& v);

	/**
	 * allows assignment of a vector
	 * @param v the second vector
	 * @return the new vector
	 */
	Vector& operator=(const Vector& v);


	//methods

	/**
	 * normalise the vector
	 */
	void normalize();


	double rms();

	/**
	 * Norm's calculation of the vector
	 * @return the norm of the vector
	 */
	double norm() const;

	void display();
	Vector clone();

	void dispose();

private:

	//attributes

	/**
	 * the length of the vector
	 */
	int _length;

	/**
	 * the values of the vector
	 */
	double* _values;


};

std::ostream& operator<<(std::ostream& os,const Vector& v);



} /* namespace isf */

#endif /* VECTOR_H_ */
