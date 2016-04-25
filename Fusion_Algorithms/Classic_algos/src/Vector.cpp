/*
 * Vector.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#include "Vector.h"

namespace om {


///////////////////////////////////////////////////////
/////        Constructors-Destructors             /////
///////////////////////////////////////////////////////

/* Default constructor */
Vector::Vector() {

   _length = -1;
   _values = 0;
}

/* initialize a null vector with a size equal to length */


Vector::Vector ( double length ){

	_length = static_cast<int>(length);

	_values = (double*)malloc(_length*sizeof(double));

   for (int i=0; i<_length; i++)
         _values[i] = 0.0;

}

Vector::Vector(const Vector& v){

	_length = v._length;
	_values = (double*)malloc(_length*sizeof(double));

    for (int i=0; i<_length; i++)
    	_values[i] = v._values[i];


}

Vector::Vector(int length, ...){

	_length = length;

	_values = (double*)malloc(_length*sizeof(double));

	va_list arguments;

	va_start ( arguments, length );

	for ( int i = 0; i < _length; i++ )
		_values[i] =  va_arg ( arguments, double );

	va_end(arguments);
}




Vector::Vector(std::vector<Vector> vec){

	_length = 0;

	 for( typename std::vector<Vector>::iterator iter = vec.begin();  iter != vec.end(); ++iter )
	       _length += (*iter)._length;

	 _values = (double*)malloc(_length*sizeof(double));

	 int index = 0;
	 for( typename std::vector<Vector>::iterator iter = vec.begin();  iter != vec.end(); ++iter ){

		 for(int i=0; i<(*iter)._length; ++i){
		  	   _values[index] = (*iter)._values[i];

		  	   index++;
		 }

	 }

}





/* free memory space */
Vector::~Vector() {

	if(_values != 0){
		free(_values);
		_values = 0;
	}

}


///////////////////////////////////////////////////////
/////              Methods                        /////
///////////////////////////////////////////////////////

/* */
void Vector::setValues(int n,...){

	va_list arguments;

	va_start ( arguments, n );

	if(_values != 0)
		dispose();

	_values = (double*)malloc(_length*sizeof(double));


	for ( int i = 0; i < _length; i++ )
		_values[i] =  va_arg ( arguments, double );

	va_end(arguments);

}

void Vector::setValues(const Vector& v){

	_length = v._length;

	if(_values != 0)
		dispose();

	_values = (double*)malloc(_length*sizeof(double));


    for (int i=0; i<_length; i++)
    	_values[i] = v._values[i];

}


/* return a clone*/
Vector Vector::clone(){

	Vector clone(static_cast<double>(_length));

	for(int i = 0; i<_length;++i)
		clone._values[i] = _values[i];

	return clone;

}

/* Norm's calculation of the vector */
double Vector::norm() const{

	double norm = 0.0;

	for(int i = 0; i<_length;++i)
		norm += _values[i]*_values[i];

	return sqrt(norm);
}

/* normalise the vector */
void Vector::normalize(){

	/* get vector's norm*/
	double n = norm();

	/* divide all values by the norm*/
	for(int i = 0; i<_length;++i)
		_values[i] /= n;

}


double Vector::rms(){

	double rms = 0.0;

	for(int i = 0; i<_length;++i)
		rms += _values[i]*_values[i];

	rms /= static_cast<double>(_length);

	return sqrt(rms);
}


/* Display vector's values in console */
void Vector::display(){

	cout.precision(5);
	cout << "[" ;
	for(int i = 0; i<_length;++i){
		if (i != _length-1)
			cout << _values[i] << ",";
		else
			cout << _values[i];
	}


	cout<< "]";


}



void Vector::dispose(){

	if(_values != 0){
		free(_values);
		_values = 0;
	}

}


///////////////////////////////////////////////////////
/////              Operators                      /////
///////////////////////////////////////////////////////


std::ostream& operator<<(std::ostream& os,const Vector& v)
{
	os.precision(10);
	os << "[" ;
	for(int i = 0; i<v.getLength();++i){
		if (i != v.getLength()-1)
			os << v.getValue(i) << ",";
		else
			os << v.getValue(i);
	}

	os<< "]";

  return os;
}


/* addition between two vector */
Vector Vector::operator +(const Vector& v){

	Vector res(static_cast<double>(_length));

	if(_length == v._length){
		for (int i = 0; i < _length; ++i){
			res._values[i] = _values[i] + v._values[i];
		}
	}else{
		 cerr << "Error !! Vector Impossible addition" << endl;
	}



	return res;
}

/* substraction between two vector */
Vector Vector::operator -(const Vector& v){

	Vector res(static_cast<double>(_length));

	if(_length == v._length){

		for (int i = 0; i < _length; ++i){
			res._values[i] = _values[i] - v._values[i];
		}
	}else{
		 cerr << "Error !! Impossible soustraction" <<endl;
	}

	return res;
}

/* Assignment operator */
Vector& Vector::operator=( const Vector& v){

    if (this != &v) {

    	// free resource;
    	if(_length != -1){
    		free(_values);
    		_values = 0;
    	}

        _length = v._length;
        _values = new double[_length];

        for (int i=0; i<_length; i++)
        	_values[i] = v._values[i];

    }
    return *this;

}


Vector Vector::operator/(const double& alpha){

	Vector res(static_cast<double>(_length));

	for (int i = 0; i < _length; ++i)
		res._values[i] = _values[i] / alpha;

	return res;


}


Vector Vector::operator*(const double& alpha){

	Vector res= Vector(static_cast<double>(_length));

	for (int i = 0; i < _length; ++i)
		res._values[i] = _values[i] * alpha;

	return res;

}

} /* namespace sdf */
