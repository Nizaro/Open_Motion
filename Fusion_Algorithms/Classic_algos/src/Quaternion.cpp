/*
 * Quaternion.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#include "Quaternion.h"

namespace om {


///////////////////////////////////////////////////////
/////        Constructors-Destructors             /////
///////////////////////////////////////////////////////


Quaternion::Quaternion() {
	_qw = -1.0f;
	_qx = -1.0f;
	_qy = -1.0f;
	_qz = -1.0f;
}


Quaternion::Quaternion(double qw,double qx,double qy,double qz){
	_qw = qw;
	_qx = qx;
	_qy = qy;
	_qz = qz;
}


Quaternion::Quaternion(double phi,double theta,double psy){


	//angle in radian
	phi /= 2.0;
	theta /= 2.0;
	psy /= 2.0;


	_qw = (cos(phi)*cos(theta)*cos(psy)) + (sin(phi)*sin(theta)*sin(psy));
	_qx = (sin(phi)*cos(theta)*cos(psy)) - (cos(phi)*sin(theta)*sin(psy));
	_qy = (cos(phi)*sin(theta)*cos(psy)) + (sin(phi)*cos(theta)*sin(psy));
	_qz= (cos(phi)*cos(theta)*sin(psy)) - (sin(phi)*sin(theta)*cos(psy));


}



Quaternion::Quaternion(double angle,Vector axis){

	// angle in degree
	double rad_angle = (angle/2.0) * DEG_TO_RAD;

	_qw = cos(rad_angle);
	_qx = axis.getValue(0)*sin( rad_angle ) ;
	_qy = axis.getValue(1)*sin( rad_angle );
	_qz = axis.getValue(2)*sin( rad_angle );


}

Quaternion::Quaternion(const Vector& v){

	if(v.getLength() == 3){

		_qw = 0.0f;
		_qx = v.getValue(0);
		_qy = v.getValue(1);
		_qz = v.getValue(2);

	}else if(v.getLength() == 4){

		_qw = v.getValue(3);
		_qx = v.getValue(0);
		_qy = v.getValue(1);
		_qz = v.getValue(2);

	}else{

		_qw = -1.0f;
		_qx = -1.0f;
		_qy = -1.0f;
		_qz = -1.0f;
	}

}



Quaternion::~Quaternion() {

}



///////////////////////////////////////////////////////
/////              Methods                        /////
///////////////////////////////////////////////////////


double Quaternion::norm(){

	return sqrtf( (_qw*_qw) + (_qx*_qx) +(_qy*_qy) + (_qz*_qz) );

}


void Quaternion::display(){

	cout << "q[" << _qw << ","<< _qx << ","<< _qy << ","<< _qz << "]";

}

void Quaternion::normalize(){

	double n = norm();
	_qw /= n;
	_qx /= n;
	_qy /= n;
	_qz /= n;

}


Quaternion Quaternion::conjugate(){

	return Quaternion(_qw,-_qx,-_qy,-_qz);
}

Quaternion Quaternion::inverse(){

	return conjugate()*(1.0f/norm());
}


///////////////////////////////////////////////////////
/////              Operators                      /////
///////////////////////////////////////////////////////

Quaternion& Quaternion::operator=(const Quaternion& p){

	   if (this != &p) {

		   _qw = p._qw;
		   _qx = p._qx;
		   _qy = p._qy;
		   _qz = p._qz;
	    }

	    return *this;



}

std::ostream& operator<<(std::ostream& os,Quaternion q)
{
	os.precision(10);
	os << "Quaternion[" << q.getQw() << "," << q.getQx() << "," << q.getQy() << "," << q.getQz() << "]";
  return os;
}



Quaternion Quaternion::operator+(const Quaternion& p){

	return Quaternion ( _qw+p._qw , _qx+p._qx , _qy+p._qy , _qz+p._qz);

}

Quaternion Quaternion::operator-(const Quaternion& p){

	return Quaternion ( _qw-p._qw , _qx-p._qx , _qy-p._qy , _qz-p._qz);

}

Quaternion Quaternion::operator*(const double& alpha){

	return Quaternion ( _qw*alpha , _qx*alpha , _qy*alpha , _qz*alpha);

}

Quaternion Quaternion::operator/(const double& alpha){

	return Quaternion ( _qw/alpha , _qx/alpha , _qy/alpha , _qz/alpha);

}



Quaternion Quaternion::operator*(const Quaternion& p){

	Quaternion product = Quaternion();

	product._qw = ( _qw*p._qw - _qx*p._qx -_qy*p._qy -_qz*p._qz);
	product._qx = ( _qw*p._qx + _qx*p._qw +_qy*p._qz -_qz*p._qy);
	product._qy = ( _qw*p._qy - _qx*p._qz +_qy*p._qw +_qz*p._qx);
	product._qz = ( _qw*p._qz + _qx*p._qy -_qy*p._qx +_qz*p._qw);

	/*
	product._qw = abs(product._qw) < EPSILON ? 0.0 : product._qw;
	product._qx = abs(product._qx) < EPSILON ? 0.0 : product._qx;
	product._qy = abs(product._qy) < EPSILON ? 0.0 : product._qy;
	product._qz = abs(product._qz) < EPSILON ? 0.0 : product._qz;
	*/


	return product;

}

} /* namespace isf */
