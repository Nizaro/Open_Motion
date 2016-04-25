/*
 * Quaternion.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#ifndef QUATERNION_H_
#define QUATERNION_H_

#include "Vector.h"

namespace om {


class Quaternion {

public:


	Quaternion();
	Quaternion(double qw,double qx,double qy,double qz);
	Quaternion(double phi,double theta,double psy);
	Quaternion(double angle,Vector axis);
	Quaternion(const Vector& v);

	virtual ~Quaternion();

	//getters

	double getQx() const {return _qx;}
	double getQw() const {return _qw;}
	double getQy() const {return _qy;}
	double getQz() const {return _qz;}

	//setters
	void setQx(double qx)  { _qx = qx;}
	void setQy(double qy)  { _qy = qy;}
	void setQz(double qz)  { _qz = qz;}
	void setQw(double qw)  { _qw = qw;}

	//operator
	Quaternion operator+(const Quaternion& p);
	Quaternion operator-(const Quaternion& p);
	Quaternion operator*(const Quaternion& p);
	Quaternion operator*(const double& alpha);
	Quaternion operator/(const double& alpha);
	Quaternion& operator=(const Quaternion& p);


	double norm();

	Quaternion conjugate();

	Quaternion inverse();

	void normalize();

	void display();

	Vector toVector(){return Vector(4,_qx,_qy,_qz,_qw);}

	Vector imaginary(){return Vector(3,_qx,_qy,_qz);}

private:
	double _qw;
	double _qx;
	double _qy;
	double _qz;
};


std::ostream& operator<< (std::ostream& os, Quaternion q);


} /* namespace isf */

#endif /* QUATERNION_H_ */
