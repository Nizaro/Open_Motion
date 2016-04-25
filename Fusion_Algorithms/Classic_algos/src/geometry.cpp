/*
 * geometry.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#include "geometry.h"

namespace om {



///////////////////////////////////////////////////////
/////               Matrix tools                  /////
///////////////////////////////////////////////////////



/* compute the jacobian matrix of a numerical function f define R^n to R^m */
Matrix jacobianMatrix(Vector x, NumericalFunction f ){

	Matrix jacobian = Matrix (f.m,x.getLength());
	double h = 0.01;

	for (int j = 0; j < jacobian.getColumns(); j++){

		Vector xhp = x.clone();
		Vector xhm = x.clone();

		xhp.setValue(j,xhp.getValue(j)+h);
		xhm.setValue(j,xhm.getValue(j)-h);

		Vector value =  (f.function(xhp) - f.function(xhm))/(2.0*h);

		for (int i = 0 ; i < jacobian.getRows(); i++){

			jacobian.setValue(i,j,value.getValue(i));

		}
	}


	return jacobian;



}

/* return the identity matrix with a size n*/
Matrix identity(int n){

	Matrix I = Matrix(n,n);

	for(int i=0;i<n;++i)
		I.setValue(i,i,1.0);

	return I;
}


/* get rotation matrix from Euler angles in convention ZYX */
Matrix rotationMatrixZYX(double phi,double theta,double psy){

	Matrix rotationPhi(3,3);
	Matrix rotationTheta(3,3);
	Matrix rotationPsy(3,3);

	rotationPhi.setValue(0,0,1.0);
	rotationPhi.setValue(0,1,0.0);
	rotationPhi.setValue(0,2,0.0);
	rotationPhi.setValue(1,0,0.0);
	rotationPhi.setValue(1,1, cos(phi));
	rotationPhi.setValue(1,2, sin(phi));
	rotationPhi.setValue(2,0,0.0);
	rotationPhi.setValue(2,1,-sin(phi));
	rotationPhi.setValue(2,2,cos(phi));

	rotationTheta.setValue(0,0,cos(theta));
	rotationTheta.setValue(0,1,0.0);
	rotationTheta.setValue(0,2,-sin(theta));
	rotationTheta.setValue(1,0,0.0);
	rotationTheta.setValue(1,1, 1.0);
	rotationTheta.setValue(1,2, 0.0);
	rotationTheta.setValue(2,0,sin(theta));
	rotationTheta.setValue(2,1,0.0);
	rotationTheta.setValue(2,2,cos(theta));

	rotationPsy.setValue(0,0,cos(psy));
	rotationPsy.setValue(0,1,sin(psy));
	rotationPsy.setValue(0,2,0.0);
	rotationPsy.setValue(1,0,-sin(psy));
	rotationPsy.setValue(1,1, cos(psy));
	rotationPsy.setValue(1,2, 0.0);
	rotationPsy.setValue(2,0,0.0);
	rotationPsy.setValue(2,1,0.0);
	rotationPsy.setValue(2,2,1.0);


	return rotationPsy*rotationTheta*rotationPhi;

}


Matrix rotationMatrixXYZ(double phi,double theta,double psy){

	return rotationMatrixZYX( phi, theta, psy).transpose();

}



Matrix skewSymetricMatrix(const Vector& x ){

	Matrix S(3,3);

	S.setValue(0,0,0.0);
	S.setValue(1,1,0.0);
	S.setValue(2,2,0.0);
	S.setValue(0,1,-x.getValue(2));
	S.setValue(0,2,x.getValue(1));
	S.setValue(1,0,x.getValue(2));
	S.setValue(1,2,-x.getValue(0));
	S.setValue(2,0,-x.getValue(1));
	S.setValue(2,1,x.getValue(0));

	return S;
}


Matrix operatorPsy(Quaternion q){

	Matrix Psy(4,3);

	Psy.setValue(0,0,q.getQw());
	Psy.setValue(0,1,q.getQz());
	Psy.setValue(0,2,-q.getQy());

	Psy.setValue(1,0,-q.getQz());
	Psy.setValue(1,1,q.getQw());
	Psy.setValue(1,2,q.getQx());

	Psy.setValue(2,0,q.getQy());
	Psy.setValue(2,1,-q.getQx());
	Psy.setValue(2,2,q.getQw());

	Psy.setValue(3,0,-q.getQx());
	Psy.setValue(3,1,-q.getQy());
	Psy.setValue(3,2,-q.getQz());

	return Psy;


}



Matrix operatorXi(Quaternion q){

	Matrix Xi(4,3);


	Xi.setValue(0,0,q.getQw());
	Xi.setValue(0,1,-q.getQz());
	Xi.setValue(0,2,q.getQy());

	Xi.setValue(1,0,q.getQz());
	Xi.setValue(1,1,q.getQw());
	Xi.setValue(1,2,-q.getQx());

	Xi.setValue(2,0,-q.getQy());
	Xi.setValue(2,1,q.getQx());
	Xi.setValue(2,2,q.getQw());

	Xi.setValue(3,0,-q.getQx());
	Xi.setValue(3,1,-q.getQy());
	Xi.setValue(3,2,-q.getQz());

	return Xi;

}

Matrix operatorGamma(const Vector& v){

	Matrix gamma(4,4);

	gamma.setValue(0,0,0.0);
	gamma.setValue(0,1,-v.getValue(2));
	gamma.setValue(0,2,v.getValue(1));
	gamma.setValue(0,3,v.getValue(0));

	gamma.setValue(1,0,v.getValue(2));
	gamma.setValue(1,1,0.0);
	gamma.setValue(1,2,-v.getValue(0));
	gamma.setValue(1,3,v.getValue(1));

	gamma.setValue(2,0,-v.getValue(1));
	gamma.setValue(2,1,v.getValue(0));
	gamma.setValue(2,2,0.0);
	gamma.setValue(2,3,v.getValue(2));

	gamma.setValue(3,0,-v.getValue(0));
	gamma.setValue(3,1,-v.getValue(1));
	gamma.setValue(3,2,-v.getValue(2));
	gamma.setValue(3,3,0.0);

	return gamma;

}

Matrix operatorOmega(const Vector& v){

	Matrix omega(4,4);

	omega.setValue(0,0,0.0);
	omega.setValue(0,1,v.getValue(2));
	omega.setValue(0,2,-v.getValue(1));
	omega.setValue(0,3,v.getValue(0));

	omega.setValue(1,0,-v.getValue(2));
	omega.setValue(1,1,0.0);
	omega.setValue(1,2,v.getValue(0));
	omega.setValue(1,3,v.getValue(1));

	omega.setValue(2,0,v.getValue(1));
	omega.setValue(2,1,-v.getValue(0));
	omega.setValue(2,2,0.0);
	omega.setValue(2,3,v.getValue(2));

	omega.setValue(3,0,-v.getValue(0));
	omega.setValue(3,1,-v.getValue(1));
	omega.setValue(3,2,-v.getValue(2));
	omega.setValue(3,3,0.0);

	/*/
	omega.setValue(0,0,0.0);
	omega.setValue(0,1,-v.getValue(0));
	omega.setValue(0,2,-v.getValue(1));
	omega.setValue(0,3,-v.getValue(2));

	omega.setValue(1,0,v.getValue(0));
	omega.setValue(1,1,0.0);
	omega.setValue(1,2,v.getValue(2));
	omega.setValue(1,3,-v.getValue(1));

	omega.setValue(2,0,v.getValue(1));
	omega.setValue(2,1,-v.getValue(2));
	omega.setValue(2,2,0.0);
	omega.setValue(2,3,v.getValue(0));

	omega.setValue(3,0,v.getValue(2));
	omega.setValue(3,1,v.getValue(1));
	omega.setValue(3,2,-v.getValue(0));
	omega.setValue(3,3,0.0);
	//*/

	return omega;
}




Matrix matSquareRoot1x1(Matrix A){

	Matrix square(1,1);

	square.setValue(0,0,sqrt(abs(A.getValue(0,0))));

	return square;

}
Matrix matSquareRoot2x2( Matrix A){

	Matrix square(2,2);
	Matrix I(identity(2));


	double detA = sqrt(abs( A.determinant()));
	double trA = sqrt(abs(A.trace() + 2.0*detA ));


	square = (A + (I*detA))/trA;


	return square;
}



Matrix matSquareRoot(Matrix A){

	Matrix squareA(A.getRows(),A.getColumns());

	if(!A.isNull()){

		if(A.getRows() == A.getColumns()){

			/*/
			squareA = A.choleskyDecomposition();
			/*/
			Matrix D = A.clone();
			int N=50;
			pair<Matrix,Matrix> QR_k = D.factorisationQR();

			Matrix P(QR_k.first);

			for(int i=0;i<N;i++){

				D = QR_k.second*QR_k.first;
				QR_k = D.factorisationQR();
				P = P*QR_k.first;

			}


			Matrix squareD(D.getRows(),D.getRows());

			for(int i=0;i<D.getRows();i++)
				squareD.setValue(i,i,sqrt(abs(D.getValue(i,i))));


			squareA = P*squareD*P.inverse();
			//*/
		}else
			cerr << "Impossible, non-square matrix" << endl;
	}

	return squareA;


}


void getEigenValuesVector(const Matrix& A,vector<double>& eigenValues,vector<Vector>& eigenVector,int N){

	Matrix D(A.getRows(),A.getColumns());

	for(int i=0;i<A.getRows();i++)
		for(int j=0;j<A.getColumns();j++)
			D.setValue(i,j,A.getValue(i,j));


	/* QR algorithm */
	pair<Matrix,Matrix> QR_k = D.factorisationQR();
	Matrix P(QR_k.first);

	for(int i=0;i<N;i++){

		D = QR_k.second*QR_k.first;
		QR_k = D.factorisationQR();
		P = P*QR_k.first;
	}

	/*/
	cout << "eigen value D" <<endl;
	D.display();
	cout << "eigen vector P"<< endl;
	P.display();
	//*/

	/* get the result */
	for(int l=0;l<A.getColumns();++l){
		eigenValues.push_back(D.getValue(l,l));
		eigenVector.push_back(P.getColumn(l).clone());
	}

}


Matrix matExponential(const Matrix& A){

	Matrix expA(A.getRows(),A.getColumns());

	if(A.getRows() == A.getColumns()){
		Matrix acc_A = A;
		double acc_n = 1.0;
		expA = identity(A.getRows()) + acc_A;

		for (int i = 2 ; i<50 ; ++i){
			acc_A = acc_A.clone() * A;
			acc_n = acc_n*static_cast<double>(i);

			expA = expA + (acc_A/acc_n);

		}

	}else
		cerr << "Impossible, non-square matrix";

	return expA;

}


Matrix transitionOmega(const Vector& gyro){

	Matrix Omega(4,4);

	double norm = gyro.norm();
	double cos_tmp = cos(0.5*norm*DELTA_T);
	double sin_tmp = sin(0.5*norm*DELTA_T);
	Vector phi(3,gyro.getValue(0)*( sin_tmp/norm ),gyro.getValue(1)*( sin_tmp/norm ),gyro.getValue(2)*( sin_tmp/norm ));

	Matrix Omega_tmp = (identity(3)*cos_tmp) - skewSymetricMatrix(phi);

	Omega.setRow(0,Vector(4,Omega_tmp.getValue(0,0),Omega_tmp.getValue(0,1),Omega_tmp.getValue(0,2),phi.getValue(0)));
	Omega.setRow(1,Vector(4,Omega_tmp.getValue(1,0),Omega_tmp.getValue(1,1),Omega_tmp.getValue(1,2),phi.getValue(1)));
	Omega.setRow(2,Vector(4,Omega_tmp.getValue(2,0),Omega_tmp.getValue(2,1),Omega_tmp.getValue(2,2),phi.getValue(2)));
	Omega.setRow(3,Vector(4,(-1.0)*phi.getValue(0),(-1.0)*phi.getValue(1),(-1.0)*phi.getValue(2),cos_tmp));

	return Omega;
}



///////////////////////////////////////////////////////
/////               Vector tools                  /////
///////////////////////////////////////////////////////


Vector rotationVector(const Vector& v,Quaternion q){

	Quaternion qv(0.0,v.getValue(0),v.getValue(1),v.getValue(2));
	Quaternion tmp = q.inverse()*qv*q;
	return tmp.imaginary();

}


Matrix vectorToMatrix(const Vector& v){

	Matrix res(v.getLength(),1);

	for(int i=0;i<v.getLength();i++){
		res.setValue(i,0,v.getValue(i));
	}

	return res;
}


Vector operatorVex(const Vector& a,const Vector& b){
	return operatorVex(skewSymetricMatrix(crossProduct(a,b)*(-1.0)));
}


Vector operatorVex(const Matrix& S){


	double x = S.getValue(2,1);
	double y = S.getValue(0,2);
	double z = S.getValue(1,0);

	Vector res (3,x,y,z);
	return res;

}


Vector projection( Vector v1, Vector v2 ){

	return v2 * dotProduct(v1,v2);
}




Vector crossProduct(const Vector& v1,const Vector& v2  ){

	Vector cross(3.0);

	if( v1.getLength() == 3 && v2.getLength() == 3){

		double cross_x = v1.getValue(1)*v2.getValue(2) - v1.getValue(2)*v2.getValue(1);
		double cross_y = v1.getValue(2)*v2.getValue(0) - v1.getValue(0)*v2.getValue(2);
		double cross_z =  v1.getValue(0)*v2.getValue(1) - v1.getValue(1)*v2.getValue(0);

		cross.setValues(3,cross_x,cross_y,cross_z);

	}


	return cross;
}


Matrix outerProduct(const Vector& v1,const Vector& v2 ){

	Matrix outer(3,3);

	outer.setValue(0,0, v1.getValue(0)* v2.getValue(0));
	outer.setValue(0,1, v1.getValue(0)* v2.getValue(1));
	outer.setValue(0,2, v1.getValue(0)* v2.getValue(2));
	outer.setValue(1,0, v1.getValue(1)* v2.getValue(0));
	outer.setValue(1,1, v1.getValue(1)* v2.getValue(1));
	outer.setValue(1,2, v1.getValue(1)* v2.getValue(2));
	outer.setValue(2,0, v1.getValue(2)* v2.getValue(0));
	outer.setValue(2,1, v1.getValue(2)* v2.getValue(1));
	outer.setValue(2,2, v1.getValue(2)* v2.getValue(2));

	return outer;



}


double dotProduct(const Vector& v1,const Vector& v2 ){

	double res = 0.0;

	if( v1.getLength() ==  v2.getLength() )
		for(int i = 0; i<v1.getLength();++i)
			res += v1.getValue(i)*v2.getValue(i);


	return res;

}

double euclidianDistance(const Vector& v1,const Vector& v2  ){

	double distance = -1.0;

	if( v1.getLength() == v2.getLength()){
		for(int i = 0; i<v1.getLength();++i)
			distance += (v1.getValue(i) - v2.getValue(i))*(v1.getValue(i) - v2.getValue(i));

		distance = sqrt(distance);
	}

	return distance;
}



///////////////////////////////////////////////////////
/////               Angles tools                  /////
///////////////////////////////////////////////////////



double angularDistance(const Vector& v1,const Vector& v2 ){

	double distance = -1.0;

	if( v1.getLength() == v2.getLength( )){

		distance = 0.0;

		for(int i=0;i<v1.getLength();++i){

			double diff = abs( v1.getValue(i) - v2.getValue(i) );
			distance += min(diff, 360.0f - diff)*min(diff, 360.0f - diff) ;
		}


		distance = sqrtf(distance);
	}



	return distance;

}



double angularSum(double alpha,double beta,double coeff_alpha,double coeff_beta){

	Vector u(2, cos(alpha),sin(alpha));
	Vector v(2, cos(beta),sin(beta));

	Vector sum = u*coeff_alpha + v*coeff_beta;
	sum.normalize();

	double tmp = sum.getValue(0);

	double angle = signeOf(asinf(sum.getValue(1))) * 2.0 * atan2f(sqrtf(1.0f - tmp*tmp),1.0f + tmp);

	if( tmp < -1.0 + EPSILON && tmp >-1.0 - EPSILON)
		angle = PI;

	return angle*RAD_TO_DEG;

}







Vector angularVelocity(Quaternion q_t,Quaternion q_tm1){

	/*/
	Quaternion r = q_t*q_tm1.inverse();

	double theta = 2.0*acos(r.getQw());
	Vector tmp = Vector(3,r.getQx(),r.getQy(),r.getQz());
	tmp.normalize();

	Vector angular_velocity = tmp*(theta/DELTA_T);
	/*/

	Quaternion dq_t = (q_t - q_tm1)/DELTA_T;
	Quaternion r = dq_t*q_t.conjugate()*2.0;

	Vector angular_velocity(3,r.getQx(),r.getQy(),r.getQz());
	//*/


	return angular_velocity;

}


///////////////////////////////////////////////////////
/////             Conversion tools                /////
///////////////////////////////////////////////////////



AxisAngle quatToAxisAngle(Quaternion q){


	AxisAngle axis_angle;

	if (  abs(q.getQw()-1.0) < EPSILON  ){

		axis_angle.angle = 0.0;
		axis_angle.axis = Vector(3.0);

	}else{

		Vector e(3,q.getQx(),q.getQy(),q.getQz());

		axis_angle.angle = 2.0*atan2( e.norm() , q.getQw() )*RAD_TO_DEG;

		/*
		if(axis_angle.angle > 180.0)
			axis_angle.angle -= 360.0;
		*/


		axis_angle.axis = e;
		axis_angle.axis.normalize();


	}

	return axis_angle;

}

Matrix quatToRotationMatrix (Quaternion q){

	Matrix rotation(3,3);
	Matrix I = identity(3);

	Vector x(3,q.getQx(),q.getQy(),q.getQz());
	Matrix S(skewSymetricMatrix(x));

	rotation = I + S*(2.0*q.getQw()) + S*S*2.0;

	return rotation;

}


Vector quatToEuler(Quaternion rotation){

	double sqx = rotation.getQx()*rotation.getQx();
	double sqy = rotation.getQy()*rotation.getQy();
	double sqz = rotation.getQz()*rotation.getQz();


	double test = 2.0* ( rotation.getQx()*rotation.getQz() - rotation.getQw()*rotation.getQy()   );

	double phi,theta,psy;


	/* north pole singularity detected */
	if ( abs(test + 1.0) < EPSILON ) {

		psy =  atan2(  (rotation.getQy()*rotation.getQx() - rotation.getQz()*rotation.getQw()) , rotation.getQx()*rotation.getQz() + rotation.getQw()*rotation.getQy()  );
		theta = PI/2.0;
		phi = 0;

	/* south pole singularity detected */
	} else if (abs(test - 1.0) < EPSILON ) {

		psy =  - atan2(  (rotation.getQy()*rotation.getQx() - rotation.getQz()*rotation.getQw()) , rotation.getQx()*rotation.getQz() + rotation.getQw()*rotation.getQy()  );
		theta = -PI/2.0;
		phi = 0;

	}else{
		phi =  atan2( (rotation.getQx()*rotation.getQw() + rotation.getQy()*rotation.getQz()) , ( 0.5 - (sqx+sqy )));
		theta = asin(- 2.0* ( rotation.getQx()*rotation.getQz() - rotation.getQw()*rotation.getQy()   ) );
		psy =  atan2( (rotation.getQz()*rotation.getQw() + rotation.getQy()*rotation.getQx()) , ( 0.5 - (sqz+sqy )));

	}

	phi = abs(phi) < EPSILON ? 0.0 : phi*RAD_TO_DEG;
	theta = abs(theta) < EPSILON ? 0.0 : theta*RAD_TO_DEG;
	psy = abs(psy) < EPSILON ? 0.0 : psy*RAD_TO_DEG;


	Vector res(3,phi,theta,psy);

	return res;

}


Matrix axisAngleToRotationMatrix (double angle,const Vector& axis){

	Matrix I = identity(3);
	Matrix S = skewSymetricMatrix(axis);

	Matrix R =  I + S*sin(angle) + S*S*(1.0 - cos(angle) );

	return R;

}



Matrix axisAngleToRotationMatrix (AxisAngle axis_angle){

	Matrix I = identity(3);
	Matrix S = skewSymetricMatrix(axis_angle.axis);

	Matrix R =  I + S*sin(axis_angle.angle) + S*S*(1.0 - cos(axis_angle.angle) );



	return R;

}


Vector rotationMatrixToEuler (const Matrix& m){

	double phi,theta,psy;

	if ( abs( m.getValue(0,2) - 1.0 ) < EPSILON ){

		theta = -PI/2.0;
		phi = atan2(- m.getValue(1,0),- m.getValue(3,1) );
		psy = 0.0;

	}else if ( abs( m.getValue(0,2) + 1.0 ) < EPSILON ){


		theta = PI/2.0;
		phi = atan2( m.getValue(1,0),m.getValue(2,0) );
		psy = 0.0;

	}else{

		phi = atan2(m.getValue(1,2),m.getValue(2,2));
		theta = asin(-m.getValue(0,2));
		psy = atan2(m.getValue(0,1),m.getValue(0,0));


	}

	phi = abs(phi) < EPSILON ? 0.0 : phi*RAD_TO_DEG;
	theta = abs(theta) < EPSILON ? 0.0 : theta*RAD_TO_DEG;
	psy = abs(psy) < EPSILON ? 0.0 : psy*RAD_TO_DEG;

	Vector res(3,phi,theta,psy);

	return res;

}





AxisAngle eulerToAxisAngle( double phi, double theta, double psy){

	double x,y,z,w,angle;

	// Assuming the angles are in radians.
	double c1 = cos(phi/2.0);
	double s1 = sin(phi/2.0);
	double c2 = cos(theta/2.0);
	double s2 = sin(theta/2.0);
	double c3 = cos(psy/2.0);
	double s3 = sin(psy/2.0);

	w =c1*c2*c3 - s1*s2*s3;
	x =c1*c2*s3 + s1*s2*c3;
	y =s1*c2*c3 + c1*s2*s3;
	z =c1*s2*c3 - s1*c2*s3;

	angle = 2.0f * acos(w);

	double norm = sqrt( x*x + y*y + z*z);

	if (norm < EPSILON) {

		x=1.0;
		y=z=0.0;

	} else {

		x = abs(x/norm) < EPSILON ? 0.0 : x/norm;
		y = abs(y/norm) < EPSILON ? 0.0 : y/norm;
		z = abs(z/norm) < EPSILON ? 0.0 : z/norm;
	}

	AxisAngle axis_angle;
	axis_angle.angle = angle;
	axis_angle.axis = Vector(3,x,y,z);

	return axis_angle;

}


Quaternion rotationMatrixToQuat(Matrix m){


	double tr = m.trace();
	double qw,qx,qy,qz;

	if (tr > 0) {

	  double S = sqrt(tr + 1.0) * 2.0;
	  qw = 0.25 * S;
	  qx = (m.getValue(2,1) - m.getValue(1,2)) / S;
	  qy = (m.getValue(0,2) - m.getValue(2,0)) / S;
	  qz = (m.getValue(1,0) - m.getValue(0,1)) / S;

	} else if ( (m.getValue(0,0) > m.getValue(1,1)) && (m.getValue(0,0) > m.getValue(2,2)) ) {

	  float S = sqrt( 1.0 + m.getValue(0,0) - m.getValue(1,1) - m.getValue(2,2) ) * 2.0; // S=4*qx
	  qw = (m.getValue(2,1) - m.getValue(1,2)) / S;
	  qx = 0.25 * S;
	  qy = (m.getValue(0,1) + m.getValue(1,0)) / S;
	  qz = (m.getValue(0,2) + m.getValue(2,0)) / S;

	} else if ( (m.getValue(1,1) > m.getValue(2,2)) ) {

	float S = sqrt(1.0 + m.getValue(1,1) - m.getValue(0,0) - m.getValue(2,2)) * 2.0; // S=4*qy
	  qw = (m.getValue(0,2) - m.getValue(2,0)) / S;
	  qx = (m.getValue(0,1) + m.getValue(1,0)) / S;
	  qy = 0.25 * S;
	  qz = (m.getValue(1,2)  + m.getValue(2,1)) / S;

	} else {

	  float S = sqrt(1.0 + m.getValue(2,2) - m.getValue(0,0) - m.getValue(1,1)) * 2.0; // S=4*qz
	  qw = (m.getValue(1,0) - m.getValue(0,1)) / S;
	  qx = (m.getValue(0,2) + m.getValue(2,0)) / S;
	  qy = (m.getValue(1,2)  + m.getValue(2,1)) / S;
	  qz = 0.25 * S;

	}

	Quaternion res(qw,qx,qy,qz);

	return res;

}







///////////////////////////////////////////////////////
/////               Quaternion tools              /////
///////////////////////////////////////////////////////

Quaternion quatExponential(Quaternion q){

	Vector e(3,q.getQx(),q.getQy(),q.getQz());
	double expo = exp(q.getQw());

	Vector tmp = (e/e.norm())*expo*sin(e.norm());

	Quaternion res( cos(e.norm())*expo, tmp.getValue(0), tmp.getValue(1), tmp.getValue(2)  );

	cout << "\n EXP res=" << res << endl;

	return res;


}

Quaternion quatLogarithm(Quaternion q){

	Vector e(3,q.getQx(),q.getQy(),q.getQz());
	double ln = log(q.norm());

	Vector tmp = (e/e.norm())*acos(q.getQw()/q.norm());

	Quaternion res(ln, tmp.getValue(0), tmp.getValue(1), tmp.getValue(2)  );

	cout << "\n LN res=" << res << endl;

	return res;
}

Quaternion quatPower(Quaternion q,double p){

	cout << "\n POWER q=" << q << endl;

	return quatExponential(quatLogarithm(q)*p);

}


double dotProduct(Quaternion q1,Quaternion q2){

	return (q1.getQw()*q2.getQw()) + (q1.getQx()*q2.getQx()) + (q1.getQy()*q2.getQy()) + (q1.getQz()*q2.getQz()) ;

}


Quaternion slerp(Quaternion q1,Quaternion q2,double t){

	return q1*quatPower(q1.inverse()*q2,t);
}

Quaternion nlerp(Quaternion q1,Quaternion q2,double t){

	q1.normalize();
	q2.normalize();

	return slerp(q1,q2,t);
}

///////////////////////////////////////////////////////
/////               Others tools                  /////
///////////////////////////////////////////////////////


double calculErrorOrientation(Quaternion q_real,Quaternion q_est){

	Quaternion dq = q_real*q_est.inverse();

	Vector dp = dq.imaginary() *
			(4.0 * (signeOf(dq.getQw())/(1.0+abs(dq.getQw()))));

	return dp.rms();


}





} /* namespace isf */
