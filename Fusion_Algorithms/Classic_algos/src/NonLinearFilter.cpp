/*
 * NonLinearFilter.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#include "NonLinearFilter.h"

namespace om {


Vector ned_gravity(3,0.0,0.0,G);
Vector ned_geographic_north(3,1.0,0.0,0.0);


///////////////////////////////////////////////////////
/////              Abstract Class                 /////
///////////////////////////////////////////////////////


NonLinearFilter::NonLinearFilter() {

	_var_gyro = _var_mag = _var_acc = _seed = -1.0;
}



void NonLinearFilter::setSensorsInfo(const Vector& bias_gyro,const Vector& bias_acc,const Vector& bias_mag,double var_gyro,double var_acc,double var_mag){

	_bias_gyro = bias_gyro;
	_bias_acc = bias_acc;
	_bias_mag = bias_mag;
	_var_gyro = var_gyro;
	_var_acc = var_acc;
	_var_mag = var_mag;

}

void NonLinearFilter::setSensorsData(const Vector& gyro,const Vector& acc,const Vector& mag){

	_gyroscope.dispose();
	_accelerometer.dispose();
	_magnetometer.dispose();

	_gyroscope = gyro;
	_accelerometer = acc;
	_magnetometer = mag;

}

void NonLinearFilter::setGyroscopeData(const Vector& data){
	_gyroscope.dispose();
	_gyroscope = data;
}

void NonLinearFilter::setAccelerometerData(const Vector& data){

	_accelerometer.dispose();
	_accelerometer = data;
}

void NonLinearFilter::setMagnetometerData(const Vector& data){

	_magnetometer.dispose();
	_magnetometer = data;

}


Vector NonLinearFilter::getGyroscopeData(){
	return _gyroscope.clone();
}

Vector NonLinearFilter::getAccelerometerData(){
	return _accelerometer.clone();
}

Vector NonLinearFilter::getMagnetometerData(){
	return _magnetometer.clone();
}



void NonLinearFilter::process(){

	preprocess();
	prediction();
	correction();

	_seed += 10.0;

}

///////////////////////////////////////////////////////
/////              MEKF 	                      /////
///////////////////////////////////////////////////////


MEKF::MEKF(){

	_a = 1.0;
	_f = 4.0;
	_h = 0.01;
	_seed = 0.0;
}

MEKF::~MEKF(){}

Vector MEKF::f_function(const Vector& x){


	Vector f_x(6.0);

	Vector b(3,x.getValue(3),x.getValue(4),x.getValue(5));

	Vector w_dp(3,_w_k.getValue(0),_w_k.getValue(1),_w_k.getValue(2));
	Vector w_b(3,_w_k.getValue(3),_w_k.getValue(4),_w_k.getValue(5));

	/*/
	Vector omega = _gyroscope - b + w_dp;

	Quaternion q( transitionOmega(omega) * _q_k_est.toVector() );
	Quaternion dq = q*_q_k_pred.inverse();

	Vector dqv (3,dq.getQx(),dq.getQy(),dq.getQz());
	Vector dp = dqv*(_f/(_a+dq.getQw()));

	for(int i=0;i<3;i++){
		f_x.setValue(i,dp.getValue(i));
		f_x.setValue(i+3,b.getValue(i)+w_b.getValue(i));
	}

	return f_x;

	/*/
	Vector dp(3,x.getValue(0),x.getValue(1),x.getValue(2));
	Vector omega = _gyroscope - b + w_dp;

	Vector delta_omega = _gyroscope - omega;


	Vector dp_var = crossProduct((omega),dp)*(-1.0) + delta_omega + crossProduct(delta_omega,dp)*(-0.5);

	for(int i=0;i<3;i++){
		f_x.setValue(i,dp_var.getValue(i));
		f_x.setValue(i+3,b.getValue(i));
	}
	return f_x + _w_k;
	//*/





}

Vector MEKF::h_function(const Vector& x){

	//*/
	Vector dp(3,x.getValue(0),x.getValue(1),x.getValue(2));

	double dqw = ( (_a*(-1))*(dp.norm()*dp.norm()) + (_f*sqrt(_f*_f + (1.0 - _a*_a)* (dp.norm()*dp.norm()) )))/((_f*_f) + (dp.norm()*dp.norm()));
	Vector dqv = dp * ( (1.0/_f) * (_a + dqw ));
	Quaternion dq(dqw,dqv.getValue(0),dqv.getValue(1),dqv.getValue(2));


	Vector tmp_res_m = quatToRotationMatrix(dq*_q_k_pred).transpose()*ned_geographic_north  + (_bias_mag );
	Vector tmp_res_a = quatToRotationMatrix(dq*_q_k_pred).transpose()*ned_gravity  + (_bias_acc );

	vector<Vector> vec;
	vec.push_back(tmp_res_a);
	vec.push_back(tmp_res_m);

	Vector h_xk(vec);

	return h_xk + _v_k;

	/*/
	Vector de(3,x.getValue(0),x.getValue(1),x.getValue(2));
	Quaternion dq(sqrt(1.0 - dotProduct(de,de)),de.getValue(0),de.getValue(1),de.getValue(2));

	Vector a(3,0.0,0.0,-G);
	Vector m(3,1.0,0.0,0.0);

	Vector tmp_res_m = quatToRotationMatrix(dq*_q_k_pred).transpose()*m  + (_bias_mag );
	Vector tmp_res_a = quatToRotationMatrix(dq*_q_k_pred).transpose()*a  + (_bias_acc );

	vector<Vector> vec;
	vec.push_back(tmp_res_a);
	vec.push_back(tmp_res_m);

	Vector h_xk(vec);

	return h_xk + _v_k;
	//*/
}

Matrix MEKF::W(const Vector& e,const Vector& v){

	//*/
	Matrix W_e_v(3,3);
	Matrix m_e = vectorToMatrix(e);
	Matrix m_v = vectorToMatrix(v);
	Matrix S(skewSymetricMatrix(v));

	double dot = dotProduct(e,e);
	dot = dot > 1.0 ? 2.0 - dot : dot;

	double alpha = sqrt(1.0 - dot);

	W_e_v = (S*(alpha*2.0)) + ( S*m_e*m_e.transpose()*(-2.0/alpha) ) + (identity(3)*(2.0*dotProduct(v,e))) + (m_e*m_v.transpose()*2.0) - (m_v*m_e.transpose()*4.0);

	return W_e_v;

	/*/
		Matrix W_e_v(3,3);
	Matrix m_e = vectorToMatrix(e);
	Matrix m_v = vectorToMatrix(v);
	Matrix S(skewSymetricMatrix(v));
	double alpha = sqrt(1.0 - dotProduct(e,e));

	W_e_v = (S*(alpha*2.0)) + ( S*m_e*m_e.transpose()*(-2.0/alpha) ) + (identity(3)*(2.0*dotProduct(v,e))) + (m_e*m_v.transpose()*2.0) - (m_v*m_e.transpose()*4.0);

	return W_e_v;
	//*/
}

void MEKF::computeJacobianF(){


	//*/

	_F_k = Matrix (6,6);
	double h = 0.01;

	for (int j = 0; j < _F_k.getColumns(); j++){

		Vector xhp = _x_k_est.clone();
		Vector xhm = _x_k_est.clone();

		xhp.setValue(j,xhp.getValue(j)+h);
		xhm.setValue(j,xhm.getValue(j)-h);

		Vector value =  (f_function(xhp) - f_function(xhm))/(2.0*h);

		for (int i = 0 ; i < _F_k.getRows(); i++){

			_F_k.setValue(i,j,value.getValue(i));

		}
	}

}

void MEKF::computeJacobianH(){


	_H_k = Matrix (6,6);
	double h = 0.01;

	for (int j = 0; j < _F_k.getColumns(); j++){

		Vector xhp = _x_k_pred.clone();
		Vector xhm = _x_k_pred.clone();

		xhp.setValue(j,xhp.getValue(j)+h);
		xhm.setValue(j,xhm.getValue(j)-h);

		Vector value =  (h_function(xhp) - h_function(xhm))/(2.0*h);

		for (int i = 0 ; i < _H_k.getRows(); i++){

			_H_k.setValue(i,j,value.getValue(i));

		}
	}

}

void MEKF::initialization(Quaternion q_0,int argc,...){


	if(argc != 0){

		double* values = new double[argc];

		va_list arguments;

		va_start ( arguments, argc );

		for ( int i = 0; i < argc; i++ )
			values[i] =  va_arg ( arguments, double );

		_a = values[0];
		_f = values[1];
		_h = values[2];

		va_end(arguments);

	}else{
		_a = 1.0;
		_f = 4.0;
		_h = 0.01;
		_seed = 0.0;
	}



	_q_k_est = Quaternion(q_0.getQw(),q_0.getQx(),q_0.getQy(),q_0.getQz());
	_x_k_est = Vector(6,0.0,0.0,0.0,_bias_gyro.getValue(0),_bias_gyro.getValue(1),_bias_gyro.getValue(2));


	/* compute Q_k (time invariant) */
	_Q_k = Matrix(6,6);

	double variance_u = 0.000031623;
	double variance_v = 0.0031623;

	// compute Q_k (time invariant)
	for(int i=0;i<3;++i){

		_Q_k.setValue(i,i,(variance_v*DELTA_T + 0.33333*variance_u*(pow(DELTA_T,3.0))));
		_Q_k.setValue(i+3,i,-(0.5*variance_u*DELTA_T*DELTA_T));
		_Q_k.setValue(i,i+3,-(0.5*variance_u*DELTA_T*DELTA_T));
		_Q_k.setValue(i+3,i+3,(variance_u*DELTA_T));

	}

	// compute R_k (time invariant)
	_R_k = Matrix  (6,6);
	for(int i=0;i<3;++i){
		_R_k.setValue(i,i,_var_acc);
		_R_k.setValue(i+3,i+3,_var_mag);
	}

	_Q_cho_k = _Q_k.choleskyDecomposition();
	for(int i=0;i<3;i++){

		double tmp= ((variance_v*variance_v)/DELTA_T) + ((variance_u*variance_u*DELTA_T)/12.0);

		_Q_cho_k.setValue(i,i,sqrt(tmp));
		_Q_cho_k.setValue(i+3,i+3,(variance_u*sqrt(DELTA_T)));
		_Q_cho_k.setValue(i+3,i,(-0.5*variance_u*sqrt(DELTA_T)));

	}


	_R_cho_k = _R_k.choleskyDecomposition();


	// init P_0 to null matrix
	_P_k = Matrix(6,6);


}

void MEKF::preprocess(){

	/* compute measurement noise vector at time k*/
	_w_k = whiteNoiseVector(0.0,_Q_cho_k,_seed);
	_seed += 6.0;

	/* compute measurement noise vector at time k*/
	_v_k = whiteNoiseVector(0.0,_R_cho_k,_seed);
	_seed += 6.0;

	/* compute measurement vector from sensor data */
	_z_k = Vector(6,_accelerometer.getValue(0),_accelerometer.getValue(1),_accelerometer.getValue(2),
			_magnetometer.getValue(0),_magnetometer.getValue(1),_magnetometer.getValue(2));

}

void MEKF::prediction(){

	// propagation of quaternion
	Vector b_pred(3,_x_k_est.getValue(3),_x_k_est.getValue(4),_x_k_est.getValue(5));
	_q_k_pred = Quaternion(  transitionOmega( _gyroscope - b_pred) * _q_k_est.toVector() );

	// compute the jacobian matrix of function f
	computeJacobianF();

	// propagation of state and covariance
	_x_k_pred = f_function(_x_k_est);
	_P_k_pred = _F_k*_P_k*_F_k.transpose() + _Q_k;



}

void MEKF::correction(){

	// compute the jacobian matrix of function h
	computeJacobianH();

	// estimate the measurement vector
	_z_k_pred = h_function(_x_k_pred);

	// compute Kalman gain
	Matrix K_k = _P_k_pred*_H_k.transpose()*( _H_k*_P_k_pred*_H_k.transpose()+_R_k ).inverse();

	// Innovation matrix
	Matrix S = (identity(6) - K_k*_H_k);

	// update state vector and covariance
	_x_k_est = _x_k_pred + K_k*(_z_k - _z_k_pred);

	//_P_k = S*_P_k_pred*S.transpose() + K_k*_R_k*K_k.transpose();
	_P_k = S*_P_k_pred;

	// update attitude quaternion
	Vector dp_est(3,_x_k_est.getValue(0),_x_k_est.getValue(1),_x_k_est.getValue(2));
	double dqw_est = ( (_a*(-1))*(dp_est.norm()*dp_est.norm()) + (_f*sqrt(_f*_f + (1.0 - _a*_a)* (dp_est.norm()*dp_est.norm()) )))/((_f*_f) + (dp_est.norm()*dp_est.norm()));
	Vector dqv_est = dp_est * ( (1.0/_f) * (_a + dqw_est ));
	Quaternion dq_est(dqw_est,dqv_est.getValue(0),dqv_est.getValue(1),dqv_est.getValue(2));

	// compute optimal estimation
	_q_k_est = 	dq_est*_q_k_pred;

	// reset error de to zero
	//_x_k_est = Vector(6,0.0,0.0,0.0,_x_k_est.getValue(3),_x_k_est.getValue(4),_x_k_est.getValue(5));

}


///////////////////////////////////////////////////////
/////              USQUE	                      /////
///////////////////////////////////////////////////////

USQUE::USQUE(){

	_f = -1.0;
	_lambda = -1.0;
	_a = -1.0;
	_variance_u = 0.000031623;
	_variance_v = 0.0031623;

	_n = 6.0;
}

USQUE::~USQUE(){}

void USQUE::initialization(Quaternion q_0,int argc,...){



	if(argc != 0){

		double* values = new double[argc];

		va_list arguments;

		va_start ( arguments, argc );

		for ( int i = 0; i < argc; i++ )
			values[i] =  va_arg ( arguments, double );

		_a = values[0];
		_f = values[1];
		_lambda = values[2];
		_variance_u = values[3];
		_variance_v = values[4];


		va_end(arguments);

	}else{

		_a = 1.0;
		_f = 4.0;
		_lambda = 1.0;
		_seed = 0.0;
		_variance_u = 0.000031623;
		_variance_v = 0.0031623;

	}


	_q_k_est = Quaternion(q_0.getQw(),q_0.getQx(),q_0.getQy(),q_0.getQz());
	_x_k_est = Vector(6,0.0,0.0,0.0,_bias_gyro.getValue(0),_bias_gyro.getValue(1),_bias_gyro.getValue(2));

	/* compute Q_k (time invariant) */
	_Q_k = Matrix(6,6);
	for(int i=0;i<3;++i){

		_Q_k.setValue(i,i,(_variance_v*DELTA_T + 0.33333*_variance_u*(pow(DELTA_T,3.0))));
		_Q_k.setValue(i+3,i,-(0.5*_variance_u*DELTA_T*DELTA_T));
		_Q_k.setValue(i,i+3,-(0.5*_variance_u*DELTA_T*DELTA_T));
		_Q_k.setValue(i+3,i+3,(_variance_u*DELTA_T));

	}

	/* compute R_k (time invariant) */
	_R_k = Matrix (6,6);
	for(int i=0;i<3;++i){
		_R_k.setValue(i,i,_var_acc);
		_R_k.setValue(i+3,i+3,_var_mag);
	}

	/* init P_0 to null matrix */
	_P_k = Matrix(6,6);


}

void USQUE::preprocess(){


	/* empty sigma points */
	if(!_sigma_points_x_k.empty())
		_sigma_points_x_k.clear();

	if(!_sigma_points_x_kp1.empty())
		_sigma_points_x_kp1.clear();

	if(!_sigma_points_z_kp1.empty())
		_sigma_points_z_kp1.clear();

	if(!_sigma_quaternion_k.empty())
		_sigma_quaternion_k.clear();


	/* compute measurement vector from sensor data */
	_z_k = Vector(6,_accelerometer.getValue(0),_accelerometer.getValue(1),_accelerometer.getValue(2),
			_magnetometer.getValue(0),_magnetometer.getValue(1),_magnetometer.getValue(2));

}

void USQUE::prediction(){

	/* compute square root */
	Matrix tmpS = (_Q_k + _P_k)*(_n+_lambda);
	Matrix S = tmpS.choleskyDecomposition();
	//Matrix S = matSquareRoot(tmpS.clone());

	/* generate sigma point at time k */
	for(int i=0;i<=(int)(_n);++i){

		if(i==0){
			_sigma_points_x_k.push_back(_x_k_est);
		}else{
			_sigma_points_x_k.push_back(_x_k_est + S.getRow(i-1));
			_sigma_points_x_k.push_back(_x_k_est - S.getRow(i-1));
		}
	}


	/* generate sigma quaternion at time k */
	for(unsigned int i=0;i<_sigma_points_x_k.size();++i){

		if(i==0){
			_sigma_quaternion_k.push_back(_q_k_est);
		}else{

			/* compute error quaternion i */
			Vector dp_i(3,_sigma_points_x_k[i].getValue(0),_sigma_points_x_k[i].getValue(1),_sigma_points_x_k[i].getValue(2));
			double dqw_i = ((_a*(-1))*(dp_i.norm()*dp_i.norm()) + (_f*sqrt(_f*_f + (1.0 - _a*_a)* (dp_i.norm()*dp_i.norm()) )))/((_f*_f) + (dp_i.norm()*dp_i.norm()));
			Vector dqv_i = dp_i * ( (1.0/_f) * (_a + dqw_i ));
			Quaternion dq_i (dqw_i,dqv_i.getValue(0),dqv_i.getValue(1),dqv_i.getValue(2));

			/* compute sigma quaternion i */
			Quaternion q_i = dq_i*_q_k_est;
			_sigma_quaternion_k.push_back(q_i);



		}
	}


	/* propagation of sigma quaternion and sigma points */
	for(unsigned int i=0;i<_sigma_points_x_k.size();++i){

		/* propagation of sigma quaternion */
		Vector b_i(3,_sigma_points_x_k[i].getValue(3),_sigma_points_x_k[i].getValue(4),_sigma_points_x_k[i].getValue(5));
		Quaternion q_pred_i(  transitionOmega(_gyroscope - b_i) * _sigma_quaternion_k[i].toVector() );

		if(i==0){

			_sigma_quaternion_k[0] = q_pred_i;
			_sigma_points_x_kp1.push_back(Vector(6,0.0,0.0,0.0,b_i.getValue(0),b_i.getValue(1),b_i.getValue(2)));

		}else{

			_sigma_quaternion_k[i] = q_pred_i;

			Quaternion dq_i = q_pred_i*_sigma_quaternion_k[0].inverse();

			Vector dqv_i (3,dq_i.getQx(),dq_i.getQy(),dq_i.getQz());
			Vector dp_i = dqv_i*(_f/(_a+dq_i.getQw()));

			_sigma_points_x_kp1.push_back(Vector(6,dp_i.getValue(0),dp_i.getValue(1),dp_i.getValue(2),b_i.getValue(0),b_i.getValue(1),b_i.getValue(2)));

		}

	}


	/* compute mean of x_k */
	Vector sum_x(6.0);
	for(unsigned int i=1;i<_sigma_points_x_kp1.size();++i){
		sum_x = sum_x + _sigma_points_x_kp1[i];
	}
	_x_k_pred = (_sigma_points_x_kp1[0]*(_lambda/(_n+_lambda))) + (sum_x*(1.0/(2.0*(_n+_lambda))));



	/* compute covarience P_k*/
	Matrix sum_P(6,6);
	for(unsigned int i=1;i<_sigma_points_x_kp1.size();++i){
		Matrix tmp_i = vectorToMatrix( _sigma_points_x_kp1[i] - _x_k_pred );
		sum_P = sum_P + tmp_i*tmp_i.transpose();
	}
	Matrix tmp_0 = vectorToMatrix( _sigma_points_x_kp1[0] - _x_k_pred );
	_P_k_pred = _Q_k + (tmp_0*tmp_0.transpose()*(_lambda/(_n+_lambda))) + (sum_P*(1.0/(2.0*(_n+_lambda))));


	/* compute mean of q_k */
	Vector b_pred(3,_x_k_est.getValue(3),_x_k_est.getValue(4),_x_k_est.getValue(5));
	_q_k_pred = Quaternion(  transitionOmega(_gyroscope - b_pred) * _q_k_est.toVector() );




}

void USQUE::correction(){





	for(unsigned int i=0;i<_sigma_points_x_kp1.size();++i){

		Vector tmp_res_m = quatToRotationMatrix(_sigma_quaternion_k[i]).transpose()*ned_geographic_north  + (_bias_mag );
		Vector tmp_res_a = quatToRotationMatrix(_sigma_quaternion_k[i]).transpose()*ned_gravity  + (_bias_acc );

		tmp_res_m.normalize();

		vector<Vector> vec;
		vec.push_back(tmp_res_a);
		vec.push_back(tmp_res_m);

		Vector z_k_i(vec);
		_sigma_points_z_kp1.push_back(z_k_i);
	}

	/* compute mean of z_k */
	Vector sum_z(6.0);
	for(unsigned int i=1;i<_sigma_points_z_kp1.size();++i){
		sum_z = sum_z + _sigma_points_z_kp1[i];
	}
	_z_k_pred = (_sigma_points_z_kp1[0]*(_lambda/(_n+_lambda))) + (sum_z*(1.0/(2.0*(_n+_lambda))));


	//cout << "\n_z_k = " << _z_k << endl;
	//cout << "_z_k_pred = " << _z_k_pred << endl;

	/* compute covariance Pxz and Pzz*/
	Matrix sum_Pxz(6,6);
	Matrix sum_Pzz(6,6);
	for(unsigned int i=1;i<_sigma_points_x_kp1.size();++i){

		Matrix tmp_z_i = vectorToMatrix( _sigma_points_z_kp1[i] - _z_k_pred );
		Matrix tmp_x_i = vectorToMatrix( _sigma_points_x_kp1[i] - _x_k_pred );

		sum_Pzz = sum_Pzz + tmp_z_i*tmp_z_i.transpose();
		sum_Pxz = sum_Pxz + tmp_x_i*tmp_z_i.transpose();

	}

	Matrix tmp_x_0 = vectorToMatrix( _sigma_points_x_kp1[0] - _x_k_pred );
	Matrix tmp_z_0 = vectorToMatrix( _sigma_points_z_kp1[0] - _z_k_pred );

	_cov_Z = _R_k + (tmp_z_0*tmp_z_0.transpose()*(_lambda/(_n+_lambda))) + (sum_Pzz*(1.0/(2.0*(_n+_lambda))));
	_cov_Z_x = (tmp_x_0*tmp_z_0.transpose()*(_lambda/(_n+_lambda))) + (sum_Pxz*(1.0/(2.0*(_n+_lambda))));


	/* Kalman gain */
	Matrix K_k = _cov_Z_x*_cov_Z.inverse();

	/* correction and update */
	_x_k_est = _x_k_pred + K_k*(_z_k - _z_k_pred);
	_P_k = _P_k_pred - (K_k*_cov_Z*K_k.transpose());


	/* update of quaternion */
	Vector dp_est(3,_x_k_est.getValue(0),_x_k_est.getValue(1),_x_k_est.getValue(2));
	double dqw_est = ( (_a*(-1))*(dp_est.norm()*dp_est.norm()) + (_f*sqrt(_f*_f + (1.0 - _a*_a)* (dp_est.norm()*dp_est.norm()) )))/((_f*_f) + (dp_est.norm()*dp_est.norm()));
	Vector dqv_est = dp_est * ( (1.0/_f) * (_a + dqw_est ));
	Quaternion dq_est (dqw_est,dqv_est.getValue(0),dqv_est.getValue(1),dqv_est.getValue(2));

	_q_k_est = 	dq_est*_q_k_pred;


}





///////////////////////////////////////////////////////
/////              CGO		                      /////
///////////////////////////////////////////////////////


CGO::CGO() {
	_k_mag = -1.0;
	_k_acc = -1.0;
	_k_I = -1.0;
	_k_P = -1.0;

}

CGO::~CGO() {
	// TODO Auto-generated destructor stub
}


void CGO::initialization(Quaternion q_0,int argc,...){


	if(argc != 0){

		double* values = new double[argc];

		va_list arguments;

		va_start ( arguments, argc );

		for ( int i = 0; i < argc; i++ )
			values[i] =  va_arg ( arguments, double );

		va_end(arguments);

		_k_I = values[0];
		_k_P = values[1];
		_k_mag = values[2];
		_k_acc = values[3];

	}else{

		_k_I = 0.3;
		_k_P = 1.0;
		_k_mag = 1.0;
		_k_acc = 1.0;

	}

	_q_k_est = q_0;
	_bias_est = _bias_gyro;




}

void CGO::preprocess(){

}

void CGO::prediction(){

	_bias_pred = _bias_est;

	//*/
	_q_k_pred = Quaternion(  transitionOmega(_gyroscope - _bias_pred) * _q_k_est.toVector() );
	/*/
	_q_k_pred = _q_k_est + Quaternion(operatorXi(_q_k_est)*(_gyroscope - _bias_pred)*(0.5*DELTA_T));
	//*/

}


void CGO::correction(){


	_v_acc_pred = quatToRotationMatrix(_q_k_pred).transpose()*ned_gravity  + (_bias_acc );
	_v_mag_pred = quatToRotationMatrix(_q_k_pred).transpose()*ned_geographic_north  + (_bias_mag );

	_accelerometer.normalize();
	_magnetometer.normalize();

	_v_acc_pred.normalize();
	_v_mag_pred.normalize();


	Matrix S_acc = skewSymetricMatrix(crossProduct(_accelerometer,_v_acc_pred)*(-1.0))*(_k_acc/2.0);
	Matrix S_mag = skewSymetricMatrix(crossProduct(_magnetometer,_v_mag_pred)*(-1.0))*(_k_mag/2.0);

	_omega = operatorVex(S_acc + S_mag)*(-1.0);

	_q_k_est = _q_k_pred + Quaternion(operatorXi(_q_k_pred)*(_omega*_k_P  )*(0.5*DELTA_T));

	_bias_est = _bias_pred - (_omega*(_k_I)*DELTA_T);


}




///////////////////////////////////////////////////////
/////                REQUEST                      /////
///////////////////////////////////////////////////////

REQUEST::REQUEST() {
	_var_gyro = _var_mag = _var_acc = -1.0;
	_lambda_m_k = -1.0;
	_mu_k = -1.0;
	_d_m_k = -1.0;
	_m_k = -1.0;
	_d_sigma_k = -1.0;
}

REQUEST::~REQUEST() {

}

void REQUEST::computeRk(){

	_R_k = Matrix(4,4);
	double n_k = 3.0;

	/* computation of R22 */
	double R22 = (2.0*_mu_k)/n_k;

	/* computation of R11 */
	Matrix R11(3,3);

	for(unsigned int i=0;i<_a.size();i++){
		Matrix m_bi = vectorToMatrix(_b[i]);
		Matrix m_ri = vectorToMatrix(_r[i]);
		Matrix S_ri = skewSymetricMatrix(_r[i]);
		Matrix I = identity(3)*(3.0 -(dotProduct(_r[i],_b[i])*dotProduct(_r[i],_b[i])));

		R11 = R11 + ( I + ( ( (m_bi*m_ri.transpose()) +(m_ri*m_bi.transpose())   )*dotProduct(_b[i],_r[i]) )  + (S_ri*(m_bi*m_bi.transpose())*S_ri.transpose())  );

	}

	R11 = R11*(_mu_k/n_k);

	_R_k.setRow(0,Vector(4,R11.getValue(0,0),R11.getValue(0,1),R11.getValue(0,2),0.0));
	_R_k.setRow(1,Vector(4,R11.getValue(1,0),R11.getValue(1,1),R11.getValue(1,2),0.0));
	_R_k.setRow(2,Vector(4,R11.getValue(2,0),R11.getValue(2,1),R11.getValue(2,2),0.0));
	_R_k.setRow(3,Vector(4,0.0,0.0,0.0,R22));



}

void REQUEST::computeQk(){


	_Q_k = Matrix(4,4);
	double eta_k = _var_gyro;

	/* computation of Q22 */
	double Q22 = eta_k*((_d_B_k*_d_B_k.transpose()).trace() + (_d_sigma_k*_d_sigma_k) + (_d_z_k.norm()*_d_z_k.norm()) );


	/* computation of Q22 */
	double tmp = ( (_d_sigma_k*_d_sigma_k) + (_d_z_k.norm()*_d_z_k.norm()) - (_d_B_k*_d_B_k.transpose()).trace());
	Matrix Q11 = ( (identity(3)*tmp) + ( (_d_B_k.transpose()*_d_B_k) - (_d_B_k*_d_B_k) - (_d_B_k.transpose()*_d_B_k.transpose()) )*2.0 );


	/* computation of Q12 et Q21 */
	Matrix M = _d_B_k*(_d_B_k - (identity(3)*_d_sigma_k));
	Vector y = operatorVex(M.transpose() - M);
	Vector Q12 = (y + (_d_B_k.transpose()*_d_z_k))*(eta_k*(-1));

	_Q_k.setRow(0,Vector(4,Q11.getValue(0,0),Q11.getValue(0,1),Q11.getValue(0,2),Q12.getValue(0)));
	_Q_k.setRow(1,Vector(4,Q11.getValue(1,0),Q11.getValue(1,1),Q11.getValue(1,2),Q12.getValue(1)));
	_Q_k.setRow(2,Vector(4,Q11.getValue(2,0),Q11.getValue(2,1),Q11.getValue(2,2),Q12.getValue(2)));
	_Q_k.setRow(3,Vector(4,Q12.getValue(0),Q12.getValue(1),Q12.getValue(2),Q22));

	_Q_k = _Q_k*(DELTA_T*DELTA_T);

}

void REQUEST::initialization(Quaternion q_0,int argc,...){


	if(argc != 0){

		double* values = new double[argc];

		va_list arguments;

		va_start ( arguments, argc );

		for ( int i = 0; i < argc; i++ )
			values[i] =  va_arg ( arguments, double );

		va_end(arguments);

		_mu_k = values[0];

	}else{

		_mu_k = 0.01;

	}

	/* init q_0 */
	_q_k_est = q_0;



	/* init weights */
	_a.push_back(0.34);
	_a.push_back(0.33);
	_a.push_back(0.33);

	/* init measuerment vector */
	_b.push_back(quatToRotationMatrix(_q_k_est).transpose()*Vector(3,1.0,0.0,0.0));
	_b.push_back(quatToRotationMatrix(_q_k_est).transpose()*Vector(3,0.0,1.0,0.0));
	_b.push_back(quatToRotationMatrix(_q_k_est).transpose()*Vector(3,0.0,0.0,1.0));

	/* init reference vector */
	_r.push_back(Vector(3,1.0,0.0,0.0));
	_r.push_back(Vector(3,0.0,1.0,0.0));
	_r.push_back(Vector(3,0.0,0.0,1.0));

	/* init m_0, sigma_0, B_0 and z_0 */
	_d_m_k = 0.0;
	_d_sigma_k = 0.0;
	_d_B_k = Matrix(3,3);
	_d_z_k = Vector(3.0);
	for(unsigned i=0;i<_a.size();++i){
		_d_m_k += _a[i];
		_d_sigma_k += _a[i]*dotProduct(_b[i],_r[i]);
		_d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
		_d_B_k = _d_B_k + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
	}

	_d_sigma_k /= _d_m_k;
	_d_B_k = _d_B_k/_d_m_k;
	_d_z_k = _d_z_k/_d_m_k;

	/* init S_0 */
	_d_S_k =_d_B_k + _d_B_k.transpose();

	/* init K_est_0 */
	_K_est_k = Matrix(4,4);
	Matrix tmp = _d_S_k - (identity(3)*_d_sigma_k);

	_K_est_k.setRow(0,Vector(4,tmp.getValue(0,0),tmp.getValue(0,1),tmp.getValue(0,2),_d_z_k.getValue(0)));
	_K_est_k.setRow(1,Vector(4,tmp.getValue(1,0),tmp.getValue(1,1),tmp.getValue(1,2),_d_z_k.getValue(1)));
	_K_est_k.setRow(2,Vector(4,tmp.getValue(2,0),tmp.getValue(2,1),tmp.getValue(2,2),_d_z_k.getValue(2)));
	_K_est_k.setRow(3,Vector(4,_d_z_k.getValue(0),_d_z_k.getValue(1),_d_z_k.getValue(2),_d_sigma_k));

	computeRk();

	_P_k = _R_k;
	_m_k = _d_m_k;
}

void REQUEST::preprocess(){

	/* update measurement vector */
	Vector acc = (_accelerometer - _bias_acc);
	Vector mag = (_magnetometer  - _bias_mag);
	acc.normalize();
	mag.normalize();

	Vector y_axis = crossProduct(acc,mag);
	y_axis.normalize();

	_b.clear();
	_b.push_back(mag.clone());
	_b.push_back(y_axis.clone());
	_b.push_back(acc.clone());


}

void REQUEST::prediction(){

	/* compute d_m_k, d_sigma_k, d_B_k, d_z_k*/
		_d_m_k = 0.0;
		_d_sigma_k = 0.0;
		_d_B_k = Matrix(3,3);
		_d_z_k = Vector(3.0);
		for(unsigned i=0;i<_a.size();++i){
			_d_m_k += _a[i];
			_d_sigma_k += _a[i]*dotProduct(_b[i],_r[i]);
			_d_z_k = _d_z_k + (crossProduct(_b[i],_r[i])*_a[i]);
			_d_B_k = _d_B_k + (vectorToMatrix(_b[i])*vectorToMatrix(_r[i]).transpose()*_a[i]);
		}

		_d_sigma_k /= _d_m_k;
		_d_B_k = _d_B_k/_d_m_k;
		_d_z_k = _d_z_k/_d_m_k;

		/* compute d_S_k */
		_d_S_k =_d_B_k + _d_B_k.transpose();

		/* compute d_K_k */
		_d_K_k = Matrix(4,4);
		Matrix tmp = _d_S_k - (identity(3)*_d_sigma_k);

		_d_K_k.setRow(0,Vector(4,tmp.getValue(0,0),tmp.getValue(0,1),tmp.getValue(0,2),_d_z_k.getValue(0)));
		_d_K_k.setRow(1,Vector(4,tmp.getValue(1,0),tmp.getValue(1,1),tmp.getValue(1,2),_d_z_k.getValue(1)));
		_d_K_k.setRow(2,Vector(4,tmp.getValue(2,0),tmp.getValue(2,1),tmp.getValue(2,2),_d_z_k.getValue(2)));
		_d_K_k.setRow(3,Vector(4,_d_z_k.getValue(0),_d_z_k.getValue(1),_d_z_k.getValue(2),_d_sigma_k));

		/* propagate Q_k */
		computeQk();

		/* propagate R_k */
		computeRk();


		/* compute state transition matrix Phy */
		Matrix Phy = transitionOmega(_gyroscope - _bias_gyro);

		/* propagate K_k */
		_K_pred_k = Phy*_K_est_k*Phy.transpose();

		/* propagate P_k */
		_P_pred_k = Phy*_P_k*Phy.transpose() + _Q_k;



}

void REQUEST::correction(){

	/* compute optimal gain */
	double rho_k_opt = ( (_m_k*_m_k)*_P_pred_k.trace())/( (_m_k*_m_k*_P_pred_k.trace()) +  (_d_m_k*_d_m_k*_R_k.trace()) );
	rho_k_opt=0.001;

	/* propagate m_k */
	double m_kp1 = ((1.0 - rho_k_opt)*_m_k) + (rho_k_opt*_d_m_k);

	/* some variable */
	double tmp_a =  (1.0 - rho_k_opt)*( _m_k/m_kp1 ) ;
	double tmp_b =  (rho_k_opt)*( _d_m_k/m_kp1 ) ;

	/* update K_k */
	_K_est_k = ( _K_pred_k * tmp_a)  +  ( _d_K_k * tmp_b);

	/* update P_k */
	_P_k = (_P_pred_k *(tmp_a*tmp_a)) +(_R_k*(tmp_b*tmp_b));

	/* update m_k */
	_m_k = m_kp1;




	/* calcul of lambda_max */
	vector<double> eigen_value;
	vector<Vector> eigen_vector;

	getEigenValuesVector(_K_est_k,eigen_value,eigen_vector,50);

	_lambda_m_k=0.0;
	int index=0;
	for(unsigned int i=0;i<eigen_vector.size();i++){

		if(_lambda_m_k < eigen_value[i]){
			_lambda_m_k = eigen_value[i];
			index=i;
		}
	}

	_q_k_est =  Quaternion(eigen_vector[index]);

	eigen_value.clear();
	eigen_vector.clear();



}


///////////////////////////////////////////////////////
/////         Attitude Estimation PF              /////
///////////////////////////////////////////////////////


PF::~PF() {

	_particles_d_q_k.clear();
	_particles_q_k.clear();
	_particles_w_k.clear();
	_particles_x_k.clear();

}

PF::PF() {

	_f = -1.0;
	_N_p = -1;
	_h = -1.0;
	seed = 0;
	_sum_w_k = 0.0;
	_sum_n_eff = 0.0;
	resample = false;

	THREESHOLD = -1.0;

	_variance_u = 0.000031623;
	_variance_v = 0.0031623;

}

void PF::initialization(Quaternion q_0,int argc,...){

	if(argc != 0){

		double* values = new double[argc];

		va_list arguments;

		va_start ( arguments, argc );

		for ( int i = 0; i < argc; i++ )
			values[i] =  va_arg ( arguments, double );

		va_end(arguments);

		_f = values[0];
		_N_p = values[1];
		THREESHOLD = values[2];
		_h = values[3];
		_variance_u = values[4];
		_variance_v = values[5];

	}else{

		_f = 4.0;
		_N_p = 2000;
		THREESHOLD = 4.0* static_cast<double>(_N_p) /7.0;
		_h = 0.1;
		_variance_u = 0.000031623;
		_variance_v = 0.0031623;

	}

	// init state
	_q_k_est = q_0;
	_x_k = Vector(6,0.0,0.0,0.0,_bias_gyro.getValue(0),_bias_gyro.getValue(1),_bias_gyro.getValue(2));

	// compute Cholesky Decomposition of matrix Q in order to generate noise
	_L_m = Matrix(6,6);
	for(int i=0;i<3;i++){

		double tmp= ((_variance_v*_variance_v)/DELTA_T) + ((_variance_u*_variance_u*DELTA_T)/12.0);

		_L_m.setValue(i,i,sqrt(tmp));
		_L_m.setValue(i+3,i+3,(_variance_u*sqrt(DELTA_T)));
		_L_m.setValue(i+3,i,(-0.5*_variance_u*sqrt(DELTA_T)));

	}

	// compute covariance matrix R
	_R = Matrix(6,6);
	for(int i=0;i<3;i++){
		_R.setValue(i,i,_var_acc);
		_R.setValue(i+3,i+3,_var_mag);
	}
	//_R = _R.inverse();


	// generation of N particles
	for(int i=0;i<_N_p;i++){


		// perturbation of x_0
		Vector noise = whiteNoiseVector(0.0,_L_m,seed);
		seed += noise.getLength();

		Vector x_0_i = _x_k + noise;

		// get d_p_0_i
		Vector d_p_0_i(3,x_0_i.getValue(0),x_0_i.getValue(1),x_0_i.getValue(2));

		// compute d_q_0_i
		double d_p_0_i_norm = d_p_0_i.norm();
		double d_p_0_i_norm_square = d_p_0_i_norm*d_p_0_i_norm;

		double d_q_w_0_i = ((_f*_f) - d_p_0_i_norm_square) / ( (_f*_f) + d_p_0_i_norm_square );
		Vector d_q_xyz_0_i  = d_p_0_i*( (1.0 + d_q_w_0_i)/_f );
		Quaternion d_q_0_i (d_q_w_0_i,d_q_xyz_0_i.getValue(0),d_q_xyz_0_i.getValue(1),d_q_xyz_0_i.getValue(2));

		// compute q_0_i
		Quaternion q_0_i = d_q_0_i*_q_k_est;
		q_0_i.normalize();

		// add to the list of particles
		_particles_x_k.push_back(x_0_i);
		_particles_d_q_k.push_back(d_q_0_i);
		_particles_q_k.push_back(q_0_i);
		_particles_wn_k.push_back(noise);

		// init weights at time 0
		_particles_w_k.push_back(1.0/static_cast<double>(_N_p));

	}

	resample = false;

}

void PF::preprocess(){

	_accelerometer.normalize();

	// compute measurement vector from sensor data
	_z_k = Vector(6,_accelerometer.getValue(0),_accelerometer.getValue(1),_accelerometer.getValue(2),
			_magnetometer.getValue(0),_magnetometer.getValue(1),_magnetometer.getValue(2));

}

void PF::prediction(){




	//*/
	//normal way

	// compute q_pred_k
	Vector angular_velocity = _gyroscope - Vector(3,_x_k.getValue(3),_x_k.getValue(4),_x_k.getValue(5));
	_q_k_pred = Quaternion( transitionOmega(angular_velocity)*_q_k_est.toVector() );

	 _sum_w_k = 0.0;

	// compute N particle of angular velocity and propagate _q_k_i to _q_(k+1)_i
	for(unsigned int i=0;i<_particles_w_k.size();i++){

		// generate normal random number
		Vector noise_w_k_i = _particles_wn_k[i].clone();

		// get bias_k_i
		Vector bias_k_i = Vector(3,_particles_x_k[i].getValue(3),_particles_x_k[i].getValue(4),_particles_x_k[i].getValue(5));

		// compute angular velocity
		Vector angular_velocity_i = _gyroscope - bias_k_i + Vector(3,noise_w_k_i.getValue(0),noise_w_k_i.getValue(1),noise_w_k_i.getValue(2));

		// propagate _q_k_i to _q_(k+1)_i
		_particles_q_k[i] = Quaternion( transitionOmega(angular_velocity_i)*_particles_q_k[i].toVector() );


		// propagate bias_k_i to bias_(k+1)_i
		bias_k_i = bias_k_i + Vector(3,noise_w_k_i.getValue(3),noise_w_k_i.getValue(4),noise_w_k_i.getValue(5));

		// compute d_q_k_i
		_particles_d_q_k[i] = _particles_q_k[i]*_q_k_pred.inverse();

		// propagate d_p_k_i to d_p_(k+1)_i
		float tmp =(_f * (signeOf(_particles_d_q_k[i].getQw())/(1.0 + abs(_particles_d_q_k[i].getQw()))));
		Vector d_p_kp_i = _particles_d_q_k[i].imaginary() * tmp;

		// update _x_k_i
		_particles_x_k[i].setValues(Vector(6,d_p_kp_i.getValue(0),d_p_kp_i.getValue(1),d_p_kp_i.getValue(2),bias_k_i.getValue(0),bias_k_i.getValue(1),bias_k_i.getValue(2)));

		// compute z_k_i
		Vector z_k_i_acc = rotationVector(Vector(3,0.0,0.0,1.0),_particles_q_k[i]);
		Vector z_k_i_mag = rotationVector(ned_geographic_north.clone(),_particles_q_k[i]);
		Vector z_k_i = Vector(6,z_k_i_acc.getValue(0),z_k_i_acc.getValue(1),z_k_i_acc.getValue(2),z_k_i_mag.getValue(0),z_k_i_mag.getValue(1),z_k_i_mag.getValue(2));

		// likehood function
		Vector s_k = _z_k - z_k_i;

		float tmp_acc = (1.0/_var_acc)* ( ( s_k.getValue(0)*s_k.getValue(0)) + ( s_k.getValue(1)*s_k.getValue(1)) +( s_k.getValue(2)*s_k.getValue(2)) );
		float tmp_mag = (1.0/_var_mag)* ( ( s_k.getValue(3)*s_k.getValue(3)) + ( s_k.getValue(4)*s_k.getValue(4)) +( s_k.getValue(5)*s_k.getValue(5)) );

		double L_k_i = exp( (tmp_acc+tmp_mag)*(-0.5) );

		// update weights
		_particles_w_k[i] *= L_k_i;

		// accumulator to normalize weights
		_sum_w_k += _particles_w_k[i];

	}

	/*/
	//revisit way

	_sum_w_k = 0.0;

	for(unsigned int i=0;i<_particles_w_k.size();i++){

		// generate normal random number
		Vector noise_w_k_i = _particles_wn_k[i].clone();

		// get bias_k_i
		Vector bias_k_i = Vector(3,_particles_x_k[i].getValue(3),_particles_x_k[i].getValue(4),_particles_x_k[i].getValue(5));

		// compute angular velocity
		Vector angular_velocity_i = _gyroscope - bias_k_i + Vector(3,noise_w_k_i.getValue(0),noise_w_k_i.getValue(1),noise_w_k_i.getValue(2));

		// propagate _q_k_i to _q_(k+1)_i

		_particles_q_k[i] = Quaternion( transitionOmega(angular_velocity_i)*_particles_q_k[i].toVector() );

		// propagate bias_k_i to bias_(k+1)_i
		bias_k_i = bias_k_i + Vector(3,noise_w_k_i.getValue(3),noise_w_k_i.getValue(4),noise_w_k_i.getValue(5));

		for(int j=0;j<3;j++)
			_particles_x_k[i].setValue(j+3,bias_k_i.getValue(j));

		// compute z_k_i
		Vector z_k_i_acc = rotationVector(ned_gravity.clone(),_particles_q_k[i]);
		Vector z_k_i_mag = rotationVector(ned_geographic_north.clone(),_particles_q_k[i]);
		Vector z_k_i = Vector(6,z_k_i_acc.getValue(0),z_k_i_acc.getValue(1),z_k_i_acc.getValue(2),z_k_i_mag.getValue(0),z_k_i_mag.getValue(1),z_k_i_mag.getValue(2));

		// likehood function
		Vector s_k = _z_k - z_k_i;

		float tmp_acc = (1.0/_var_acc)* ( ( s_k.getValue(0)*s_k.getValue(0)) + ( s_k.getValue(1)*s_k.getValue(1)) +( s_k.getValue(2)*s_k.getValue(2)) );
		float tmp_mag = (1.0/_var_mag)* ( ( s_k.getValue(3)*s_k.getValue(3)) + ( s_k.getValue(4)*s_k.getValue(4)) +( s_k.getValue(5)*s_k.getValue(5)) );

		double L_k_i = exp( (tmp_acc+tmp_mag)*(-0.5) );

		// update weights
		_particles_w_k[i] *= L_k_i;

		// accumulator to normalize weights
		_sum_w_k += _particles_w_k[i];

	}

	Matrix M(4,4);

	_sum_n_eff = 0.0;
	_q_k_pred = Quaternion(0.0,0.0,0.0,0.0);
	for(unsigned int i=0;i<_particles_w_k.size();i++){

		if(_sum_w_k > 0.0)
			_particles_w_k[i] /= _sum_w_k;
		else
			_particles_w_k[i] = 1.0/static_cast<double>(_N_p);

		_sum_n_eff += (_particles_w_k[i]*_particles_w_k[i]);

		Matrix tmp = vectorToMatrix(_particles_q_k[i].toVector());
		M = M + ((tmp*tmp.transpose())*_particles_w_k[i]);

		//_q_k_pred = _q_k_pred + (_particles_q_k[i]*_particles_w_k[i]);
	}

	vector<double> eigen_value;
	vector<Vector> eigen_vector;

	getEigenValuesVector(M,eigen_value,eigen_vector,50);

	int index_max=0;
	double max_eigen=0.0;

	for(unsigned int i=0;i<eigen_vector.size();i++){

		if(max_eigen<eigen_value[i]){
			index_max=i;
			max_eigen = eigen_value[i];
		}
	}

	_q_k_pred = Quaternion(eigen_vector[index_max]);

	eigen_value.clear();
	eigen_vector.clear();

	cout << _q_k_pred << endl;
	//*/

}

void PF::swap(int i,int j){

	//swap weight
	double temp = _particles_w_k[(unsigned int)i];
	_particles_w_k[(unsigned int)i] = _particles_w_k[(unsigned int)j];
	_particles_w_k[(unsigned int)j] = temp;

	//swap particle
	Vector vec_tmp = _particles_x_k[(unsigned int)i];
	_particles_x_k[(unsigned int)i].setValues(_particles_x_k[(unsigned int)j]);
	_particles_x_k[(unsigned int)j].setValues(vec_tmp);

}

void PF::quicksort( int left, int right){
	int min = (left+right)/2;

	int i = left;
	int j = right;
    double pivot = _particles_w_k[(unsigned int)min];

    while(left<j || i<right)
    {
        while(_particles_w_k[(unsigned int)i]<pivot)
        i++;
        while(_particles_w_k[(unsigned int)j]>pivot)
        j--;

        if(i<=j){
            swap(i,j);
            i++;
            j--;
        }
        else{
            if(left<j)
                quicksort(left, j);
            if(i<right)
                quicksort( i,right);
            return;
        }
    }
}

void PF::resampling2(){

	// set bool to false
	resample = false;

	double N_eff;

	N_eff = floor(1.0/_sum_n_eff);

	//cout << "N_eff = " << N_eff << endl;

	if( N_eff < THREESHOLD){

		resample = true;
		quicksort(0,_N_p-1);

		int n = (int)N_eff/3;

		vector<Vector> particle_chosen;

		for (int l = 0; l < n; ++l) {
			particle_chosen.push_back(_particles_x_k[(unsigned int)(_N_p-l-1)]);

		}

		int index=0;

		for (int k = 0; k < (_N_p - n); ++k) {
				index = k%n;
				_particles_x_k[(unsigned int)(k)].setValues(particle_chosen[index]);
				_particles_w_k[(unsigned int)(k)] = (1.0/static_cast<double>(_N_p));


		}

		particle_chosen.clear();


	}


}

void PF::correction(){

	/*/
	//revisit way


	for(unsigned int i=0;i<_particles_w_k.size();i++){

		// compute _d_q_(k+1)_i from _q_(k+1)_i
		_particles_d_q_k[i] = _particles_q_k[i]*_q_k_pred.inverse();

		// propagate d_p_k_i to d_p_(k+1)_i
		Vector d_p_kp1_i = _particles_d_q_k[i].imaginary() *
				(_f * (signeOf(_particles_d_q_k[i].getQw())/(1.0+abs(_particles_d_q_k[i].getQw()))));

		// update _x_k_i
		for(int j = 0;j<3;++j)
			_particles_x_k[i].setValue(j,d_p_kp1_i.getValue(j));


	}

	_x_k = Vector(6.0);
	Vector x_k_sigma(6.0);

	for(unsigned int i=0;i<_particles_w_k.size();i++){
		_x_k = _x_k + (_particles_x_k[i]*_particles_w_k[i]);
		x_k_sigma = x_k_sigma + _particles_x_k[i];

	}

	/*/

	// compute mean x_k_f
	_x_k = Vector(6.0);
	Vector x_k_sigma(6.0);

	_sum_n_eff = 0.0;

	// update weights
	for(unsigned int i=0;i<_particles_w_k.size();i++){

		_particles_w_k[i] /= _sum_w_k;

		_sum_n_eff += (_particles_w_k[i]*_particles_w_k[i]);

		_x_k = _x_k + (_particles_x_k[i]*_particles_w_k[i]);
		x_k_sigma = x_k_sigma + _particles_x_k[i];

	}

	//*/

	x_k_sigma = x_k_sigma *(1.0/(static_cast<double>(_N_p)));


	/* compute _q_est_k */
	Vector d_p_k(3,_x_k.getValue(0),_x_k.getValue(1),_x_k.getValue(2));
	double d_p_k_norm = d_p_k.norm();
	double d_p_k_norm_square = d_p_k_norm*d_p_k_norm;

	double d_q_w = (_f*_f - d_p_k_norm_square)/(_f*_f + d_p_k_norm_square);
	Vector d_q_xyz  = d_p_k*((1.0/_f)*(1.0+d_q_w));

	Quaternion d_q_k_f (d_q_w,d_q_xyz.getValue(0),d_q_xyz.getValue(1),d_q_xyz.getValue(2));
	_q_k_est = d_q_k_f*_q_k_pred;

	/* resampling step */
	resampling2();

	if(resample){

		/* compute covariance Sigma_k*/
		_Sigma_k = Matrix(6,6);
		for(unsigned int i=0;i<_particles_x_k.size();i++){

			Matrix mat_x_k_i_s = vectorToMatrix(_particles_x_k[i] - x_k_sigma );
			_Sigma_k = _Sigma_k + (mat_x_k_i_s*mat_x_k_i_s.transpose());
		}
		_Sigma_k = _Sigma_k*(1.0/(static_cast<double>(_N_p-1)));


		/* perturb x_k with some noise */
		Matrix Sigma_k_square_root = (_Sigma_k*(_h*_h)).choleskyDecomposition();

		if(Sigma_k_square_root.containsNan()){
			Sigma_k_square_root = matSquareRoot(_Sigma_k.clone()*(_h*_h));
		}


		for(unsigned int i=0;i<_particles_x_k.size();i++){

			// generate normal random number
			Vector noise_i = whiteNoiseVector(0.0,Sigma_k_square_root,seed);
			seed += noise_i.getLength();

			// perturb x_k */
			_particles_x_k[i].setValues(_particles_x_k[i] + noise_i);

			// "adding" pertupation to _particles_q_k
			Vector d_p_k_i(3,_particles_x_k[i].getValue(0),_particles_x_k[i].getValue(1),_particles_x_k[i].getValue(2));

			double norm_d_p_k_i = d_p_k_i.norm();
			double norm_d_p_k_i_square = norm_d_p_k_i*norm_d_p_k_i;

			double d_q_w_i = (_f*_f - norm_d_p_k_i_square)/(_f*_f + norm_d_p_k_i_square);
			Vector d_q_xyz_i  = d_p_k_i*((1.0+d_q_w_i)/_f);

			Quaternion d_q_k_f_i (d_q_w_i,d_q_xyz_i.getValue(0),d_q_xyz_i.getValue(1),d_q_xyz_i.getValue(2));

			_particles_d_q_k[i] = d_q_k_f_i;
			_particles_q_k[i] = d_q_k_f_i*_q_k_pred;


		}

 	}

}




} /* namespace isf */

