/*
 * NonLinearFilter.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#ifndef NONLINEARFILTER_H_
#define NONLINEARFILTER_H_


#include "geometry.h"

namespace om {


extern Vector ned_gravity;
extern Vector ned_geographic_north;

class NonLinearFilter {

public:

	NonLinearFilter();
	virtual ~NonLinearFilter(){}


	virtual void prediction() = 0;
	virtual void correction() = 0;
	virtual void preprocess() = 0;
	virtual void initialization(Quaternion q_0,int argc=0,...) = 0;

	void process();

	Quaternion getOrientation(){ return _q_k_est; }
	void setSensorsInfo(const Vector& bias_gyro,const Vector& bias_acc,const Vector& bias_mag,double var_gyro,double var_acc,double var_mag);
	void setSensorsData(const Vector& gyro,const Vector& acc,const Vector& mag);

	void setGyroscopeData(const Vector& data);
	void setAccelerometerData(const Vector& data);
	void setMagnetometerData(const Vector& data);

	Vector getGyroscopeData();
	Vector getAccelerometerData();
	Vector getMagnetometerData();




protected:

	Vector _gyroscope;
	Vector _accelerometer;
	Vector _magnetometer;

	Vector _bias_gyro;
	Vector _bias_acc;
	Vector _bias_mag;

	double _var_gyro;
	double _var_acc;
	double _var_mag;

	double _seed;

	Quaternion _q_k_est;
	Quaternion _q_k_pred;

};

class MEKF: public NonLinearFilter{

public:

	MEKF();
	virtual ~MEKF();

	virtual void prediction();
	virtual void correction();
	virtual void preprocess();
	virtual void initialization(Quaternion q_0,int argc=0,...);




private:

	void computeJacobianF();
	void computeJacobianH();

	Vector f_function(const Vector& x);
	Vector h_function(const Vector& x);

	Matrix W(const Vector& e,const Vector& v);

	Vector _x_k_est;
	Matrix _P_k;
	Vector _z_k;

	Vector _v_k;
	Vector _w_k;

	Vector _z_k_pred;
	Vector _x_k_pred;
	Matrix _P_k_pred;

	Matrix _Q_k;
	Matrix _R_k;

	Matrix _Q_cho_k;
	Matrix _R_cho_k;

	Matrix _F_k;
	Matrix _H_k;

	double _h;
	double _a;
	double _f;

};

class USQUE: public NonLinearFilter{

public:

	USQUE();
	virtual ~USQUE();

	 void prediction();
	 void correction();
	 void preprocess();
	 void initialization(Quaternion q_0,int argc=0,...);


private:

	Vector _x_k_est;
	Matrix _P_k;
	Vector _z_k;

	Vector _z_k_pred;
	Vector _x_k_pred;
	Matrix _P_k_pred;

	Matrix _Q_k;
	Matrix _R_k;

	Matrix _cov_Z;
	Matrix _cov_Z_x;

	vector<Vector> _sigma_points_x_k;
	vector<Vector> _sigma_points_x_kp1;
	vector<Vector> _sigma_points_z_kp1;
	vector<Quaternion> _sigma_quaternion_k;

	double _f;
	double _a;
	double _lambda;
	double _n;
	double _variance_u;
	double _variance_v;

};

class CGO : public NonLinearFilter{

public:

	CGO();
	virtual ~CGO();

	 void prediction();
	 void correction();
	 void preprocess();
	 void initialization(Quaternion q_0,int argc=0,...);

private:

	Vector _bias_est;
	Vector _bias_pred;

	Vector _omega;

	double _k_mag;
	double _k_acc;

	double _k_I;
	double _k_P;

	Vector _v_acc_pred;
	Vector _v_mag_pred;

};

class REQUEST  : public NonLinearFilter{
public:
	REQUEST();
	virtual ~REQUEST();

	 void prediction();
	 void correction();
	 void preprocess();
	 void initialization(Quaternion q_0,int argc=0,...);




private:

	void computeQk();
	void computeRk();

	Matrix _P_pred_k;
	Matrix _P_k;

	Matrix _R_k;
	Matrix _Q_k;


	Matrix _K_pred_k;
	Matrix _K_est_k;
	Matrix _d_K_k;

	Matrix _d_B_k;
	Matrix _d_S_k;
	Vector _d_z_k;
	double _d_m_k;
	double _d_sigma_k;

	double _lambda_m_k;
	double _mu_k;
	double _m_k;

	vector<double> _a;
	vector<Vector> _b;
	vector<Vector> _r;
};


class PF: public NonLinearFilter{
public:
	PF();
	virtual ~PF();

	 void prediction();
	 void correction();
	 void preprocess();
	 void initialization(Quaternion q_0,int argc=0,...);

	void resampling2();

	void quicksort(int left, int right);
	void swap(int i,int j);


private:

    double THREESHOLD;
	int _N_p;
	double _f;
	double _h;
	double _variance_u;
	double _variance_v;
	double _sum_w_k;
	double _sum_n_eff;

	int seed;

	Vector _x_k;
	Vector _z_k;

	vector<Vector> _particles_x_k;
	vector<Quaternion> _particles_q_k;
	vector<Quaternion> _particles_d_q_k;
	vector<Vector> _particles_wn_k;
	vector<double> _particles_w_k;

	Matrix _P_k;
	Matrix _P_k_f;

	Matrix _Sigma_k;

	Matrix _L_m;
	Matrix _R;

	bool resample;

};





} /* namespace isf */

#endif /* NONLINEARFILTER_H_ */
