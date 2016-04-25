/*
 * isf.cpp
 *
 *  Created on: 30 Nov, 2015
 *      Author: thomas
 */

#include "om.h"

namespace om {


SensorFusionManager::SensorFusionManager(){

	param_1 = param_2 = param_3 = param_4 = param_5 = -1.0;

	methodName = M_HYOUKF;
	IOAType = UnitQuaternion;

	/* default param */
	Vector vector_bias_gyroscope = Vector(3,0.0031623,0.0031623,0.0031623);
	Vector vector_bias_accelerometer = Vector(3,0.03,0.03,0.03);
	Vector vector_bias_magnetometer = Vector(3,0.03,0.01,-0.06);

	double double_variance_gyro = 0.0031623;
	double double_variance_acc = 0.02;
	double double_variance_mag = 0.05;

	filter = new HyOKF();
	filter->setSensorsInfo(vector_bias_gyroscope,vector_bias_accelerometer,vector_bias_magnetometer,double_variance_gyro,double_variance_acc,double_variance_mag);
	filter->initialization(Quaternion(1.0,0.0,0.0,0.0));

}

void SensorFusionManager::setAccelerometerData(const Vector& data){

	filter->setAccelerometerData(data);

}

void SensorFusionManager::setGyroscopeData(const Vector& data){

	filter->setGyroscopeData(data);

}

void SensorFusionManager::setMagnetometerData(const Vector& data){

	filter->setMagnetometerData(data);

}


void SensorFusionManager::setGyroscopeData(double x_value,double y_value,double z_value){

	filter->setGyroscopeData(Vector(3,x_value,y_value,z_value));

}


void SensorFusionManager::setAccelerometerData(double x_value,double y_value,double z_value){

	filter->setAccelerometerData(Vector(3,x_value,y_value,z_value));

}

void SensorFusionManager::setMagnetometerData(double x_value,double y_value,double z_value){

	filter->setMagnetometerData(Vector(3,x_value,y_value,z_value));

}

Vector SensorFusionManager::getGyroscopeData(){
	return filter->getGyroscopeData();
}

Vector SensorFusionManager::getAccelerometerData(){
	return filter->getAccelerometerData();
}

Vector SensorFusionManager::getMagnetometerData(){
	return filter->getMagnetometerData();
}




void SensorFusionManager::readXMLFile(char* file_name){
	cout << "Not ready yet! Cannot read :" << file_name << endl;
}

template<typename T>
T SensorFusionManager::getAttitude(){
	cerr << "Output Type no specified !" << endl;
	return NULL;
}

template<>
Quaternion SensorFusionManager::getAttitude<Quaternion>(void){

	filter->process();
	Quaternion q_est = filter->getOrientation();

	return q_est;
}


template<>
Matrix SensorFusionManager::getAttitude<Matrix>(void){

	filter->process();
	Quaternion q_est = filter->getOrientation();

	return quatToRotationMatrix(q_est);

}

template<>
AxisAngle SensorFusionManager::getAttitude<AxisAngle>(void){

	filter->process();
	Quaternion q_est = filter->getOrientation();

	return quatToAxisAngle(q_est);

}






} /* namespace isf */
