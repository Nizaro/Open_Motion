/*
 * isf.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */



#ifndef OM_H_
#define OM_H_


#include "NonLinearFilter.h"

namespace om {

enum MethodName{
	M_MEKF,M_USQUE,M_PF,M_CGO,M_REQUEST,M_HYOUKF
};

enum InputOutputAttitudeType{
	RotationMatrix,UnitQuaternion,AngleAxis,EulerAngle,Rodriguez
};



class SensorFusionManager {

public:

	SensorFusionManager();
	virtual ~SensorFusionManager(){}

	void setGyroscopeData(const Vector& data);
	void setGyroscopeData(double x_value,double y_value,double z_value);

	void setAccelerometerData(const Vector& data);
	void setAccelerometerData(double x_value,double y_value,double z_value);

	void setMagnetometerData(const Vector& data);
	void setMagnetometerData(double x_value,double y_value,double z_value);


	Vector getGyroscopeData();
	Vector getAccelerometerData();
	Vector getMagnetometerData();

	void readXMLFile(char* file_name);

	template<typename T>
	T getAttitude(void);


private:

	MethodName methodName;
	InputOutputAttitudeType IOAType;

	NonLinearFilter* filter;

	double param_1;
	double param_2;
	double param_3;
	double param_4;
	double param_5;


};

template<>
Quaternion SensorFusionManager::getAttitude<Quaternion>(void);

template<>
Matrix SensorFusionManager::getAttitude<Matrix>(void);

template<>
AxisAngle SensorFusionManager::getAttitude<AxisAngle>(void);




} /* namespace isf */

#endif /* ISF_H_ */

