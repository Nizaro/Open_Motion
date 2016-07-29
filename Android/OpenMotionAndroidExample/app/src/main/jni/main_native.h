//
// Created by thomas on 6/6/16.
//

#ifndef OPENMOTIONANDROIDEXAMPLE_MAIN_NATIVE_H
#define OPENMOTIONANDROIDEXAMPLE_MAIN_NATIVE_H


#include <android/log.h>
#include <android/sensor.h>
#include <android/looper.h>
#include <time.h>

#include <jni.h>
#include <om.h>


// sensor fusion manager
omSensorFusionManager manager;

// nonlinear filter USQUE
omNonLinearFilter_USQUE filter_usque;

// transform the magnetic field to the geographic north (depend on the location)
void modification_magneto(struct omVector *input,struct omVector* output,struct omQuaternion *q_est);

// initialization of the sensor fusion manager
void init_manager(omSensorFusionManager *manager);


///////////////////////////////////////////////////////
/////               JNI Functions                 /////
///////////////////////////////////////////////////////


JNIEXPORT void JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeSetGyroData
        (JNIEnv * jenv, jobject thiz, jfloat gyroX,jfloat gyroY,jfloat gyroZ);

JNIEXPORT void JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeSetAccData
        (JNIEnv * jenv, jobject thiz, jfloat accX,jfloat accY,jfloat accZ);

JNIEXPORT void JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeSetMagData
        (JNIEnv * jenv, jobject thiz, jfloat magX,jfloat magY,jfloat magZ);

JNIEXPORT jlong JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeCreateObject
        (JNIEnv * jenv, jobject thiz);

JNIEXPORT jlong JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeDestroyObject
        (JNIEnv * jenv, jobject thiz);

JNIEXPORT jfloatArray JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeGetQuatAttitude
        (JNIEnv * jenv, jobject thiz);

JNIEXPORT jfloatArray JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeGetEulerAngle
        (JNIEnv * jenv, jobject thiz);


#endif //OPENMOTIONANDROIDEXAMPLE_MAIN_NATIVE_H
