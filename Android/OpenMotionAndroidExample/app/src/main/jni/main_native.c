//
// Created by thomas on 6/6/16.
//

#include "main_native.h"

// transform the magnetic field to the geographic north (depend on the location)
void modification_magneto(struct omVector *input,struct omVector* output,struct omQuaternion *q_est){


    omQuaternion q_y;
    omQuaternion q_z;
    omQuaternion q_y_inv;
    omQuaternion q_z_inv;
    omMatrix M;
    omVector tmp;
    omVector tmp2;
    omVector tmp3;
    omAxisAngle aa_y;
    omAxisAngle aa_z;

    // this values works in singapore area
    aa_y._angle = -14.36;
    aa_z._angle = 0.23;

    om_vector_create(&aa_y._axis,3,0.0,1.0,0.0);
    om_vector_create(&aa_z._axis,3,0.0,0.0,1.0);

    om_convert_axisAngle2quaternion(&aa_y,&q_y);
    om_convert_axisAngle2quaternion(&aa_z,&q_z);


    om_quat_inverse(&q_y,&q_y_inv);
    om_quat_inverse(&q_z,&q_z_inv);

    om_matrix_create(&M,3,3);

    om_convert_quaternion2matrix(q_est,&M);

    om_vector_create(&tmp,3,0.0,0.0,0.0);
    om_vector_create(&tmp2,3,0.0,0.0,0.0);
    om_vector_create(&tmp3,3,0.0,0.0,0.0);
    om_vector_create(output,3,0.0,0.0,0.0);

    om_operator_matrix_vector_mul(&M,input,&tmp);

    om_rotate_vector_quaternion(&q_z_inv,&tmp,&tmp2);
    om_rotate_vector_quaternion(&q_y_inv,&tmp2,&tmp3);

    om_rotate_vector_quaternion(q_est,&tmp3,output);

    om_matrix_free(&M);
    om_vector_free(&tmp);
    om_vector_free(&tmp2);
    om_vector_free(&tmp3);
    om_vector_free(&aa_y._axis);
    om_vector_free(&aa_z._axis);


}


void init_manager(omSensorFusionManager *manager){

    // allocation for sensor bias
    om_vector_create(&manager->imu_params.bias_accelerometer,3,0.0,0.0,0.0);
    om_vector_create(&manager->imu_params.bias_magnetometer,3,0.0,0.0,0.0);
    om_vector_create(&manager->imu_params.bias_gyroscope,3,0.00031623,0.00031623,0.00031623);

    // temporary values for sensor variance (this won't work on every device)
    manager->imu_params.variance_accelerometer = 0.05;
    manager->imu_params.variance_magnetometer = 0.023;
    manager->imu_params.variance_gyroscope = 0.0031623;

    // quaternion representation has been chosen here
    manager->type = Quarternion;

    // allocation for sensor output
    om_vector_create(&manager->imu_data.data_accelerometer,3);
    om_vector_create(&manager->imu_data.data_magnetometer,3);
    om_vector_create(&manager->imu_data.data_gyroscope,3);

    // set the appropriate function
    manager->initialization_filter = &om_usque_initialization;
    manager->process_filter = &om_usque_process;
    manager->free_filter = &om_usque_free;

    // initialization of the quaternion
    om_quat_create(&manager->output.quaternion,1.0,0.0,0.0,0.0);

}



/**
 * Native function createObject
 */

JNIEXPORT jlong JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeCreateObject(JNIEnv * jenv, jobject thiz)
{


    jlong result = 0;

    // initialization of the sensor fusion manager
    init_manager(&manager);

    // initialization of the nonlinear filter
    manager.initialization_filter(&manager,&filter_usque);


    return result;
}

JNIEXPORT jlong JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeDestroyObject(JNIEnv * jenv, jobject thiz)
{

    jlong result = 0;

    // release nonlinear filter
    manager.free_filter(&filter_usque);

    return result;
}


JNIEXPORT jfloatArray JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeGetQuatAttitude(JNIEnv * jenv, jobject thiz) {

    // process the method
     manager.process_filter(&manager,&filter_usque);

    // get the sign of the real part of the quaternion
    float sign = manager.output.quaternion._qw > 0.0 ? 1.0 : -1.0;

    // get values
    float* tmp = (float*)malloc(sizeof(float)*4);
    tmp[0]=(float)manager.output.quaternion._qx * sign;
    tmp[1]=(float)manager.output.quaternion._qy * sign;
    tmp[2]=(float)manager.output.quaternion._qz * sign;
    tmp[3]=(float)manager.output.quaternion._qw * sign;

    // allocate jfloatArray
    jfloatArray res = (*jenv)->NewFloatArray(jenv,4);
    if (res == NULL) {
        return NULL; /* out of memory error thrown */
    }

    (*jenv)->SetFloatArrayRegion(jenv,res,0,4,tmp);

    free(tmp);


    return res;

}


JNIEXPORT jfloatArray JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeGetEulerAngle
        (JNIEnv * jenv, jobject thiz){


    float* tmp = (float*)malloc(sizeof(float)*3);


    if(!isnan(manager.output.quaternion._qw) &&!isnan(manager.output.quaternion._qx)
       &&!isnan(manager.output.quaternion._qy) && !isnan(manager.output.quaternion._qz)){

        omEulerAngle eulerAngle;
        om_convert_quaternion2euler(&manager.output.quaternion,&eulerAngle);

        tmp[0]=eulerAngle._pitch;
        tmp[1]=eulerAngle._roll;
        tmp[2]=eulerAngle._yaw;

    }else{
        tmp[0]=0.0f;
        tmp[1]=0.0f;
        tmp[2]=0.0f;

    }


    jfloatArray res = (*jenv)->NewFloatArray(jenv,3);
    if (res == NULL) {
        return NULL; /* out of memory error thrown */
    }

    (*jenv)->SetFloatArrayRegion(jenv,res,0,3,tmp);

    free(tmp);


    return res;





}






/**
 * Native function setGyroData
 */
JNIEXPORT void JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeSetGyroData
        (JNIEnv * jenv, jobject thiz, jfloat gyroX,jfloat gyroY,jfloat gyroZ)
{

    // send values to the sensor fusion manager
    om_vector_setValues(&manager.imu_data.data_gyroscope,3,(double)gyroX, (double)gyroY, (double)gyroZ);

}




/**
 * Native function setAccData
 */

JNIEXPORT void JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeSetAccData
        (JNIEnv * jenv, jobject thiz, jfloat accX,jfloat accY,jfloat accZ)
{

    // send values to the sensor fusion manager
    om_vector_setValues(&manager.imu_data.data_accelerometer,3,(double)accX, (double)accY, (double)accZ);

}


/**
 * Native function setMagData
 */
JNIEXPORT void JNICALL Java_fr_cnrs_ipal_openmotionandroidexample_NativeManager_nativeSetMagData
        (JNIEnv * jenv, jobject thiz, jfloat magX,jfloat magY,jfloat magZ){


    omVector mag;
    omVector mag_calib;

    // allocation
    om_vector_create(&mag,3,(double)magX, (double)magY, (double)magZ);

    // transform the magnetic field to the geographic north (depend on the location)
    modification_magneto(&mag,&mag_calib,&manager.output.quaternion);

    // normalize the vector
    om_vector_normalize(&mag_calib);

    // send values to the sensor fusion manager
    om_vector_setValues(&manager.imu_data.data_magnetometer,3,om_vector_getValue(&mag_calib,0),
                                                              om_vector_getValue(&mag_calib,1),
                                                              om_vector_getValue(&mag_calib,2));

    // release memory
    om_vector_free(&mag);
    om_vector_free(&mag_calib);


}

