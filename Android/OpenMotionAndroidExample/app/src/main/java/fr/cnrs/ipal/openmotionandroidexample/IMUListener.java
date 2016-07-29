package fr.cnrs.ipal.openmotionandroidexample;

import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;

/**
 * Created by thomas on 4/12/16.
 */
public class IMUListener implements SensorEventListener {

    // output of the accelerometer
    private float[] mAccelerometerData;

    // output of the gyroscope
    private float[] mGyroscopeData;

    // output of the magnetometer
    private float[] mMagnetometerData;

    // output of OpenMotion (quaternion representation has been chosen)
    private float[] quaternion_openmotion;

    // object who manage the native part of the app
    private NativeManager mNativeManager;

    // sensor manger
    private SensorManager sensorManager;

    // array of boolean use for IMU synchronization
    private boolean[] init_imu;

    // constructor
    public IMUListener(Context context) {

        // initialization of boolean array
        init_imu = new boolean[3];
        for(int i=0;i<3;++i)
            init_imu[i] = false;

        // initialization of native manager
        mNativeManager = new NativeManager();

        // initialization of IMU output
        mAccelerometerData = new float[3];
        mGyroscopeData = new float[3];
        mMagnetometerData = new float[3];

        // initialization of OpenMotion output
        quaternion_openmotion = new float[4];

        // initialization of sensor manager
        sensorManager = (SensorManager) context.getSystemService(Context.SENSOR_SERVICE);


    }


    public void onPause() {

        sensorManager.unregisterListener(this);

    }

    public void onDestroy() {

        sensorManager.unregisterListener(this);
        mNativeManager.destroyObject();

    }


    public void onResume() {

        sensorManager.registerListener(this, sensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER), SensorManager.SENSOR_DELAY_GAME);
        sensorManager.registerListener(this, sensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE), SensorManager.SENSOR_DELAY_GAME);
        sensorManager.registerListener(this, sensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD), SensorManager.SENSOR_DELAY_GAME);

    }

    @Override
    public void onAccuracyChanged(Sensor arg0, int arg1) {
    }

    @Override
    public void onSensorChanged(SensorEvent event) {

        Sensor sensor = event.sensor;

        // get accelerometer output
        // OpenMotion works only if all component from the IMU belong to the same reference frame
        if (sensor.getType() == Sensor.TYPE_ACCELEROMETER){

            mAccelerometerData[0] = event.values[0];
            mAccelerometerData[1] = event.values[1];
            mAccelerometerData[2] = event.values[2]*(-1.0f);

            mNativeManager.setAccData(mAccelerometerData[0],mAccelerometerData[1],mAccelerometerData[2]);

            init_imu[0] = true;
        }

        // get gyroscope output
        // OpenMotion works only if all component from the IMU belong to the same reference frame
        else if (sensor.getType() == Sensor.TYPE_GYROSCOPE){
            float rad2deg = 15.0f;

            mGyroscopeData[0] = Math.abs(event.values[0]) > 0.01f ? event.values[0]*rad2deg : 0.0f;
            mGyroscopeData[1] = Math.abs(event.values[1]) > 0.01f ? event.values[1]*rad2deg : 0.0f;
            mGyroscopeData[2] = Math.abs(event.values[2]) > 0.01f ? event.values[2]*rad2deg : 0.0f;

            mNativeManager.setGyroData(mGyroscopeData[0],mGyroscopeData[1],mGyroscopeData[2]);
            init_imu[1] = true;
        }

        // get magnetometer output
        // OpenMotion works only if all component from the IMU belong to the same reference frame
        else if (sensor.getType() == Sensor.TYPE_MAGNETIC_FIELD){
            mMagnetometerData[0] = event.values[0]*(-1.0f);
            mMagnetometerData[1] = event.values[1]*(-1.0f);
            mMagnetometerData[2] = event.values[2];
            mNativeManager.setMagData(mMagnetometerData[0],mMagnetometerData[1],mMagnetometerData[2]);
            init_imu[2] = true;
        }

        // if all values from the IMU has been get
        if( init_imu[0] && init_imu[1] && init_imu[2] ){

            //  estimation of the attitude from OpenMotion
            quaternion_openmotion = mNativeManager.getQuatAttitude();

            // re-set boolean to false
            for(int i=0;i<3;++i)
                init_imu[i] = false;
        }

    }


    // get the euler angle
    public float[] getEulerAngle(){
        return mNativeManager.getEulerAngle();
    }

    // get the gyroscope output
    public float[] getGyroData(){
        return mGyroscopeData;
    }

    // get the accelerometer output
    public float[] getAccData(){
        return mAccelerometerData;
    }

    // get the magnetometer output
    public float[] getMagData(){
        return mMagnetometerData;
    }

    // get the attitude estimation from OpenMotion
    public float[] getOpenMotionQuat(){
        return quaternion_openmotion;
    }


}
