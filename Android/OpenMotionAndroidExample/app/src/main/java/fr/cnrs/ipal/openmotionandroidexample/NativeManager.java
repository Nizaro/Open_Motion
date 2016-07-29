package fr.cnrs.ipal.openmotionandroidexample;

/**
 * Created by thomas on 4/12/16.
 */
public class NativeManager {


    private long mNativeObj = 0;

    // load native library
    static {
        System.loadLibrary("native_lib");
    }

    // constructor
    public NativeManager() {

        // initialization of native objects
        mNativeObj = nativeCreateObject();

    }


    public void destroyObject(){

        // release all native objects
        mNativeObj = nativeDestroyObject();
    }


    // send gyroscope data to native objects
    public void setGyroData(float gyroX,float gyroY,float gyroZ){
        nativeSetGyroData(gyroX, gyroY, gyroZ);
    }

    // send accelerometer data to native objects
    public void setAccData(float accX,float accY,float accZ){
        nativeSetAccData( accX, accY, accZ);
    }

    // send magnetometer data to native objects
    public void setMagData(float magX,float magY,float magZ){

        nativeSetMagData(magX, magY, magZ);
    }

    // get estimated attitude from native library
    public float[] getQuatAttitude() {

        float[] res = nativeGetQuatAttitude();

        if (res == null) {
            res = new float[4];
            res[3] = 1.0f;
        }
        return res;
    }

    // get euler angle from native library
    public float[] getEulerAngle(){
        float[] res = nativeGetEulerAngle();
        if (res == null) {
            res = new float[3];
        }

        return res;
    }

    public void release() {

        mNativeObj = 0;
    }


    //all native methods
    private static native long nativeCreateObject();
    private static native long nativeDestroyObject();
    private static native float[] nativeGetQuatAttitude();
    private static native float[] nativeGetEulerAngle();
    private static native void nativeSetGyroData( float gyroX,float gyroY,float gyroZ);
    private static native void nativeSetAccData( float accX,float accY,float accZ);
    private static native void nativeSetMagData( float accX,float magY,float magZ);

}
