//
//  om_Matlab.c
//  
//
//  Created by Nizar Ouarti on 13/06/17.
//
//

#include "om_Matlab.h"
#include "utils_.h"
#include <stdio.h>
#include <mex.h>
#include <stdint.h>
typedef uint16_t uint16;

//#include "matrix.h"

void buildStruct(double DELTA_T,double var_accel, double var_gyro, double var_magn,double *local_magn,double *bias_accel,double *bias_gyro,double *bias_magn,double *init_quat, Init_Params *I_P){
    
    
    // time gap between 2 sensor recording
    I_P->DELTA_T = DELTA_T;
    
    
    // this value depends of your location on earth
    // magnetic field in micro Tesla ( according to https://www.ngdc.noaa.gov/geomag-web/ )
    I_P->local_magn[0]=local_magn[0];I_P->local_magn[1]=local_magn[1];I_P->local_magn[2]=local_magn[2];
    
    // these values depends of your hardware: variance and bias
    I_P->var_accel=var_accel;
    I_P->var_magn=var_magn;
    I_P->var_gyro=var_gyro;
    
    
    I_P->bias_accel[0]=bias_accel[0];I_P->bias_accel[1]=bias_accel[1];I_P->bias_accel[2]=bias_accel[2];
    I_P->bias_magn[0]=bias_magn[0];I_P->bias_magn[1]=bias_magn[1];I_P->bias_magn[2]=bias_magn[2];
    I_P->bias_gyro[0]=bias_gyro[0];I_P->bias_gyro[1]=bias_gyro[1];I_P->bias_gyro[2]=bias_gyro[2];
    
    
    // this value depends of your initial orientation
    I_P->init_quat[0]=init_quat[0];I_P->init_quat[1]=init_quat[1];I_P->init_quat[2]=init_quat[2];I_P->init_quat[3]=init_quat[3];
    
    

}

void Compute(double *data,Init_Params *I_P, char *method, int cols,double *quats){

//    if(strcmp(method,"USQUE")){
        // variable for data fusion
        omSensorFusionManager manager;
        omNonLinearFilter_USQUE filter;//HAVE TO BE MODIFIED!!!!!!! if define here cannot be seen outside these brackets
    
        // initialization of the sensor fusion manager and the nonlinear filter
        init_manager(&manager,&filter,USQUE,I_P);
//    }
    
    
    double gyroX;double gyroY;double gyroZ;
    double accX;double accY;double accZ;
    double magX;double magY;double magZ;
    int i;
    for (i=0;i<cols;i++){
        
        //extract data from the matrix
        accX=data[9*i];accY=data[9*i+1];accZ=data[9*i+2];
        magX=data[9*i+3];magY=data[9*i+4];magZ=data[9*i+5];
        gyroX=data[9*i+6];gyroY=data[9*i+7];gyroZ=data[9*i+8];
        
        // set gyro data
        om_vector_setValues(&manager.imu_data.data_gyroscope,3,(double)gyroX, (double)gyroY, (double)gyroZ);
        
        // set gyroelerometer data
        om_vector_setValues(&manager.imu_data.data_accelerometer,3,(double)accX,(double)accY,(double)accZ);
        
        // set magnetometer data
        om_vector_setValues(&manager.imu_data.data_magnetometer,3,(double)magX,(double)magY,(double)magZ);
        
        // process filtering
        manager.process_filter(&manager,&filter);
        
    
        quats[4*i]=manager.output.quaternion._qw;
        quats[4*i+1]=manager.output.quaternion._qx;
        quats[4*i+2]=manager.output.quaternion._qy;
        quats[4*i+3]=manager.output.quaternion._qz;
        //printf("%d: %lf \n",i,quats[4*i]);
    }

    
}

void mexFunction(int nout, mxArray *out[],
                 int nin, const mxArray *in[])
{
    const int *dimsdata;
    int ndims;
    uint16 rows, cols;
    
    if (nin != 3)
    {
        mexErrMsgTxt("Usage: [quats] = omCompute(data,Init_Params,method);");
    }
    if (nout > 1)
    {
        mexErrMsgTxt("At most one output arguments allowed");
    }
    
    if (!(mxIsClass(in[0], "double")))
    {
        mexErrMsgTxt("arg 1 of the function, data must be of type double");
    }
    if (!(mxIsClass(in[1], "struct")))
    {
        mexErrMsgTxt("arg 2 of the function, Init_Params is a structure create by omCreateInitParams");
    }
    if (!(mxIsClass(in[2], "char")))
    {
        mexErrMsgTxt("arg 3 of the function, method is a string");
    }
    
    ndims = mxGetNumberOfDimensions(in[0]);
    dimsdata = mxGetDimensions(in[0]);
    if (ndims != 2)
    {
        mexErrMsgTxt("data is a matrix with each line representing a sample from different sensors at a given time t");
    }
    
    rows = (uint16)dimsdata[0];
    cols = (uint16)dimsdata[1];
    
    if (rows != 9)
    {
        mexErrMsgTxt("The number of lines have to be 9 (3 sensors each x,y,z) in the order accelerometer, magnetometer and gyro");
    }
    
    data= (double *)mxGetData(in[0]);// data pointer

    //access to fields of the structure
    
    mxArray *DEL=mxGetField(in[1],0,"DELTA_T");
    double DELTA_T=(double)mxGetScalar(DEL) ;
    mxArray *VA=mxGetField(in[1],0,"var_accel");
    double var_accel=(double)mxGetScalar(VA) ;
    mxArray *VM=mxGetField(in[1],0,"var_magn");
    double var_magn=(double)mxGetScalar(VM) ;
    mxArray *VG=mxGetField(in[1],0,"var_gyro");
    double var_gyro=(double)mxGetScalar(VG) ;
    
    mxArray *LM=mxGetField(in[1],0,"local_magn");
    double *local_magn=(double*)mxGetData(LM) ;
    mxArray *BA=mxGetField(in[1],0,"bias_accel");
    double *bias_accel=(double*)mxGetData(BA) ;
    mxArray *BM=mxGetField(in[1],0,"bias_magn");
    double *bias_magn=(double*)mxGetData(BM) ;
    mxArray *BG=mxGetField(in[1],0,"bias_gyro");
    double *bias_gyro=(double*)mxGetData(BG) ;
    mxArray *IQ=mxGetField(in[1],0,"init_quat");
    double *init_quat=(double*)mxGetData(IQ) ;
   
    //access to the string for the method
    size_t buflen;
    int status;
    buflen = mxGetN(in[2]) + 1;
    method = mxMalloc(buflen);
    
    status = mxGetString(in[2], method, (mwSize)buflen);
    if (~status) printf("The fusion method chosen is:  %s\n", method);
    
    //Output of the function
    out[0] = mxCreateNumericMatrix(4, cols, mxDOUBLE_CLASS, mxREAL);// quaternion 4 components
    double *quats = (double*)mxGetData(out[0]);// output pointer
    
    Init_Params I_P;//set the parameters of the structure Init_Params
    buildStruct(DELTA_T,var_accel,var_gyro,var_magn,local_magn,bias_accel,bias_gyro,bias_magn,init_quat, &I_P);
    
    Compute(data, &I_P, method, cols, quats);// compute the fusion
    
    // free memory
    mxFree(method);
    
    
    
}
