function quats=om_Matlab(data,O_P,method)
% quat=om_Matlab(data,O_P,'USQUE')
% + data : have to be a 9XN array with lines 1:3 are accelerometer
% X,Y,Z, lines 4:6 are the magnetometer X<Y<Z and lines 7:9 are the gyro X,Y,Z 
% Init_Params:  is a structure create by omCreateInitParams
% + method: the different methods accessible are :
%USQUE
%CGO
%MEKF
%
%
% + quats: is the quaternion q(t) that is the estimation of the attitude over
% the time.
%
%
% Usage:
% OMPATH='/Users/scchia/Desktop/Open_Motion/'
% cd ([OMPATH,'Matlab'])
% DELTA_T=0.03;
% var_accel=0.75;
% var_magn=0.075;
% var_gyro=0.031623;
% local_magn=[40.810,10.1609,10.2636];
% bias_accel=[0,0,0];
% bias_magn=[0,0,0];
% bias_gyro=[0.000031623,0.0000316230,0.00003162];
% init_quat=[1,0,0,0];
% [ I_P ] = omCreateInitParams( DELTA_T, var_accel,var_magn,var_gyro,local_magn, bias_accel, bias_magn, bias_gyro, init_quat);
% tmp = csvread('../Fusion_Algorithms/Samples/imu_data.csv',1,0); 
% data=[tmp(:,[12,13,14])';tmp(:,[15,16,17])';tmp(:,[9,10,11])'];%ensure that data is 9XN 
% quats=om_Matlab(data,I_P,'USQUE');
% figure, plot(quats')

%end

