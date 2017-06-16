function [ I_P ] = omCreateInitParams( DELTA_T, var_accel,var_magn,var_gyro, local_magn, bias_accel, bias_magn, bias_gyro, init_quat)
%[ I_P ] = omCreateInitParams( DELTA_T, var_accel,var_magn,var_gyro, local_magn, bias_accel, bias_magn, bias_gyro, init_quat)
%This function will create the structure for the OpenMotion library
% example :
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


if (nargin~=9)
    error('This function needs 9 arguments');
end

if (~(isscalar(DELTA_T) || isscalar(var_accel)|| isscalar(var_magn) || isscalar(var_gyro)))
    error('DELTA_T, var_accel,var_magn and var_gyro are scalars, check your input');
    
end

if ( (length(local_magn(:))~= 3) || (length(bias_accel(:))~= 3) || (length(bias_magn(:))~= 3) || (length(bias_gyro(:))~= 3) || (length(init_quat(:))~=4) )
    error('local_magn, bias_accel, bias_magn, bias_gyro are vectors of size 3 and init_quat is a vector of size 4, check your input');
end

I_P.DELTA_T=DELTA_T
I_P.var_accel=var_accel
I_P.var_magn=var_magn
I_P.var_gyro=var_gyro
I_P.local_magn=local_magn
I_P.bias_accel=bias_accel
I_P.bias_magn=bias_magn
I_P.bias_gyro=bias_gyro
I_P.init_quat=init_quat
