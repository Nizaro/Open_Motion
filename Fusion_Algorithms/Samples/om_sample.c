/*
 * om_example.c
 *
 *  Created on: 1 Aug, 2016
 *      Author: Thomas BRAUD, Nizar OUARTI
 *      Laboratory: Image Pervasive Access Lab
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <err.h>
#include <sys/types.h>
#include <openmotion/om.h>
#include "utils_.h"

/////////////////////////////
//   OPENMOTION EXAMPLE    //
/////////////////////////////




/*
 * main function
 */
int main(int argc,char** argv){

	//variables file I/O
	char csvfile_input[100];
	struct line_reader lr;
	FILE *f;
	size_t len;
	char *line;

	// variable performance
	struct timeval tbegin_tot,tend_tot;
	struct Groundtruth gt;
	int index = 0;
	double mean_error = 0.0;
	double mean_time = 0.0;

	// variable for data fusion
	omSensorFusionManager manager;
	omNonLinearFilter_USQUE filter;
	
	// get the input csv file
	strcpy(csvfile_input,argv[2]);
    
    
    
    
    /***** Initialisation of the parameters *****/
    Init_Params I_P;
    // time gap between 2 sensor recording (here the sensor rate is 33hz)
    I_P.DELTA_T = 0.03;// time between 2 samples of sensor recording (here the sensor rate is 33hz)
    
    
    // this value depends of your location on earth
    // magnetic field in singapore in micro Tesla ( according to https://www.ngdc.noaa.gov/geomag-web/ )
    I_P.local_magn[0]=40.8101;I_P.local_magn[1]=0.1609;I_P.local_magn[2]=10.2636;
    
    // these values depends of your hardware: variance and bias
    I_P.var_accel=0.5;
    I_P.var_magn=0.15;
    I_P.var_gyro=0.00031623;
    
    
    I_P.bias_accel[0]=0.0;I_P.bias_accel[1]=0.0;I_P.bias_accel[2]=0.0;
    I_P.bias_magn[0]=0.0;I_P.bias_magn[1]=0.0;I_P.bias_magn[2]=0.0;
    I_P.bias_gyro[0]=0.000031623;I_P.bias_gyro[1]=0.0000316230;I_P.bias_gyro[2]=0.00003162;
    
    
    // this value depends of your initial orientation
    I_P.init_quat[0]=1.0;I_P.init_quat[1]=0.0;I_P.init_quat[2]=0.0;I_P.init_quat[3]=0.0;
    
    /***** End Initialisation of the parameters *****/

	// initialization of the sensor fusion manager and the nonlinear filter
	init_manager(&manager,&filter,USQUE, &I_P);
	
	
	// open the input file
	f = fopen(csvfile_input, "r");
	if (f == NULL) {
		perror(csvfile_input);
		exit(1);
	}

	// initialization of the reader
	lr_init(&lr, f);

	// allocation
	om_vector_create(&gt.position,3);

	// Start timer
	gettimeofday(&tbegin_tot,NULL);

	// for each line of the csv file
    // index is beginning at 1
    while ((line = next_line(&lr, &len))) {

		// get the values
		char** tokens;
	    tokens = str_split(line, ',');
        
	   
	    if ((index>0) && tokens)//first line with the text
	    {

            // get ground truth
            om_vector_setValues(&gt.position,3,atof(tokens[1]),atof(tokens[2]),atof(tokens[3]));
            om_quat_create(&gt.q_true,atof(tokens[4]),atof(tokens[5]),atof(tokens[6]),atof(tokens[7]));

            // get magnetometer values from csv file
            double accX = atof(tokens[11]);
            double accY = atof(tokens[12]);
            double accZ = atof(tokens[13]);

            // set accelerometer data
            om_vector_setValues(&manager.imu_data.data_accelerometer,3,(double)accX,(double)accY,(double)accZ);

            // get magnetometer values from csv file
            double magX = atof(tokens[14]);
            double magY = atof(tokens[15]);
            double magZ = atof(tokens[16]);

            // set magnetometer data
            om_vector_setValues(&manager.imu_data.data_magnetometer,3,(double)magX,(double)magY,(double)magZ);
                
                
            // get gyroscope values from csv file
            double gyroX = atof(tokens[8]);
            double gyroY = atof(tokens[9]);
            double gyroZ = atof(tokens[10]);
            
            // send values to the sensor fusion manager
            om_vector_setValues(&manager.imu_data.data_gyroscope,3,(double)gyroX, (double)gyroY, (double)gyroZ);

            // Chronometer variables
            double texec=0.0;
            struct timeval tbegin,tend;

            // Start timer
            gettimeofday(&tbegin,NULL);

            // process filtering
            manager.process_filter(&manager,&filter);
        
            // Stop timer
            gettimeofday(&tend,NULL);

            // Compute execution time
            texec=((double)(1000.0*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000.0)))/1000.0;

            // Compute error attitude
            double error = calculErrorOrientation(&gt.q_true,&manager.output.quaternion);
                
            // Compute mean error and mean time
            mean_error += error;
            mean_time += texec;

	    }

	    // free tokens
        for (int i = 0; *(tokens + i); i++){
	    	free(*(tokens + i));
        }
	    free(tokens);

	    // increase the index
	    index++;

	    // show the advancements
	    displayLoadingBar(index);

	}

	// compute performances
	mean_error /= (double)(index +1);
	mean_time /= (double)(index +1);

	// Stop timer
	gettimeofday(&tend_tot,NULL);

	// Compute execution time
	double texec_tot=((double)(1000.0*(tend_tot.tv_sec-tbegin_tot.tv_sec)+((tend_tot.tv_usec-tbegin_tot.tv_usec)/1000.0)))/1000.0;

	// print results
	printf("\n");
	printf("DONE\n");
	printf("Total Time : %f sec\n\n",texec_tot);
	printf("PERFORMANCE METHOD:\n");
	printf("\t mean attitude error = %.*f\n",10,mean_error);
	printf("\t mean computation time per sample = %.*f\n\n",10,mean_time);


	if (!feof(f)) {
		perror("next_line");
		exit(1);
	}

	//free file reader
	lr_free(&lr);

	// end program
	return EXIT_SUCCESS;
}






