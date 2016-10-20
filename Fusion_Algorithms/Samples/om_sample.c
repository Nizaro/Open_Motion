/*
 * om_example.c
 *
 *  Created on: 1 Aug, 2016
 *      Author: Thomas BRAUD, Nizar OUARTI, Vivien BILLAUD
 *      Company: Image Pervasive Access Lab
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <err.h>
#include <sys/types.h>
#include <openmotion/om.h>



/////////////////////////////
//   FILE READER Manager   //
/////////////////////////////



/*
 * File reader
 */
struct line_reader {

	// All members are private.
	FILE	*f;
	char	*buf;
	size_t	 siz;

};

/*
 * Initializes a line reader _lr_ for the stream _f_.
 */
void lr_init(struct line_reader *lr, FILE *f){
	lr->f = f;
	lr->buf = NULL;
	lr->siz = 0;
}

/*
 * Reads the next line. If successful, returns a pointer to the line,
 * and sets *len to the number of characters, at least 1. The result is
 * _not_ a C string; it has no terminating '\0'. The returned pointer
 * remains valid until the next call to next_line() or lr_free() with
 * the same _lr_.
 *
 * next_line() returns NULL at end of file, or if there is an error (on
 * the stream, or with memory allocation).
 */
char* next_line(struct line_reader *lr, size_t *len){

	size_t newsiz;
	int c;
	char *newbuf;

	*len = 0;			/* Start with empty line. */
	for (;;) {
		c = fgetc(lr->f);	/* Read next character. */
		if (ferror(lr->f))
			return NULL;

		if (c == EOF) {
			/*
			 * End of file is also end of last line,
		`	 * unless this last line would be empty.
			 */
			if (*len == 0)
				return NULL;
			else
				return lr->buf;
		} else {
			/* Append c to the buffer. */
			if (*len == lr->siz) {
				/* Need a bigger buffer! */
				newsiz = lr->siz + 4096;
				newbuf = realloc(lr->buf, newsiz);
				if (newbuf == NULL)
					return NULL;
				lr->buf = newbuf;
				lr->siz = newsiz;
			}
			lr->buf[(*len)++] = c;

			/* '\n' is end of line. */
			if (c == '\n')
				return lr->buf;
		}
	}
}

/*
 * Frees internal memory used by _lr_.
 */
void lr_free(struct line_reader *lr){
	free(lr->buf);
	lr->buf = NULL;
	lr->siz = 0;
}



/*
 * Split a string into a array of string according to a delimitor
 */
char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    // Count how many elements will be extracted.
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    // Add space for trailing token.
    count += last_comma < (a_str + strlen(a_str) - 1);

    // Add space for terminating null string so caller
    // knows where the list of returned strings ends.
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}


/////////////////////////////
//   OPENMOTION EXAMPLE    //
/////////////////////////////


/*
 * The ground truth composed by the true position in a 3D scene and the true orientation represented by a quaternion
 */
struct Groundtruth{

	omVector position;
	omQuaternion q_true;

};


/*
 * allows to display a loading bar
 */
void displayLoadingBar(int index){

	double percentage = ((double)index / 60000.0) * 100;

   printf("LOADING %.2f%%\r", percentage);
   fflush(stdout);
}

/*
 * give the sign of a real. It will return -1.0 if the number is negative
 */
double signeOf(double d){
	return d > 0.0 ? 1.0 : -1.0;
}


/*
 * compute the rms attitude error between two quaternion
 */
double calculErrorOrientation(omQuaternion *q_real,omQuaternion *q_est){

	// some variables
	omQuaternion q_est_inv;
	omQuaternion dq;
	omVector dp;

	// compute dq = q_real*q_est.inverse();
	om_quat_inverse(q_est,&q_est_inv);
	om_operator_quat_mul(q_real,&q_est_inv,&dq);

	// convert dq into modified rodriguez vector representation
	double tmp = (4.0 * (signeOf(dq._qw)))/(1.0 + fabs(dq._qw));
	om_quat_imaginary(&dq,&dp);
	om_operator_vector_scal_mul(&dp,tmp,&dp);

	// return the rms attitude error
	return om_vector_rms(&dp);

}



/*
 * Transform the magnetic field into an estimate of the geographic north
 */
void modification_magneto(struct omVector *input,struct omVector* output,struct omQuaternion *q_est){

	// variables
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

	// represents the declinaison and the inclinaison of the magnetic field.
	// this values works in Singapore
	// go see http://www.ngdc.noaa.gov/geomag-web/ for more details
	aa_y._angle = -14.36*DEG_TO_RAD;
	aa_z._angle = 0.23*DEG_TO_RAD;

	//allocation
	om_vector_create(&aa_y._axis,3,0.0,1.0,0.0);
	om_vector_create(&aa_z._axis,3,0.0,0.0,1.0);
	om_matrix_create(&M,3,3);
	om_convert_quaternion2matrix(q_est,&M);
	om_vector_create(&tmp,3,0.0,0.0,0.0);
	om_vector_create(&tmp2,3,0.0,0.0,0.0);
	om_vector_create(&tmp3,3,0.0,0.0,0.0);
	om_vector_create(output,3,0.0,0.0,0.0);

	// convert the declinaison and inclinaison into quaternion
	om_convert_axisAngle2quaternion(&aa_y,&q_y);
	om_convert_axisAngle2quaternion(&aa_z,&q_z);

	// get the inverse
	om_quat_inverse(&q_y,&q_y_inv);
	om_quat_inverse(&q_z,&q_z_inv);

	// transform the magnetic field defined in the reference frame relative to the IMU into the North-East-Down(NED) frame
	om_operator_matrix_vector_mul(&M,input,&tmp);

	// correct the magnetic field to get the geographic north defined in the NED frame
	om_rotate_vector_quaternion(&q_z_inv,&tmp,&tmp2);
	om_rotate_vector_quaternion(&q_y_inv,&tmp2,&tmp3);

	// transform the geographic north defined in in the NED frame into the reference frame relative to the IMU
	om_rotate_vector_quaternion(q_est,&tmp3,output);

	// free memory
	om_matrix_free(&M);
	om_vector_free(&tmp);
	om_vector_free(&tmp2);
	om_vector_free(&tmp3);
	om_vector_free(&aa_y._axis);
	om_vector_free(&aa_z._axis);


}



/*
 * initialization of the manager components and the nonlinear filter
 */
void init_manager(omSensorFusionManager *manager,void* filter){

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

    // initialization of the nonlinear filter
    manager->initialization_filter(manager,(omNonLinearFilter_USQUE*)filter);
}


/*
 * main function
 */
int main(int argc,char** argv){

	//variables
	char csvfile_input[100];
	struct line_reader lr;
	FILE *f;
	size_t len;
	char *line;
	struct timeval tbegin_tot,tend_tot;
	struct Groundtruth gt;
	int index = 0;
	double mean_error = 0.0;
	double mean_time = 0.0;
	omSensorFusionManager manager;
	omNonLinearFilter_USQUE filter_usque;

	// get the input csv file
	strcpy(csvfile_input,argv[2]);

	// initialization of the sensor fusion manager and the nonlinear filter
	init_manager(&manager,&filter_usque);

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
	while ( (line = next_line(&lr, &len))) {

		// get the values
		char** tokens;
	    tokens = str_split(line, ',');

	    // if t=0 we set the gyroscope values
	    if(index == 1){

	    	// get gyroscope values from csv file
	    	double gyroX = atof(tokens[8]);
	    	double gyroY = atof(tokens[9]);
	    	double gyroZ = atof(tokens[10]);

	    	// send values to the sensor fusion manager
	        om_vector_setValues(&manager.imu_data.data_gyroscope,3,(double)gyroX, (double)gyroY, (double)gyroZ);

	    }
	    // for the rest
	    else if (index > 1 && tokens)
	    {

	    	// get ground truth
	    	om_vector_setValues(&gt.position,3,atof(tokens[1]),atof(tokens[2]),atof(tokens[3]));
	    	om_quat_create(&gt.q_true,atof(tokens[4]),atof(tokens[5]),atof(tokens[6]),atof(tokens[7]));

	    	// get magnetometer values from csv file
	    	double accX = atof(tokens[11]);
	    	double accY = atof(tokens[12]);
	    	double accZ = atof(tokens[13]);

	    	// set accelerometer data
	    	om_vector_setValues(&manager.imu_data.data_accelerometer,3,accX,accY,accZ);

	    	// get magnetometer values from csv file
	    	double magX = atof(tokens[14]);
	    	double magY = atof(tokens[15]);
	    	double magZ = atof(tokens[16]);

			// set magnetometer data
			omVector mag;
			om_vector_create(&mag,3,magX,magY,magZ);


			// rectification magnetometer
			modification_magneto(&mag,&manager.imu_data.data_magnetometer,&manager.output.quaternion);

			// Chronometer variables
			double texec=0.0;
			struct timeval tbegin,tend;

			// Start timer
			gettimeofday(&tbegin,NULL);

			// process filtering
			manager.process_filter(&manager,&filter_usque);

			// Stop timer
			gettimeofday(&tend,NULL);


			// Compute execution time
			texec=((double)(1000.0*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000.0)))/1000.0;

			// Compute error attitude
			double error = calculErrorOrientation(&gt.q_true,&manager.output.quaternion);

			// Update
			mean_error += error;
			mean_time += texec;

		    	// get gyroscope values from csv file
		    	double gyroX = atof(tokens[8]);
	    		double gyroY = atof(tokens[9]);
		    	double gyroZ = atof(tokens[10]);


		    	// send values to the sensor fusion manager
		        om_vector_setValues(&manager.imu_data.data_gyroscope,3,(double)gyroX, (double)gyroY, (double)gyroZ);

	    }

	    // free tokens
	    for (int i = 0; *(tokens + i); i++)
	    	free(*(tokens + i));
	    free(tokens);

	    // increase the index
	    index++;

	    // show the advancements
		displayLoadingBar(index);

	}

	// compute performances
	mean_error /= (double)(index - 2);
	mean_time /= (double)(index - 2);

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
	printf("\t mean computation time = %.*f\n\n",10,mean_time);


	if (!feof(f)) {
		perror("next_line");
		exit(1);
	}

	//free file reader
	lr_free(&lr);


}






