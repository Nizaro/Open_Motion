/*
 * random.c
 *
 *  Created on: 17 May, 2016
 *      Author: thomas
 */


#include "random.h"



void om_random_generateWhiteNoise(int n,double mu,double sigma,double seed,struct omVector *out){

	om_vector_create(out,n);

	for(int i=0;i<n;i++){

		double normal_value = om_random_normalDistribution(mu,sigma,seed+((double)(i)*100.0));
		om_vector_setValue(out,i,normal_value);

	}

}

void om_random_generateWhiteNoiseFromCovarianceMatrix(double mu,struct omMatrix *cov_L,double seed,struct omVector *out){

	int n = cov_L->_rows;

	if(!out->_length)
		om_vector_create(out,n);

	omVector z;
	omVector mu_vec;

	om_random_generateWhiteNoise(n,0.0,1.0,seed,&z);
	om_vector_create(&mu_vec,n);

	for(int i=0;i<n;i++)
		om_vector_setValue(&mu_vec,i,mu);

	om_operator_matrix_vector_mul(cov_L,&z,out);
	om_operator_vector_add(&mu_vec,out,out);

	om_vector_free(&z);
	om_vector_free(&mu_vec);


}




double om_random_weibullDistribution(double a, double lambda,double seed){

	srand(seed);

	double u = (double)(rand())/(double)(RAND_MAX);
	double x = pow(((-1.0)*log(u)),(1.0/a))  / (lambda);

	return x;
}


double om_random_poissonDistribution(double lambda, int seed){

    srand(seed);

    double k=0.0;
    double L = exp(-lambda);
    double p = 1.0;

    do{
    	k++;
    	//double u = uniformDistribution(seed);
    	double u = (double)(rand())/(double)(RAND_MAX);

    	p *= u;
    }while(p > L);

    return k;
}



/* Generate random value with exponential distribution */
double om_random_exponentialDistribution(double lambda,double seed){

	srand(seed);

	double u = (double)(rand())/(double)(RAND_MAX);
	double x = -log(u)/(lambda);

	return x;

}

/* Generate random value with geometric distribution */
double om_random_geometricDistribution(double p,double seed){

	srand(seed);

	double u = (double)(rand())/(double)(RAND_MAX);
	double x = floor(log(u)/log(1.0-p));

	return x;
}

/* Generate random value with bernouilli distribution */
double om_random_bernouilliDistribution(double p,double seed){

	srand(seed);
	double u = (double)(rand())/(double)(RAND_MAX);

	double x;

	if(u <= p){
		x=1.0;
	}else{
		x=0.0;
	}

	return x;
}

/* Generate random value with uniform distribution */
double om_random_uniformDistribution(double seed){

	srand(seed);
	double u = (double)(rand()) /((double)(RAND_MAX));
	return u;
}




double om_random_gammaDistribution(double alpha,double beta,double seed){

	srand(seed);
	double x;

	if(alpha > 1.0){

		int stop=0;

		do{

			double u1 = (double)(rand())/(double)(RAND_MAX);
			double u2 = (double)(rand())/(double)(RAND_MAX);

			double v = ( u1*(alpha - (1.0/(6.0*alpha))))/(u2*(alpha-1.0));

			x = (alpha-1.0)*v;

			stop |= (  (2.0*(u2-1.0) /(alpha-1.0) ) + v + (1.0/v)  ) <= 2.0;
			stop |= (  (2.0*log10(u2)/(alpha-1.0) ) - log10(v) + v  ) <= 1.0;

		}while(!stop);


	}else{

		double t = 0.07 + (0.75*sqrt(1.0-alpha));
		double b= 1.0 + (exp((-1.0)*t)*alpha/t);
		int stop=0;

		do{

			double u1 = (double)(rand())/(double)(RAND_MAX);
			double u2 = (double)(rand())/(double)(RAND_MAX);

			double v = b*u1;

			if(v<=1.0){
				x = t*pow(v,1.0/alpha);
				stop |= u2 <= (2.0 - x)/(2.0+x);
				stop |= u2 <= exp(x*(-1.0));
			}else{
				x = -log10((t*(b-v)/alpha));
				double y=x/t;
				stop |= (u2*(alpha + y*(1.0-alpha))) <= 1.0;
				stop |= u2 <= pow(y,alpha-1.0);
			}


		}while(!stop);



	}


	return x*beta;

}


double om_random_normalDistribution(double mean,double variance,double seed){


	const double epsilon = -DBL_MAX;

	srand(seed);

	static double z0, z1;

	double u1, u2;
	do{
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 =  rand() * (1.0 / RAND_MAX);
	}while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(2.0*PI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(2.0*PI * u2);

	return (z0 * variance) + mean;


}


double om_random_brownianMotion(double mean, double a,double seed){

	  double z1, z2, r, d;
	  double x = 0.0;

	  srand(seed);

	  do{
	      z1 = 2.0*( rand() * (1.0 / RAND_MAX)) - 1.0;
	      z2 = 2.0*( rand() * (1.0 / RAND_MAX)) - 1.0;
	      r = (z1 * z1) + (z2 * z2);
	  }while(r >= 1.0);

	  d = a * sqrt(-2.0 * log(r) / r);
	  x += mean + (d * z1);

	  return x;

}



