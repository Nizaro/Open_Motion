/*
 * random.cpp
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */


#include "random.h"


namespace om
{

///////////////////////////////////////////////////////
/////          Random number generators           /////
///////////////////////////////////////////////////////

Vector whiteNoiseVector(int n,double mu,double sigma,double seed){


	Vector noise(static_cast<double>(n));

	for(int i=0;i<n;i++){
		double normal_value = whiteNoise(mu,sigma,seed+(static_cast<double>(i)*100.0));
		noise.setValue(i,normal_value);
	}

	return noise;
}


Vector whiteNoiseVector(double mu,const Matrix& covariance_L,double seed){

	Vector noise(static_cast<double>(covariance_L.getRows()));

	Vector z = whiteNoiseVector(covariance_L.getRows(),0.0,1.0,seed);

	Vector mu_vec(static_cast<double>(covariance_L.getRows()));
	for(int i=0;i<covariance_L.getRows();i++)
		mu_vec.setValue(i,mu);

	noise = mu_vec + covariance_L*z;

	return noise;
}



double weibullDistribution(double a, double lambda,double seed){

	srand(seed);

	double u = double(rand())/double(RAND_MAX);
	double x = pow(((-1.0)*log(u)),(1.0/a))  / (lambda);

	return x;


}


double poissonDistribution(double lambda, int seed){

    srand(seed);

    double k=0.0;
    double L = exp(-lambda);
    double p = 1.0;

    do{
    	k++;
    	//double u = uniformDistribution(seed);
    	double u = double(rand())/double(RAND_MAX);

    	p *= u;
    }while(p > L);

    return k;
}



/* Generate random value with exponential distribution */
double exponentialDistribution(double lambda,double seed){

	srand(seed);

	double u = double(rand())/double(RAND_MAX);
	double x = -log(u)/(lambda);

	return x;

}

/* Generate random value with geometric distribution */
double geometricDistribution(double p,double seed){

	srand(seed);

	double u = double(rand())/double(RAND_MAX);
	double x = floor(log(u)/log(1.0-p));

	return x;
}

/* Generate random value with bernouilli distribution */
double bernouilliDistribution(double p,double seed){

	srand(seed);
	double u = double(rand())/double(RAND_MAX);

	double x;

	if(u <= p){
		x=1.0;
	}else{
		x=0.0;
	}

	return x;
}

/* Generate random value with uniform distribution */
double uniformDistribution(double seed){

	srand(seed);
	double u = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX));
	return u;
}


/* Random number generator between 2 values (min and max). It's base of the function rand() from c */
double randomDouble(double min,double max){
	return min + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(max-min)));
}


double gammaDistribution(double alpha,double beta,double seed){

	srand(seed);
	double x;

	if(alpha > 1.0){

		bool stop=false;

		do{

			double u1 = double(rand())/double(RAND_MAX);
			double u2 = double(rand())/double(RAND_MAX);

			double v = ( u1*(alpha - (1.0/(6.0*alpha))))/(u2*(alpha-1.0));

			x = (alpha-1.0)*v;

			stop |= (  (2.0*(u2-1.0) /(alpha-1.0) ) + v + (1.0/v)  ) <= 2.0;
			stop |= (  (2.0*log10(u2)/(alpha-1.0) ) - log10(v) + v  ) <= 1.0;

		}while(!stop);


	}else{

		double t = 0.07 + (0.75*sqrt(1.0-alpha));
		double b= 1.0 + (exp((-1.0)*t)*alpha/t);
		bool stop=false;

		do{

			double u1 = double(rand())/double(RAND_MAX);
			double u2 = double(rand())/double(RAND_MAX);

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


double whiteNoise(double mean,double variance,double seed){


	const double epsilon = std::numeric_limits<double>::min();

	srand(seed);

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * variance + mean;

	double u1, u2;
	do{
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 =  rand() * (1.0 / RAND_MAX);
	}while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(2.0*PI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(2.0*PI * u2);

	return (z0 * variance) + mean;


}


double brownianNoise(double mean, double a,double seed){

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




void poisson(bool** pos_imp,int N,double lambda, int newseed){
    // lambda is the mean of the occurence of an impulse in a poisson distribution
    // lamdba=10 meaning that 10 occurrences are needed to have an impulse in average
    int i;

    //initialisation
    for (i=0;i<N;i++){
        (*pos_imp)[i]=false;
    }

    srand(newseed);
    int pos=0;
    double unif;

    for (i=0;i<N;i++){
        unif=double(rand())/double(RAND_MAX);

        pos = pos + (int)(-log(unif)*lambda);

        if (pos>=N) break;

        (*pos_imp)[pos]=true;
        i=pos;



    }


}





}/* namespace isf */
