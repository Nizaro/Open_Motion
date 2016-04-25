/*
 * random.h
 *
 *  Created on: 16 Nov, 2015
 *      Author: thomas
 */

#ifndef RANDOM_H_
#define RANDOM_H_


#include "Quaternion.h"
#include "Matrix.h"


namespace om
{



///////////////////////////////////////////////////////
/////          Random number generators           /////
///////////////////////////////////////////////////////


Vector whiteNoiseVector(int n,double mu,double sigma,double seed);

Vector whiteNoiseVector(double mu,const Matrix& covariance_L,double seed);



/**
 * Random number generator between 2 values (min and max). It's base of the function rand() from c
 *
 * @return a number in the interval min and max
 * @param min : minimal value
 * @param max : minimal value
 */
double randomDouble(double min,double max);

/**
 * Gaussian white noise generator
 *
 * @param mu : minimal value
 * @param max : minimal value
 */
double whiteNoise(double mu,double sigma,double seed);

/**
 * Red (brownian) noise generator
 *
 * @return a number in the interval min and max
 * @param mu : minimal value
 * @param max : minimal value
 */
double brownianNoise(double mean, double a,double seed);

/**
 * Generate random value with bernouilli distribution
 * @param p: a real between 0 and 1
 */
double bernouilliDistribution(double p,double seed);

/**
 * Generate random value with uniform distribution
 *
 */
double weibullDistribution(double a, double lambda,double seed);


/**
 * Generate random value with uniform distribution
 *
 */
double uniformDistribution(double seed);

/**
 * Generate random value with uniform distribution
 *
 */
double geometricDistribution(double p,double seed);


double gammaDistribution(double alpha,double beta,double seed);

/**
 * Generate random value with exponential distribution
 *
 */
double exponentialDistribution(double lambda,double seed);

/**
 * Generate random value with poisson distribution
 *
 */
double poissonDistribution(double lambda, int seed);



void poisson(bool** pos_imp,int N,double lambda, int newseed);





} /* namespace sdf */




#endif /* RANDOM_H_ */
