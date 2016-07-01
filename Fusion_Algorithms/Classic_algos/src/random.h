
/*
 * random.h
 *
 *  Created on: 17 May, 2016
 *      Author: thomas
 */

#ifndef RANDOM_H_
#define RANDOM_H_


#include "algebra.h"


void om_random_generateWhiteNoise(int n,double mu,double sigma,double seed,struct omVector* out);

void om_random_generateWhiteNoise(int n,double mu,double sigma,double seed,struct omVector* out);
/**
 * Gaussian white noise generator
 *
 * @param mu : minimal value
 * @param max : minimal value
 */
double om_random_normalDistribution(double mu,double sigma,double seed);

/**
 * Red (brownian) noise generator
 *
 * @return a number in the interval min and max
 * @param mu : minimal value
 * @param max : minimal value
 */
double om_random_brownianMotion(double mean, double a,double seed);

/**
 * Generate random value with bernouilli distribution
 * @param p: a real between 0 and 1
 */
double om_random_bernouilliDistribution(double p,double seed);

/**
 * Generate random value with weibull distribution
 *
 */
double om_random_weibullDistribution(double a, double lambda,double seed);

/**
 * Generate random value with uniform distribution
 */
double om_random_uniformDistribution(double seed);

/**
 * Generate random value with geometric distribution
 */
double om_random_geometricDistribution(double p,double seed);

/**
 * Generate random value with gamma distribution
 */
double om_random_gammaDistribution(double alpha,double beta,double seed);

/**
 * Generate random value with exponential distribution
 */
double om_random_exponentialDistribution(double lambda,double seed);

/**
 * Generate random value with poisson distribution
 *
 */
double om_random_poissonDistribution(double lambda, int seed);



#endif /* RANDOM_H_ */
