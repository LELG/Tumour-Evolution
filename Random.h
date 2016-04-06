/*
	Thisis a header file for the Clone Structure
*/

#ifndef RANDOM_H
#define RANDOM_H

#include "Random.h"
#include <string>
#include <map>
#include <memory>           // std::unique_ptr
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <random> 			// std::generator
#include <tuple>
#include <ctime>
#include <unistd.h>
#include <time.h>
#include <ctime>
#include <random>
#include <cmath>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>

class Random
{

	/**
	    Variable Definitions
	**************************/
public:
	Random();
    Random(std::mt19937::result_type seed);
    ~Random();
    unsigned int Recovery_After_Replication();

private:        
    std::mt19937 eng{std::random_device{}()};


};

#endif