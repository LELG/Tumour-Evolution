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

Random::Random()
{}

Random::Random(std::mt19937::result_type seed) : eng(seed) {}

Random::~Random() 
{}

unsigned int Random::Recovery_After_Replication()
{
    return std::uniform_int_distribution<unsigned int>{12, 15}(eng);
}

