//Config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <string>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
// Definitions 
extern const double MUT_RATE; // = 0.0000001;//0.0002002
extern const double DEATH_RATE; //= 0.02;
extern const double PROLIFERATION_RATE; //= 0.03;

extern const unsigned long long MAXIMUM_POPULATION_SIZE; //=   4500000000; //1000000000;
extern const unsigned long long MAXIMUM_POPULATION_SIZE_TO_STOP ;;//=   MAXIMUM_POPULATION_SIZE + 2000000000; //1000000000;
extern const double  PS ; //= (double) MAXIMUM_POPULATION_SIZE;


extern const std::map<std::string, bool> CLONE_VARIABLES;
extern gsl_rng *r_global ;	
extern const size_t K; //= 4;

extern const unsigned long long DETECTABLE_POPULATION_SIZE ;//= 4500000000; //1000000000;
extern const unsigned int STOP_GROWTH_AFTER_DIAGNOSIS ;//= 1; //[hours]

extern const bool SAVE_CLONAL_EVOLUTION ;//= false;

extern const bool REMOVE_DEATH_CLONES ;//= false;

extern const double DIFF ;//= PROLIFERATION_RATE - DEATH_RATE ;
extern const unsigned int dt ;//= 3600;

extern const bool MULTIPLE_TUMOURS ;//= false;
extern const int PENALTY; //= 3;

extern const double KILLER_PROBABILITY ;//= 0.025;
extern const double DRIVER_PROBABILITY ;//= 0.025;
extern const double PASSENGER_PROBAILIBITY ;//= 0.70;
extern const double DELETERIOUS_PROBABILITY ;//= 0.20;
extern const double BENEFICIAL_PROBABILITY ;//= 0.050;

extern const bool MULTIPLE_MUTATIONS;// = true;

extern const int MUTATION_MODEL; //= 1;

extern const unsigned long long MAX_CLONE_SIZE;

  const unsigned int P_STAYING  = 3;
  const unsigned int P_DYING  = 4;
  const unsigned int STAYING_CELLS   = 0;
  const unsigned int DYING_CELLS     = 1;
  const unsigned int EXITING_CELLS   = 2;


#endif