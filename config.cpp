#include "config.h"
#include <map>
#include <string>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <limits>
#include <random>


//-1 MR for testing
const double MUT_RATE           = 0.00000000001; //0.0000001
const double DEATH_RATE         = 0.001;
const double PROLIFERATION_RATE = 0.002;

const unsigned long long MAXIMUM_POPULATION_SIZE =  17088669 ; //100000; //1000000000;
const unsigned long long MAXIMUM_POPULATION_SIZE_TO_STOP =   MAXIMUM_POPULATION_SIZE + 1000; // + 2000000000;
const unsigned int STOP_AFTER_DIAGNOSIS_COUNTER = 100;
const double  PS = (double) MAXIMUM_POPULATION_SIZE;


const unsigned long long DETECTABLE_POPULATION_SIZE = 100000; //1000000000;
const unsigned int STOP_GROWTH_AFTER_DIAGNOSIS = 1; //[hours]

const bool SAVE_CLONAL_EVOLUTION = false;

const bool REMOVE_DEATH_CLONES = false;

const unsigned long long MAX_CLONE_SIZE = (unsigned long long) std::numeric_limits<unsigned int>::max();


const double P_DRUG_RESISTANCE = 0.0000001;

// We can overide thios in a config file, this are default
const std::map<std::string, bool> CLONE_VARIABLES = {	{"Mutation_Rate", true}, \
														{"Proliferation_Rate", true}, \
														{"Death_Rate", true}, \
														{"Number_of_Memebers_to_Start_Heterogeneity", true},\
														{"Generation_ID_Counter", true}, \
														{"Driver_10_fold_accumulation", true}, \
														{"clone_extinct", true}, \
														{"Number_of_Mutations", true}, \
														{"Clone_Size", true}, \
														{"Initiall_Expasion_Period", true}, \
														{"In_G0_Phase", true}, \
														{"In_G1_Phase", true}, \
														{"In_S_Phase", true}, \
														{"In_G2_Phase", true}, \
														{"In_M_Phase", true}, \
														{"Remaining_Time_in_G1_Phase", true}, \
														{"Remaining_Time_in_S_Phase", true}, \
														{"Remaining_Time_in_G2_Phase", true}, \
														{"Remaining_Time_in_M_Phase", true}, \
														{"init_PR", true}, \
														{"max_PR", true}, \
														{"final_PR", true}, \
														{"Generation_ID", true}\
													};	

const size_t K = 4;	

const double DIFF = PROLIFERATION_RATE - DEATH_RATE ;
const unsigned int dt = 3600;

const bool MULTIPLE_TUMOURS = false;
const int PENALTY = 1;

const double KILLER_PROBABILITY = 0.025;
const double DRIVER_PROBABILITY = 0.025;
const double PASSENGER_PROBAILIBITY = 0.80;
const double DELETERIOUS_PROBABILITY = 0.10;
const double BENEFICIAL_PROBABILITY = 0.050;

const bool MULTIPLE_MUTATIONS = true;

const int MUTATION_MODEL = 1;

gsl_rng *r_global ;		



 // const unsigned int P_STAYING = 3;
 // const unsigned int P_DYING   = 4;
 // const unsigned int STAYING_CELLS = 0;
 // const unsigned int DYING_CELLS = 1;
 // const unsigned int EXITING_CELLS = 2;									

// const gsl_rng_type * T2;
// gsl_rng * r2;
// srand(time(NULL)); // srand(time(NULL));
// //unsigned int Seed2 = 1234567; // rand();
// gsl_rng_env_setup();

// T2 = gsl_rng_default;
// r2 = gsl_rng_alloc (T2);
// gsl_rng_set (r2, 1234567);								

