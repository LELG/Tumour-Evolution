/*
	Thisis a header file for the Clone Structure
*/

#ifndef RANDOM_H
#define RANDOM_H

#include "config.h"
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
#include <math.h>
#include <iostream>         // std::cout,std::endl
#include <vector>
#include <fstream>

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
    unsigned int G1();
    unsigned int S();
    unsigned int G2();
    unsigned int Poisson();
    double Z();
    double Binomial_dying(unsigned int Clone_Size, double Death_Rate);
    double Binomial_newborn(unsigned int Clone_Size, double Adjusted_Proliferation_Rate);
    unsigned int Binomial_Mutants(unsigned int NewBorn_Size, double Mutation_Rate);

    double Uniform_Mutation_Rate(double mu_rate);
    double Update_Proliferation_Rate(double Proliferation_Rate);
    double Uniform_Mutation_Rate_2(double Parent_mu_rate);
    unsigned int * Newborn_G0_and_G1(unsigned long long int &Newborn_cells, double & p_enter_G0);
    unsigned int * Update_G0_Phase(unsigned int & G0_cells, double & p_staying_G0, double & p_dying_in_G0, double & p_exiting_G0);
    unsigned int * Update_G1_Phase(unsigned int & G1_cells, double & p_staying_G1, double & p_dying_in_G1, double & p_exiting_G1);
    unsigned int * Update_G2_Phase(unsigned int & G2_cells, double & p_staying_G2, double & p_dying_in_G2, double & p_exiting_G2);
    unsigned int * Update_S_Phase(unsigned int & S_cells, double & p_staying_S, double & p_dying_S, double & p_exiting_S);
    unsigned int * Update_M_Phase(unsigned int & M_cells, double & p_staying_M, double & p_dying_M, double & p_exiting_M);
    unsigned int * Clonal_Functions(unsigned int & Available_cells, double & p_of_going_G0, double & p_entering_mitosis, double & p_dying );
    unsigned int * Newborn( unsigned int & Newborn_cells, double & p_idle, double & p_mutate, double & p_go_to_G0 );
    unsigned int * Mutational_Effects();
    void Valid_Limits(double & effect);
    void Laplace(double & effect);
    void Laplace(double & effect, unsigned int & Clone_size);
    unsigned int * Mutational_Proportions( double & p_go_to_G0, double & p_go_to_G1, unsigned int & Mutant_Cells );

private:        
    std::mt19937 eng{std::random_device{}()};


};

#endif