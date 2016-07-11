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
    void Multiple_Mutations_Poisson(unsigned int & number_of_mutations);
    double Z();
    double Binomial_dying(unsigned int Clone_Size, double Death_Rate);
    double Binomial_newborn(unsigned int Clone_Size, double Adjusted_Proliferation_Rate);
    unsigned int Binomial_Mutants(unsigned int NewBorn_Size, double Mutation_Rate);

    void Uniform_Gain_Update(const double & Parent_mu_rate, const double & gain, double & updated_MR);
    void Uniform_PR_Update(const double & Parent_Proliferation_Rate, const double & Variance, double & updated_PR);
    double Uniform_Mutation_Rate(double mu_rate);
    double Update_Proliferation_Rate(double Proliferation_Rate);
    void Update_Proliferation_Rate_V2(const double & Parent_Proliferation_Rate, double & new_PR );
    double Uniform_Mutation_Rate_2(double Parent_mu_rate);

    void Basic_Clonal_Expansion_Sampling_V1( const unsigned long long int & Clone_Size, const  double & p_dr, const  double & p_pr, const double & p_idle, std::vector<unsigned int> & NewBorn_Cells);

    void Basic_Clonal_Expansion_Sampling_V1_HPC(const unsigned long long int & Clone_Size, const double & p_dr, const double & p_pr, std::vector<unsigned long long int> & NewBorn_Cells);
    void Basic_Clonal_Expansion_Sampling_V2_HPC(const unsigned long long int & Clone_Size, const double & p_dr, const double & p_pr, const double & p_mutate, std::vector<unsigned long long int> & NewBorn_Cells);

    void Basic_Clonal_Expansion_Sampling_V2R_HPC(const unsigned long long int & Clone_Size, const double & p_dr, const double & p_pr, const double & p_mutate, std::vector<unsigned long long int> & NewBorn_Cells);

    void Binomial_Mutant( std::tuple<unsigned long long int, 
                                     unsigned long long int, 
                                     unsigned long long int, 
                                     bool, bool> & Dying_and_Newborn,
                                     const double & mutation_rate);

    unsigned int * Newborn_G0_and_G1(unsigned long long int &Newborn_cells, double & p_enter_G0);
    unsigned int * Update_G0_Phase(const unsigned int & G0_cells, const double & p_staying_G0, const double & p_dying_in_G0, const double & p_exiting_G0);
    unsigned int * Update_G1_Phase(const unsigned int & G1_cells, const double & p_staying_G1, const double & p_dying_in_G1, const double & p_exiting_G1);
    unsigned int * Update_G2_Phase(const unsigned int & G2_cells, const double & p_staying_G2, const double & p_dying_in_G2, const double & p_exiting_G2);
    unsigned int * Update_S_Phase(const unsigned int & S_cells, const  double & p_staying_S, const double & p_dying_S, const double & p_exiting_S);
    unsigned int * Update_M_Phase(const unsigned int & M_cells, const double & p_staying_M, const double & p_dying_M, const double & p_exiting_M);
    unsigned int * Clonal_Functions(const unsigned int & Available_cells, const double & p_of_going_G0, const double & p_entering_mitosis, const double & p_dying );
    void Newborn( const unsigned int & Newborn_cells, const  double & p_idle, const  double & p_mutate, const  double & p_go_to_G0, std::vector<unsigned int> & Division_Model);
    unsigned int * Mutational_Effects();
    void Valid_Limits(double & effect);
    void Laplace(double & effect);
    void Laplace(double & effect, unsigned int & Clone_size);
    void Mutational_Proportions( const double & p_go_to_G0,const  double & p_go_to_G1, const unsigned int & Mutant_Cells, std::vector<unsigned int> & Mutations );

    void Drug_Resistance_Strength( double & strenght );

private: 
    //std::mt19937 rng;
    //rng.seed(std::random_device()());       
    std::mt19937 eng{std::random_device{}()};


};

#endif