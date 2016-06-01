#include "Random.h"
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
#include <fstream>
#include <vector>

Random::Random()
{}

Random::Random(std::mt19937::result_type seed) : eng(seed) {}

Random::~Random() 
{}

unsigned int Random::Recovery_After_Replication()
{
    return std::uniform_int_distribution<unsigned int>{12, 15}(eng);
}

unsigned int Random::G1()
{
	return std::uniform_int_distribution<unsigned int>{12, 15}(eng);
}

 unsigned int Random::S()
{
	return std::uniform_int_distribution<unsigned int>{5, 10}(eng);
}

unsigned int Random::G2()
{
	return std::uniform_int_distribution<unsigned int>{3, 5}(eng);
}

unsigned int Random::Poisson()
{
	return std::poisson_distribution<unsigned int>{1}(eng);
}

void Random::Multiple_Mutations_Poisson(unsigned int & number_of_mutations)
{
	number_of_mutations = std::poisson_distribution<unsigned int>{1}(eng);
	if(number_of_mutations == 0)
		number_of_mutations = 1;

}

double Random::Z()
{
	return std::normal_distribution<double>{0.0,1.0}(eng);
}

double Random::Binomial_dying(unsigned int Clone_Size, double Death_Rate)
{
	return std::binomial_distribution<unsigned int>{Clone_Size, Death_Rate}(eng);
}
	
double Random::Binomial_newborn(unsigned int Clone_Size, double Adjusted_Proliferation_Rate)
{
	return std::binomial_distribution<unsigned int>{Clone_Size, Adjusted_Proliferation_Rate }(eng);
}

unsigned int Random::Binomial_Mutants(unsigned int NewBorn_Size, double Mutation_Rate)
{
	return std::binomial_distribution<unsigned int>{NewBorn_Size, Mutation_Rate}(eng);
}

double Random::Uniform_Mutation_Rate(double mu_rate)
{
	return std::uniform_real_distribution<double>{mu_rate/2.0, mu_rate*2.0}(eng);
}

// double Random::Update_Proliferation_Rate(double Parent_Proliferation_Rate)
// {
// 	double U = normal_distribution<double>{Parent_Proliferation_Rate, 0.000001 }(eng);//0.000001
// 	if(U > 1.0)
// 	{
// 		U = 1.0;
// 	}
// 	if (U < 0.0)
// 	{
// 		U = 0.0;
// 	}
		
// 	return U;
// }

// 
// New Models
//
void Random::Uniform_PR_Update(const double & Parent_Proliferation_Rate, const double & Variance, double & updated_PR)
{
	double U_sample = std::normal_distribution<double>{Parent_Proliferation_Rate, Variance }(eng);
	if(U_sample > 1.0)
	{
		U_sample = 1.0;
	}
	if (U_sample < 0.0)
	{
		U_sample = 0.0;
	}
	updated_PR = U_sample;
}

void Random::Uniform_Gain_Update(const double & Parent_mu_rate, const double gain, double & updated_MR)
{
	double U = 0.0;
	U = fabs(std::normal_distribution<double>{Parent_mu_rate * gain, Parent_mu_rate/gain  }(eng)) ;
	updated_MR = U;

}

double Random::Update_Proliferation_Rate(double Parent_Proliferation_Rate)
{
	/*
	double Proliferation_Gain = 0.05;//gsl_ran_beta(r_global, 1.0, (double) Accumulated_Drivers );
	//cout << "B: " << Proliferation_Gain << " Par: " << Parent_Proliferation_Rate << " AC " << Accumulated_Drivers;
	Proliferation_Gain += Parent_Proliferation_Rate;
	//cout << " New: " << Proliferation_Gain << endl;
	//getchar();
		
	if(Proliferation_Gain > 1.0)
		Proliferation_Gain = 1.0;
	
	if (Proliferation_Gain < 0.0)
		Proliferation_Gain = 0.0;

	return Proliferation_Gain;
	*/

	//0.000001
	double U = std::normal_distribution<double>{Parent_Proliferation_Rate, 0.001 }(eng);//0.000001
	if(U > 1.0)
	{
		U = 1.0;
	}
	if (U < 0.0)
	{
		U = 0.0;
	}
		
	return U;

	
}

void Random::Update_Proliferation_Rate_V2(const double & Parent_Proliferation_Rate, double & new_PR )
{
	

	new_PR = Parent_Proliferation_Rate + std::normal_distribution<double>{0.01, 0.05 }(eng);
	if(new_PR >= 0.2)
		new_PR = 0.2;	

	if(new_PR < 0.0)
		new_PR = 0.0;




}
//TODO this is weird
double Random::Uniform_Mutation_Rate_2(double Parent_mu_rate)
{

	double U;// = normal_distribution<double>{Parent_mu_rate, mut_rate }(eng);
	//U = gsl_ran_beta (r_global, 1.0, 1000000);// 5000000
	U = fabs(std::normal_distribution<double>{Parent_mu_rate*0.25, Parent_mu_rate/0.25  }(eng)) ;
	//U = gsl_cdf_beta_P (Parent_mu_rate, 2.0,  2.0);
	//cout << " parent MR " << Parent_mu_rate << " Child  "<< U << endl;
	//getchar();
		
	return U;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////// HIGH PERFOMANCE DECLARATION OF FUNCTIONS //////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Random::Basic_Clonal_Expansion_Sampling_V1_HPC(const unsigned long long int & Clone_Size, const double & p_dr, const double & p_pr, std::vector<unsigned long long int> & NewBorn_Cells)
{
	double                     p_vector[] = { p_dr, p_pr, 1.0 - (p_dr + p_pr) };
	unsigned long long int cumulative_N   = 0;
	double                 cumulative_p   = 0.0;
    
	for(unsigned long long int i = 0; i < 2; i++)
	{
		if(p_vector[i] == 0.0)
		{
			NewBorn_Cells.push_back(0);
		}
		else
		{
			std::binomial_distribution<unsigned long long int> distribution( (Clone_Size - cumulative_N), p_vector[i] / (1.0 -cumulative_p) );
			NewBorn_Cells.push_back( distribution(eng) );
		}
		cumulative_N += NewBorn_Cells.at(i);
		cumulative_p += p_vector[i];
	
	}
	NewBorn_Cells.push_back( Clone_Size - cumulative_N );
}


void Random::Binomial_Mutant( std::tuple<unsigned long long int, 
										 unsigned long long int, 
										 unsigned long long int, 
										 bool, bool> & Dying_and_Newborn,
										 const double & mutation_rate)
{

    std::binomial_distribution<unsigned long long int> distribution( std::get<1>(Dying_and_Newborn), mutation_rate );
    std::get<2>(Dying_and_Newborn) = distribution(eng);

}


// For V1 smapling model//fix
void Random::Basic_Clonal_Expansion_Sampling_V1( const unsigned long long int & Clone_Size, const  double & p_dr, const  double & p_pr, const double & p_idle, std::vector<unsigned int> & NewBorn_Cells)
{


	NewBorn_Cells.push_back(0);
	NewBorn_Cells.push_back(0);
	NewBorn_Cells.push_back(0);

	static unsigned int buffer[3] = {0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = {p_dr, p_pr, 1.0 - (p_dr + p_pr)};
	const size_t K_S = 3;


	gsl_ran_multinomial ( r_global,
						  K_S,
						  static_cast<unsigned int>(Clone_Size),
						  p_vector,
						  buffer
						);

	NewBorn_Cells[0] = buffer[0];
	NewBorn_Cells[1] = buffer[1];
	NewBorn_Cells[2] = buffer[2];
	
}


/**
	Mitosis Stochastic Network
	I) Once Newborn cells are estimated this is the first step of the network

	[B_{G0,i}(t), B_{G1,i}(t)] ~ Mult(Bi(t),[p_{q,i}(t), 1 - p_{q,i}(t) ])
	Where:
	Bi(t) is the Newborn cells at time t
	B_{G0,i}(t) is the number of cells that enter G0 phase
	B_{G1,i}(t) is the number of cells that enter G1 phase

****************************************************************************/
unsigned int *  Random::Newborn_G0_and_G1( unsigned long long int &Newborn_cells, double & p_enter_G0)
{
	static unsigned int buffer[2] = {0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { p_enter_G0, 1.0 - p_enter_G0 };
	const size_t K_S = 2;	

	if(Newborn_cells > MAX_CLONE_SIZE)
		std::cout << "Clone is too Big... [OVERFLOW]" << std::endl; 
	else
		gsl_ran_multinomial ( r_global, // Random Seed Generator 
							  K_S, 		// Number of otput elements
							  static_cast<unsigned int>(Newborn_cells), //Bi(t)
							  p_vector, //[p_{q,i}, 1 - p_{q,i} ]
							  buffer
							);

	return buffer;
} 

/**
	Mitosis Stochastic Network
	II) Updating G0 Phase

	[G0{S,i}(t), G0{Dy,i}(t), G0{G0->.,i}(t)] ~ Mult(G0_i(t),[p^G0_{S,i}(t), p^G0_{Dy,i}(t), p^G0_{Ex,i}(t) ])
	Where:
	G0_i(t) is the number of cells in G0 phase at time t
	G0{S,i}(t) is the number of cells that stay in G0 Phase
	G0{Dy,i}(t) is the number of cells that die in G0 Phase
	G0{G0->.,i}(t) is the number of cells that are exiting G0 Phase

****************************************************************************/
unsigned int * Random::Update_G0_Phase(const unsigned int &G0_cells, const double & p_staying_G0, const double & p_dying_in_G0, const double & p_exiting_G0)
{
	static unsigned int buffer[3] = {0, 0, 0};
	double p_vector[] = {p_staying_G0, p_dying_in_G0, p_exiting_G0};
	const size_t K_S = 3;

	gsl_ran_multinomial ( r_global, // Random Seed Generator 
						  K_S, 		// Number of output elements
						  G0_cells, // G0_i(t)
						  p_vector, // [p^G0_{S,i}(t), p^G0_{Dy,i}(t), p^G0_{Ex,i}(t) ]
						  buffer
						);

	return buffer;
}

/**
	Mitosis Stochastic Network
	III) Updating G1 Phase

	[G1{S,i}(t), G1{G1->.,i}(t)] ~ Mult(G1_i(t),[p^G1_{S,i}(t), p^G1_{G1->G2}(t) ])
	Where:
	G1_i(t) is the number of cells in G1 phase at time t
	G1{S,i}(t) is the number of cells that stay in G1 Phase
	G1{G1->.,i}(t) is the number of cells exiting G1 phase

****************************************************************************/
unsigned int * Random::Update_G1_Phase(const unsigned int & G1_cells, const double & p_staying_G1, const double & p_dying_in_G1, const double & p_exiting_G1)
{
	static unsigned int buffer[4] = {0, 0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { p_staying_G1, p_dying_in_G1, p_exiting_G1, 1.0 - (p_staying_G1 + p_dying_in_G1 + p_exiting_G1) };
	const size_t K_S = 4;	

	gsl_ran_multinomial ( r_global, // Random Seed Generator 
						  K_S, 		// Number of output elements
						  G1_cells, // G1_i(t)
						  p_vector, //[p_{q,i}, 1 - p_{q,i} ]
						  buffer
						);
	return buffer;
}

/**
	Mitosis Stochastic Network
	IV) Updating G2 Phase

	[G2{S,i}(t), G2{G2->.,i}(t)] ~ Mult(G2_i(t),[p^G2_{S,i}(t), p^G2_{G2->S}(t) ])
	Where:
	G2_i(t) is the number of cells in G2 phase at time t
	G2{S,i}(t) is the number of cells that stay in G2 Phase
	G2{G2->.,i}(t) is the number of cells exiting G2 phase

****************************************************************************/
unsigned int * Random::Update_G2_Phase(const unsigned int & G2_cells, const double & p_staying_G2, const double & p_dying_in_G2, const double & p_exiting_G2)
{
	static unsigned int buffer[3] = {0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { p_staying_G2, p_dying_in_G2, p_exiting_G2 };
	const size_t K_S = 3;	

	gsl_ran_multinomial ( r_global, // Random Seed Generator 
						  K_S, 		// Number of output elements
						  G2_cells, // G2_i(t)
						  p_vector, //[p_{q,i}, 1 - p_{q,i} ]
						  buffer
						);
	return buffer;
}

/**
	Mitosis Stochastic Network
	V) Updating S Phase

	[S{S,i}(t), S{M->.,i}(t)] ~ Mult(S_i(t),[p^S_{S,i}(t), p^S_{S->M}(t) ])
	Where:
	S_i(t) is the number of cells in S phase at time t
	S{S,i}(t) is the number of cells that stay in S Phase
	S{M->.,i}(t) is the number of cells exiting S phase

****************************************************************************/
unsigned int * Random::Update_S_Phase(const unsigned int & S_cells, const double & p_staying_S, const double & p_dying_S, const double & p_exiting_S)
{
	static unsigned int buffer[3] = {0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { p_staying_S, p_dying_S, p_exiting_S };
	const size_t K_S = 3;	

	gsl_ran_multinomial ( r_global, // Random Seed Generator 
						  K_S, 		// Number of output elements
						  S_cells, // G2_i(t)
						  p_vector, //[p_{q,i}, 1 - p_{q,i} ]
						  buffer
						);
	return buffer;
}

/**
	Mitosis Stochastic Network
	VI) Updating M Phase

	[M{S,i}(t), M{EXIT->.,i}(t)] ~ Mult(M_i(t),[p^M_{S,i}(t), p^M_{M->EXIT}(t) ])
	Where:
	M_i(t) is the number of cells in M phase at time t
	M{S,i}(t) is the number of cells that stay in M Phase
	M{EXIT->.,i}(t) is the number of cells exiting M phase

****************************************************************************/
unsigned int * Random::Update_M_Phase(const unsigned int & M_cells, const double & p_staying_M, const double & p_dying_M, const double & p_exiting_M)
{
	static unsigned int buffer[3] = {0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { p_staying_M, p_dying_M, p_exiting_M };
	const size_t K_S = 3;	

	gsl_ran_multinomial ( r_global, // Random Seed Generator 
						  K_S, 		// Number of output elements
						  M_cells, // G2_i(t)
						  p_vector, //[p_{q,i}, 1 - p_{q,i} ]
						  buffer
						);
	return buffer;

}

unsigned int * Random::Clonal_Functions(const unsigned int & Available_cells, const double & p_of_going_G0, const double & p_entering_mitosis, const double & p_dying )
{
	static unsigned int buffer[4] = {0, 0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { 1.0 - (p_of_going_G0 + p_entering_mitosis + p_dying ), p_dying, p_of_going_G0, p_entering_mitosis };
	const size_t K_S = 4;	

	gsl_ran_multinomial ( r_global, // Random Seed Generator 
						  K_S, 		// Number of otput elements
						  Available_cells, //Bi(t)
						  p_vector, //[p_{q,i}, 1 - p_{q,i} ]
						  buffer
						);

	return buffer;
}

void Random::Newborn( const unsigned int & Newborn_cells, const  double & p_idle, const  double & p_mutate, const  double & p_go_to_G0, std::vector<unsigned int> & Division_Model)
{

	Division_Model.push_back(0);
	Division_Model.push_back(0);
	Division_Model.push_back(0);
	Division_Model.push_back(0);

	static unsigned int buffer[4] = {0, 0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = {p_idle, p_mutate, p_go_to_G0, 1.0 - (p_idle + p_mutate + p_go_to_G0) };
	const size_t K_S = 4;


	gsl_ran_multinomial ( r_global,
						  K_S,
						  2 * Newborn_cells,
						  p_vector,
						  buffer
						);

	Division_Model[0] = buffer[0];
	Division_Model[1] = buffer[1];
	Division_Model[2] = buffer[2];
	Division_Model[3] = buffer[3];



}

unsigned int * Random::Mutational_Effects()
{
	double p_vector[] = {KILLER_PROBABILITY, DRIVER_PROBABILITY, 1.0 - (KILLER_PROBABILITY + DRIVER_PROBABILITY) };
	static unsigned int buffer[3] = {0, 0, 0};
	const size_t K_S = 3;

	gsl_ran_multinomial ( r_global,
						  K_S,
						  1,
						  p_vector,
						  buffer
						);
	return buffer;
}

void Random::Laplace(double & effect)
{
	effect = gsl_ran_laplace (r_global, 0.001);

}

void Random::Valid_Limits(double & effect)
{
	if(effect < 0.0)
		effect = 0.0;
	else if( effect > 1.0 )
		effect = 1.0;

}

void Random::Laplace(double & effect, unsigned int & Clone_size)
{
	double a = (double) 1/Clone_size;
	//TODO debuging with this a = 0.001 
	a = 0.0001;
	effect += gsl_ran_laplace (r_global, a);
	Valid_Limits( effect );
	//std::cout << "Clone_Size " << Clone_size << " a " << a <<std::endl;
	
}

void Random::Mutational_Proportions( const double & p_go_to_G0,const  double & p_go_to_G1, const unsigned int & Mutant_Cells, std::vector<unsigned int> & Mutations )
{
	Mutations.push_back(0);
	Mutations.push_back(0);
	Mutations.push_back(0);

	static unsigned int buffer[3] = {0, 0, 0} ;// This will contain B_{G0,i}(t) and B_{G1,i}(t)
	double p_vector[] = { p_go_to_G0, p_go_to_G1, 1.0 - ( p_go_to_G0 + p_go_to_G1 ) };
	const size_t K_S = 3;

	gsl_ran_multinomial ( r_global,
						  K_S,
						  Mutant_Cells,
						  p_vector,
						  buffer
						);

	Mutations[0] = buffer[0];
	Mutations[1] = buffer[1];
	Mutations[2] = buffer[2];


	
}



