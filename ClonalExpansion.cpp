#include "ClonalExpansion.h"
#include "Clone.h"
#include "config.h"
#include "clonalFun.h"
#include "inout_funs.h"
#include "Random.h"
#include <map>
#include <memory>           // std::unique_ptr
#include <string>
#include <cstring>
#include <vector>
#include <iostream>         // std::cout,std::endl
#include <tuple>
#include <random>
#include <fstream>
#include <libgen.h>
#include <time.h>
#include <stdio.h>
#include <algorithm>
#include <sys/stat.h>
#include <mpi.h>
Clonal_Expansion::Clonal_Expansion()
 :
 	Population_Size( 0 ),
 	feedback( 0.0 ),
 	Mitosis( Clonal_Expansion::Delayed_V1 ),
 	Version( Clonal_Expansion::Tester ),
 	SD_Penalty( Clonal_Expansion::standard ),
 	PR_Sampling( Clonal_Expansion::Uniform ),
 	MR_Sampling( Clonal_Expansion::uniform),
 	Mutational_Effects( Clonal_Expansion::Normal_gradient ),
 	Passenger_Distribution(Clonal_Expansion::Passenger_Beta),
 	Deleterious_Distribution(Clonal_Expansion::Deleterious_Beta),
 	Beneficial_Distribution(Clonal_Expansion::Beneficial_Beta),
 	Penalty_Max_Population_Size(PS),
	PR_DR_Difference(DIFF),
	Multiple_Mutations(false),
	Newborn_Clones_at_t(0),
	Population_at_prev_t(0),
	Clonality_at_prev_t(0),
	Minimum_Clone_Size_Wrtite_TE(100),
	Delta_Population_Write_TEM(600),
	Delta_Clone_Write_TEM(10)
 	{       }

Clonal_Expansion::~Clonal_Expansion()
 {       }

 unsigned int Clonal_Expansion::getSize()
 {
 	return Tumour -> size();
 }

 void Clonal_Expansion::printParameters(void)
 {
 	std::cout << "\n\n##################################\n\n CLONAL EXPANSION DS \n{ " << std::endl;
 	std::cout << "\tSelective Pressure: " << feedback << std::endl;
 	std::cout << "\tPopulation Size: " << Population_Size << std::endl;
 	std::cout << "\tNumber of Clones: " << Tumour -> size() << std::endl;
 	std::cout << "}\n#################################" << std::endl; 
 }



 void Clonal_Expansion::Init_Random(void)
 {
 	long seed;
	r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  	seed = time (NULL) * getpid();    
  	gsl_rng_set (r_global, seed);  
 }

 void Clonal_Expansion::printValues (std::unique_ptr<Clone> const & ith_clone)
 {
 	std::cout << "\n\n################# C L O N E   V A L U E S  #################" << "\n\n";
 	std::cout << "################# G E N E R A L  V A L U E S  #################" << "\n\n";
 	std::cout << "Number of Members to Start Seeding Heterogeneity: " << ith_clone -> Number_of_Memebers_to_Start_Heterogeneity << "\n";
 	std::cout << "Generation Branch ID: " << ith_clone ->  Generation_ID_Counter << "\n";
 	std::cout << "Log(10) Penalty Accumulation: " << ith_clone ->  Driver_10_fold_accumulation << "\n";
 	std::cout << "Is the Clone extinct: " << ith_clone ->  clone_extinct << "\n";
 	std::cout << "Initial Proliferation Rate: " << ith_clone ->  init_PR << "\n";
 	std::cout << "Maximum Proliferation Rate Achieved: " << ith_clone ->  max_PR << "\n";
 	std::cout << "Final Proliferation Rate: " << ith_clone -> final_PR << "\n";
 	std::cout << "\n";
 	std::cout << "################# M I T O S I S #################" << "\n\n";
 	std::cout << "Is In Initial Expansion Period: " <<ith_clone ->  Initiall_Expasion_Period << "\n";
 	std::cout << "In G0: " << ith_clone -> In_G0_Phase << "\n";
 	std::cout << "In G1: " << ith_clone -> In_G1_Phase << "\n";
 	std::cout << "In S: " << ith_clone -> In_S_Phase << "\n";
 	std::cout << "In G2: " << ith_clone -> In_G2_Phase << "\n";
 	std::cout << "In M: " << ith_clone -> In_M_Phase << "\n";
 	std::cout << "Remaining Time in G1: " << ith_clone -> Remaining_Time_in_G1_Phase << "\n";
 	std::cout << "Remaining Time in G2: " << ith_clone -> Remaining_Time_in_G2_Phase << "\n";
 	std::cout << "Remaining Time in S: " << ith_clone -> Remaining_Time_in_S_Phase << "\n";
	std::cout << "Remaining Time in M: " << ith_clone -> Remaining_Time_in_M_Phase << "\n";
 	std::cout << "\n";
 	std::cout << "################# E X P A N S I O N  #################" << "\n\n";
 	std::cout << "Clone Size: " << ith_clone -> Clone_Size << "\n";
 	std::cout << "Death Rate : " << ith_clone ->  P_Expansion[0] << "\n";
 	std::cout << "Proliferation Rate : " << ith_clone -> P_Expansion[1] << "\n";
 	std::cout << "Mutation Rate : " << ith_clone -> P_Expansion[2] << "\n";
 	std::cout << "\n";
 	std::cout << "############################################" << "\n\n";

 }


 void Clonal_Expansion::add_Clone_to_Tumour(void)
 {
 	Tumour -> push_back ( get_Clone_DS() );
 }

std::string Clonal_Expansion::getMitosis_type(void)
 {
 	std::string types[] = {"Delayed_V1", "Delayed_V2", "Network"};
 	return types[Mitosis]; 
 }

std::string Clonal_Expansion::getVersion_type(void)
 {
 	std::string types[] = {"Tester", "V1", "V2", "V2R", "V3", "VM"};
 	return types[Version]; 
 }

 std::string Clonal_Expansion::getSD_Penalty(void)
 {
 	std::string types[] = {"standard", "mean", "quartile"};
 	return types[SD_Penalty]; 
 }

 std::string Clonal_Expansion::getMutational_Effect_Dist(void)
 {
 	std::string types[] = {"Normal_gradient","Exponential","Double_Beta"};
 	return types[SD_Penalty]; 
 }

 void Clonal_Expansion::setMutational_Effect_Dist(const std::string & str_mutational_effect_dist)
 {
 	if( str_mutational_effect_dist == "Normal_gradient" )
 	{
 		Mutational_Effects = Clonal_Expansion::Normal_gradient;
 	}
 	else if(str_mutational_effect_dist == "Exponential")
 	{
 		Mutational_Effects = Clonal_Expansion::Exponential;
 	}
 	else if( str_mutational_effect_dist == "Double_Beta" )
 	{
 		Mutational_Effects = Clonal_Expansion::Double_Beta;
 	}
 	else
 	{
 		Mutational_Effects = Clonal_Expansion::Normal_gradient;
 	}
 }




void Clonal_Expansion::setMR_Sampling_type(const std::string & str_mr_sampling)
{

	if(str_mr_sampling == "uniform")
	{
		MR_Sampling = Clonal_Expansion::uniform;
	}
	else if (str_mr_sampling == "uniform_Gain")
	{
		MR_Sampling = Clonal_Expansion::uniform_Gain;
	}
	else if (str_mr_sampling == "beta")
	{
		MR_Sampling = Clonal_Expansion::beta;
	}
	else 
	{
		MR_Sampling = Clonal_Expansion::uniform;
	}
}

void Clonal_Expansion::setPR_Sampling_type(const std::string & str_sampling)
{
	if(str_sampling == "Uniform")
	{
		PR_Sampling = Clonal_Expansion::Uniform;
	}
	else if(str_sampling == "Uniform_V1")
	{
		PR_Sampling = Clonal_Expansion::Uniform_V1;
	}
	else if(str_sampling == "Beta")
	{
		PR_Sampling = Clonal_Expansion::Beta;
	}
	else if(str_sampling == "Pareto" )
	{
		PR_Sampling = Clonal_Expansion::Pareto;
	}
	else if(str_sampling == "Gamma")
	{
		PR_Sampling = Clonal_Expansion::Gamma;
	}
	else if(str_sampling == "Laplace")
	{
		PR_Sampling = Clonal_Expansion::Laplace;
	}
	else
	{
		PR_Sampling = Clonal_Expansion::Uniform;
	}
}

 void Clonal_Expansion::setMistosis_type(const std::string & str_mitosis)
 {
 	if(str_mitosis == "Delayed_V1")
 	{
 		Mitosis = Clonal_Expansion::Delayed_V1;
 	}
 	else if(str_mitosis == "Delayed_V2")
 	{
 		Mitosis = Clonal_Expansion::Delayed_V2;
 	}
 	else if(str_mitosis == "Network")
 	{
 		Mitosis = Clonal_Expansion::Network;
 	}
 	else
 	{
 		Mitosis = Clonal_Expansion::Delayed_V1;
 	}
 }

  void Clonal_Expansion::setVersion_type(const std::string & str_version)
 {

 	if(str_version == "Tester")
 	{
 		Version = Clonal_Expansion::Tester;
 	}
 	else if(str_version == "V1")
 	{
 		Version = Clonal_Expansion::V1;

 	}
 	else if(str_version == "V2")
 	{
 		Version = Clonal_Expansion::V2;

 	}
 	else if(str_version == "V2R")
 	{
 		Version = Clonal_Expansion::V2R;

 	}
 	else if(str_version == "V3")
 	{
 		Version = Clonal_Expansion::V3;

 	}
 	else if(str_version == "VM")
 	{
 		Version = Clonal_Expansion::VM;

 	}
 	else
 	{
 		Version = Clonal_Expansion::VM;
 	}

 }

void Clonal_Expansion::setSD_Penalty(const std::string & str_Penalty)
{
	if(str_Penalty == "standard")
	{
		SD_Penalty = Clonal_Expansion::standard;
	}
	else if(str_Penalty == "mean")
	{
		SD_Penalty = Clonal_Expansion::mean;
	}
	else if(str_Penalty == "quartile")
	{
		SD_Penalty = Clonal_Expansion::quartile;
	}
	else
	{
		SD_Penalty = Clonal_Expansion::standard;
	}

}

void Clonal_Expansion::setPassenger_Distribution( const std::string & str_Passenger_Distribution )
{
	if( str_Passenger_Distribution == "Passenger_Beta" )
	{
		Passenger_Distribution = Clonal_Expansion::Passenger_Beta;
	}
	else if( str_Passenger_Distribution == "Passenger_Pareto" )
	{
		Passenger_Distribution = Clonal_Expansion::Passenger_Pareto;
	}
	else if(str_Passenger_Distribution == "Passenger_Laplace" )
	{
		Passenger_Distribution = Clonal_Expansion::Passenger_Laplace;
	}
	else if(str_Passenger_Distribution == "Passenger_Gamma" )
	{
		Passenger_Distribution = Clonal_Expansion::Passenger_Gamma;
	}
	else
	{
		Passenger_Distribution = Clonal_Expansion::Passenger_Beta;
	}
}

void Clonal_Expansion::setDeleterious_Distribution( const std::string & str_Deleterious_Distribution )
{
	if( str_Deleterious_Distribution == "Deleterious_Beta" )
	{
		Deleterious_Distribution = Clonal_Expansion::Deleterious_Beta;
	}
	else if(str_Deleterious_Distribution == "Deleterious_Pareto")
	{
		Deleterious_Distribution = Clonal_Expansion::Deleterious_Pareto;
	}
	else if( str_Deleterious_Distribution == "Deleterious_Laplace" )
	{
		Deleterious_Distribution = Clonal_Expansion::Deleterious_Laplace;
	}
	else if(str_Deleterious_Distribution == "Deleterious_Gamma")
	{
		Deleterious_Distribution = Clonal_Expansion::Deleterious_Gamma;
	}
	else
	{
		Deleterious_Distribution = Clonal_Expansion::Deleterious_Beta;
	}
}

void Clonal_Expansion::setBeneficial_Distribution(const std::string & str_Beneficial_Distribution )
{
	if(str_Beneficial_Distribution == "Beneficial_Beta" )
	{
		Beneficial_Distribution = Clonal_Expansion::Beneficial_Beta;
	}
	else if( str_Beneficial_Distribution == "Beneficial_Pareto" )
	{
		Beneficial_Distribution = Clonal_Expansion::Beneficial_Pareto;
	}
	else if (str_Beneficial_Distribution == "Beneficial_Laplace" )
	{
		Beneficial_Distribution = Clonal_Expansion::Beneficial_Laplace;
	}
	else if ( str_Beneficial_Distribution == "Beneficial_Gamma")
	{
		Beneficial_Distribution = Clonal_Expansion::Beneficial_Gamma;
	} 
	else 
	{
		Beneficial_Distribution = Clonal_Expansion::Beneficial_Beta;
	}
}


void Clonal_Expansion::Passenger_Beta_Sampling( double & Penalty)
{
	Penalty = gsl_ran_beta (r_global, Passenger_Distribution_Parameters.at(0), Passenger_Distribution_Parameters.at(1));
}

void Clonal_Expansion::Select_Passenger_Distribution(double & Penalty)
{
	switch(Passenger_Distribution)
	{
		case 0:
			//std::cout << "Passenger_Beta " << std::endl;
			Passenger_Beta_Sampling( Penalty );
		break;
		case 1:
			std::cout << "Passenger_Pareto " << std::endl;
		break;
		case 2:
			std::cout << " Passenger_Laplace " << std::endl;
		break;
		case 3: 
			std::cout << "Passenger_Gamma " << std::endl;
		break;
		default:
			std::cout << "Passenger_Beta " << std::endl;
		break;


	}
}

void Clonal_Expansion::Deleterious_Beta_Sampling( double & Penalty )
{
	Penalty = gsl_ran_beta (r_global, Deleterious_Distribution_Parameters.at(0), Deleterious_Distribution_Parameters.at(1));
}

void Clonal_Expansion::Select_Deleterious_distribution(double & Penalty )
{
	switch(Deleterious_Distribution)
	{
		case 0:
			//std::cout << "Deleterious_Beta " << std::endl;
			Deleterious_Beta_Sampling( Penalty );
		break;
		case 1:
			std::cout << "Deleterious_Pareto " << std::endl;
		break;
		case 2:
			std::cout << "Deleterious_Laplace " << std::endl;
		break;
		case 3: 
			std::cout << "Deleterious_Gamma " << std::endl;
		break;
		default:
			std::cout << "Deleterious_Beta " << std::endl;
		break;
	}
}

void Clonal_Expansion::Beneficial_Beta_Sampling( double & Penalty )
{
	Penalty = gsl_ran_beta (r_global, Beneficial_Distribution_Parameters.at(0), Beneficial_Distribution_Parameters.at(1));
}

void Clonal_Expansion::Select_Beneficial_distribution(double & Penalty)
{
		switch(Beneficial_Distribution)
	{
		case 0:
			//std::cout << "Beneficial_Beta " << std::endl;
			Beneficial_Beta_Sampling( Penalty );
		break;
		case 1:
			std::cout << "Beneficial_Pareto " << std::endl;
		break;
		case 2:
			std::cout << "Beneficial_Laplace " << std::endl;
		break;
		case 3: 
			std::cout << "Beneficial_Gamma " << std::endl;
		break;
		default:
			std::cout << "Beneficial_Beta " << std::endl;
		break;
	}
}



void Clonal_Expansion::Map_Feedback_Penalty()
{
	feedback =  (PR_DR_Difference - 0.0) * (( (double) Population_Size - 0.0) / ( Penalty_Max_Population_Size - 0.0));
	//std::cout << "Feedback " << feedback << std::endl;
}

void Clonal_Expansion::Induce_Multiple_Mutations(unsigned int & number_of_mutations)
{
	if(Multiple_Mutations)
	{
		r.Multiple_Mutations_Poisson(number_of_mutations);
	}
	else
	{
		number_of_mutations = 1;
	}
}


//Mutational effects distribution
void Clonal_Expansion::Normal_Gradient_Mutational_Effect(const unsigned int & ith_clone, std::tuple<unsigned long long int, 
																									  unsigned long long int, 
																									  unsigned long long int, 
																									  bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years)
{
	unsigned int buffer[5];
	unsigned long long int Mutants_in_clone = 0;
	bool Death_or_Driver = false;
	
	double p[] = { 
				   MEFF_Sampling_Parameters.at(0), // Killer P
				   MEFF_Sampling_Parameters.at(1), // Driver P
				   MEFF_Sampling_Parameters.at(2), // Passenger P
				   MEFF_Sampling_Parameters.at(3), // Deleterious P
				   MEFF_Sampling_Parameters.at(4)  // Beneficial P
				};

	// std::cout 
	// << MEFF_Sampling_Parameters.at(0) 
	// << " " <<  MEFF_Sampling_Parameters.at(1) 
	// << " " <<  MEFF_Sampling_Parameters.at(2) 
	// << " " <<  MEFF_Sampling_Parameters.at(3)
	// << " " <<  MEFF_Sampling_Parameters.at(4)  
	// <<std::endl;
	// getchar();

	unsigned int number_of_mutations;

	for(unsigned int j = 0; j < std::get<2>(Dying_and_Newborn) ; j++)
	{
		Induce_Multiple_Mutations(number_of_mutations);
		for(unsigned int i = 0; i < number_of_mutations; i++)
		{
			gsl_ran_multinomial ( r_global, 5, 1, p, buffer);
			if(buffer[0] > 0) // Killr mutation
			{
				Death_or_Driver = true;
				Tumour -> at(ith_clone) -> Clonal_Mutations++;
				break;
			}
			else if(buffer[1] > 0) //Driver mutation
			{
				Death_or_Driver = true;
				Tumour -> at(ith_clone) -> Clonal_Mutations++;
				Carcinogenesis_From_Driver_V2( ith_clone, hours, years );
				break;
			}
			else if(buffer[2] > 0)//Passenger
			{
				double Penalty;
				Select_Passenger_Distribution( Penalty );
				Tumour -> at(ith_clone) -> Clonal_Mutations++;
				Tumour -> at(ith_clone) -> Clonal_Mutational_Burden -= Penalty;
				
				Tumour -> at(ith_clone) -> P_Expansion[1] -= Penalty;//apply penlaty

				if(Tumour -> at(ith_clone) -> P_Expansion[1] < 0.0)
					Tumour -> at(ith_clone) -> P_Expansion[1] = 0.0;


			}
			else if(buffer[3] > 0)//Deleterious
			{
				double Penalty;
				Select_Deleterious_distribution( Penalty );
				Tumour -> at(ith_clone) -> Clonal_Mutations++;
				Tumour -> at(ith_clone) -> Clonal_Mutational_Burden -= Penalty;
				
				Tumour -> at(ith_clone) -> P_Expansion[1] -= Penalty; // apply penalty to clone

				if(Tumour -> at(ith_clone) -> P_Expansion[1] < 0.0)
					Tumour -> at(ith_clone) -> P_Expansion[1] = 0.0;
			}
			else if (buffer[4] > 0)// Beneficial Mutation
			{
				double last_recorded_max = Tumour -> at(ith_clone) -> max_PR;
				double Penalty;
				Select_Beneficial_distribution( Penalty );
				Tumour -> at(ith_clone) -> Clonal_Mutations++;
				Tumour -> at(ith_clone) -> Clonal_Mutational_Burden += Penalty;
				
				Tumour -> at(ith_clone) -> P_Expansion[1] += Penalty; // apply Penalty



				if(Tumour -> at(ith_clone) -> P_Expansion[1] > 0.2)
					Tumour -> at(ith_clone) -> P_Expansion[1] = 0.2;

				if(last_recorded_max < Tumour -> at(ith_clone) -> P_Expansion[1] )
					Tumour -> at(ith_clone) -> max_PR = Tumour -> at(ith_clone) -> P_Expansion[1];

			}
			else
			{
				std::cout << "NOT FOUND " << std::endl;
			}
		}// for mutants
		// if killer or passenger mutation, adjust number of mutants
		if(Death_or_Driver)
			Mutants_in_clone++;

	}
	//substract clones that died or got driver mutations
	if(Mutants_in_clone > std::get<2>(Dying_and_Newborn) )
		std::get<2>(Dying_and_Newborn) = 0;
	else
		std::get<2>(Dying_and_Newborn) -=  Mutants_in_clone;

}


void Clonal_Expansion::Select_Mutational_Effects_Distribution(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years, 
															   std::tuple<unsigned long long int, 
																			   unsigned long long int, 
																			   unsigned long long int, 
																			   bool, bool> & Dying_and_Newborn)
{
	switch(Mutational_Effects)
	{
		case 0:
			//std::cout << "Send to Normal_gradient mutational effects " << std::endl;
			Normal_Gradient_Mutational_Effect( ith_clone, Dying_and_Newborn, hours, years );
			break;
		case 1:
			std::cout << "Send to Exponenetial " << std::endl;
			break;
		case 2:
			std::cout << " Send to Double Beta " << std::endl;
			break;
		default:
			std::cout << "Defualt case " << std::endl;
			break;
	}
}


void Clonal_Expansion::Select_Mitosis_Type(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	switch(Mitosis)
	{
		case 0:
			Dealyed_Mitosis_V1(ith_clone, hours, years);
			break;
		case 1:
			std::cout << "Binomial " << std::endl;
			break;
		case 2:
			std::cout << "Network " << std::endl;
			break;
		default:
			Dealyed_Mitosis_V1(ith_clone, hours, years);
			break;
	}
}

void Clonal_Expansion::Select_Size_Dependant_Penalty(void)
{
	switch(SD_Penalty)
	{
		case 0:
			Map_Feedback_Penalty();
			break;
		case 1:
			std::cout << " Mean Penalty " << std::endl;
			break;
		case 2:
			std::cout << " Quartile " << std::endl;
			break;
		default:
			std::cout << "STD default " << std::endl;
			break;

	}//switch
}


 void Clonal_Expansion::Select_Carcinogenesis(void)
 {
 	switch(Version)
 	{
 		case 0:
            std::cout << "0 Version: " << Version << std::endl;
            break;
        case 1:
            std::cout << "1 Version: " << Version << std::endl;
            carcinogenesis_V1();
            break;
        case 2:
            std::cout << "2 Version: " << Version << std::endl;
            break;
        case 3:
            std::cout << "3 Version: " << Version << std::endl;
            break;
        case 4:
            std::cout << "4 Version: " << Version << std::endl;
            break;
        case 5:
            std::cout << "5 Version: " << Version << std::endl;
            break;
        default:
        	 std::cout << "Def Version: " << Version << std::endl;
        	 break;
 	}//switch	
 
 }

void Clonal_Expansion::Select_PR_Samplig(const double & Parent_Proliferation_Rate, double & updated_PR)
{

	switch(PR_Sampling)
	{
		case 0:
			r.Uniform_PR_Update( Parent_Proliferation_Rate, PR_Sampling_Parameters.at(0), updated_PR);
			break;
		case 1: 
			std::cout << " Uniform V1 " << std::endl;
			break;
		case 2:
			std::cout << " Beta " << std::endl;
			break;
		case 3:
			std::cout << " Pareto " << std::endl;
			break;
		case 4:
			std::cout << " Gamma " << std::endl;
			break;
		case 5:
			std::cout << " Laplace " << std::endl;
			break;
		default:
			std::cout << " Uniform " << std::endl;
	}
}

void Clonal_Expansion::Select_MR_Sampling(const double & Parent_mu_rate, double & updated_MR)
{
	switch(MR_Sampling)
	{
		case 0:
			std::cout << "Uniform " << std::endl;
			break;
		case 1:
			//std::cout << " unfirom Gain " << std::endl;
			r.Uniform_Gain_Update( Parent_mu_rate,  MR_Sampling_Parameters.at(0), updated_MR);
			break;

		case 2:
			std::cout << " beta " << std::endl;
			break;
		default:
			std::cout << " default " << std::endl; 
	}

}

void Clonal_Expansion::Select_Penalty_Type(void)
{

	std::cout << Clonal_Expansion::SD_Penalty << "  penalty " << std::endl;

	switch(SD_Penalty)
	{
		case 0:
			std::cout << "Standard " << std::endl;
		break;
		case 1:
		 	std::cout << "Mean " << std::endl;
		 break;
		 case 2:
		 	std::cout<< " quartile " << std::endl;
		 	break;
		 default:
		 	std::cout << "default " << std::endl;
		 	break;
	}
}



/**
	Adds a parent clone with size 1 in the tumour.
	Steps:
	1) Ask for a new Clone DS and append to vector Tumour
	2) From the recently added clone, initialise values

*******************************************/
void Clonal_Expansion::carcinogenesis(void)
{

	// 1) Add a Clone
	// Adding a new clone
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;
	Tumour -> back() -> clone_extinct = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();

	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[P_Expansion_PR];
	Tumour -> back() -> init_PR  = Tumour -> back() -> P_Expansion[P_Expansion_PR];


	Population_Size++;
	Mitosis = Clonal_Expansion::Delayed_V1; // Requires modification from input parameter
	
	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
}

/*
	Initialisation of Carcinogenesis of Version 1
*/
void Clonal_Expansion::carcinogenesis_V1(void)
{
	//std::cout << "Carcinogenesis Version 1 " << std::endl;
	// 1) Add a Clone
	//Random r;
	// Adding a new clone
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	

	Population_Size = Tumour -> back() -> Clone_Size;

}


void Clonal_Expansion::carcinogenesis_V2(void)
{
	//std::cout << "Carcinogenesis Version 1 " << std::endl;
	// 1) Add a Clone
	//Random r;
	// Adding a new clone
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> Clonal_Mutations = 1;
	Tumour -> back() -> Clonal_Mutational_Burden = 0.0;
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[1];

	Population_Size = Tumour -> back() -> Clone_Size;
}

void Clonal_Expansion::carcinogenesis_V1_File_Input(const double & mutation_rate)
{
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	

	Population_Size = Tumour -> back() -> Clone_Size;
	Tumour -> back() -> P_Expansion[2] = mutation_rate;
}

void Clonal_Expansion::carcinogenesis_V2_File_Input(const double & mutation_rate, const double & proliferation_rate)
{
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> Clonal_Mutations = 1;
	Tumour -> back() -> Clonal_Mutational_Burden = 0.0;
	Tumour -> back() -> max_PR = proliferation_rate;
	Tumour -> back() -> init_PR = proliferation_rate;

	Population_Size = Tumour -> back() -> Clone_Size;

	Tumour -> back() -> P_Expansion[2] = mutation_rate;
	Tumour -> back() -> Mutation_Rate = mutation_rate;
	Tumour -> back() -> P_Expansion[1] = proliferation_rate;
	//std::cout << "MR of parent clone " << Tumour -> back() -> P_Expansion[2] << std::endl; 
	//std::cout << "PR of parent clone " << Tumour -> back() -> P_Expansion[1] << std::endl; 
	

}




void Clonal_Expansion::Generate_Clone_Generation_ID( std::string & newGeneration_ID,
								  const int & Parent_Generation_ID_Counter, 
								  const std::string & Parent_ID, 
								  const unsigned int & years, 
								  const unsigned int & hours)
{

	std::string heredity_pattern = Parent_ID.substr( 0, Parent_ID.find("-") );
	std::string Clone_ID ("");
		//cout << "I GET " << heredity_pattern << endl;

		if(heredity_pattern.compare("P") == 0)
		{
			//cout << "I GET 1 " << heredity_pattern << endl;
			std::string Clone_ID_tag ("PC");
			Clone_ID_tag.append(std::to_string(Parent_Generation_ID_Counter));Clone_ID_tag.append("-");
			Clone_ID_tag.append(std::to_string(years));Clone_ID_tag.append(":");Clone_ID_tag.append(std::to_string(hours));
			Clone_ID = Clone_ID_tag;
		}
		else
		{
			
			std::string Clone_ID_tag ("");
			Clone_ID_tag.append(heredity_pattern); Clone_ID_tag.append(","); Clone_ID_tag.append(std::to_string(Parent_Generation_ID_Counter));
			Clone_ID_tag.append("-");
			Clone_ID_tag.append(std::to_string(years));Clone_ID_tag.append(":");Clone_ID_tag.append(std::to_string(hours));
			Clone_ID = Clone_ID_tag;

		}

	//std::cout << "ID " << Parent_ID << std::endl;
//	std::cout << "HP " << heredity_pattern << std::endl;

	newGeneration_ID = Clone_ID;

}




void Clonal_Expansion::carcinogenesis_from_driver(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours)
{
	

	double mr = Tumour -> at( ith_clone ) -> Mutation_Rate;
	unsigned int AC = Tumour -> at( ith_clone ) -> Driver_10_fold_accumulation;
	//unsigned int Accumulated_Drivers = Tumour -> at( ith_clone ) -> Driver_10_fold_accumulation;
	//unsigned long long int NOM = Tumour -> at( ith_clone ) -> Number_of_Mutations;
	double pr = Tumour -> at( ith_clone ) -> P_Expansion[1];

	double p_st_g0 =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status );
	double p_dy_g0 =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status );

	double p_st_g1 =  std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status );
	double p_dy_g1 =  std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status );

	double p_st_g2 =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status );
	double p_dy_g2 =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status );

	double p_st_s =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status );
	double p_dy_s =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status );

	double p_st_m =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status );
	double p_dy_m =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status );

	std::string cloneName = "";	
	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) -> Generation_ID_Counter;
	Generate_Clone_Generation_ID(cloneName, 
								 Parent_Generation_ID_Counter, 
								 Tumour -> at( ith_clone ) -> Generation_ID,
								 years, 
								 hours);

	Tumour -> push_back( get_Clone_DS() );

	Tumour -> back() -> Generation_ID = cloneName;
	Tumour -> back() -> Initiall_Expasion_Period = false;
	Tumour -> back() -> Clone_Size = 1;
	//Tumour -> back() -> Number_of_Mutations = NOM + 1;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Generation_ID_Counter = 0;
	Tumour -> back() -> Mutation_Rate = r.Uniform_Mutation_Rate_2(mr);
	Tumour -> back() -> P_Expansion[2] = Tumour -> back() -> Mutation_Rate;
	double PR = r.Update_Proliferation_Rate(pr);
	Tumour -> back() -> Driver_10_fold_accumulation = AC + 5;
	Tumour -> back() -> P_Expansion[1] =  PR;
	Tumour -> back() -> init_PR = PR;
	Tumour -> back() -> max_PR = PR;
	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;

	Tumour -> back() -> G0_cells = 0;
	std::get<P_STAYING>( Tumour -> back() -> G0_status ) = p_st_g0;
	std::get<P_DYING>( Tumour -> back() -> G0_status ) = p_dy_g0;

	Tumour -> back() -> G1_cells = 1;
	std::get<P_STAYING_G1>( Tumour -> back() -> G1_status ) = p_st_g1;
	std::get<P_DYING_G1>( Tumour -> back() -> G1_status ) = p_dy_g1;

	 Tumour -> back() -> G2_cells = 0;
	 std::get<P_STAYING>( Tumour -> back() -> G2_status ) = p_st_g2;
	 std::get<P_DYING>( Tumour -> back() -> G2_status ) = p_dy_g2;


	 Tumour -> back() -> S_cells = 0;
	 std::get<P_STAYING>( Tumour -> back() -> S_status ) = p_st_s;
	 std::get<P_DYING>( Tumour -> back() -> S_status ) = p_dy_s;

	 Tumour -> back() -> M_cells = 0;
	 std::get<P_STAYING>( Tumour -> back() -> M_status ) = p_st_m;
	 std::get<P_DYING>( Tumour -> back() -> M_status ) = p_dy_m;
	

	Population_Size++;

	//Tumour -> back() -> printValues();


	//std::cout << "\n\n \t\t\t VS \n\n";

	//Tumour -> at( ith_clone ) -> printValues();
	
	//ith_clone = Tumour -> at(Tumour -> size() - 2 ) const;
	//printValues(Tumour -> at (Tumour -> size() - 2 ) );


}

/**
	Adds a parent clone with size neoplastic_cells in the tumour.
	Steps:
	1) Ask for a new Clone DS and append to vector Tumour
	2) From the recently added clone, initialise values

*******************************************/
void Clonal_Expansion::carcinogenesis(unsigned long long int &neoplastic_cells)
{
	

	Tumour -> push_back( get_Clone_DS() );
	Tumour -> back() -> Clone_Size = neoplastic_cells;

	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[P_Expansion_PR];

	Population_Size = Tumour -> back() -> Clone_Size;
	Mitosis = Clonal_Expansion::Delayed_V1;
		
	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
}

/**
	Check if the ith_clone is in G1 mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_G1_phase( const unsigned int & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G1_Phase && Tumour -> at (ith_clone) -> Remaining_Time_in_G1_Phase !=0);
}


/**
	Check if the ith_clone is exiting the G1 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_G1_phase( const unsigned int  & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G1_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase == 0);
}

/**
	Check if the ith_clone is in S mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_S_phase( const unsigned int & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_S_Phase && Tumour -> at (ith_clone) -> Remaining_Time_in_S_Phase !=0);
}

/**
	Check if the ith_clone is exiting the G1 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_S_phase( const unsigned int  & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_S_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase == 0);
}

/**
	Check if the ith_clone is in G2 mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_G2_phase( const unsigned int  & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G2_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase !=0);
}

/**
	Check if the ith_clone is exiting the G2 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_G2_phase( const unsigned int & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G2_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase == 0);
}

/**
	ith_Clone from G1 -> S phase
******************************************************/
void Clonal_Expansion::Transition_From_G1_S(const unsigned int & ith_clone)
{
	Tumour -> at(ith_clone) -> In_G1_Phase  = false; 
	Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> at(ith_clone) -> In_S_Phase = true;
}

/**
	ith_Clone from S -> G2 phase
******************************************************/
void Clonal_Expansion::Transition_From_S_G2( const unsigned int  & ith_clone)
{
	Tumour -> at (ith_clone) -> In_S_Phase  = false; 
	Tumour -> at (ith_clone) -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> at (ith_clone) -> In_G2_Phase = true;
}

/**
	ith_Clone from G2 -> M phase
******************************************************/
void Clonal_Expansion::Transition_From_G2_M(const unsigned int  & ith_clone)
{
	Tumour -> at(ith_clone) -> In_G2_Phase  = false; 
	Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> at(ith_clone) -> In_M_Phase = true;
}



// TODO Asymetric division
void Clonal_Expansion::Update_Population_After_Division_V1(const unsigned int  & ith_clone)
{

	unsigned long long int cells_after_division = Tumour -> at (ith_clone) -> Clone_Size * 2;
	Population_Size =  (Population_Size - Tumour -> at(ith_clone) -> Clone_Size ) + cells_after_division;
	Tumour -> at(ith_clone) -> Clone_Size = cells_after_division;
}

void Clonal_Expansion::Estimate_Mutational_Effects(unsigned int & Mutant_Cells, unsigned long long int & Clone_Size) //years and hours
{

	

	unsigned int cs = static_cast<unsigned int>(Clone_Size);
	unsigned int * mutational_effect_buffer;
	double effect = 0.0;

	unsigned int kills = 0;
	unsigned int drivers = 0;
	double mutational_burden = 0.0;

	for(unsigned int i = 0; i < Mutant_Cells  ; i ++)
	{
		for(unsigned int j = 0; j < add_multiple_mutations() ; j++)
		{
			//std::cout << add_multiple_mutations() << " mutations to add for cell " << i << std::endl; 
			mutational_effect_buffer = r.Mutational_Effects();
			std::cout << "Killer Mutation " << mutational_effect_buffer[0]
					  << " Driver Mutation " <<  mutational_effect_buffer[1]
					  << " Gradient Effect "<<   mutational_effect_buffer[2]
					  << std::endl;


			if(mutational_effect_buffer[0] > 0)
			{
				std::cout << "Highly deleterious mutation... inducing cell death" << std::endl;
				kills++;
			}
			else if(mutational_effect_buffer[1] > 0)
			{
				std::cout << "Driver mutation... generating a new clone " << std::endl;
				drivers++;
			}
			else
			{
				r.Laplace(effect, cs );
				mutational_burden += effect;
				std::cout << " Mutation with effect of " << effect << std::endl;
			}

		}
	}

	std::cout<< " TOTAL KILLS " << kills << " DRIVERS " << drivers << " mutational burden " << mutational_burden << std::endl;  
	
}


void Clonal_Expansion::Check_for_Clonal_Extintion_at_min_size(const unsigned int & ith_clone)
{

	if( Tumour -> at(ith_clone) -> Clone_Size == 1 )
	{
		Tumour -> at(ith_clone) -> Clone_Size = 0 ;
		Tumour -> at(ith_clone) -> clone_extinct = true; 
	}	

}



void Clonal_Expansion::Update_Clonal_Mutational_Burden(const unsigned int & ith_clone, const unsigned int & Mutant_cells, const unsigned int & years, const unsigned int & hours, std::vector<unsigned int> & Mutations )
{
	
	unsigned int * mutational_effect_buffer;
	unsigned int counter = 0;
	unsigned int death = 0;
	double effect_gain = static_cast<double>( Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] );
	//bool mut_effect_flag = false;
	

	unsigned int j = 0;


	while( (j < Mutant_cells  ) && !(Tumour -> at(ith_clone) -> clone_extinct) )
	{
		unsigned int i = 0;
		while( i < add_multiple_mutations() )
		{
			// check for mutations
			mutational_effect_buffer = r.Mutational_Effects();
			// buffer is size 3 --> need to check kills drivers and effects
			std::cout << "P(Ki) = " << KILLER_PROBABILITY << "; P(Dr) = " << DRIVER_PROBABILITY << " ; P(P) =  " << 1.0 -(KILLER_PROBABILITY+DRIVER_PROBABILITY) << " [ " << i << " cells]:  "<< " KILL = " <<  mutational_effect_buffer[0] << " DR = " <<  mutational_effect_buffer[1] << " PASS= " << mutational_effect_buffer[2] << std::endl;
			if ( mutational_effect_buffer[PASSENGER_CELLS] > 0 )// Passenger
			{
			 	std::cout << " MUTATIONAL EFFECT " << std::endl; 
			 	unsigned int CS = static_cast<unsigned int>(Tumour -> at(ith_clone) -> Clone_Size);
			 	// std::cout << "Before Expansion " 
			 	// 		  << Tumour -> at(ith_clone) -> P_Expansion[1] << " " 
			 	// 		  << Tumour -> at(ith_clone) -> Clone_Size << std::endl;

			 	r.Laplace( Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR], CS );

			 	// std::cout << "After Expansion " 
			 	// 		  << Tumour -> at(ith_clone) -> P_Expansion[1] << " " 
			 	// 		  << Tumour -> at(ith_clone) -> Clone_Size << std::endl;

			 	counter++;
			 	//getchar();

			 	//mut_effect_flag = true;

			}
			else if(mutational_effect_buffer[DRIVER_CELLS] > 0)// Driver
			{
				std::cout << " DRIVER " << std::endl;
				carcinogenesis_from_driver( ith_clone,  years,  hours);
				
				//getchar();
			}
			else if(mutational_effect_buffer[0] > 0)
			{
				std::cout << " DEATH " << std::endl;
				Check_for_Clonal_Extintion_at_min_size( ith_clone );
				death++;
				//getchar();
				// if(Tumour -> at(ith_clone) -> clone_extinct)
				// 	break;
			}
			i++;
		}//innner for
		j++;
	}// outer loop

	// VERY CAREFUL WITH THIS
	
	effect_gain = effect_gain -  Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] ;
	std::cout << "AFT Gain " << effect_gain << " counter " << counter << std::endl;
	std::cout << " Final PR " << Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] << std::endl;
	effect_gain = fabs(effect_gain);


	//unsigned int MCC = Mutant_cells ;
	
	//Mutations[0] = 0; Mutations[1] = 0; Mutations[2] = 0;

	//if (!Tumour -> at(ith_clone) -> clone_extinct )
	//{

		//Mutations = 
		r.Mutational_Proportions( effect_gain, Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] , Mutant_cells, Mutations );
	
	//}
	

	std::cout << "[ Muts  =  " << Mutant_cells<< " ] "
							 << "P(G0) = "  << effect_gain
							 << " ; P(G1) = " << Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR]
							 << " ; P(I) = " << 1.0 -(effect_gain + Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR])
							 << " ; to G0 = " << Mutations[0]
							 << " ; to G1 = " << Mutations[1]
							 << " ; IDLE  = " << Mutations[2]
							 << " ; Deaths = " << death << std::endl;

	// std::cout << "Values Muts: " << "Muts(G0) " << Mutations[MUTANTS_to_G0]
	// 			  << "  Muts(G1) " << Mutations[MUTANTS_to_G1]
	// 			  << " Muts(IDLE) " << Mutations[MUTANTS_to_IDLE]
	// 			  << " Muts (Death) " << death 
	// 			  << " Mutant cells "<< Mutant_cells 
	// 			  << " effect gain " <<  effect_gain << std::endl;



	

}

void Clonal_Expansion::Probabilities_of_Cell_Division (const unsigned int & ith_clone, double & p_idle, double & p_go_to_G0 )
{
	p_go_to_G0 = fabs( Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR] - Tumour -> at( ith_clone ) -> max_PR);
	p_idle = 0.62 * (1.0 - ( Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR] + p_go_to_G0 ) );
	

}


void Clonal_Expansion::Complementary_Probability( double & comp_probability_G0, 
												  double & comp_probability_G1, 
												  double & comp_probability_G2,
												  double & comp_probability_S,
												  double & comp_probability_M,
												  const unsigned int & ith_clone )
{

	comp_probability_G0 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) ); 
	comp_probability_G1 = 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) ) - 0.1;
	comp_probability_G2 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status ) );
	comp_probability_S  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status ) );
	comp_probability_M  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status ) );

}


// Update G0
void Clonal_Expansion::Update_G0_Phase( const double & comp_probability_G0, const unsigned int & ith_clone)
{
	unsigned int * G0_buffer = r.Update_G0_Phase( Tumour -> at( ith_clone ) -> G0_cells, 
	 							std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ), 
	 							std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G0_status ),
	 							comp_probability_G0);

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[DIE];// Die 
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[NEXT_STAGE];// to G1
}

//Update G1
void Clonal_Expansion::Update_G1_Phase(const double & comp_probability_G1, const unsigned int & ith_clone)
{
	
	unsigned int * G1_buffer = r.Update_G1_Phase( Tumour -> at( ith_clone ) -> G1_cells,
	 							std::get<P_STAYING_G1> ( Tumour -> at( ith_clone ) -> G1_status ),
	 							std::get<P_DYING_G1>   ( Tumour -> at( ith_clone ) -> G1_status),
	 							comp_probability_G1);

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[DIE];// to Dying
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[NEXT_STAGE];// to G2
	std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[REVERT_TO_G0];// to G0
}

// Update G2
void Clonal_Expansion::Update_G2_Phase( const double & comp_probability_G2, const unsigned int & ith_clone)
{
	unsigned int * G2_buffer = r.Update_G2_Phase( Tumour -> at( ith_clone ) -> G2_cells,
	  							   std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ),
	  							   std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G2_status ),
	  							   comp_probability_G2 );

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[DIE];// to Dying
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[NEXT_STAGE];
}

//Update S
void Clonal_Expansion::Update_S_Phase( const double & comp_probability_S, const unsigned int & ith_clone)
{
	unsigned int *  S_buffer = r.Update_S_Phase( Tumour -> at( ith_clone ) -> S_cells,
	  							std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ),
	  							std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> S_status ),
	  							comp_probability_S );

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[DIE];// die
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[NEXT_STAGE];// to M

}

//Update M

void Clonal_Expansion::Update_M_Phase( const double & comp_probability_M, const unsigned int & ith_clone)
{
	unsigned int *  M_buffer = r.Update_M_Phase( Tumour -> at( ith_clone ) -> M_cells,
	  							 std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status),
	  							 std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> M_status ),
	  							 comp_probability_M );

	 std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[STAY];// Stay
	 std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[DIE];// Dy
	 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[NEXT_STAGE];// to Division Model
}



//This can run in parallel 
// Mitosis
// https://www.youtube.com/watch?v=koudmJdil60
void Clonal_Expansion::Check_Mitosis_Network_Status(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours )
{
	// All of this runs
	//Random r;
	


	double comp_probability_G0 = 0;
	double comp_probability_G1 = 0;
	double comp_probability_G2 = 0;
	double comp_probability_S = 0;
	double comp_probability_M = 0;
	double p_idle = 0.0;
	double p_go_to_G0 = 0.0;

	std::vector<unsigned int> Division_Model;

	//std::vector<unsigned int> Division_Model = 

	//unsigned int * Mutations;

	std::vector<unsigned int> Mutations;

	Complementary_Probability( comp_probability_G0, 
							   comp_probability_G1, 
							   comp_probability_G2,
							   comp_probability_S,
							   comp_probability_M,
							   ith_clone );

	// // // for G0
	// // double comp_probability_G0 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) ); 
	// // double comp_probability_G1 = 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) ) - 0.1;
	// // double comp_probability_G2 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status ) );
	// // double comp_probability_S  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status ) );
	// // double comp_probability_M  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status ) );
	
	// std::cout << "STATS G0 >  " <<  Tumour -> at( ith_clone ) -> G0_cells << std::endl;

	 Update_G0_Phase( comp_probability_G0, ith_clone);
	 Update_G1_Phase( comp_probability_G1, ith_clone);
	 Update_G2_Phase( comp_probability_G2, ith_clone);
	 Update_S_Phase(  comp_probability_S,  ith_clone);
	 Update_M_Phase(  comp_probability_M, ith_clone);

	 unsigned int AC_at_t = Tumour -> at( ith_clone ) -> Available_cells;
	 unsigned int * Clonal_update_buffer = r.Clonal_Functions( Tumour -> at( ith_clone ) -> Available_cells,
	   									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR],
	   									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR] - feedback, 
	   									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR] ) ;
	
	 Probabilities_of_Cell_Division ( ith_clone, p_idle, p_go_to_G0 );

	 r.Newborn( std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ), 
	 										   p_idle, 
	 										   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR],
	 										   p_go_to_G0,
	 										   Division_Model
	 										  );

	 

	


	//  // unsigned int * G0_buffer = r.Update_G0_Phase( Tumour -> at( ith_clone ) -> G0_cells, 
	//  // 							  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ), 
	//  // 							  std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G0_status ),
	//  // 							  comp_probability_G0);

	//  // unsigned int * G1_buffer = r.Update_G1_Phase( Tumour -> at( ith_clone ) -> G1_cells,
	//  // 							   std::get<P_STAYING_G1> ( Tumour -> at( ith_clone ) -> G1_status ),
	//  // 							   std::get<P_DYING_G1>   ( Tumour -> at( ith_clone ) -> G1_status),
	//  // 							   comp_probability_G1);

	//   // unsigned int * G2_buffer = r.Update_G2_Phase( Tumour -> at( ith_clone ) -> G2_cells,
	//   // 							   std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ),
	//   // 							   std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G2_status ),
	//   // 							   comp_probability_G2 );

	//   // unsigned int *  S_buffer = r.Update_S_Phase( Tumour -> at( ith_clone ) -> S_cells,
	//   // 							  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ),
	//   // 							  std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> S_status ),
	//   // 							  comp_probability_S );

	//  // unsigned int *  M_buffer = r.Update_M_Phase( Tumour -> at( ith_clone ) -> M_cells,
	//  //  							 std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status),
	//  //  							 std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> M_status ),
	//  //  							 comp_probability_M );


	//  // TODO Change this
	//  //double to_G0 = 0.001; // probability of expoansion 
	//  // unsigned int AC_at_t = Tumour -> at( ith_clone ) -> Available_cells;
	//  // unsigned int * Clonal_update_buffer = r.Clonal_Functions( Tumour -> at( ith_clone ) -> Available_cells,
	//  //  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR],
	//  //  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR], 
	//  //  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR] ) ;

	

	// // //for G0
	// // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[STAY];// Stay
	// // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[DIE];// Die 
	// // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[NEXT_STAGE];// to G1

	// std::cout << "ASSIGN G0 >  " <<  Tumour -> at( ith_clone ) -> G0_cells << std::endl;

	// //for G1
	// // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[STAY];// Stay
	// // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[DIE];// to Dying
	// // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[NEXT_STAGE];// to G2
	// // std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[REVERT_TO_G0];// to G0
	// //for G2
	//  // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[STAY];// Stay
	//  // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[DIE];// to Dying
	//  // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[NEXT_STAGE];
	// // //for S
	//  // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[STAY];// Stay
	//  // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[DIE];// die
	//  // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[NEXT_STAGE];// to M
	// // //for M
	//  // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[STAY];// Stay
	//  // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[DIE];// Dy
	//  // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[NEXT_STAGE];// to Division Model

	
	std::cout << "[ G0 =   " << Tumour -> at( ith_clone ) -> G0_cells << " ] "
							 << "P(st) = " << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; P(Dy) = " << std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; P(G1) = " << comp_probability_G0  
							 << " ; G0(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; G0(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; G0(Trans) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) << std::endl;


	std::cout << "[ G1  =  " << Tumour -> at ( ith_clone ) -> G1_cells << " ] "
							 << "P(st) = "  << std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) ->G1_status ) 
							 << " ; P(Dy) = " << std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status )
							 << " ; P(tra) = "<< comp_probability_G1
							 << " ; P(toG0) = " << 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) ->G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + comp_probability_G1)
							 << " ; G1(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) 
							 << " ; G1(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) 
							 << " ; G1(Tr) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) 
							 << " ; G1(G0) = " << std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) << std::endl;

	std::cout << "[ G2 = "   << Tumour -> at ( ith_clone ) -> G2_cells << " ] "
							 << " P(St) = "  << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) 
							 << " ; P(Dy) = "  << std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status )
							 << " ; P(Tra) = " << comp_probability_G2
							 << " ; G2(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) 
							 << " ; G2(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) 
							 << " ; G2(Tr) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) << std::endl;

	std::cout << "[ S =  " << Tumour -> at ( ith_clone ) -> S_cells << " ] "
						   << " P(St) = "  << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ) 
						   << " ; P(Dy) = " << std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status )
						   << " ; P(Tra) = "<< comp_probability_S
						   << " ; S(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) 
						   << " ; S(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) 
						   << " ; S(Tr) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) << std::endl;

	std::cout << "[ M = " << Tumour -> at( ith_clone ) -> M_cells << " ] " 
	 					  << " P(St) = "  << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status ) 
						  << " ; P(Dy) = " << std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status )
						  << " ; P(Tra) = "<< comp_probability_M
						  << " ; M(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) ->M_status ) 
						  << " ; M(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) 
						  << " ; M(Tr) =  " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) << std::endl;


	std::cout << "[ DV = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) << " ] "
			  << "; PDiv(Idle) = "   << p_idle 
			  << "; P(G0) = " << p_go_to_G0
			  << " ; P(MR) = "       <<  Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR]
			  << " ; P(G1) = " 		<< 1.0 - ( p_idle + p_go_to_G0 + Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR])
			  << " : sum = " << ( 1.0 - ( p_idle + p_go_to_G0 + Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR])) + ( p_idle + p_go_to_G0 + Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR])
			  << " ; Idling = "   <<  Division_Model[DIV_CELLS_IDLING]
			  << " ; Muts = "     << Division_Model[DIV_CELLS_MUTATING]
			  << " ; To G0 = "    << Division_Model[DIV_CELLS_to_G0]
			  << " ; To G1 = "    << Division_Model[DIV_CELLS_to_G1]
			  << std::endl;


	std::cout << "[ AC(t) = " << AC_at_t << " ] => "
							   << " ; PR = " << Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR]
							   << " ; penalty = " << feedback
							   << " ; DR = " << Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR]
							   << " ; QR = " << Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR]
	 						   << " ; Idling = " <<  Clonal_update_buffer[AV_CELLS_IDLING]
							   << " ; Dying = " << Clonal_update_buffer[AV_CELLS_DYING]
							   << " ; To G0 = " << Clonal_update_buffer[AV_CELLS_to_G0]
							   << " ; To G1 = " << Clonal_update_buffer[AV_CELLS_to_G1] << std::endl;

	

	
	
	// //Send this as part of the model and then soft  max 
	

	

	//std::cout << "[ DV = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status << " ] "
					    //    << "PDiv(Idle): "   << p_idle 
					    //    << " ; PDiv(G0) = " << p_go_to_G0 
					    //    << " ; MR = "       <<  Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR] 
						   // << " ; Idling = "   <<  Division_Model[DIV_CELLS_IDLING]
						   // << " ; Muts = "     << Division_Model[DIV_CELLS_MUTATING]
						   // << " ; To G0 = "    << Division_Model[DIV_CELLS_to_G0]
						   // << " ; To G1 = "    << Division_Model[DIV_CELLS_to_G1] 
						   //<< std::endl;


	// // //update mutations
	// // unsigned int Mutant_cells = 10; 
	// // //Estimate_Mutational_Effects( Mutant_cells, ith_clone -> Clone_Size );



	
	//if( Division_Model[DIV_CELLS_MUTATING] >  0)
	//{
	 	Update_Clonal_Mutational_Burden( ith_clone, Division_Model[DIV_CELLS_MUTATING], years, hours, Mutations);
	//}

	// //
	// // UPDATE Values
	// //

	std::cout << " G0(t) = " <<  Tumour -> at( ith_clone ) -> G0_cells << " to ";

	Tumour -> at( ith_clone ) -> G0_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) + 
											std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) + 
											Division_Model[DIV_CELLS_to_G0] + 
											Mutations[MUTANTS_to_G0] + 
											Clonal_update_buffer[AV_CELLS_to_G0];

	std::cout << " G0(t + dt) = " <<  Tumour -> at( ith_clone ) -> G0_cells  << std::endl;
											
	std::cout << " G1(t) = " <<  Tumour -> at( ith_clone ) -> G1_cells << " to "  ;

	Tumour -> at ( ith_clone ) -> G1_cells =  std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) + 
											  std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) + 
											  Division_Model[DIV_CELLS_to_G1] + 
											  Clonal_update_buffer[AV_CELLS_to_G1] + 
											  Mutations[MUTANTS_to_G1] ;
	
	std::cout << " G1(t + dt) = " <<  Tumour -> at( ith_clone ) -> G1_cells  << std::endl;

	std::cout << " S(t) = " <<  Tumour -> at( ith_clone ) -> S_cells  << " to ";

	Tumour -> at ( ith_clone ) -> S_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status )  +
											std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) ;

	std::cout << " S(t + dt) = " <<  Tumour -> at( ith_clone ) -> S_cells  << std::endl;

	std::cout << " G2(t) = " <<  Tumour -> at( ith_clone ) -> G2_cells  << " to " ;

	Tumour -> at ( ith_clone ) -> G2_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) +
											 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) ;

	std::cout << " G2(t + dt) = " <<  Tumour -> at( ith_clone ) -> G2_cells  << std::endl;

	std::cout << " M(t) = " <<  Tumour -> at( ith_clone ) -> M_cells  << " to " ;

	Tumour -> at ( ith_clone ) -> M_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) +
											std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) ;

	std::cout << " M(t + dt) = " <<  Tumour -> at( ith_clone ) -> M_cells  << std::endl;

	std::cout << " AC(t ) = " <<  Tumour -> at( ith_clone ) -> Available_cells  << " to ";

	Tumour -> at( ith_clone ) -> Available_cells = Clonal_update_buffer[AV_CELLS_IDLING] +
												   Division_Model[DIV_CELLS_IDLING] +
												   Mutations[MUTANTS_to_IDLE] ;

	std::cout << " AC(t + dt) = " <<  Tumour -> at( ith_clone ) -> Available_cells  << std::endl;

	Tumour -> at( ith_clone ) -> Clone_Size = static_cast<unsigned long long int>(Tumour -> at(ith_clone) -> Available_cells + 
	 																			  Tumour -> at(ith_clone) -> G0_cells        + 
	 																			  Tumour -> at(ith_clone) -> G1_cells        + 
	 																			  Tumour -> at(ith_clone) -> G2_cells        + 
	 																			  Tumour -> at(ith_clone) -> S_cells         + 
	 																			  Tumour -> at(ith_clone) -> M_cells);

	std::cout << " Clone_S( " <<  ith_clone << " ) = "<< Tumour -> at( ith_clone ) -> Clone_Size  << std::endl;
	
}




void Clonal_Expansion::map_Feedback(  )
{

	feedback =  0.0 + (DIFF - 0.0) * (( (double) Population_Size - 0.0) / ((double) PS - 0.0));	
	//std::cout << "Population Feedback " << feedback << " = " << DIFF << " * " << Population_Size << " / " << PS << std::endl;
}


void Clonal_Expansion::Update_Tumour_Size()
{
	unsigned int ith_clone = 0;
	unsigned long long int population = 0;
	for (ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++) 
	{
		population += Tumour -> at(ith_clone) -> Clone_Size;
	}
	Population_Size = population;

	//std::cout << "Tumour SIZE = " <<  Population_Size << std::endl;
}


// void Clonal_Expansion::Non_Mutagenic_Mitosis_Standard(std::unique_ptr<Clone> & ith_clone)
// {

// 	if( check_G1_phase(ith_clone) )
// 		ith_clone -> Remaining_Time_in_G1_Phase--;
// 	else if( exiting_G1_phase(ith_clone) )
// 		Transition_From_G1_S(ith_clone);
// 	else if( check_G1_phase(ith_clone) )
// 		ith_clone -> Remaining_Time_in_S_Phase--;
// 	else if( exiting_S_phase(ith_clone) )
// 		Transition_From_S_G2(ith_clone);
// 	else if( check_G2_phase(ith_clone) )
// 		ith_clone -> Remaining_Time_in_G2_Phase--;
// 	else if( exiting_G2_phase(ith_clone) )
// 		Transition_From_G2_M(ith_clone);
// 	else if( ith_clone -> In_M_Phase)
// 	{
// 		std::cout << "SEND TO DIVISION MODEL " << std::endl;
// 		Update_Population_After_Division(ith_clone);
// 	}
// }

void Clonal_Expansion::Reset_to_G1(const unsigned int & ith_clone)
{
	Tumour -> at (ith_clone) -> In_M_Phase = false;
	Tumour -> at (ith_clone) -> In_G1_Phase = true;
}

void Clonal_Expansion::Valid_Size_Heterogeneity_V1(const unsigned int & ith_clone)
{
	if( Tumour -> at(ith_clone) -> Clone_Size >= Tumour -> at(ith_clone) -> Number_of_Memebers_to_Start_Heterogeneity)
	{
		Tumour -> at(ith_clone) ->  Initiall_Expasion_Period = false;		
	}
}


// V1 Final version

void Clonal_Expansion::Non_Mutagenic_Mitosis_Standard_V1(const unsigned int & ith_clone)
{
	if( check_G1_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone);
	else if( check_S_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase--;	
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone);
	else if( check_G2_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone);
	else if(Tumour -> at (ith_clone) -> In_M_Phase)
	{
		//std::cout << "SEND TO DIVISION MODEL " << std::endl;
		Update_Population_After_Division_V1(ith_clone);
		Valid_Size_Heterogeneity_V1(ith_clone);
	}
	else
	{
		std::cout << "ERROR " << std::endl;
		getchar();
	}



}

void Clonal_Expansion::Checking_for_Mutants( std::tuple<unsigned long long int, // Dying
													   unsigned long long int, // Newborn
													   unsigned long long int, // Mutants
													   bool, bool> & Dying_and_Newborn, 
													   const unsigned int & ith_clone )
{
	std::get<2>(Dying_and_Newborn) = 0;

	if(std::get<1>(Dying_and_Newborn) > 0)
	{	
		r.Binomial_Mutant( Dying_and_Newborn, Tumour -> at(ith_clone) -> Mutation_Rate );
		
		std::get<3>(Dying_and_Newborn) =  true;// flag that indicates that there are mutant cells
	
		if(std::get<2>(Dying_and_Newborn) > std::get<1>(Dying_and_Newborn))
		{
			std::cout << "More mutant cells than newborn? " << " MC: " << std::get<2>(Dying_and_Newborn)  << " NB: " << std::get<1>(Dying_and_Newborn) << std::endl;
		}
		else
		{
			std::get<1>(Dying_and_Newborn) -= std::get<2>(Dying_and_Newborn);
		}
	}

}

void Clonal_Expansion::Update_Newborn_Parameters_V1(const std::vector<unsigned long long int> & NewBorn_Cells, std::tuple<unsigned long long int, 
																														  unsigned long long int, 
																														  unsigned long long int, 
																														  bool, bool> & Dying_and_Newborn)
{

	std::get<0>(Dying_and_Newborn) = NewBorn_Cells.at(0);// For Dying Clones 		
	std::get<1>(Dying_and_Newborn) = 2 * NewBorn_Cells.at(1); // For Newborn Clones 

}

void Clonal_Expansion::Update_Newborn_Parameters_V2( const std::vector<unsigned long long int> & NewBorn_Cells, std::tuple<unsigned long long int, 
																															unsigned long long int,
																															unsigned long long int,
																															bool, bool> & Dying_and_Newborn)
{
	std::get<0>(Dying_and_Newborn) = NewBorn_Cells.at(0);// For Dying Clones 		
	std::get<1>(Dying_and_Newborn) = 2 * NewBorn_Cells.at(1); // For Newborn Clones 
	std::get<2>(Dying_and_Newborn) = 0;

	if( NewBorn_Cells.at(2) > 0 )// Check for mutants
	{
		std::get<2>(Dying_and_Newborn) = NewBorn_Cells.at(2);
		std::get<3>(Dying_and_Newborn) = true;
	
	}//end if

}

// Ad dtom ain structure
void Clonal_Expansion::Update_Newborn_Parameters_V2R( const std::vector<unsigned long long int> & NewBorn_Cells, std::tuple<unsigned long long int, 
																															unsigned long long int,
																															unsigned long long int,
																															bool, bool> & Dying_and_Newborn)
{
	std::get<0>(Dying_and_Newborn) = NewBorn_Cells.at(0);// For Dying Clones 		
	std::get<1>(Dying_and_Newborn) = 2 * NewBorn_Cells.at(1); // For Newborn Clones 
	std::get<2>(Dying_and_Newborn) = 0;

	if( NewBorn_Cells.at(2) > 0 )// Check for mutants
	{
		std::get<2>(Dying_and_Newborn) = NewBorn_Cells.at(2);
		std::get<3>(Dying_and_Newborn) = true;
	
	}//end if


}

void Clonal_Expansion::Adjust_Feedback_PR(double & P_NB)
{
	if(feedback >= P_NB)
		P_NB = 0.0;
	else
		P_NB -= feedback;
}


void Clonal_Expansion::Check_Clonal_Extinction_BCE_V1(const unsigned int & ith_clone, 
													  std::tuple<unsigned long long int, // Dying 
													  			 unsigned long long int, // Newborn
													  			 unsigned long long int, // Mutant
													  			 bool, bool> & Dying_and_Newborn)
{
	bool extiction =  ( Tumour -> at(ith_clone) -> Clone_Size + std::get<1>(Dying_and_Newborn) ) < ( std::get<0>(Dying_and_Newborn) + std::get<2>(Dying_and_Newborn) ) ;
	if( extiction )
	{
		Tumour -> at(ith_clone) -> clone_extinct = true;
		std::get<4>(Dying_and_Newborn) =  false;
		Population_Size -=  Tumour  -> at(ith_clone) -> Clone_Size ;
		Tumour -> at(ith_clone) -> Clone_Size = 0;
	}

}

void Clonal_Expansion::Check_Clonal_Extinction_BCE_V2(const unsigned int & ith_clone, 
													  std::tuple<unsigned long long int, // Dying 
													  			 unsigned long long int, // Newborn
													  			 unsigned long long int, // Mutant
													  			 bool, bool> & Dying_and_Newborn)
{
	bool extiction =  ( Tumour -> at(ith_clone) -> Clone_Size + std::get<1>(Dying_and_Newborn) ) < ( std::get<0>(Dying_and_Newborn) + std::get<2>(Dying_and_Newborn) ) ;
	if( extiction )
	{
		Tumour -> at(ith_clone) -> clone_extinct = true;
		std::get<4>(Dying_and_Newborn) =  false;
		Population_Size -=  Tumour  -> at(ith_clone) -> Clone_Size ;
		Tumour -> at(ith_clone) -> Clone_Size = 0;
	}

}







void Clonal_Expansion::Print_Debug_Carcinogeneisis_From_Driver(const unsigned int & ith_clone)
{
	std::cout << " TUMOUR STATUS " << std::endl;
	printParameters();

	std::cout << "ID: " << Tumour -> at(ith_clone) -> Generation_ID << " to " <<Tumour -> back() -> Generation_ID << std::endl;
	std::cout << "CS: " << Tumour -> at(ith_clone) -> Clone_Size << " to " <<Tumour -> back() -> Clone_Size << std::endl;
	//std::cout << "NOM: " << Tumour -> at(ith_clone) -> Number_of_Mutations << " to " <<Tumour -> back() -> Number_of_Mutations << std::endl;
	std::cout << "MR: " <<  Tumour -> at(ith_clone) -> Mutation_Rate << " to " <<Tumour -> back() -> Mutation_Rate << std::endl;
	std::cout << "PR: " << Tumour -> at(ith_clone) -> P_Expansion[1] << " to " <<Tumour -> back() -> P_Expansion[1] << std::endl;

	//getchar();


}


void Clonal_Expansion::Carcinogenesis_From_Driver_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	//double mr = Tumour -> at(ith_clone) -> Mutation_Rate;
	//unsigned long long int NOM = Tumour -> at (ith_clone) -> Number_of_Mutations;
	double Updated_pr = 0.0; 
	double Updated_mr = 0.0;
	//double pr = Tumour -> at(ith_clone) -> P_Expansion[1]; // with back gaves good results
	
	std::string cloneName = "";	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) -> Generation_ID_Counter++;
	Generate_Clone_Generation_ID(cloneName, 
								 Parent_Generation_ID_Counter, 
								 Tumour -> at( ith_clone ) -> Generation_ID,
								 years, 
								 hours);
	//Tumour -> at(ith_clone) -> Generation_ID_Counter++;
	//std::cout << "NID: " << cloneName << std::endl;

	Tumour -> push_back( get_Clone_DS() );
	Tumour -> back() -> Generation_ID = cloneName;
	
	Tumour -> back() -> Initiall_Expasion_Period = false;//true
	Tumour -> back() -> Clone_Size = 1;
	//Tumour -> back() -> Number_of_Mutations = NOM + 1;
	
	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	
	Select_MR_Sampling(Tumour -> at(ith_clone) -> Mutation_Rate, Updated_mr);
	Tumour -> back() -> Mutation_Rate = Updated_mr;
	//r.Uniform_Mutation_Rate_2(mr);
	
	Select_PR_Samplig(Tumour -> at(ith_clone) -> P_Expansion[1], Updated_pr);
	Tumour -> back() -> P_Expansion[1] =  Updated_pr;

	//r.Update_Proliferation_Rate(pr);
	
	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
	Population_Size++;
	
	//Print_Debug_Carcinogeneisis_From_Driver( ith_clone );

}



void Clonal_Expansion::Carcinogenesis_From_Driver_V2(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	//double mr = Tumour -> at(ith_clone) -> Mutation_Rate;
	//unsigned long long int NOM = Tumour -> at (ith_clone) -> Number_of_Mutations;
	double Updated_pr = 0.0; 
	double Updated_mr = 0.0;
	//double pr = Tumour -> at(ith_clone) -> P_Expansion[1]; // with back gaves good results
	unsigned long long int number_of_mutations = Tumour -> at(ith_clone) -> Clonal_Mutations;
	double parten_mut_burden = Tumour -> at(ith_clone) -> Clonal_Mutational_Burden;
	std::string cloneName = "";	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) -> Generation_ID_Counter;
	Generate_Clone_Generation_ID(cloneName, 
								 Parent_Generation_ID_Counter, 
								 Tumour -> at( ith_clone ) -> Generation_ID,
								 years, 
								 hours);
	//std::cout << "NID: " << cloneName << std::endl;
	Tumour -> at(ith_clone) -> Generation_ID_Counter++;
	Tumour -> push_back( get_Clone_DS() );
	Tumour -> back() -> Generation_ID = cloneName;
	Tumour -> back() -> Initiall_Expasion_Period = false;//true
	Tumour -> back() -> Clone_Size = 1;
	//Tumour -> back() -> Number_of_Mutations = NOM + 1;
	
	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	
	Select_MR_Sampling(Tumour -> at(ith_clone) -> Mutation_Rate, Updated_mr);
	Tumour -> back() -> Mutation_Rate = Updated_mr;
	Tumour -> back() -> P_Expansion[2] = Updated_mr;
	Tumour -> back() -> Clonal_Mutations = number_of_mutations;
	
	//r.Uniform_Mutation_Rate_2(mr);
	
	//Select_PR_Samplig(Tumour -> at(ith_clone) -> P_Expansion[1], Updated_pr);
	r.Update_Proliferation_Rate_V2( Tumour -> at(ith_clone) -> P_Expansion[1], Updated_pr ); // not subject to variance chsange
	Tumour -> back() -> P_Expansion[1] =  Updated_pr;
	Tumour -> back() -> Clonal_Mutational_Burden = parten_mut_burden + Updated_pr;
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[1];
	Tumour -> back() -> init_PR = Tumour -> back() -> P_Expansion[1];

	//r.Update_Proliferation_Rate(pr);
	

	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
	Population_Size++;
	
	//Print_Debug_Carcinogeneisis_From_Driver( ith_clone );

}

void Clonal_Expansion::Carcinogenesis_From_Drug_Resistance_Growth(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	//double mr = Tumour -> at(ith_clone) -> Mutation_Rate;
	//unsigned long long int NOM = Tumour -> at (ith_clone) -> Number_of_Mutations;
	double Parent_PR = Tumour -> at(ith_clone) -> P_Expansion[1]; 
	double Parent_MR = Tumour -> at(ith_clone) -> P_Expansion[2];
	//double pr = Tumour -> at(ith_clone) -> P_Expansion[1]; // with back gaves good results
	double strength = 0.0;
	unsigned long long int number_of_mutations = Tumour -> at(ith_clone) -> Clonal_Mutations;
	double mut_burden = Tumour -> at(ith_clone) -> Clonal_Mutational_Burden;
	std::string cloneName = "";	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) -> Generation_ID_Counter;
	Generate_Clone_Generation_ID(cloneName, 
								 Parent_Generation_ID_Counter, 
								 Tumour -> at( ith_clone ) -> Generation_ID,
								 years, 
								 hours);
	//std::cout << "NID: " << cloneName << std::endl;
	Tumour -> at(ith_clone) -> Generation_ID_Counter++;
	Tumour -> push_back( get_Clone_DS() );
	Tumour -> back() -> Generation_ID = cloneName;
	Tumour -> back() -> Initiall_Expasion_Period = false;//true
	Tumour -> back() -> Clone_Size = 1;
	//Tumour -> back() -> Number_of_Mutations = NOM + 1;
	
	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Generation_ID_Counter = 0; 

	Tumour -> back() -> Drug_Resistant = true;
	r.Drug_Resistance_Strength(strength);
	Tumour -> back() -> drug_resistance_strength = strength;//
	Tumour -> back() -> resistance_in_treatment = false;
	Tumour -> back() -> year_of_resistance = years;
	Tumour -> back() -> hour_of_resistance = hours;
	
	//Select_MR_Sampling(Tumour -> at(ith_clone) -> Mutation_Rate, Updated_mr);
	Tumour -> back() -> Mutation_Rate = Parent_MR;
	Tumour -> back() -> P_Expansion[2] = Parent_MR;
	Tumour -> back() -> Clonal_Mutations = number_of_mutations;
	Tumour -> back() -> Clonal_Mutational_Burden = mut_burden;
	//r.Uniform_Mutation_Rate_2(mr);
	
	//Select_PR_Samplig(Tumour -> at(ith_clone) -> P_Expansion[1], Updated_pr);
	//r.Update_Proliferation_Rate_V2( Tumour -> at(ith_clone) -> P_Expansion[1], Updated_pr );
	Tumour -> back() -> P_Expansion[1] =  Parent_PR;
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[1];

	//r.Update_Proliferation_Rate(pr);
	

	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
	Population_Size++;
	
	//Print_Debug_Carcinogeneisis_From_Driver( ith_clone );

}


void Clonal_Expansion::Mutant_Effects_V1(const unsigned int & ith_clone, const unsigned long long int & Mutant_Cells,  const unsigned int & hours, const unsigned int & years)
{
	
	unsigned int i =0;
	for( i = 0; i < Mutant_Cells; i++)
		Carcinogenesis_From_Driver_V1( ith_clone, hours, years);
}

void Clonal_Expansion::Apply_Penalties_to_Mutants_V1(const unsigned int & ith_clone, const std::tuple<unsigned long long int, 
																									  unsigned long long int, 
																									  unsigned long long int, 
																									  bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years)
{
	
	if( std::get<3>(Dying_and_Newborn) && (std::get<2>(Dying_and_Newborn) > 0) )
	{
		Mutant_Effects_V1(ith_clone, std::get<2>(Dying_and_Newborn), hours, years);
	}

	unsigned long long int New_poulation_Size = ( Tumour -> at(ith_clone) -> Clone_Size  - std::get<0>(Dying_and_Newborn) ) 
										      + 	std::get<1>(Dying_and_Newborn) ;

		// Update population size by
		// Pop size =Pop_Size + (CS[t] - CS[t-1]) 
		Population_Size +=  New_poulation_Size - Tumour  -> at(ith_clone) -> Clone_Size ;
		Tumour -> at(ith_clone) -> Clone_Size = New_poulation_Size ;
}


//Mutant effects
// void Clonal_Expansion::Mutational_Effects_V2(const unsigned int & ith_clone, const std::tuple<unsigned long long int,
// 																							  unsigned long long int,
// 																							  unsigned long long int,
// 																							  bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years)
// {

// }




void Clonal_Expansion::Apply_Penalties_to_Mutants_V2(const unsigned int & ith_clone,  std::tuple<unsigned long long int, 
																								 unsigned long long int, 
																								 unsigned long long int, 
																								 bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years)
{


	if( std::get<3>(Dying_and_Newborn) && (std::get<2>(Dying_and_Newborn) > 0) )
	{
		//std::cout << "To mutational Effects We have " << std::get<2>(Dying_and_Newborn) << " mutant cells " << std::endl;
		Select_Mutational_Effects_Distribution(ith_clone, hours, years, Dying_and_Newborn);
		
		
	}

	unsigned long long int New_poulation_Size = ( Tumour -> at(ith_clone) -> Clone_Size  - std::get<0>(Dying_and_Newborn) ) 
										      + 	(std::get<1>(Dying_and_Newborn) +  std::get<2>(Dying_and_Newborn) ) ; // MR is reduced by drivers and killers

	// Update population size by
	// Pop size =Pop_Size + (CS[t] - CS[t-1]) 
	Population_Size +=  New_poulation_Size - Tumour  -> at(ith_clone) -> Clone_Size ;
	Tumour -> at(ith_clone) -> Clone_Size = New_poulation_Size ;



}


void Clonal_Expansion::Apply_Penalties_to_Mutants_V2R(const unsigned int & ith_clone,  std::tuple<unsigned long long int, 
																								 unsigned long long int, 
																								 unsigned long long int, 
																								 bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years)
{


	if( std::get<3>(Dying_and_Newborn) && (std::get<2>(Dying_and_Newborn) > 0) )
	{
		//std::cout << "To mutational Effects We have " << std::get<2>(Dying_and_Newborn) << " mutant cells " << std::endl;
		Select_Mutational_Effects_Distribution(ith_clone, hours, years, Dying_and_Newborn);
	}

	unsigned long long int New_poulation_Size = ( Tumour -> at(ith_clone) -> Clone_Size  - std::get<0>(Dying_and_Newborn) ) 
										      + 	(std::get<1>(Dying_and_Newborn) +   std::get<2>(Dying_and_Newborn) ) ;

	// Update population size by
	// Pop size =Pop_Size + (CS[t] - CS[t-1]) 
	Population_Size +=  New_poulation_Size - Tumour  -> at(ith_clone) -> Clone_Size ;
	Tumour -> at(ith_clone) -> Clone_Size = New_poulation_Size ;



}

void Clonal_Expansion::Update_Mutant_Mutational_Effects_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years, const std::tuple<unsigned long long int, 
																																								 unsigned long long int, 
																																								 unsigned long long int, 
																																								 bool, bool> & Dying_and_Newborn)
{

	if(std::get<4>(Dying_and_Newborn))
		Apply_Penalties_to_Mutants_V1( ith_clone, Dying_and_Newborn, hours,  years);
		
}


void Clonal_Expansion::Update_Mutant_Mutational_Effects_V2(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years,  std::tuple<unsigned long long int, 
																																								 unsigned long long int, 
																																								 unsigned long long int, 
																																								 bool, bool> & Dying_and_Newborn)
{

	if(std::get<4>(Dying_and_Newborn))
	{
		Apply_Penalties_to_Mutants_V2( ith_clone,  Dying_and_Newborn,  hours,  years );	
	}
		
}

void Clonal_Expansion::Update_Mutant_Mutational_Effects_V2R(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years,  std::tuple<unsigned long long int, 
																																								 unsigned long long int, 
																																								 unsigned long long int, 
																																								 bool, bool> & Dying_and_Newborn)
{

	if(std::get<4>(Dying_and_Newborn))
	{
		//Apply_Penalties_to_Mutants_V2( ith_clone,  Dying_and_Newborn,  hours,  years );	
		Apply_Penalties_to_Mutants_V2R( ith_clone, Dying_and_Newborn, hours, years );
		//std::cout << "Drug resistance model -> to mutational effects  " << std::endl;
		//getchar();
	}
		
}

void Clonal_Expansion::Basic_Clonal_Expansion_V1(const unsigned int & ith_clone, std::tuple<unsigned long long int, 
																							unsigned long long int, 
																							unsigned long long int, bool, bool> & Dying_and_Newborn)
{
	std::vector<unsigned long long int> NewBorn_Cells;
	
	//1) Let's get basline parameters
	double P_DR =  Tumour -> at(ith_clone) -> P_Expansion[0];
	double P_NB =  Tumour -> at(ith_clone) -> P_Expansion[1];
	//2) Apply enviromental penalty to growth variables
	Adjust_Feedback_PR(P_NB);
	
	r.Basic_Clonal_Expansion_Sampling_V1_HPC(Tumour -> at (ith_clone) -> Clone_Size, P_DR, P_NB, NewBorn_Cells);

	//std::cout<< NewBorn_Cells.at(0) << " " << NewBorn_Cells.at(1) << " " << NewBorn_Cells.at(2) << " T: " << NewBorn_Cells.at(0) + NewBorn_Cells.at(1)+ NewBorn_Cells.at(2)  <<std::endl;

	Update_Newborn_Parameters_V1(NewBorn_Cells, Dying_and_Newborn);
	Checking_for_Mutants( Dying_and_Newborn, ith_clone );
	Check_Clonal_Extinction_BCE_V1( ith_clone, Dying_and_Newborn );
	
}

void Clonal_Expansion::Basic_Clonal_Expansion_V2(const unsigned int & ith_clone, std::tuple<unsigned long long int, 
																							unsigned long long int, 
																							unsigned long long int, bool, bool> & Dying_and_Newborn)
{
	std::vector<unsigned long long int> NewBorn_Cells;
	
	//1) Let's get basline parameters
	double P_DR =  Tumour -> at(ith_clone) -> P_Expansion[0];
	double P_NB =  Tumour -> at(ith_clone) -> P_Expansion[1];
	double P_MR =  Tumour -> at(ith_clone) -> P_Expansion[2];
	//2) Apply enviromental penalty to growth variables
	Adjust_Feedback_PR(P_NB);
	
	r.Basic_Clonal_Expansion_Sampling_V2_HPC(Tumour -> at (ith_clone) -> Clone_Size, P_DR, P_NB, P_MR, NewBorn_Cells);

	//std::cout<< "V2 .. " << NewBorn_Cells.at(0) << " " << NewBorn_Cells.at(1) << " " << NewBorn_Cells.at(2) << " " << NewBorn_Cells.at(3) <<" D P M O |" << " T: " << NewBorn_Cells.at(0) + NewBorn_Cells.at(1)+ NewBorn_Cells.at(2)  <<std::endl;
	Update_Newborn_Parameters_V2(NewBorn_Cells, Dying_and_Newborn);
	Check_Clonal_Extinction_BCE_V2( ith_clone, Dying_and_Newborn );

	// if(NewBorn_Cells.at(2) > 0)
	// 	getchar();
	
		
}

//dying sand newbporn
void Clonal_Expansion::Create_Drug_Resiatance_Clones(const std::vector<unsigned long long int> & NewBorn_Cells, 
													 const unsigned int & ith_clone, 
													 const unsigned int & hours, 
													 const unsigned int & years,
													 std::tuple<unsigned long long int, 
																unsigned long long int,
																unsigned long long int,
																bool, bool> & Dying_and_Newborn)
{
	if(  NewBorn_Cells.at(3) > 0  )
	{
		for(unsigned int i = 0; i < NewBorn_Cells.at(3) ; i++)// Create DR Clones
		{
			Carcinogenesis_From_Drug_Resistance_Growth( ith_clone, hours, years );
		}
		if (Tumour -> at(ith_clone) -> Clone_Size - NewBorn_Cells.at(3) <= 0)
		{
			Population_Size -= NewBorn_Cells.at(3);
			Tumour -> at(ith_clone) -> clone_extinct = true;
			Tumour -> at(ith_clone) -> final_PR = Tumour -> at (ith_clone) -> P_Expansion[1];
			std::get<4>(Dying_and_Newborn) = false;
			Tumour -> at(ith_clone) -> Clone_Size = 0;
		}
		else
		{
			Tumour -> at(ith_clone) -> Clone_Size -= NewBorn_Cells.at(3);
		}

		//Adjust Pop size
	}
}

void Clonal_Expansion::Basic_Clonal_Expansion_V2R(const unsigned int & ith_clone, std::tuple<unsigned long long int, 
																							unsigned long long int, 
																							unsigned long long int, bool, bool> & Dying_and_Newborn, 
																							const unsigned int & hours, const unsigned int & years)
{
	std::vector<unsigned long long int> NewBorn_Cells;
	
	//1) Let's get basline parameters
	double P_DR =  Tumour -> at(ith_clone) -> P_Expansion[0];
	double P_NB =  Tumour -> at(ith_clone) -> P_Expansion[1];
	double P_MR =  Tumour -> at(ith_clone) -> P_Expansion[2];
	//2) Apply enviromental penalty to growth variables
	Adjust_Feedback_PR(P_NB);
	r.Basic_Clonal_Expansion_Sampling_V2R_HPC(Tumour -> at (ith_clone) -> Clone_Size, P_DR, P_NB, P_MR, NewBorn_Cells);
	// if DR generate new clone
	//std::cout<< "DR .. " << NewBorn_Cells.at(0) << " " << NewBorn_Cells.at(1) << " " << NewBorn_Cells.at(2) << " " << NewBorn_Cells.at(3) << " " << NewBorn_Cells.at(4) <<" D P M R O |" << " T: " << NewBorn_Cells.at(0) + NewBorn_Cells.at(1)+ NewBorn_Cells.at(2)  <<std::endl;

	// if(NewBorn_Cells.at(3) >= 1)
	// {
	// 	std::cout << " Mutant " << NewBorn_Cells.at(3) << " Mutant " << std::endl;
	// 	getchar();

	// }
	Update_Newborn_Parameters_V2R(NewBorn_Cells, Dying_and_Newborn); 
	Create_Drug_Resiatance_Clones(  NewBorn_Cells, ith_clone, hours, years, Dying_and_Newborn );
	Check_Clonal_Extinction_BCE_V1( ith_clone, Dying_and_Newborn ); 
	

}



void Clonal_Expansion::Clonal_Extinction_From_Size(const unsigned int & ith_clone)
{
	if(Tumour -> at(ith_clone) -> Clone_Size == 0)
		Tumour-> at(ith_clone) -> clone_extinct = true;

}


void Clonal_Expansion::Dealyed_Mitosis_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	
	std::tuple<unsigned long long int, 
			   unsigned long long int, 
			   unsigned long long int, bool, bool> Dying_and_Newborn (0, 0, 0, false, true);

	if( check_G1_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone);
	else if( check_S_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase--;	
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone);
	else if( check_G2_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone);
	else if(Tumour -> at (ith_clone) -> In_M_Phase)
	{
		Basic_Clonal_Expansion_V1(ith_clone, Dying_and_Newborn);
		Update_Mutant_Mutational_Effects_V1(ith_clone,  hours, years, Dying_and_Newborn);
		Clonal_Extinction_From_Size(ith_clone);
		Reset_to_G1(ith_clone);
	}

}

//modify this ^ for version 2

void Clonal_Expansion::Dealyed_Mitosis_V2(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	std::tuple<unsigned long long int, 
			   unsigned long long int, 
			   unsigned long long int, bool, bool> Dying_and_Newborn (0, 0, 0, false, true);

	if( check_G1_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone);
	else if( check_S_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase--;	
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone);
	else if( check_G2_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone);
	else if(Tumour -> at (ith_clone) -> In_M_Phase)
	{
		Basic_Clonal_Expansion_V2(ith_clone, Dying_and_Newborn);
		Update_Mutant_Mutational_Effects_V2(ith_clone, hours, years, Dying_and_Newborn);
		Clonal_Extinction_From_Size(ith_clone);
		Reset_to_G1(ith_clone);
	}

}

void Clonal_Expansion::Dealyed_Mitosis_V2R(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years)
{
	std::tuple<unsigned long long int, 
			   unsigned long long int, 
			   unsigned long long int, bool, bool> Dying_and_Newborn (0, 0, 0, false, true);

	if( check_G1_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone);
	else if( check_S_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase--;	
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone);
	else if( check_G2_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone);
	else if(Tumour -> at (ith_clone) -> In_M_Phase)
	{
		Basic_Clonal_Expansion_V2R(ith_clone, Dying_and_Newborn, hours, years);
		Update_Mutant_Mutational_Effects_V2R(ith_clone, hours, years, Dying_and_Newborn);
		Clonal_Extinction_From_Size(ith_clone);
		Reset_to_G1(ith_clone);
	}

}


//select mitosis type
void Clonal_Expansion::Update_Population_V2(const unsigned int & hours, const unsigned int & years)
{
	unsigned int ith_clone = 0;
	unsigned int number_of_clones = Tumour -> size();

	for(ith_clone = 0; ith_clone < number_of_clones ; ith_clone++)
	{
		if(!Tumour -> at (ith_clone) -> clone_extinct)// if not extinct
		{
			if( Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
			{
				Non_Mutagenic_Mitosis_Standard_V1(ith_clone);
			}
			else
			{
				Dealyed_Mitosis_V2( ith_clone, hours, years );
			}
		}
	}


}

void Clonal_Expansion::Update_Population_V2R(const unsigned int & hours, const unsigned int & years)
{
		unsigned int ith_clone = 0;
	unsigned int number_of_clones = Tumour -> size();

	for(ith_clone = 0; ith_clone < number_of_clones ; ith_clone++)
	{
		if(!Tumour -> at (ith_clone) -> clone_extinct)// if not extinct
		{
			if( Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
			{
				Non_Mutagenic_Mitosis_Standard_V1(ith_clone);
			}
			else
			{
				//std::cout << " DRUG RESISTANCE KERNEL " << std::endl;
				Dealyed_Mitosis_V2R( ith_clone, hours, years );
			}
		}
	}
}


void Clonal_Expansion::Update_Population_V1(const unsigned int & hours, const unsigned int & years)
{
	unsigned int ith_clone = 0;
	unsigned int number_of_clones = Tumour -> size();

	for(ith_clone = 0; ith_clone < number_of_clones ; ith_clone++)
	{
		if(!Tumour -> at (ith_clone) -> clone_extinct)// if not extinct
		{
			if( Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
			{
				Non_Mutagenic_Mitosis_Standard_V1(ith_clone); // Change name
			}
			else
			{
				Select_Mitosis_Type( ith_clone, hours, years );
				//Dealyed_Mitosis_V1(ith_clone, hours, years);	// choose mitosis component 
			}
		}
		// else
		// {
		// 	std::cout << "Skiping clone " << ith_clone << std::endl;
		// }
	
	}//for
}

void Clonal_Expansion::print_Final_Status(const unsigned int & hours, const unsigned int & years, int & myID)
{
	usleep( (myID*1) + 10);
		std::cout << "ID: " << myID << " [ALIVE CELLS]: " <<  Population_Size  
			 << " [CLONES]: " << Tumour -> size() 
			 << "   H: " << hours  
			 << " Y: " << years 
			 << " FD: " << feedback  
			 << std::endl;
}

void Clonal_Expansion::print_Final_Status_DR(const unsigned int & hours, const unsigned int & years, int & myID)
{
	unsigned int ith_clone = 0;
	unsigned int number_of_clones = Tumour -> size();
	unsigned long long int Reistant_clones = 0;
	unsigned long long int Sensitive_Clones = 0; 
	for(ith_clone = 0; ith_clone < number_of_clones ; ith_clone++)
	{
		if(!Tumour -> at (ith_clone) -> clone_extinct)
		{
			if(Tumour -> at(ith_clone) -> Drug_Resistant)
			{
				Reistant_clones++;
			}
			else
			{
				Sensitive_Clones++;
			}
		}
	}


	usleep( (myID*1) + 10);
		std::cout << "ID: " << myID << " [ALIVE CELLS]: " <<  Population_Size  
			 << " [CLONES]: " << Tumour -> size() 
			 << "   H: " << hours  
			 << " Y: " << years 
			 << " FD: " << feedback 
			 << " DRC: " << Reistant_clones
			 << " SC: " << Sensitive_Clones
			 << std::endl;
}



void Clonal_Expansion::print_Status( const unsigned int & hours, const unsigned int & years,  const bool & each_100, int & myID)
{
	if(each_100 && (hours % 100 == 0))
	{	
		usleep( (myID*1) + 10);
		std::cout << "ID: " << myID << " [ALIVE CELLS]: " <<  Population_Size  
			 << " [CLONES]: " << Tumour -> size() 
			 << "   H: " << hours  
			 << " Y: " << years 
			 << " FD: " << feedback  
			 << std::endl;
	}
	//else
	//{
	//	usleep( (myID*1) + 10);
	//	std::cout<< "ID: " <<myID << " [ALIVE CELLS]: " <<  Population_Size  
	//		 << " [CLONES]: " << Tumour -> size() 
	//		 << "   H: " << hours  
	//		 << " Y: " << years 
	//		 << " FD: " << feedback  
	//		 << std::endl;
	//}

}

void Clonal_Expansion::Get_Competitive_Clonal_Size(unsigned int & Effective_Tumour_Size )
{
	if(Tumour -> size() > 1 )
	{
		for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++)
		{
			if(Tumour -> at (ith_clone) -> Clone_Size > Minimum_Clone_Size_Wrtite_TE )
				Effective_Tumour_Size++;
		}
	}
	else if(Tumour -> size() == 1)
	{
		Effective_Tumour_Size = 1;
	}
	else
	{
		Effective_Tumour_Size = 0;
	}
}

void Clonal_Expansion::Write_To_File( std::ofstream & te_file, unsigned long long int & elapsed_hours)
{

	unsigned int Effective_Tumour_Size = 0;
	Get_Competitive_Clonal_Size( Effective_Tumour_Size );
	te_file << Population_Size << "\t" << elapsed_hours << "\t"<< Effective_Tumour_Size   << "\n";
	elapsed_hours++;

}


void Clonal_Expansion::Update_Hours(unsigned int & hours, unsigned int & years)
{
	if(hours == 8764)//8764
	{
		hours = 0; years ++;
	}
}

void Clonal_Expansion::Write_Final_Population_Status_V1( const std::string & data_path)
{
	std::ofstream clonal_data;
	std::string fname =  data_path + "V1_Final_Population.txt";
	clonal_data.open(fname);
	clonal_data << "Generation_ID" << "\tClone_Size" << "\tMutation_Rate" << "\tProliferation_Rate" << std::endl;

	for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++)
	{
		clonal_data 
		<< Tumour -> at( ith_clone ) -> Generation_ID << "\t" 
		<< Tumour -> at (ith_clone) -> Clone_Size  << "\t" 
		<< Tumour -> at (ith_clone) -> Mutation_Rate << "\t" 
		<< Tumour -> at (ith_clone) -> P_Expansion[1] << std::endl;
		
	}
	clonal_data.close();

}

void Clonal_Expansion::Write_Final_Population_Status_V2( const std::string & data_path)
{
	std::ofstream clonal_data;
	std::string fname =  data_path + "V2_Final_Population.txt";
	clonal_data.open(fname);
	clonal_data << "Generation_ID" << "\tClone_Size" << "\tMutation_Rate" << "\tProliferation_Rate" << "\tMutations" << "\tBurden" << "\tMax_PR" <<std::endl;

	for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++)
	{
		clonal_data 
		<< Tumour -> at( ith_clone ) -> Generation_ID << "\t" 
		<< Tumour -> at (ith_clone) -> Clone_Size  << "\t" 
		<< Tumour -> at (ith_clone) -> Mutation_Rate << "\t" 
		<< Tumour -> at (ith_clone) -> P_Expansion[1] << "\t" 
		<< Tumour -> at (ith_clone) -> Clonal_Mutations << "\t"
		<< Tumour -> at (ith_clone) -> Clonal_Mutational_Burden << "\t"
		<< Tumour -> at (ith_clone) -> max_PR
		<< std::endl;
		
	}
	clonal_data.close();
}


void Clonal_Expansion::Valid_Configurtation(bool & valid)
{
	for(unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++)
	{
		if(Tumour -> at(ith_clone) -> Clone_Size > 100)
		{
			valid = true;
			break;
		}
	}
}

void Clonal_Expansion::Append_Tumour_Evolution_File(const unsigned long long int & elapsed_hours, const std::string & data_path)
{
	std::ofstream clonal_data;
	std::string fname =  data_path + "Tumour_Evolution.txt";
	clonal_data.open(fname, std::ofstream::out | std::ofstream::app);
	
	clonal_data << "T " <<  std::to_string( elapsed_hours ) 
				<< "\tN " <<  std::to_string( Newborn_Clones_at_t )
				<< "\tC " <<  std::to_string( Tumour -> size() )
				<< "\tS " <<  std::to_string( feedback)
				<< std::endl;
				Newborn_Clones_at_t = 0;


	for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++)
	{
		if( !Tumour -> at(ith_clone) -> clone_extinct && (Tumour -> at(ith_clone) -> Clone_Size > 100) )
		{
			clonal_data 
				<< Tumour -> at( ith_clone ) -> Generation_ID << "\t" 
				<< Tumour -> at (ith_clone) -> Clone_Size  << "\t" 
				<< Tumour -> at (ith_clone) -> Mutation_Rate << "\t" 
				<< Tumour -> at (ith_clone) -> P_Expansion[1] << std::endl;
		}
	}
	clonal_data.close();
}


void Clonal_Expansion::Population_Difference(bool & delta )
{
	double delta_pop = fabs(static_cast< double >( Population_Size ) - static_cast< double >( Population_at_prev_t ));
	
	if( delta_pop >= Delta_Population_Write_TEM )
		delta = true;

}



void Clonal_Expansion::Clonal_Difference(bool & delta_cl )
{
	double delta_clone = fabs(static_cast< double >( Clonality_at_prev_t ) - static_cast< double >( Tumour -> size() ));
	
	if( delta_clone >= Delta_Clone_Write_TEM )
		delta_cl = true;

}

void Clonal_Expansion::Write_Tumour_Evolution_Single(const unsigned long long int & elapsed_hours, const std::string & data_path)
{
	bool valid = false;
	bool delta = false;
	bool delta_cl = false;

	Valid_Configurtation(valid);
	Population_Difference( delta );
	Clonal_Difference( delta_cl );

	//if( elapsed_hours % 100 == 0 || delta )
	if( delta || delta_cl )
		if(valid)
			Append_Tumour_Evolution_File( elapsed_hours,  data_path );

}

void Clonal_Expansion::Write_Tumour_Evolution_Split(const unsigned long long int & elapsed_hours, const std::string & data_path)
{
	std::ofstream clonal_data;
	std::string fname =  data_path + "t_" + std::to_string( elapsed_hours ) + ".txt";
	clonal_data.open(fname);
	//clonal_data << "Generation_ID" << "\tClone_Size" << "\tMutation_Rate" << "\tProliferation_Rate" << std::endl;
	for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++)
	{
		if( !Tumour -> at(ith_clone) -> clone_extinct && (Tumour -> at(ith_clone) -> Clone_Size > 100) )
		{
			clonal_data 
			<< Tumour -> at( ith_clone ) -> Generation_ID << "\t" 
			<< Tumour -> at (ith_clone) -> Clone_Size  << "\t" 
			<< Tumour -> at (ith_clone) -> Mutation_Rate << "\t" 
			<< Tumour -> at (ith_clone) -> P_Expansion[1] << std::endl;
		}
	}
	clonal_data.close();
}

// Load all properties function
void Clonal_Expansion::Write_Tumour_Evolution(const unsigned long long int & elapsed_hours, const std::string & data_path, const bool single)
{
	if(single)
		Write_Tumour_Evolution_Single( elapsed_hours,  data_path);
	else
		Write_Tumour_Evolution_Split( elapsed_hours,  data_path);

}

void Clonal_Expansion::Stop_Gowth_Condition_Counter(unsigned int & stop_growth_counter, const unsigned int & _STOP_AFTER_DIAGNOSIS_COUNTER, const unsigned long long int & _DETECTABLE_POPULATION_SIZE )
{
	if(Population_Size > _DETECTABLE_POPULATION_SIZE)
	{
		stop_growth_counter++;

	}

}

void Clonal_Expansion::set_Penalty_Max_Pop_Size(const unsigned long long int _MAXIMUM_POPULATION_SIZE_GROWTH)
{
	Penalty_Max_Population_Size = (double) _MAXIMUM_POPULATION_SIZE_GROWTH;
}

void Clonal_Expansion::set_PR_DR_Difference(void)
{
	PR_DR_Difference =  Tumour -> at(0) -> P_Expansion[1] - Tumour -> at(0) -> P_Expansion[0];

}

void Clonal_Expansion::Create_Folder_TE(const std::string & data_path)
{
	//std::cout << "NEED TO CREATE FOLDER " << std::endl;
    mkdir(data_path.c_str(), ACCESSPERMS);
   // std::cout << "Folder Exists " << FolderExists(data_path) << std::endl;
}

void Clonal_Expansion::Check_Data_Path_Folder_Creation(const std::string & data_path)
{
	if( FolderExists(data_path) )
	{
    	//std::cout << "FOLDER " << data_path << " EXISTS AND VALID" << std::endl;
	}
    else
    {
    	Create_Folder_TE(data_path);
    }
    
}

void Clonal_Expansion::Create_TimeStamp_Path( std::string & data_path )
{
	std::time_t result = std::time(nullptr);
	std::string timestamp = std::asctime(std::localtime(&result));
	std::cout << timestamp << std::endl;
	data_path = data_path + "/" + ReplaceAll(timestamp, " ", "_") ;
	data_path = ReplaceAll(data_path, "\n", "");
	data_path = ReplaceAll(data_path, ":", "_");
	trim(data_path);
	data_path = data_path + "/";
   // std::cout <<" New path: " << data_path  << std::endl;

}

void Clonal_Expansion::setVersion_Folder( std::string & data_path, const std::map<std::string, std::string> &logic )
{

	std::cout << " Folder Data exists " << std::endl;
	data_path = data_path +  "/" + logic.at("Version");
	std::cout << "Ver path " << data_path << std::endl;
	if( FolderExists(data_path) )
	{
		Create_TimeStamp_Path( data_path );
    	Check_Data_Path_Folder_Creation( data_path );
    			
	}
	else
	{
		std::cout << "ERROR: " << data_path << " DOES NOT EXIST (Validate within the function)" << std::endl;
	}

}

void Clonal_Expansion::setTumour_Evolution_Folder(std::string & data_path, const std::map<std::string, std::string> &logic )
{
	char CWD[1024];
 	const std::string core_path = std::string(
               		            		 dirname(
                    	            	    dirname(
                        	            	    getcwd( CWD, sizeof(CWD) )
                            	            	    )
                                				)
                          					);

	//std::cout << "PATH: " << core_path << std::endl;
	std::vector<std::string> path_tokens =  split (core_path, '/');
	if( path_tokens.at(path_tokens.size()-1) =="core" )
	{
		//std::cout << "TRUE " << std::endl;
		data_path = core_path + "/Data";
		//std::cout << data_path << std::endl; 

		if(FolderExists( data_path ) )
			setVersion_Folder( data_path, logic );
		else
			std::cout << "ERROR: " << data_path << " INVALID PATH" << std::endl;

	}

}

//////////////////////// **** CONFIGURATION SECTION **** /////////////////////////

void Clonal_Expansion::Logic_File_Kernel_Configurations_V2(const std::map<std::string, std::string> &logic, bool & print_time_step, bool & write_tumour_evolution, bool & each_100)
{
	setSD_Penalty(logic.at("Penalty"));
	setMistosis_type(logic.at("Mitosis"));

	if(logic.at("Print_TS") == "true")
		print_time_step = true;

	if( logic.at("Print_TS.Each_100") == "true" )
		each_100 = true;

	if(logic.at("Write_Tumour_Evolution") == "true" )
		write_tumour_evolution = true;

	//std::cout << "Mutiple Mutations " << logic.at("Multiple_Mutations") << " " << Multiple_Mutations<< std::endl;
	if(logic.at("Multiple_Mutations") == "true")
		Multiple_Mutations = true;

	//std::cout << "Mutiple Mutations " << Multiple_Mutations << std::endl;

	//std::cout << "LOGIC " << logic.at("Write_TEM_Min_CS") << std::endl;
	Minimum_Clone_Size_Wrtite_TE = std::stoull(logic.at("Write_TEM_Min_CS"));
	//std::cout << "MEMEBER " << Minimum_Clone_Size_Wrtite_TE << std::endl;

	//std::cout << "LOGIC POP " << logic.at("Write_Delta_Population") << std::endl;
	//std::cout << "LOGIC CLO " << logic.at("Write_Delta_Clonality") << std::endl;

	Delta_Population_Write_TEM = std::stod(logic.at("Write_Delta_Population")) ;
	Delta_Clone_Write_TEM = std::stod(logic.at("Write_Delta_Clonality"));

	//std::cout << "LOGIC POP " << Delta_Population_Write_TEM << std::endl;
	//std::cout << "LOGIC CLO " << Delta_Clone_Write_TEM << std::endl;
	
}

void Clonal_Expansion::Logic_File_Kernel_Configurations_V1(const std::map<std::string, std::string> &logic, bool & print_time_step, bool & write_tumour_evolution, bool & each_100)
{
	setSD_Penalty(logic.at("Penalty"));
	setMistosis_type(logic.at("Mitosis"));

	if(logic.at("Print_TS") == "true")
		print_time_step = true;

	if( logic.at("Print_TS.Each_100") == "true" )
		each_100 = true;

	if(logic.at("Write_Tumour_Evolution") == "true" )
		write_tumour_evolution = true;	
}


void Clonal_Expansion::Setting_Population_Parameters(const std::map<std::string, std::string> &logic, 
													 unsigned long long int & _MAXIMUM_POPULATION_SIZE_GROWTH, 
													 unsigned long long int & _DETECTABLE_POPULATION_SIZE,
													 unsigned int & _STOP_AFTER_DIAGNOSIS_COUNTER )
{
	if ( logic.find("MAXIMUM_POPULATION_SIZE_GROWTH") == logic.end() ) 
	{ // not found
 		// std::cout << "NOT FOUND " << std::endl;
 		 _MAXIMUM_POPULATION_SIZE_GROWTH = MAXIMUM_POPULATION_SIZE;
	} 
	else 
	{ // found
  		//std::cout << "POP GROWTH " << logic.at("MAXIMUM_POPULATION_SIZE_GROWTH") << std::endl;
  		_MAXIMUM_POPULATION_SIZE_GROWTH  = std::stoull(logic.at("MAXIMUM_POPULATION_SIZE_GROWTH"));
  		//std::cout << "POP G UL: " << _MAXIMUM_POPULATION_SIZE_GROWTH << std::endl;
	}

	if( logic.find("MAXIMUM_POPULATION_SIZE_GROWTH") == logic.end() )
	{
		 //std::cout << "NOT FOUND " << std::endl;
		 _DETECTABLE_POPULATION_SIZE = DETECTABLE_POPULATION_SIZE;
	}
	else
	{
		_DETECTABLE_POPULATION_SIZE = std::stoull(logic.at("DETECTABLE_POPULATION_SIZE"));
		//std::cout << "POP G UL: " << _DETECTABLE_POPULATION_SIZE << std::endl;
	}

	if(logic.find("STOP_AFTER_DIAGNOSIS_COUNTER") == logic.end())
	{
		std::cout << "NOT FOUND " << std::endl;
		//_STOP_AFTER_DIAGNOSIS_COUNTER = STOP_AFTER_DIAGNOSIS_COUNTER;
	}
	else
	{
		_STOP_AFTER_DIAGNOSIS_COUNTER = std::stoull(logic.at("STOP_AFTER_DIAGNOSIS_COUNTER"));
		//std::cout << "POP G UL: " << _STOP_AFTER_DIAGNOSIS_COUNTER << std::endl;
	}
}

void Clonal_Expansion::Configure_Sampling_Distribution_Structures_V2(const std::map<std::string, std::string> &logic)
{
	
	std::vector<std::string> pr_vals = split(logic.at("Dist_Params"), ','); 
	if(pr_vals.size() > 0)
	{
		for(unsigned int i = 0; i < pr_vals.size(); i++)
		{
			PR_Sampling_Parameters.push_back(std::stod (pr_vals.at(i)));
			//std::cout << " val [" << i << "] " << trimmed(pr_vals.at(i)) << " vals " << PR_Sampling_Parameters.at(i) <<std::endl;
		}
	}// end of if

	setPR_Sampling_type(logic.at("PR_Dist"));
	setMR_Sampling_type(logic.at("MR_Dist"));

	std::vector<std::string> mr_vals = split(logic.at("MR_Dist_Params"), ',');
	if(mr_vals.size() > 0)
	{
		for(unsigned int i = 0; i < mr_vals.size(); i++)
		{	
			MR_Sampling_Parameters.push_back(std::stod (mr_vals.at(i)  ));
			//std::cout << " val [" << i << "] " << trimmed(mr_vals.at(i)) << " vals " << MR_Sampling_Parameters.at(i) <<std::endl;
		}
	}

	setMutational_Effect_Dist(logic.at("Mutational_Effect"));
	//std::cout << "Mut effect dist " << Mutational_Effects << std::endl;
	
	std::vector<std::string> meff_vals = split(logic.at("Mutational_Effect_Params"), ',');
	if(meff_vals.size() > 0)
	{
		for(unsigned int i = 0; i < meff_vals.size(); i++)
		{	
			MEFF_Sampling_Parameters.push_back(std::stod (meff_vals.at(i)  ));
			//std::cout << " val [" << i << "] " << trimmed(meff_vals.at(i)) << " vals " << MEFF_Sampling_Parameters.at(i) <<std::endl;
		}
	}

	std::vector<std::string> pass_vals = split(logic.at("Passenger_Distribution_Parameters"), ',');
	if(pass_vals.size() > 0)
	{
		for(unsigned int i = 0; i < pass_vals.size(); i++)
		{	
			Passenger_Distribution_Parameters.push_back(std::stod (pass_vals.at(i)  ));
			//std::cout << " val [" << i << "] " << trimmed(pass_vals.at(i)) << " vals " << Passenger_Distribution_Parameters.at(i) <<std::endl;
		}
	}

	/////// Deleterious 
	std::vector<std::string> dele_vals = split(logic.at("Deleterious_Distribution_Parameters"), ',');
	if(dele_vals.size() > 0)
	{
		for(unsigned int i = 0; i < dele_vals.size(); i++)
		{	
			Deleterious_Distribution_Parameters.push_back(std::stod (dele_vals.at(i)  ));
			//std::cout << " val [" << i << "] " << trimmed(dele_vals.at(i)) << " vals " << Deleterious_Distribution_Parameters.at(i) <<std::endl;
		}
	}

	/////// Beneficial 
	std::vector<std::string> bene_vals = split(logic.at("Beneficial_Distribution_Parameters"), ',');
	if(bene_vals.size() > 0)
	{
		for(unsigned int i = 0; i < bene_vals.size(); i++)
		{	
			Beneficial_Distribution_Parameters.push_back(std::stod (bene_vals.at(i)  ));
			//std::cout << " val [" << i << "] " << trimmed(bene_vals.at(i)) << " vals " << Beneficial_Distribution_Parameters.at(i) <<std::endl;
		}
	}

	setPassenger_Distribution( logic.at("Passenger_Distribution") );
	//std::cout << "Pass Dist " << Passenger_Distribution << std::endl;

	setDeleterious_Distribution( logic.at("Deleterious_Distribution")  );
	//std::cout << "Del Dist " << Deleterious_Distribution << std::endl;

	setBeneficial_Distribution( logic.at("Beneficial_Distribution") );
	//std::cout << "Ben Dist " << Beneficial_Distribution << std::endl;

}

void Clonal_Expansion::Configure_Sampling_Distribution_Structures_V1(const std::map<std::string, std::string> &logic)
{
	std::vector<std::string> pr_vals = split(logic.at("Dist_Params"), ','); 
	if(pr_vals.size() > 0)
	{
		for(unsigned int i = 0; i < pr_vals.size(); i++)
		{
			PR_Sampling_Parameters.push_back(std::stod (pr_vals.at(i)));
			std::cout << " val [" << i << "] " << trimmed(pr_vals.at(i)) << " vals " << PR_Sampling_Parameters.at(i) <<std::endl;
		}
	}// end of if

	setPR_Sampling_type(logic.at("PR_Dist"));
	setMR_Sampling_type(logic.at("MR_Dist"));

	std::vector<std::string> mr_vals = split(logic.at("MR_Dist_Params"), ',');
	if(mr_vals.size() > 0)
	{
		for(unsigned int i = 0; i < mr_vals.size(); i++)
		{	
			MR_Sampling_Parameters.push_back(std::stod (mr_vals.at(i)  ));
			std::cout << " val [" << i << "] " << trimmed(mr_vals.at(i)) << " vals " << MR_Sampling_Parameters.at(i) <<std::endl;
		}
	}
}


void Clonal_Expansion::Set_Data_Storage_Folders(const std::string & data_path, std::string & growth_path, std::string & evolution_path, int & myID, unsigned int & replicates)
{
	std::string iteration ="";
	if(myID == 0)
	{
		iteration = data_path + "IT_" + std::to_string(replicates) +"/";
		//std::cout << "Iteration path:" << iteration << std::endl;//<< "/" << "ID_" << std::to_string(myID)<< std::endl;
		mkdir(iteration.c_str(), ACCESSPERMS);
		Path_Bcast_From_Master( iteration );
	}
	else
	{
		Path_Bcast_From_Salves( iteration, myID);
	}
	

	std::string Core_ID_path = iteration + "ID_" + std::to_string(myID) +"/";
	//std::cout << "Core_ID_path " << Core_ID_path << std::endl;//<< "/" << "ID_" << std::to_string(myID)<< std::endl;
	mkdir(Core_ID_path.c_str(), ACCESSPERMS);
	growth_path = Core_ID_path + "Growth/";
	//std::cout << "Growth " << growth_path << std::endl;//<< "/" << "ID_" << std::to_string(myID)<< std::endl;
	mkdir(growth_path.c_str(), ACCESSPERMS);
	growth_path = growth_path + "Tumour_Growth.txt";
	 evolution_path = Core_ID_path + "Evolution/";
	//std::cout << "Evolution " << evolution_path << std::endl;
	mkdir(evolution_path.c_str(), ACCESSPERMS);
}


void Clonal_Expansion::Path_Bcast_From_Master(const std::string & data_path)
{
	char path[1024];
	strncpy(path, data_path.c_str(), sizeof(path));
	path[sizeof(path) - 1] = 0;

	int pathLength = sizeof(path);
	MPI_Bcast (&pathLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (path, pathLength, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void Clonal_Expansion::Path_Bcast_From_Salves(std::string & data_path, int & myID)
{
	unsigned pathLength;
	char * Path;
	MPI_Bcast (&pathLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
   	Path = (char *) malloc (pathLength);
   	MPI_Bcast (Path, pathLength, MPI_CHAR, 0, MPI_COMM_WORLD);
    printf ("Process %d: %s\n", myID, Path);
    data_path = std::string(Path);
    //std::cout << "Data recieved " << myID << " is " << data_path << std::endl;
}

//this will become the main
void Clonal_Expansion::Compute_Tumour_Growth_V1(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path)
{
	
	bool print_time_step = false;
	bool write_tumour_evolution = false;
	bool each_100 = true;
	bool single = false;

	std::cout << "LOGIC " << logic.at("TUMOUR_EVOLUTION_FILE") << std::endl;
	if(logic.at("TUMOUR_EVOLUTION_FILE") == "single")
	{
		std::cout << "SINGLE " << std::endl;
		single = true;
	}

	std::cout << "Mutation rate from file " << logic.at("MUTATION_RATE") << std::endl;

	Logic_File_Kernel_Configurations_V1(logic, print_time_step, write_tumour_evolution, each_100);

	std::cout << "LOGIC " << logic.at("Write_Tumour_Evolution_Hours") <<std::endl;
	

	unsigned long long int _MAXIMUM_POPULATION_SIZE_GROWTH = 0;
	unsigned long long int _DETECTABLE_POPULATION_SIZE = 0;
	unsigned int _STOP_AFTER_DIAGNOSIS_COUNTER = 0;

	Setting_Population_Parameters(logic,  _MAXIMUM_POPULATION_SIZE_GROWTH, _DETECTABLE_POPULATION_SIZE, _STOP_AFTER_DIAGNOSIS_COUNTER );

	Init_Random();
	Configure_Sampling_Distribution_Structures_V1( logic );

	//std::string data_path = "";
	if( data_path == "")
	{
		if(myID == 0)
		{
			setTumour_Evolution_Folder( data_path, logic );
			Path_Bcast_From_Master( data_path );
		}
		else
		{
			Path_Bcast_From_Salves( data_path, myID);
		}

		std::cout << "ID: " << myID << " path: " << data_path << std::endl;
		
	}
	

	unsigned int seconds = 0;
	unsigned int hours = 0;
	unsigned int years = 0;
	unsigned long long int elapsed_hours = 0;


	/* File stream */
	carcinogenesis_V1();
	carcinogenesis_V1_File_Input( std::stod( logic.at("MUTATION_RATE") ) );
	set_PR_DR_Difference();
	set_Penalty_Max_Pop_Size( _MAXIMUM_POPULATION_SIZE_GROWTH );
	if(myID==0) {printValues (Tumour -> at(0));}

	// uncomment form here

	 std::string growth_path = "";
	 std::string evolution_path = "";


	 Set_Data_Storage_Folders( data_path,  growth_path, evolution_path, myID, replicates);

	
	// getchar();

	std::ofstream te_file;
	te_file.open (growth_path);
	unsigned int stop_growth_counter = 0;

	// std::cout << "Max pop " << _MAXIMUM_POPULATION_SIZE_GROWTH << std::endl;
	// std::cout << "DETECT " << _DETECTABLE_POPULATION_SIZE << std::endl;
	// std::cout << "STOP aft diag " << _STOP_AFTER_DIAGNOSIS_COUNTER << std::endl;

	//getchar();

	while( ( stop_growth_counter < _STOP_AFTER_DIAGNOSIS_COUNTER) && !( Population_Size == 0) )
	{
		seconds += dt;
		if(seconds == 3600)
		{
			seconds = 0; hours ++;
			Update_Population_V1(hours, years);
			Select_Size_Dependant_Penalty();
			
			if( print_time_step )
				print_Status( hours, years, each_100, myID );
		}
		
		Update_Hours( hours, years );

		if(write_tumour_evolution )
			Write_Tumour_Evolution( elapsed_hours, evolution_path, single );

		Write_To_File( te_file, elapsed_hours );
		Stop_Gowth_Condition_Counter( stop_growth_counter, _STOP_AFTER_DIAGNOSIS_COUNTER, _DETECTABLE_POPULATION_SIZE );
	}

	te_file.close();
	//Print_Clones();
	print_Final_Status( hours, years, myID );

	//Print final sttus here

	Write_Final_Population_Status_V1( evolution_path );

	//print clone values
	// for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++) 			//(1) Is my clone not extinct?	tmr -> Tumour -> size()
	// {	
	// 	std::cout << "C_S[ " << ith_clone << " ] = " << Tumour -> at( ith_clone ) -> Clone_Size <<  "\t\t PR: " <<  Tumour -> at (ith_clone) -> P_Expansion[1] << " MR: " << Tumour -> at (ith_clone) -> Mutation_Rate  <<std::endl;
	// }

}



void Clonal_Expansion::Print_Clones(void)
{
	for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++) 			//(1) Is my clone not extinct?	tmr -> Tumour -> size()
	{	
		if(Tumour -> at(ith_clone) -> Clone_Size > 100 )
			std::cout << "C_S[ " << Tumour -> at( ith_clone ) -> Generation_ID << " ] = " << Tumour -> at( ith_clone ) -> Clone_Size <<  "\t\t PR: " <<  Tumour -> at (ith_clone) -> P_Expansion[1] << " MR: " << Tumour -> at (ith_clone) -> Mutation_Rate  <<std::endl;
	}
}



void Clonal_Expansion::Print_Clone_V2(void)
{
	for (unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++) 			//(1) Is my clone not extinct?	tmr -> Tumour -> size()
	{	
		if(Tumour -> at(ith_clone) -> Clone_Size > 100 )
			std::cout << "C_S[ " << Tumour -> at( ith_clone ) -> Generation_ID << " ] = " 
					  << Tumour -> at( ith_clone ) -> Clone_Size <<  "\tPR: " 
					  <<  Tumour -> at (ith_clone) -> P_Expansion[1] 
					  << " MR: " << Tumour -> at (ith_clone) -> Mutation_Rate  
					  << " \tMuts " << Tumour -> at (ith_clone) -> Clonal_Mutations
					  << " \t burden " <<   Tumour -> at (ith_clone) -> Clonal_Mutational_Burden
					  << " \t Max PR " << Tumour -> at (ith_clone) -> max_PR
					  << std::endl;
	}
}

///////////////////////////////////
/////////////////////////////////
// Quantile computation move later to Rsndom
///////
////////////

double Clonal_Expansion::Lerp(double v0, double v1, double t)
{
    return (1 - t)*v0 + t*v1;
}

std::vector<double> Clonal_Expansion::Quantile(const std::vector<double>& inData, const std::vector<double>& probs)
{
    if (inData.size() <= 2 || probs.empty())
    {
        throw std::runtime_error("Invalid input");
    }

    std::vector<double> data = inData;
    std::sort(data.begin(), data.end());
    std::vector<double> quantiles;

    for (size_t i = 0; i < probs.size(); ++i)
    {
        double center = Lerp(-0.5, data.size() - 0.5, probs[i]);

        size_t left = std::max(int64_t(std::floor(center)), int64_t(0));
        size_t right = std::min(int64_t(std::ceil(center)), int64_t(data.size() - 1));

        double datLeft = data.at(left);
        double datRight = data.at(right);

        double quantile = Lerp(datLeft, datRight, center - left);

        quantiles.push_back(quantile);
    }

    return quantiles;
}


// std::vector<std::vector<double>> Clonal_Expansion::Data( const std::vector<std::vector<double>> & Data)
// {

// }


/////////////////////////
//  Reshape the time series as a matrix with headers:
//  PR MR Mutations Clonal_Burden init_PR max_PR Clone_Size
/////////////////////////
void Clonal_Expansion::Quartile_Computation_Growth_Simulation(const std::string & data_path)
{
	//std::cout << data_path << std::endl;
	if( (Population_Size > 0) && (Tumour -> size() > 2 ) )
	{

		std::vector<std::vector<double>> Simulation_Data( 8, std::vector<double>( Tumour -> size() ) ); 
		std::vector<std::vector<double>> Quartile_Data( 5 , std::vector<double>( 8 ) );
	
		for(unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone ++)
		{
			Simulation_Data[0][ith_clone] = Tumour -> at(ith_clone) -> P_Expansion[1] ;
			Simulation_Data[1][ith_clone] = Tumour -> at(ith_clone) -> Mutation_Rate;
			Simulation_Data[2][ith_clone] = Tumour -> at(ith_clone) -> Clonal_Mutations ;
			Simulation_Data[3][ith_clone] = Tumour -> at(ith_clone) -> Clonal_Mutational_Burden ;
			Simulation_Data[4][ith_clone] = Tumour -> at(ith_clone) -> init_PR ;
			Simulation_Data[5][ith_clone] = Tumour -> at(ith_clone) -> max_PR ;
			Simulation_Data[6][ith_clone] = Tumour -> at(ith_clone) -> Clone_Size;

			//std::vector<std::string>  cid = split(split( (Tumour -> at(ith_clone) -> Generation_ID).c_str(),'-')[1], ':');
			//std::cout << cid[0] << " " <<cid[1] << std::endl;

			std::vector<std::string>  cid = split(split( (Tumour -> at(ith_clone) -> Generation_ID).c_str(),'-')[1], ':');
			
			double hrs = stod(cid[0]) * 8760.0 + stod(cid[1]);
			//std::cout << " hrs " << hrs << std::endl;

			Simulation_Data[7][ith_clone] = hrs;
		}
		
		for( unsigned int column = 0; column < 8; column++)
		{
			std::vector<double> quartiles = Quantile(Simulation_Data[column], { 0.0, 0.25, 0.5, 0.75, 1.0 });
			for (unsigned int i = 0; i < quartiles.size(); i++)
			{
				Quartile_Data[i][column] = quartiles[i];	
			}
		}

		// this si how to save into the file
		std::ofstream summary_data;
		std::string fname =  data_path + "Population_Summary.txt";
		summary_data.open(fname);
		summary_data << "PR" << "\tMR" << "\tMutations" << "\tBurden" << "\tiPR" << "\tmPR" << "\tCS" << "\thrs" << std::endl;

		//std::cout << "Quartile data DATA " << std::endl;
		for(unsigned int col = 0; col < 5; col ++)
		{
			for(unsigned int row = 0; row < 8; row++)
			{
				//std::cout << Quartile_Data[col][row] << "\t";
				if(row == 7)
				{
					summary_data << Quartile_Data[col][row] << std::endl;
				}
				else
				{
					summary_data << Quartile_Data[col][row] << "\t";
				}
			}
			//std::cout<<std::endl;
		}
		//std::cout << std::endl;
		summary_data.close();

	}
	else if( (Population_Size == 0) && ( Tumour -> size() > 2 ) )
	{
		std::vector<std::vector<double>> Simulation_Data( 7, std::vector<double>( Tumour -> size() ) ); 
		std::vector<std::vector<double>> Quartile_Data( 5 , std::vector<double>( 7 ) );

		for(unsigned int ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone ++)
		{
			Simulation_Data[0][ith_clone] = Tumour -> at(ith_clone) -> P_Expansion[1] ;
			Simulation_Data[1][ith_clone] = Tumour -> at(ith_clone) -> Mutation_Rate;
			Simulation_Data[2][ith_clone] = Tumour -> at(ith_clone) -> Clonal_Mutations ;
			Simulation_Data[3][ith_clone] = Tumour -> at(ith_clone) -> Clonal_Mutational_Burden ;
			Simulation_Data[4][ith_clone] = Tumour -> at(ith_clone) -> init_PR ;
			Simulation_Data[5][ith_clone] = Tumour -> at(ith_clone) -> max_PR ;

			std::vector<std::string>  cid = split(split( (Tumour -> at(ith_clone) -> Generation_ID).c_str(),'-')[1], ':');
			
			double hrs = stod(cid[0]) * 8760.0 + stod(cid[1]);
			//std::cout << " hrs " << hrs << std::endl;

			Simulation_Data[6][ith_clone] = hrs;
		}
		for( unsigned int column = 0; column < 7; column++)
		{
			std::vector<double> quartiles = Quantile(Simulation_Data[column], { 0.0, 0.25, 0.5, 0.75, 1.0 });
			for (unsigned int i = 0; i < quartiles.size(); i++)
			{
				Quartile_Data[i][column] = quartiles[i];	
			}
		}
		// this si how to save into the file
		//std::cout << "Quartile data DATA " << std::endl;
		std::ofstream summary_data;
		std::string fname =  data_path + "Population_Summary.txt";
		summary_data.open(fname);
		summary_data << "PR" << "\tMR" << "\tMutations" << "\tBurden" << "\tiPR" << "\tmPR" << "\thrs" << std::endl;
		for(unsigned int col = 0; col < 5; col ++)
		{
			for(unsigned int row = 0; row < 7; row++)
			{
				//std::cout << Quartile_Data[col][row] << "\t\t\t";
				if(row == 6)
				{
					summary_data << Quartile_Data[col][row] << std::endl;
				}
				else
				{
					summary_data << Quartile_Data[col][row] << "\t";
				}
			}
			//std::cout<<std::endl;
		}
		//std::cout << std::endl;
	}//else if
	else
	{
		//std::cout << " PARWNT CLONE DIED WITHOUT OFFSPRING " << std::endl;
	}


}


void Clonal_Expansion::Compute_Tumour_Growth_V2(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path)
{
	//std::cout << "In VERSION 2 " << std::endl;
	bool print_time_step = false;
	bool write_tumour_evolution = false;
	bool each_100 = false;
	bool single = false;

	//This should be in a separate function
	//std::cout << "LOGIC " << logic.at("TUMOUR_EVOLUTION_FILE") << std::endl;
	if(logic.at("TUMOUR_EVOLUTION_FILE") == "single")
	{
		// if(myID == 0)
		// {
		// 	std::cout << "SINGLE " << std::endl;	
		// }
		single = true;
	}

	// This too
	// std::cout << "mutation rate from file " << logic.at("MUTATION_RATE") << std::endl;
	 	
	Logic_File_Kernel_Configurations_V2(logic, print_time_step, write_tumour_evolution, each_100);

	unsigned long long int _MAXIMUM_POPULATION_SIZE_GROWTH = 0;
	unsigned long long int _DETECTABLE_POPULATION_SIZE = 0;
	unsigned int _STOP_AFTER_DIAGNOSIS_COUNTER = 0;

	Setting_Population_Parameters(logic, _MAXIMUM_POPULATION_SIZE_GROWTH, _DETECTABLE_POPULATION_SIZE, _STOP_AFTER_DIAGNOSIS_COUNTER );

	Init_Random();
	Configure_Sampling_Distribution_Structures_V2(logic);

	//std::string data_path = "";
	if(data_path == "")
	{
		if(myID == 0)
		{
			setTumour_Evolution_Folder( data_path, logic );
			Path_Bcast_From_Master( data_path );
		}
		else
		{
			Path_Bcast_From_Salves( data_path, myID);
		}
		std::cout << "ID: " << myID << " path: " << data_path << std::endl;
	}
	

	unsigned int seconds = 0;
	unsigned int hours = 0;
	unsigned int years = 0;
	unsigned long long int elapsed_hours = 0;


	/* File stream */
	//carcinogenesis_V2();// Change this to validate :: quick patch
	//std::cout << "MR: " << std::stod( logic.at("MUTATION_RATE") ) << std::endl;
	//std::cout << "PR: " << std::stod( logic.at("PROLIFERATION_RATE") ) << std::endl;
	

	carcinogenesis_V2_File_Input( std::stod( logic.at("MUTATION_RATE") ), std::stod( logic.at("PROLIFERATION_RATE") )  );
	//getchar();

	set_PR_DR_Difference();
	set_Penalty_Max_Pop_Size( _MAXIMUM_POPULATION_SIZE_GROWTH );
	//printValues (Tumour -> at(0));

	 std::string growth_path = "";
	 std::string evolution_path = "";

	 Set_Data_Storage_Folders( data_path,  growth_path, evolution_path, myID, replicates);
		
	//te_file.open ("./Mitosis_Model/te_file_33_V2.txt");
	std::ofstream te_file;
	te_file.open (growth_path);
	unsigned int stop_growth_counter = 0;

	
	//getchar();
	// We should add a years limit but those dynamics are stilk interesting
	while( ( stop_growth_counter < _STOP_AFTER_DIAGNOSIS_COUNTER) && !( Population_Size == 0) )
	{
		seconds += dt;
		
		if(seconds == 3600)
		{
			seconds = 0; hours ++;
			Population_at_prev_t = Population_Size;
			Clonality_at_prev_t = Tumour -> size();
			//std::cout << "P " << Population_Size << " \t p(-1) " << Population_at_prev_t << std::endl;
			Update_Population_V2( hours, years );
			Select_Size_Dependant_Penalty();
			
			 if(print_time_step)
			 	print_Status( hours, years, each_100, myID );
		}
		Update_Hours( hours, years );

		if(write_tumour_evolution)
			Write_Tumour_Evolution( elapsed_hours, evolution_path, single );


		Write_To_File( te_file, elapsed_hours );
		Stop_Gowth_Condition_Counter( stop_growth_counter, _STOP_AFTER_DIAGNOSIS_COUNTER, _DETECTABLE_POPULATION_SIZE );
		//Newborn_Clones_at_t = 0;

	}// while

	te_file.close();
	//Print_Clone_V2();
	
	print_Final_Status( hours, years, myID );

	//std::cout << "elapsed_hours " << elapsed_hours << std::endl;

	Write_Final_Population_Status_V2( evolution_path );

	
	//std::vector<double> PR_vec;
	Quartile_Computation_Growth_Simulation(evolution_path);


}


void Clonal_Expansion::Compute_Tumour_Growth_V2_DR(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path)
{
	std::cout << "In drug Resistance Verison " << std::endl;

	bool print_time_step = false;
	bool write_tumour_evolution = false;
	bool each_100 = false;
	bool single = false;

	std::cout << "LOGIC " << logic.at("TUMOUR_EVOLUTION_FILE") << std::endl;
	if(logic.at("TUMOUR_EVOLUTION_FILE") == "single")
	{
		std::cout << "SINGLE " << std::endl;
		single = true;
	}

	std::cout << "Mutation rate from file " << logic.at("MUTATION_RATE") << std::endl;
	getchar();

	//getchar();
	
	Logic_File_Kernel_Configurations_V2(logic, print_time_step, write_tumour_evolution, each_100);

	unsigned long long int _MAXIMUM_POPULATION_SIZE_GROWTH = 0;
	unsigned long long int _DETECTABLE_POPULATION_SIZE = 0;
	unsigned int _STOP_AFTER_DIAGNOSIS_COUNTER = 0;

	Setting_Population_Parameters(logic, _MAXIMUM_POPULATION_SIZE_GROWTH, _DETECTABLE_POPULATION_SIZE, _STOP_AFTER_DIAGNOSIS_COUNTER );

	Init_Random();
	Configure_Sampling_Distribution_Structures_V2(logic);

		if(data_path == "")
	{
		if(myID == 0)
		{
			setTumour_Evolution_Folder( data_path, logic );
			Path_Bcast_From_Master( data_path );
		}
		else
		{
			Path_Bcast_From_Salves( data_path, myID);
		}
		//std::cout << "ID: " << myID << " path: " << data_path << std::endl;
	}
	

	unsigned int seconds = 0;
	unsigned int hours = 0;
	unsigned int years = 0;
	unsigned long long int elapsed_hours = 0;


	/* File stream */
	//carcinogenesis_V2();// chaneg this
	std::cout << "MR: " << std::stod( logic.at("MUTATION_RATE") ) << std::endl;
	std::cout << "PR: " << std::stod( logic.at("PROLIFERATION_RATE") ) << std::endl;
	getchar();


	carcinogenesis_V2_File_Input( std::stod( logic.at("MUTATION_RATE") ), std::stod( logic.at("PROLIFERATION_RATE") )  );

	set_PR_DR_Difference();
	set_Penalty_Max_Pop_Size( _MAXIMUM_POPULATION_SIZE_GROWTH );
	//printValues (Tumour -> at(0));

	 std::string growth_path = "";
	 std::string evolution_path = "";

	 Set_Data_Storage_Folders( data_path,  growth_path, evolution_path, myID, replicates);
		
	//te_file.open ("./Mitosis_Model/te_file_33_V2.txt");
	std::ofstream te_file;
	te_file.open (growth_path);
	unsigned int stop_growth_counter = 0;

	//getchar();

	// TODO, TAKE THIS OUT !!
	while( ( stop_growth_counter < _STOP_AFTER_DIAGNOSIS_COUNTER) && !( Population_Size == 0) )
	{
		seconds += dt;
		
		if(seconds == 3600)
		{
			seconds = 0; hours ++;
			Population_at_prev_t = Population_Size;
			Clonality_at_prev_t = Tumour -> size();
			// DR KERNEL
			Update_Population_V2R( hours, years);
			Select_Size_Dependant_Penalty(); 

			 if(print_time_step)
			 	print_Status( hours, years, each_100, myID );

		}
		Update_Hours( hours, years );
		if(write_tumour_evolution)
			Write_Tumour_Evolution( elapsed_hours, evolution_path, single );


		Write_To_File( te_file, elapsed_hours );
		Stop_Gowth_Condition_Counter( stop_growth_counter, _STOP_AFTER_DIAGNOSIS_COUNTER, _DETECTABLE_POPULATION_SIZE );
		//Newborn_Clones_at_t = 0;

	}//while 


	te_file.close();
	//Print_Clone_V2();
	print_Final_Status_DR( hours, years, myID );
	std::cout << "elapsed_hours " << elapsed_hours << std::endl;

	Write_Final_Population_Status_V2( evolution_path );

	if(Population_Size == 0)
	{
		// Modify this to allow growth file and summary file
		std::string failed_to_grow = data_path + "/IT_"+ std::to_string(replicates) + "/ID_" + std::to_string(myID) + "/";
		std::string remove_folder = std::string("rm -rf ") + "\"" + failed_to_grow + "\"";
		std::cout << "Removing Folder " << failed_to_grow << std::endl;
		system( remove_folder.c_str() );
	}
	else
	{
		std::cout << "Simulation in process " << myID << " successfully grew" << std::endl;
	}



}

void Clonal_Expansion::Compute_Tumour_Growth(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path)
{
	switch(Version)
 	{
 		case 0:
            std::cout << "0 Version: " << Version << std::endl;
            break;
        case 1:
            std::cout << "1 Version: " << Version << std::endl;
            Compute_Tumour_Growth_V1(logic, replicates, myID,  data_path);
            break;
        case 2:
        	// if(myID == 0)
        	// {
        	// 	std::cout << "2 Version: " << Version << std::endl;	
        	// }
            Compute_Tumour_Growth_V2(logic, replicates, myID, data_path);
            MPI_Barrier(MPI_COMM_WORLD);
            break;
        case 3:
            std::cout << "3) DR Version: " << Version << std::endl;
            Compute_Tumour_Growth_V2_DR(logic, replicates, myID, data_path);
            MPI_Barrier(MPI_COMM_WORLD);
            break;
        case 4:
            std::cout << "4 Version: " << Version << std::endl;
            break;
        case 5:
            std::cout << "5 Version: " << Version << std::endl;
            break;
        default:
        	 std::cout << "Def Version: " << Version << std::endl;
        	 break;
 	}//switch	
}



