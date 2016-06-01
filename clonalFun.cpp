#include "clonalFun.h"
#include "inout_funs.h"
#include "Random.h"
#include "config.h"
#include "ClonalExpansion.h"
#include "Clone.h"
#include <iostream>     // std::cout
#include <sstream>
#include <memory>           // std::unique_ptr
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <map>
#include <functional>

std::unique_ptr<Clone> get_Clone_DS()
{
    return std::unique_ptr<Clone>( new Clone( ) );
} // end function


inline bool FileExist (const std::string& name) 
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

std::unique_ptr<Clonal_Expansion> get_Clonal_Expasion_DS()
{
	return std::unique_ptr<Clonal_Expansion>( new Clonal_Expansion() );
}

unsigned long long int Division_Model(unsigned long long int Clone_Size) 
{
		//Clone_Size +=  Clone_Size * 2;
	return  (Clone_Size =  Clone_Size * 2);
}

// Call this function as Division_Model(&Clone_Size)
void Division_Model(unsigned long long int *Clone_Size)
{
	*Clone_Size = *Clone_Size * 2;
}


bool check_Clonal_Extinction(unsigned long long int clone_size, unsigned int death_cells, unsigned int newborne_cells)
{
	bool extinct = false;
	signed long long int _newborne = (signed long long int) newborne_cells;
	signed long long int _death = (signed long long int) death_cells;
	signed long long int _clone_size = (signed long long int) clone_size;
	signed long long int o = ( _clone_size + _newborne ) -  _death;
	if(o <= 0)
		extinct = true;
			

	return extinct;
}

//By ptr
bool check_Clonal_Extinction(unsigned long long int *clone_size, unsigned int *death_cells, unsigned int *newborne_cells)
{
	bool extinct = false;
	signed long long int _newborne = (signed long long int) *newborne_cells;
	signed long long int _death = (signed long long int) *death_cells;
	signed long long int _clone_size = (signed long long int) *clone_size;
	signed long long int o = ( _clone_size + _newborne ) -  _death;
	if(o <= 0)
		extinct = true;

	return extinct;
}

void Enviromental_Penalty(double *feedback, double *NewBorn_probability)
{
	if( *feedback >= *NewBorn_probability )
	{
		*NewBorn_probability = 0.0;
	}
	else
	{
		*NewBorn_probability -= *feedback;
	}
}


std::tuple<unsigned long long int, unsigned long long int, unsigned long long int, bool, bool> Basic_Clonal_Expansion(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int hours)
{

	std::tuple<unsigned long long int, // Dying Cells  1
					unsigned long long int, // Newborne cells 2
					unsigned long long int, // Mutant cells 3
									bool, // if there are mutant cells 4
									bool> Dying_and_Newborn (0, 0, 0, false, true);

	std::get<3>(Dying_and_Newborn) =  false;
	std::get<4>(Dying_and_Newborn) =  true;
	unsigned int buffer[4];
	//unsigned int mutant_cells = 0;
	/*
		1) Let's get basline parameters
	*/
	double P_DR =  tmr -> Tumour -> at(Generation_ID) -> P_Expansion[0];
	double P_NB =  tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];
	double P_MR =  tmr -> Tumour -> at(Generation_ID) -> Mutation_Rate;

	Enviromental_Penalty(&tmr -> feedback, &P_NB );

	double P_NT = 1.0 - (P_DR + P_NB + P_MR);
	double p_vector[] = {P_DR, P_NB, P_MR, P_NT};

	gsl_ran_multinomial ( r_global, K, tmr -> Tumour -> at (Generation_ID) -> Clone_Size, p_vector, buffer);

	std::get<0>(Dying_and_Newborn) = (unsigned long long int) buffer[0];  // For Dying Clones 		
	std::get<1>(Dying_and_Newborn) = (unsigned long long int) buffer[1]; // For Newborn Clones 
	std::get<2>(Dying_and_Newborn) = (unsigned long long int) buffer[2]; // For Mutants

	if(buffer[2] > 0)
	{
		std::get<3>(Dying_and_Newborn) =  true;
	}
	if( (std::get<2>(Dying_and_Newborn) == 0) && 
				( check_Clonal_Extinction(  &tmr -> Tumour -> at(Generation_ID) -> Clone_Size,   &buffer[0],  &buffer[1] ) ) 
			)
		{
			tmr -> Tumour -> at(Generation_ID) -> clone_extinct = true;
			tmr -> Tumour -> at(Generation_ID) -> final_PR = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];
			std::get<4>(Dying_and_Newborn) =  false;
			tmr -> Population_Size -= tmr -> Tumour  -> at(Generation_ID) -> Clone_Size ;
			tmr -> Tumour -> at(Generation_ID) -> Clone_Size = 0;
		}

	return Dying_and_Newborn;

}




// This is the pointer version
// Reassess this function
void Basic_Clonal_Expansion(std::unique_ptr<Clonal_Expansion>  & tmr, 
							int *clone_indx, 
							unsigned int *hours,
							unsigned long long int *dying_cells,
							unsigned long long int *newborn_cells,
							unsigned long long int *mutant_cells,
							bool *mutant_flag,
							bool *alive)// send data
{
	*mutant_flag = false;
	*alive = true;

	// Stat buffer
	unsigned int buffer[4];

	double prob_dying   = tmr -> Tumour -> at(*clone_indx) -> P_Expansion[0];
	double prob_newborn = tmr -> Tumour -> at(*clone_indx) -> P_Expansion[1];
	double prob_mutant  = tmr -> Tumour -> at(*clone_indx) -> Mutation_Rate;

	Enviromental_Penalty( &tmr ->  feedback, &prob_newborn  );

	double prob_idle = 1.0 - (prob_dying + prob_newborn + prob_mutant);
	double p_vec[] = {prob_dying, prob_newborn, prob_mutant, prob_idle};

	gsl_ran_multinomial ( r_global, K, tmr -> Tumour -> at (*clone_indx) -> Clone_Size, p_vec, buffer);

	*dying_cells    = (unsigned long long int) buffer[0];
	*newborn_cells = (unsigned long long int) buffer[1];
	*mutant_cells   = (unsigned long long int) buffer[2];

	std::cout << "p(nb) " << prob_newborn << std::endl;
	std::cout << "p(exp) " << tmr -> Tumour -> at(*clone_indx) -> P_Expansion[1] << std::endl;
	std::cout << "DY " << *dying_cells << std::endl;
	std::cout << "NB " << *newborn_cells << std::endl;
	std::cout << "IDLe " << buffer[3] << std::endl;
	std::cout << "MC " << *mutant_cells << std::endl;
	std::cout << "indx " << *clone_indx << std::endl;

	//Check for mutants
	if( *mutant_cells > 0)
	{
		*mutant_flag = true;
	}
	//Check for extintion if tehre are no mutants
	if( (*mutant_cells == 0) && check_Clonal_Extinction( &tmr -> Tumour -> at (*clone_indx) -> Clone_Size, &buffer[0], &buffer[1]) )
	{
		tmr -> Tumour -> at (*clone_indx) -> clone_extinct = true;
		tmr -> Tumour -> at (*clone_indx) -> final_PR = tmr -> Tumour -> at(*clone_indx) -> P_Expansion[1];
		*alive = false;
		tmr -> Population_Size -= tmr -> Tumour -> at(*clone_indx) -> Clone_Size;
		tmr -> Tumour -> at (*clone_indx) -> Clone_Size = 0;
	}
	
	
}

void Generate_Clone_Generation_ID(std::string & newGeneration_ID,
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

//Generate clone name Id 
//pointer version
// this is a tricky function
void Generate_Clone_Generation_ID(std::string *newGeneration_ID,
								  int *Parent_Generation_ID_Counter, 
								  const std::string &Parent_ID, 
								  unsigned int *years, 
								  unsigned int *hours)
{

	std::string heredity_pattern = Parent_ID.substr( 0, Parent_ID.find("-") );
	std::string Clone_ID ("");
		//cout << "I GET " << heredity_pattern << endl;

		if(heredity_pattern.compare("P") == 0)
		{
			//cout << "I GET 1 " << heredity_pattern << endl;
			std::string Clone_ID_tag ("PC");
			Clone_ID_tag.append(std::to_string(*Parent_Generation_ID_Counter));Clone_ID_tag.append("-");
			Clone_ID_tag.append(std::to_string(*years));Clone_ID_tag.append(":");Clone_ID_tag.append(std::to_string(*hours));
			Clone_ID = Clone_ID_tag;
		}
		else
		{
			
			std::string Clone_ID_tag ("");
			Clone_ID_tag.append(heredity_pattern); Clone_ID_tag.append(","); Clone_ID_tag.append(std::to_string(*Parent_Generation_ID_Counter));
			Clone_ID_tag.append("-");
			Clone_ID_tag.append(std::to_string(*years));Clone_ID_tag.append(":");Clone_ID_tag.append(std::to_string(*hours));
			Clone_ID = Clone_ID_tag;

		}

	//std::cout << "ID " << Parent_ID << std::endl;
//	std::cout << "HP " << heredity_pattern << std::endl;

	*newGeneration_ID = Clone_ID;

}


/**
	Adds a parent clone with size 1 in the tumour.
	Steps:
	1) Ask for a new Clone DS and append to vector Tumour
	2) From the recently added clone, initialise values

*******************************************/
void carcinogenesis(std::unique_ptr<Clonal_Expansion>  &tmr)
{
	Random r;

	tmr -> Tumour -> push_back( get_Clone_DS() );

	tmr -> Tumour -> back() -> Initiall_Expasion_Period = true;
	tmr -> Tumour -> back() -> Clone_Size = 1;

	tmr -> Tumour -> back() -> In_G0_Phase = false;
	tmr -> Tumour -> back() -> In_G1_Phase = true;
	tmr -> Tumour -> back() -> In_S_Phase = false;
	tmr -> Tumour -> back() -> In_G2_Phase = false;
	tmr -> Tumour -> back() -> In_M_Phase = false;

	tmr -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	tmr -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	tmr -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	tmr -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	tmr -> Tumour -> back() -> Generation_ID_Counter = 0; 
	tmr -> Tumour -> back() -> Generation_ID = "P-0:0";
	tmr -> Tumour -> back() -> max_PR = tmr -> Tumour -> back() -> P_Expansion[1];

	tmr -> Population_Size++;

	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
}


/**
	Adds a parent clone with size neoplastic_cells in the tumour.
	Steps:
	1) Ask for a new Clone DS and append to vector Tumour
	2) From the recently added clone, initialise values

*******************************************/
void carcinogenesis(std::unique_ptr<Clonal_Expansion>  &tmr, unsigned long long int &neoplastic_cells)
{
	Random r;

	tmr -> Tumour -> push_back( get_Clone_DS() );
	tmr -> Tumour -> back() -> Clone_Size = neoplastic_cells;

	tmr -> Tumour -> back() -> In_G0_Phase = false;
	tmr -> Tumour -> back() -> In_G1_Phase = true;
	tmr -> Tumour -> back() -> In_S_Phase = false;
	tmr -> Tumour -> back() -> In_G2_Phase = false;
	tmr -> Tumour -> back() -> In_M_Phase = false;

	tmr -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	tmr -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	tmr -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	tmr -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	tmr -> Tumour -> back() -> Generation_ID_Counter = 0; 
	tmr -> Tumour -> back() -> Generation_ID = "P-0:0";
	tmr -> Tumour -> back() -> max_PR = tmr -> Tumour -> back() -> P_Expansion[1];

	tmr -> Population_Size = tmr -> Tumour -> back() -> Clone_Size;
	
	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
}


void carcinogensis_from_Driver(std::unique_ptr<Clonal_Expansion>  & tmr,
								int *clone_indx,
								unsigned int *years,
								unsigned int *hours)
{

	Random r;
	 double mr = tmr -> Tumour -> at(*clone_indx) -> Mutation_Rate;
	 unsigned int AC = tmr -> Tumour -> at(*clone_indx) -> Driver_10_fold_accumulation;
	// unsigned int Accumulated_Drivers = tmr -> Tumour -> at(*clone_indx) -> Driver_10_fold_accumulation;
	// unsigned long long int NOM = tmr -> Tumour -> at (*clone_indx) -> Number_of_Mutations;
	 double pr = tmr -> Tumour -> at(*clone_indx) -> P_Expansion[1];

	
	std::string cloneName = "";	
	//std::cout << "cloneName Before: " << cloneName << std::endl;
	int Parent_Generation_ID_Counter = (int)tmr -> Tumour -> at(*clone_indx) -> Generation_ID_Counter;

	Generate_Clone_Generation_ID(&cloneName, 
								 &Parent_Generation_ID_Counter, 
								 tmr -> Tumour -> at(*clone_indx) -> Generation_ID,
								 years, 
								 hours);

	//std::cout << "cloneName After: " << cloneName << std::endl;

	tmr -> Tumour -> push_back( get_Clone_DS() ); // This will put a new clone in the vector list

	tmr -> Tumour -> back() -> Generation_ID = cloneName;
	tmr -> Tumour -> back() -> Initiall_Expasion_Period = false;
	tmr -> Tumour -> back() -> Clone_Size = 1;
	//tmr -> Tumour -> back() -> Number_of_Mutations = NOM + 1;

	tmr -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	tmr -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	tmr -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	tmr -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	tmr -> Tumour -> back() -> In_G0_Phase = false;
	tmr -> Tumour -> back() -> In_G1_Phase = true;
	tmr -> Tumour -> back() -> In_S_Phase = false;
	tmr -> Tumour -> back() -> In_G2_Phase = false;
	tmr -> Tumour -> back() -> In_M_Phase = false;

	tmr -> Tumour -> back() -> clone_extinct = false;
	tmr -> Tumour -> back() -> Generation_ID_Counter = 0;
	tmr -> Tumour -> back() -> Mutation_Rate = r.Uniform_Mutation_Rate_2(mr);
	double PR = r.Update_Proliferation_Rate(pr);
	tmr -> Tumour -> back() -> Driver_10_fold_accumulation = AC + 5;
	tmr -> Tumour -> back() -> P_Expansion[1] =  PR;
	tmr -> Tumour -> back() -> init_PR = PR;
	tmr -> Tumour -> back() -> max_PR = PR;
	tmr -> Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
	tmr -> Population_Size++;


}

double map_Model( double Current_Population_Size)
{
	return 0.0 + (1.0 - 0.0) * ((1.0 - 0.0) / (Current_Population_Size - 0.0));	
}

void map_Model( double *Current_Population_Size, double *Mapped)
{
	*Mapped = 0.0 + (1.0 - 0.0) * ((1.0 - 0.0) / (*Current_Population_Size - 0.0));	
}

void map_Feedback_Standard(double *Feedback, double *Population_Size )
{
	*Feedback = 0.0 + (DIFF - 0.0) * ((  *Population_Size - 0.0) / ((double) PS - 0.0));	
}

void map_Feedback( std::unique_ptr<Clonal_Expansion>  & tmr, double scale  )
{
	unsigned int ith_clone = 0;
	unsigned int size = tmr -> Tumour  -> size();
	double avg = 0.0;
	unsigned long long int k = 0;
	
	for (ith_clone = 0; ith_clone < size ; ith_clone++)
		if( !(tmr -> Tumour -> at (ith_clone) -> clone_extinct) )
		{
			avg += tmr -> Tumour -> at(ith_clone) -> P_Expansion[1] ;
			k++;
		}

	double diff = avg/((double) k * scale);

	tmr -> feedback = 0.0 + (diff - 0.0) * (( (double) tmr -> Population_Size - 0.0) / ((double) PS - 0.0));	
}//map_feedback

//This function is for deciding wich feedback function to use
void Size_Dependent_Penalty(std::unique_ptr<Clonal_Expansion>  & tmr)
{
		
	switch (PENALTY)
	{
		case 1:
		{
			map_Feedback(tmr, 1.0);
			break;
		}

		case 2:
		{
			map_Feedback(tmr, 2.0);
			break;
		}

		case 3:
		{
			map_Feedback_Standard( &tmr -> feedback, reinterpret_cast<double*>(tmr->Population_Size) ); /// LOOK THIS CAST
			break;
		}

	}	
			
}//funtion


//Check if can be optimised 

unsigned int add_multiple_mutations()
{
	Random r;
	unsigned int number_of_mutations = 1;
	if(MULTIPLE_MUTATIONS)
		number_of_mutations = r.Poisson();

	if(number_of_mutations == 0)
		number_of_mutations = 1;

	return number_of_mutations;
}

void mutations_in_cell_division(unsigned int *number_of_mutations)
{
	Random r;
	if(MULTIPLE_MUTATIONS)
		*number_of_mutations = r.Poisson();

	if(*number_of_mutations == 0)
		*number_of_mutations = 1;


}


unsigned long long int Mutational_Effects_Normal_Kernell(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int Number_of_Mutants, unsigned int years, unsigned int hours)
{
	unsigned long long int Mutants_in_clone = Number_of_Mutants;
	double p_vector[] = {KILLER_PROBABILITY, DRIVER_PROBABILITY, PASSENGER_PROBAILIBITY, DELETERIOUS_PROBABILITY, BENEFICIAL_PROBABILITY}; 

	unsigned int buffer[5];

	//unsigned int number_of_mutations = 1;

	/*
		(1) This is a mutant cell we need to know how
			many mutations picked at mitosis
 		(2) Determine the number of existing mustations
 	*/
	for(unsigned int j = 0; j < Number_of_Mutants; j++)
	{
		for(unsigned int i = 0; i < add_multiple_mutations() ; i++)
		{
			gsl_ran_multinomial ( r_global, 5, 1, p_vector, buffer);


			//Define mutational outcomes
			if(buffer[0] > 0)
			{
				Mutants_in_clone--;
			}
			if(buffer[1] > 0)
			{
				carcinogensis_from_Driver(tmr, &Generation_ID, &years, &hours);
				Mutants_in_clone--;
			}
			if(buffer[2] > 0)
			{
				//tmr -> Tumour-> at(Generation_ID) -> Number_of_Mutations++;
				double beta_penalty = gsl_ran_beta( r_global, 0.1, 200.0);
				double PR_t = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];
				if(beta_penalty >= PR_t)
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
				else
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= beta_penalty;
			}

			if(buffer[3] > 0)
			{
				//tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;
				double beta_penalty = gsl_ran_beta(r_global, 0.1, 25.0);
				double PR_t = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];
				if(beta_penalty >= PR_t)
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
				else
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= beta_penalty;
			}

			if(buffer[4]> 0)
			{
			//	tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;
				double beta_gain = gsl_ran_beta (r_global, 0.1, 40.0);

				tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] += beta_gain;
			}

		}//inner for
	}

	return Mutants_in_clone;


}


unsigned long long int Basic_Tester(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int Number_of_Mutants, unsigned int years, unsigned int hours)
{
	unsigned long long int Mutants_in_clone = Number_of_Mutants;
	double p_vector[] = {KILLER_PROBABILITY, DRIVER_PROBABILITY, PASSENGER_PROBAILIBITY, DELETERIOUS_PROBABILITY, BENEFICIAL_PROBABILITY};
	unsigned int buffer[5];

	for(unsigned int j = 0; j <Number_of_Mutants; j ++ )
	{
		for(unsigned int i = 0; i < add_multiple_mutations(); i++ )
		{
			gsl_ran_multinomial(r_global, 5, 1, p_vector, buffer);

			if(buffer[0] > 0 ) //Killer Mutation
			{
				Mutants_in_clone--;
			}
			if(buffer[1] > 0)
			{
				carcinogensis_from_Driver(tmr, &Generation_ID, &years, &hours);
				Mutants_in_clone--;
			}
			if(buffer[2] > 0 )
			{
				//tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;
				//if(tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations % 700 == 0)
				//	tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= 0.002;
			}
			if(buffer[3]>0)
			{
				if( (tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] - 0.002) <= 0.0)
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
				else
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= 0.002;
			}
			if(buffer[4] >0)
			{
				tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] += 0.001;
			}


		}//inner for
	}

	return Mutants_in_clone;


}

unsigned long long int Size_Mapped(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned long long int Number_of_Mutants, unsigned int years, unsigned int hours)
{
	unsigned long long int Mutants_in_clone = Number_of_Mutants;
	double p_vector[] = {KILLER_PROBABILITY, DRIVER_PROBABILITY, PASSENGER_PROBAILIBITY, DELETERIOUS_PROBABILITY, BENEFICIAL_PROBABILITY};
	unsigned int buffer[5];

	for(unsigned int j = 0; j <Number_of_Mutants; j ++ )
	{
		for(unsigned int i = 0; i < add_multiple_mutations(); i++)
		{
			gsl_ran_multinomial (r_global, 5, 1, p_vector, buffer);

			if(buffer[0] > 0)
			{
				Mutants_in_clone--;
			}
			if(buffer[1] > 0)
			{
				carcinogensis_from_Driver(tmr, &Generation_ID, &years, &hours);
				Mutants_in_clone--;
			}
			if(buffer[2] > 0 )
			{
				//tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;
				double beta_penalty = gsl_ran_beta(r_global, 0.1, 500.0);
				double PR_t = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];
				if(beta_penalty >= PR_t)
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
				else
					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= beta_penalty;
			}

			if(buffer[3] > 0)
			{
				double N = ( (double) tmr -> Tumour -> at(Generation_ID) -> Clone_Size );
				double Inv_N = 1.0/((double) tmr -> Tumour -> at(Generation_ID) -> Clone_Size ); 
				double Stochastic_Penalty = gsl_ran_beta(r_global, 1.0, 30.0 + N); //20
				double Total_Penalty = Inv_N + Stochastic_Penalty;
				double PR_t = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];

				if( Total_Penalty > PR_t )
 					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
 				else
 					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= Total_Penalty;

 				//tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;

			}

			if(buffer[4] > 0)
			{
				double N = ( (double) tmr -> Tumour -> at(Generation_ID) -> Clone_Size );
				double Inv_N = 1.0/((double) tmr -> Tumour -> at(Generation_ID) -> Clone_Size ); 
				double Stochastic_Gain = gsl_ran_beta(r_global, 1, 40 + N); //40
				double Total_Gain = Inv_N + Stochastic_Gain;
 				tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] += Total_Gain;
 				//tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;
			}
		}//inner for
	}
	return Mutants_in_clone;
}

unsigned long long int Laplace_Effect(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned long long int Number_of_Mutants, unsigned int years, unsigned int hours)
{
	unsigned long long int Mutants_in_clone = Number_of_Mutants;
	double p_vector[] = {KILLER_PROBABILITY, DRIVER_PROBABILITY, 1.0 - (KILLER_PROBABILITY + DRIVER_PROBABILITY) }; 
	unsigned int buffer[3];

	for(unsigned int j = 0; j <Number_of_Mutants; j ++ )
	{
		for(unsigned int i = 0; i < add_multiple_mutations()  ; i++) /* ADD MULTIPLE MUTATIONS */
 		{
 			gsl_ran_multinomial ( r_global, 3, 1, p_vector, buffer);

 			if(buffer[0] > 0)
 			{
 				Mutants_in_clone--;
 			}
 			if(buffer[1] >0)
 			{
 				carcinogensis_from_Driver(tmr, &Generation_ID, &years, &hours);
 				Mutants_in_clone--;
 			}
 			if(buffer[2]>0)
 			{
 				//tmr -> Tumour -> at(Generation_ID) -> Number_of_Mutations++;
 				double penalty = (gsl_ran_beta(r_global, 2.0, 9.0)/1000.0)-0.0004 ; 
 				double PR_t = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];

 				if(penalty <= 0.0)
 				{
 					if( (PR_t - penalty) <= 0.0 )
 					{
 						tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
 					}
 					else
 					{
 						tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] += penalty;
 					}
 				}
 				else
 				{
 					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] += penalty;
 				}

 				if(tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] < 0.0)
 					tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1] = 0.0;
 			}


 		}//innner for
	}

	return Mutants_in_clone;

}

unsigned long long int Mutational_Effect_From_Mutants(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned long long int Number_of_Mutants, unsigned int years, unsigned int hours  )
{
	unsigned long long int Mutants_in_clone = 0;
		switch (MUTATION_MODEL)
		{
			case 1:
			{
				Mutants_in_clone = Mutational_Effects_Normal_Kernell(tmr, Generation_ID, Number_of_Mutants, years,  hours);
				break;
			}

			case 2:
			{
				Mutants_in_clone = Size_Mapped(tmr, Generation_ID, Number_of_Mutants, years,  hours);
				break;
			}

			case 3:
			{
				Mutants_in_clone = Laplace_Effect(tmr, Generation_ID, Number_of_Mutants, years,  hours);
				break;
			}
			case 4:
			{
				Mutants_in_clone = 0;
				break;
			}
			case 5:
			{
				Mutants_in_clone = Basic_Tester( tmr, Generation_ID, Number_of_Mutants, years,  hours );
				break;
			}

			return Mutants_in_clone;

		}//SWITCH	
		return Mutants_in_clone;
}


void Check_for_size(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, signed long long int clone_size, signed long long int death_cells, signed long long int newborn_cells, signed long long int clonal_mutants, signed long long int total_mutants)
{
		
	signed long long int New_Size = (clone_size - death_cells) + (newborn_cells + (total_mutants - clonal_mutants) );
	unsigned long long int  _new_size = 0;
	if(New_Size <= 0)// clone is dead
	{
		tmr -> Tumour -> at(Generation_ID) -> clone_extinct = true;
		tmr -> Tumour -> at(Generation_ID) -> final_PR = tmr -> Tumour -> at(Generation_ID) -> P_Expansion[1];
	}
	else
	{
		_new_size = (unsigned long long int) New_Size;
	}

	tmr -> Population_Size += _new_size - tmr -> Tumour  -> at(Generation_ID) -> Clone_Size;
	tmr -> Tumour -> at(Generation_ID) -> Clone_Size = _new_size;

} 

void compute_Mutations(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID ,std::tuple<unsigned long long  int, unsigned long long  int, unsigned long long int, bool, bool> Dying_and_Newborn, unsigned int years, unsigned int hours )
{
		/*
			Version 1 => requires that all new mutants are considered new clones
		*/
			unsigned long long int non_driver_mutants = 0;
		if( std::get<3>(Dying_and_Newborn) ) // Check the mutants if they are new drivers
		{
			non_driver_mutants = Mutational_Effect_From_Mutants(tmr, Generation_ID, std::get<2>(Dying_and_Newborn), years, hours );

		}
		 Check_for_size( tmr, Generation_ID, (signed long long int) tmr -> Tumour -> at(Generation_ID) -> Clone_Size,
		 																		(signed long long int) std::get<0>(Dying_and_Newborn), 
		 																		(signed long long int) std::get<1>(Dying_and_Newborn), 
		 																		(signed long long int) non_driver_mutants, 
		 																		(signed long long int) std::get<2>(Dying_and_Newborn) );
				
}


void Initial_Expasion_Mitosis(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int hours, unsigned int years)
{

	Random r;


	if(
		tmr -> Tumour -> at (Generation_ID) -> In_G1_Phase &&
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase != 0
		)
	{
		tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_G1_Phase--;
	}
	else if (
				tmr -> Tumour -> at (Generation_ID) -> In_G1_Phase &&
				tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase == 0
		     )
	{
		tmr -> Tumour -> at (Generation_ID) -> In_G1_Phase = false;
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase = r.G1();
		tmr -> Tumour -> at (Generation_ID) -> In_S_Phase = true;
	}
	else if (
				tmr -> Tumour -> at (Generation_ID) -> In_S_Phase &&
				tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase != 0
		     )
	{
		tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase--;
	}
	else if (
				tmr -> Tumour -> at(Generation_ID) -> In_S_Phase &&
				tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase ==0
			)
	{
		tmr -> Tumour -> at(Generation_ID) -> In_S_Phase = false;
		tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase = r.S();
		tmr -> Tumour -> at(Generation_ID) -> In_G2_Phase = true;
	}
	else if (
				tmr -> Tumour -> at(Generation_ID) -> In_G2_Phase &&
				tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_G2_Phase != 0
			)
	{
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase--;
	}
	else if (
				tmr -> Tumour -> at(Generation_ID) -> In_G2_Phase &&
				tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase == 0
			)
	{
		tmr -> Tumour -> at (Generation_ID) -> In_G2_Phase = false;
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase = r.G2();
		tmr -> Tumour -> at (Generation_ID) -> In_M_Phase = true;
	}
	else if(
				tmr -> Tumour -> at (Generation_ID) -> In_M_Phase
			)
	{
		unsigned long long int cells_after_division = Division_Model(tmr -> Tumour -> at (Generation_ID) ->  Clone_Size);
		tmr -> Population_Size = (tmr -> Population_Size - tmr -> Tumour -> at(Generation_ID) -> Clone_Size ) + cells_after_division;
		tmr -> Tumour -> at (Generation_ID) -> Clone_Size = cells_after_division;

		if(tmr -> Tumour -> at (Generation_ID) -> Clone_Size >= tmr -> Tumour -> at (Generation_ID) -> Number_of_Memebers_to_Start_Heterogeneity )
		{
			tmr -> Tumour -> at(Generation_ID) -> Initiall_Expasion_Period = false;
		}

		tmr -> Tumour -> at(Generation_ID) -> In_M_Phase = false;
		tmr -> Tumour -> at(Generation_ID) -> In_G1_Phase = true;

	}


}//

void Mitosis(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int hours, unsigned int years)
{
	Random r;
	if(
		tmr -> Tumour -> at(Generation_ID) -> In_G1_Phase &&
		tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_G1_Phase != 0
	  )
	{
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase--;
	}
	else if (
				tmr -> Tumour -> at (Generation_ID) -> In_G1_Phase &&
				tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_G1_Phase == 0
			)
	{
		tmr -> Tumour -> at (Generation_ID) -> In_G1_Phase = false;
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase = r.G1();
		tmr -> Tumour -> at (Generation_ID) -> In_S_Phase = true;
	}
	else if (
				tmr -> Tumour -> at (Generation_ID) -> In_S_Phase &&
				tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase != 0
			)
	{
		tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase--;
	}
	else if ( 
				tmr -> Tumour -> at(Generation_ID) -> In_S_Phase &&
				tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase == 0
			)
	{
		tmr -> Tumour -> at (Generation_ID) -> In_S_Phase = false;
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase = r.S();
		tmr -> Tumour -> at (Generation_ID) -> In_G2_Phase = true;
	}
	else if (
				tmr -> Tumour -> at (Generation_ID) -> In_G2_Phase &&
				tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase != 0
			)
	{
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase--;
	}
	else if (
				tmr -> Tumour -> at(Generation_ID) -> In_G2_Phase &&
				tmr -> Tumour -> at(Generation_ID) -> Remaining_Time_in_G2_Phase == 0
			)
	{
		tmr -> Tumour -> at (Generation_ID) -> In_G2_Phase = false;
		tmr -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase = r.G2();
		tmr -> Tumour -> at (Generation_ID) -> In_M_Phase = true;
	}
	else if (
				tmr -> Tumour -> at (Generation_ID) -> In_M_Phase
			)
	{
		std::tuple<unsigned long long int, unsigned long long int, unsigned long long int, bool, bool> Dying_and_Newborn = Basic_Clonal_Expansion( tmr, Generation_ID, hours ); 
		if(std::get<4>(Dying_and_Newborn))
		{
			compute_Mutations(tmr, Generation_ID, Dying_and_Newborn, years, hours );
		}
		tmr -> Tumour -> at (Generation_ID) -> In_M_Phase = false;
		tmr -> Tumour -> at (Generation_ID) -> In_G1_Phase = true;

	}



}//Mitosis 


void update_population(std::unique_ptr<Clonal_Expansion>  & tmr, unsigned int hours, unsigned int years)
{
	unsigned int size = tmr -> Tumour -> size ();
	unsigned int ith_clone = 0;

	for (ith_clone = 0; ith_clone < size ; ith_clone++) 			//(1) Is my clone not extinct?	
	{
		if( !(tmr -> Tumour -> at (ith_clone) -> clone_extinct) )	// (1.a) if not, is my clone in free division?
		{
			if(tmr -> Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
			{
				Initial_Expasion_Mitosis(tmr, ith_clone, hours, years);
			}
			else
			{
				Mitosis( tmr, ith_clone, hours, years );
			}
		}// if
	}// if

}

void print_Status(std::unique_ptr<Clonal_Expansion>  & tmr, unsigned int hours, unsigned int years, int myID,  bool each_100 /*=false*/ )
{
		
	if(each_100 && (hours % 100 == 0))
	{	
		usleep( (myID*1) + 10);
		std::cout << " \n\n ACTIVE CELLS " <<  tmr -> Population_Size  
			 << " CLONES " << tmr -> Tumour -> size() 
				<< "   H: " << hours  
				<< " Y: " << years 
			 	<< " FD: " << tmr -> feedback 
			 	<< " ID: " << myID  
			 	<< std::endl;
		 }
		 else if(each_100 == false)
			std::cout << " \n\n ACTIVE CELLS " <<  tmr -> Population_Size  
				 << " CLONES " << tmr -> Tumour -> size() 
				 << "   H: " << hours  
				 << " Y: " << years 
			 	<< " FD: " << tmr -> feedback
			 	<< " ID: " << myID  
			 	<< std::endl;
		
			
}// print command

unsigned int Abort_Condition(unsigned long long int Population_Size, unsigned int times_to_wait)
{
	if(Population_Size >  DETECTABLE_POPULATION_SIZE )
	{
		times_to_wait++;
		std::cout <<"times to wait " <<times_to_wait<<std::endl;
		
	}

	return times_to_wait;
}

// Modify this model to support a comptrenhensive framework of paralell and estimation
void Compute_Tumour_Growth(std::unique_ptr<Clonal_Expansion>  & tmr, std::string BasePath, int myID)
{
	unsigned int seconds = 0;
	unsigned int hours = 0;
	unsigned int years = 0;
	//unsigned int elapsed_hours = 0;
	unsigned int times_to_wait = 0;
	//bool print = false;

	carcinogenesis(tmr);
	tmr -> printParameters();

	while( ((times_to_wait < STOP_GROWTH_AFTER_DIAGNOSIS) && (years < 60) ) && (tmr -> Population_Size > 0) )
	{
		seconds += dt;
		if(seconds == 3600)
		{
			seconds = 0; hours ++;
			update_population(tmr, hours, years);
			//std::cout << "H" ;
			Size_Dependent_Penalty( tmr );
			print_Status(tmr, hours, years, myID, true );
			times_to_wait = Abort_Condition(tmr -> Population_Size, times_to_wait);
		}
		if(hours == 8764)//8764
		{
			hours = 0; years ++;
			//print = true;
		}

	}//while
	
}

void trim_path_to_vector(std::vector<std::string> & tokens, const std::string & s_path)
{
	tokens = split(s_path, '/');
}

void build_core_path(std::vector<std::string> & tokens, std::string & path_to_core)
{
	if(tokens.size() != 0)
	{
		for (auto i = tokens.begin(); i != tokens.end(); ++i)
		{
			if(*i == "core")
			{
				path_to_core += *i + "/";
				break;
			}
			else
			{
				path_to_core += *i + "/";
			}
		}
	}

}

bool check_valid_folders(const std::string & path_to_core, const std::string & path_to_settings)
{
	return FolderExists(path_to_core) & FolderExists(path_to_settings) ;
}


void build_settings_path(const std::string & path_to_core, std::string & path_to_settings)
{
	path_to_settings = path_to_core + "settings/";
	if( check_valid_folders( path_to_core, path_to_settings) )
	{
		path_to_settings +=  "logic_flow.lgic";
		if(!FileExists (path_to_settings))
		{
			std::cout << "ERROR: FILE NOT FOUND \n" ;
		}
	}
	else
	{
		std::cout << "ERROR: INVALID PATH \n" ;
	}
}


void Load_Logic_File(std::map<std::string, std::string> & logic )
{
	//using namespace std;
	char buffer[1000];
	char *path = getcwd(buffer, sizeof(buffer));
	std::vector<std::string> tokens; 
	std::string s_path;

	std::string path_to_core = "";
	std::string path_to_logic_file = "";

	if (path)
	{
    	s_path = path;
    	std::cout << "DIR: " << s_path << std::endl; 
    	trim_path_to_vector(tokens, s_path);
    	build_core_path(tokens, path_to_core);
    	build_settings_path(path_to_core,  path_to_logic_file);
	
	}
	std::cout << path_to_logic_file << std::endl;
	logic = ReadAndParse( path_to_logic_file );
	
}


void print_logic_file(const std::map<std::string, std::string> & logic)
{
	std::cout << "///// RUNNING PARAMETERS ////////" << std::endl;
	for(auto it = logic.cbegin(); it != logic.cend(); ++it)
	{
    	std::cout << it->first << " = " << it->second << "\n";
	}
}


// 	////////////////////////////////////////////////////////////////
// 	std::cout << "DATA TEST VALIDATION " << std::endl;
// 	std::cout << "Gen ID: " << tmr -> Tumour -> back() -> Generation_ID << std::endl;
// 	std::cout << "Initiall_Expasion_Period " << tmr -> Tumour -> back() -> Initiall_Expasion_Period << std::endl;
// 	std::cout << "Clone_Size " << tmr -> Tumour -> back() -> Clone_Size << std::endl;
// 	std::cout << "Number_of_Mutations " << tmr -> Tumour -> back() -> Number_of_Mutations << std::endl;
// 	std::cout << "Remaining_Time_in_G1_Phase "<< tmr -> Tumour -> back() -> Remaining_Time_in_G1_Phase << std::endl;
// 	std::cout << "Remaining_Time_in_S_Phase " << tmr -> Tumour -> back() -> Remaining_Time_in_S_Phase << std::endl;
// 	std::cout << "Remaining_Time_in_G2_Phase " << tmr -> Tumour -> back() -> Remaining_Time_in_G2_Phase << std::endl;
// 	std::cout << "Remaining_Time_in_M_Phase " << tmr -> Tumour -> back() -> Remaining_Time_in_M_Phase << std::endl;
// 	std::cout << "In_G0_Phase " << tmr -> Tumour -> back() -> In_G0_Phase << std::endl;
// 	std::cout << "In_G1_Phase " << tmr -> Tumour -> back() -> In_G1_Phase << std::endl;
// 	std::cout << "In_G2_Phase " << tmr -> Tumour -> back() -> In_G2_Phase << std::endl;
// 	std::cout << "In_M_Phase " << tmr -> Tumour -> back() -> In_M_Phase << std::endl;
// 	std::cout << "clone_extinct " << tmr -> Tumour -> back() -> clone_extinct << std::endl;
// 	std::cout << "Generation_ID_Counter " << tmr -> Tumour -> back() -> Generation_ID_Counter << std::endl;
// 	std::cout << "Mutation_Rate " << tmr -> Tumour -> back() -> Mutation_Rate << std::endl;
// 	std::cout << "Driver_10_fold_accumulation " << tmr -> Tumour -> back() -> Driver_10_fold_accumulation << std::endl;
// 	std::cout << "P_Expansion[1]  " << tmr -> Tumour -> back() -> P_Expansion[1]  << std::endl;
// 	std::cout << "init_PR " << tmr -> Tumour -> back() -> init_PR << std::endl;
// 	std::cout << "max_PR " << tmr -> Tumour -> back() -> max_PR << std::endl;
// 	std::cout << "Number_of_Memebers_to_Start_Heterogeneity " << tmr -> Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity << std::endl;
// 	std::cout << "Population " << tmr -> Population_Size << std::endl;


// //std::cout << "Clone name in new clone struct " << tmr -> Tumour -> back() -> Generation_ID << std::endl;


// 	////////////////////////////////////////////////////////////////
// 	std::cout << "DATA TEST VALIDATION AFTER" << std::endl;
// 	std::cout << "Gen ID: " << tmr -> Tumour -> back() -> Generation_ID << std::endl;
// 	std::cout << "Initiall_Expasion_Period " << tmr -> Tumour -> back() -> Initiall_Expasion_Period << std::endl;
// 	std::cout << "Clone_Size " << tmr -> Tumour -> back() -> Clone_Size << std::endl;
// 	std::cout << "Number_of_Mutations " << tmr -> Tumour -> back() -> Number_of_Mutations << std::endl;
// 	std::cout << "Remaining_Time_in_G1_Phase "<< tmr -> Tumour -> back() -> Remaining_Time_in_G1_Phase << std::endl;
// 	std::cout << "Remaining_Time_in_S_Phase " << tmr -> Tumour -> back() -> Remaining_Time_in_S_Phase << std::endl;
// 	std::cout << "Remaining_Time_in_G2_Phase " << tmr -> Tumour -> back() -> Remaining_Time_in_G2_Phase << std::endl;
// 	std::cout << "Remaining_Time_in_M_Phase " << tmr -> Tumour -> back() -> Remaining_Time_in_M_Phase << std::endl;
// 	std::cout << "In_G0_Phase " << tmr -> Tumour -> back() -> In_G0_Phase << std::endl;
// 	std::cout << "In_G1_Phase " << tmr -> Tumour -> back() -> In_G1_Phase << std::endl;
// 	std::cout << "In_G2_Phase " << tmr -> Tumour -> back() -> In_G2_Phase << std::endl;
// 	std::cout << "In_M_Phase " << tmr -> Tumour -> back() -> In_M_Phase << std::endl;
// 	std::cout << "clone_extinct " << tmr -> Tumour -> back() -> clone_extinct << std::endl;
// 	std::cout << "Generation_ID_Counter " << tmr -> Tumour -> back() -> Generation_ID_Counter << std::endl;
// 	std::cout << "Mutation_Rate " << tmr -> Tumour -> back() -> Mutation_Rate << std::endl;
// 	std::cout << "Driver_10_fold_accumulation " << tmr -> Tumour -> back() -> Driver_10_fold_accumulation << std::endl;
// 	std::cout << "P_Expansion[1]  " << tmr -> Tumour -> back() -> P_Expansion[1]  << std::endl;
// 	std::cout << "init_PR " << tmr -> Tumour -> back() -> init_PR << std::endl;
// 	std::cout << "max_PR " << tmr -> Tumour -> back() -> max_PR << std::endl;
// 	std::cout << "Number_of_Memebers_to_Start_Heterogeneity " << tmr -> Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity << std::endl;
// 	std::cout << "Population " << tmr -> Population_Size << std::endl;




















	

