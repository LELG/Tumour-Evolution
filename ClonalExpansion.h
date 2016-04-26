/*
	Thisis a header file for the Clone Structure
*/

#ifndef CLONAL_EXPANSION_H
#define CLONAL_EXPANSION_H

#include "Clone.h"
#include <string>
#include <map>
#include <memory>           // std::unique_ptr
#include <vector>
#include "config.h"
class Clonal_Expansion
{

	/**
	    Variable Definitions
	**************************/
public:
	unsigned long long int Population_Size;

	double feedback; // This is the output of a proportional control system
	enum Mitosis 
	{
		Standard = 0, 
		Binomial, 
		Network
	} mitosis;
		
		/**
			This is the main Data Structure to store all clones, this will
			reflect the heterogeneity.

			NOTE: This is a map of of unorder maps, this will allow us to model
			micro enviroments. The layer can be computed in MPI
		*/
		//unordered_map<string, unique_ptr<Clone> > Tumour;
		//unordered_map<int, struct Clone*> Tumour;
		/* We can create as many constructures as we want */

		 std::vector<std::unique_ptr<Clone> > *Tumour = new std::vector<std::unique_ptr<Clone> >;

public:
	Clonal_Expansion();
	~Clonal_Expansion();
	unsigned int getSize();
	void printParameters();
	void printValues (std::unique_ptr<Clone> const & ith_clone);
	void add_Clone_to_Tumour(void);
	void carcinogenesis(void);
	void carcinogenesis(unsigned long long int &neoplastic_cells);
	std::string getMitosis_type(void);
	bool check_G1_phase(std::unique_ptr<Clone> const & ith_clone);
	bool exiting_G1_phase(std::unique_ptr<Clone> const & ith_clone);
	bool check_S_phase(std::unique_ptr<Clone> const & ith_clone);
	bool exiting_S_phase(std::unique_ptr<Clone> const & ith_clone);
	bool check_G2_phase(std::unique_ptr<Clone> const & ith_clone);
	bool exiting_G2_phase(std::unique_ptr<Clone> const & ith_clone);
	void Transition_From_G1_S(std::unique_ptr<Clone> const & ith_clone);
	void Non_Mutagenic_Mitosis_Standard(std::unique_ptr<Clone> const & ith_clone);
	void Transition_From_S_G2(std::unique_ptr<Clone> const & ith_clone);
	void Transition_From_G2_M(std::unique_ptr<Clone> const & ith_clone);
	void Update_Population_After_Division(std::unique_ptr<Clone> const & ith_clone);
	void Estimate_Mutational_Effects(unsigned int & Mutant_Cells, unsigned long long int & Clone_Size);
	unsigned int * Update_Clonal_Mutational_Burden(unsigned int & ith_clone, unsigned int & Mutant_cells, unsigned int & years, unsigned int & hours );
	void carcinogenesis_from_driver(unsigned int & ith_clone, unsigned int & years, unsigned int & hours);
	void Check_for_Clonal_Extintion_at_min_size(unsigned int & ith_clone);
	void Probabilities_of_Cell_Division (unsigned int & ith_clone, double & p_idle, double & p_go_to_G0 );
	void Check_Mitosis_Network_Status(unsigned int & ith_clone, unsigned int & years, unsigned int & hours);
	void Grow_A_Tumour();


	
};

#endif