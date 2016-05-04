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
#include "Random.h"
#include <tuple>
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

	enum Version
	{
		Tester = 0,
		V1     = 1,
		V2 	   = 2,
		V2R    = 3,
		V3 	   = 4,
		VM     = 5
	} Version;
		
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
	std::string getVersion_type(void);
	//void getMitosis_type(void);

 	void setVersion_type(const std::string & str_version);
 	void Select_Carcinogenesis();
 	void carcinogenesis_V1(void);

	void Generate_Clone_Generation_ID( std::string & newGeneration_ID,
								  const int & Parent_Generation_ID_Counter, 
								  const std::string & Parent_ID, 
								  const unsigned int & years, 
								  const unsigned int & hours);


	bool check_G1_phase( const unsigned int & ith_clone);
	bool exiting_G1_phase( const unsigned int & ith_clone);
	bool check_S_phase(const unsigned int & ith_clone);
	bool exiting_S_phase(const unsigned int & ith_clone);
	bool check_G2_phase(const unsigned int & ith_clone);
	bool exiting_G2_phase( const unsigned int & ith_clone);
	void Transition_From_G1_S(const unsigned int  & ith_clone, Random & r);
	//void Non_Mutagenic_Mitosis_Standard(std::unique_ptr<Clone>  & ith_clone);
	void Non_Mutagenic_Mitosis_Standard_V1(const unsigned int & ith_clone, Random & r);
	void Transition_From_S_G2(const unsigned int  & ith_clone, Random & r);
	void Transition_From_G2_M(const unsigned int  & ith_clone, Random & r);
	void Update_Population_After_Division_V1(const unsigned int  & ith_clone);
	void Dealyed_Mitosis_V1(const unsigned int & ith_clone, Random & r, unsigned int & hours, unsigned int & years);
	void Estimate_Mutational_Effects(unsigned int & Mutant_Cells, unsigned long long int & Clone_Size);
	void Update_Clonal_Mutational_Burden(const unsigned int & ith_clone, const unsigned int & Mutant_cells, const unsigned int & years, const unsigned int & hours, std::vector<unsigned int> & Mutations, Random & r );
	void carcinogenesis_from_driver(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours, Random & r);
	void Check_for_Clonal_Extintion_at_min_size(const unsigned int & ith_clone);
	void Probabilities_of_Cell_Division (const unsigned int & ith_clone, double & p_idle, double & p_go_to_G0 );

	void Complementary_Probability( double & comp_probability_G0, 
									double & comp_probability_G1, 
									double & comp_probability_G2,
									double & comp_probability_S,
									double & comp_probability_M,
									const unsigned int & ith_clone );

	void Update_G0_Phase( const double & comp_probability_G0, const unsigned int & ith_clone, Random & r);
	void Update_G1_Phase( const double & comp_probability_G1, const unsigned int & ith_clone, Random & r);
	void Update_G2_Phase( const double & comp_probability_G2, const unsigned int & ith_clone, Random & r);
	void Update_S_Phase( const double & comp_probability_S, const unsigned int  & ith_clone, Random & r);
	void Update_M_Phase( const double & comp_probability_M, const unsigned int & ith_clone, Random & r);


	void Check_Mitosis_Network_Status(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours, Random & r);
	
	void Update_Tumour_Size();
	void map_Feedback(  );
	void Basic_Clonal_Expansion_V1(const unsigned int & ith_clone, std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn, Random & r);
	void Reset_to_G1(const unsigned int & ith_clone);
	void Valid_Size_Heterogeneity_V1(const unsigned int & ith_clone);
	void Update_Population_V1(unsigned int & hours, unsigned int & years, Random & r);
	void Compute_Tumour_Growth(const std::map<std::string, std::string> &logic);
	void Compute_Tumour_Growth_V1(const std::map<std::string, std::string> &logic);

	
};

#endif