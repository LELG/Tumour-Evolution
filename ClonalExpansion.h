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
#include <mpi.h>
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
		Delayed_V1 = 0, 
		Delayed_V2 = 1, 
		Network  = 2
	} Mitosis;

	enum Version
	{
		Tester = 0,
		V1     = 1,
		V2 	   = 2,
		V2R    = 3,
		V3 	   = 4,
		VM     = 5
	} Version;

	enum SD_Penalty
	{
		standard = 0,
		mean     = 1,
		quartile  = 2

	} SD_Penalty;

	enum PR_Sampling
	{
		Uniform     = 0,
		Uniform_V1  = 1,
		Beta        = 2,
		Pareto      = 3,
		Gamma       = 4,
		Laplace     = 5
	} PR_Sampling;

	enum MR_Sampling
	{
		uniform = 0,
		uniform_Gain =1,
		beta 	=2

	}MR_Sampling;

	enum Mutational_Effects
	{
		Normal_gradient = 0,
		Exponential     = 1,
		Double_Beta     = 2
	} Mutational_Effects;

	enum Passenger_Distribution
	{
		Passenger_Beta    = 0,
		Passenger_Pareto  = 1,
		Passenger_Laplace = 2,
		Passenger_Gamma   = 3
	} Passenger_Distribution;

	enum Deleterious_Distribution
	{
		Deleterious_Beta    = 0,
		Deleterious_Pareto  = 1,
		Deleterious_Laplace = 2,
		Deleterious_Gamma   = 3

	} Deleterious_Distribution;

	enum Beneficial_Distribution
	{
		Beneficial_Beta    = 0,
		Beneficial_Pareto  = 1,
		Beneficial_Laplace = 2,
		Beneficial_Gamma   = 3

	} Beneficial_Distribution;

	Random r;

	std::vector<double> PR_Sampling_Parameters;
	std::vector<double> MR_Sampling_Parameters;
	std::vector<double> MEFF_Sampling_Parameters;
	std::vector<double> Passenger_Distribution_Parameters;
	std::vector<double> Deleterious_Distribution_Parameters;
	std::vector<double> Beneficial_Distribution_Parameters;

	double Penalty_Max_Population_Size;
	double PR_DR_Difference;
	bool Multiple_Mutations;
	
	//std::vector<std::string> distribution_values;
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
	void Init_Random(void);
	void printValues (std::unique_ptr<Clone> const & ith_clone);
	void add_Clone_to_Tumour(void);
	void carcinogenesis(void);
	void carcinogenesis(unsigned long long int &neoplastic_cells);
	std::string getMitosis_type(void);
	std::string getVersion_type(void);
	std::string getSD_Penalty(void);
	std::string getMutational_Effect_Dist(void);
	//void getMitosis_type(void);

	void setMutational_Effect_Dist(const std::string & str_mutational_effect_dist);
	void setMR_Sampling_type(const std::string & str_mr_sampling);
	void setPR_Sampling_type(const std::string & str_sampling);
	void setMistosis_type(const std::string & str_mitosis);
 	void setVersion_type(const std::string & str_version);
 	void setSD_Penalty(const std::string & str_Penalty);
 	void setDeleterious_Distribution( const std::string & str_Deleterious_Distribution );
 	void setPassenger_Distribution( const std::string & str_Passenger_Distribution );
 	void setBeneficial_Distribution(const std::string & str_Beneficial_Distribution );
 	void Map_Feedback_Penalty(void);
 	void Select_MR_Sampling(const double & Parent_mu_rate, double & updated_MR);
 	void Select_Mutational_Effects_Distribution(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years, 
												 std::tuple<unsigned long long int, 
																 unsigned long long int, 
																 unsigned long long int, 
																 bool, bool> & Dying_and_Newborn
																 );

 	void Select_Mitosis_Type(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years);
 	void Select_Size_Dependant_Penalty(void);
 	void Select_Carcinogenesis(void);
 	void Select_PR_Samplig(const double & Parent_Proliferation_Rate, double & updated_PR);
 	void Select_Penalty_Type(void);
 	void Passenger_Beta_Sampling( double & Penalty );
 	void Select_Passenger_Distribution( double & Penalty );
 	void Deleterious_Beta_Sampling( double & Penalty );
 	void Select_Deleterious_distribution(double & Penalty);
 	void Beneficial_Beta_Sampling( double & Penalty );
 	void Select_Beneficial_distribution(double & Penalty);
 	void carcinogenesis_V1(void);
 	void carcinogenesis_V2(void);

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
	void Transition_From_G1_S(const unsigned int  & ith_clone);
	//void Non_Mutagenic_Mitosis_Standard(std::unique_ptr<Clone>  & ith_clone);
	void Non_Mutagenic_Mitosis_Standard_V1(const unsigned int & ith_clone);
	void Transition_From_S_G2(const unsigned int  & ith_clone);
	void Transition_From_G2_M(const unsigned int  & ith_clone);
	void Update_Population_After_Division_V1(const unsigned int  & ith_clone);
	void Clonal_Extinction_From_Size(const unsigned int & ith_clone);
	void Dealyed_Mitosis_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years);
	void Dealyed_Mitosis_V2(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years);
	void Estimate_Mutational_Effects(unsigned int & Mutant_Cells, unsigned long long int & Clone_Size);
	void Update_Clonal_Mutational_Burden(const unsigned int & ith_clone, const unsigned int & Mutant_cells, const unsigned int & years, const unsigned int & hours, std::vector<unsigned int> & Mutations );
	void carcinogenesis_from_driver(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours);
	void Check_for_Clonal_Extintion_at_min_size(const unsigned int & ith_clone);
	void Probabilities_of_Cell_Division (const unsigned int & ith_clone, double & p_idle, double & p_go_to_G0 );

	void Complementary_Probability( double & comp_probability_G0, 
									double & comp_probability_G1, 
									double & comp_probability_G2,
									double & comp_probability_S,
									double & comp_probability_M,
									const unsigned int & ith_clone );

	void Update_G0_Phase( const double & comp_probability_G0, const unsigned int & ith_clone);
	void Update_G1_Phase( const double & comp_probability_G1, const unsigned int & ith_clone);
	void Update_G2_Phase( const double & comp_probability_G2, const unsigned int & ith_clone);
	void Update_S_Phase( const double & comp_probability_S, const unsigned int  & ith_clone);
	void Update_M_Phase( const double & comp_probability_M, const unsigned int & ith_clone);


	void Check_Mitosis_Network_Status(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours);
	
	void Update_Tumour_Size();
	void map_Feedback(  );
	void Adjust_Feedback_PR(double & P_NB);
	void Check_Clonal_Extinction_BCE_V1(const unsigned int & ith_clone, std::tuple<unsigned long long int, unsigned long long int, unsigned long long int, bool, bool> & Dying_and_Newborn);
	void Print_Debug_Carcinogeneisis_From_Driver(const unsigned int & ith_clone);
	void Carcinogenesis_From_Driver_V2(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years);
	void Carcinogenesis_From_Driver_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years);
	void Checking_for_Mutants(std::tuple<unsigned long long int, 
													   unsigned long long int, 
													   unsigned long long int, 
													   bool, bool> & Dying_and_Newborn, 
													   const unsigned int & ith_clone );

	void Mutant_Effects_V1(const unsigned int & ith_clone, const unsigned long long int & Mutant_Cells,  const unsigned int & hours, const unsigned int & years);

	void Apply_Penalties_to_Mutants_V1(const unsigned int & ith_clone, const std::tuple<unsigned long long int, 
																						unsigned long long int, 
																						unsigned long long int, 
																						bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years);

	void Induce_Multiple_Mutations(unsigned int & number_of_mutations);
	
	void Normal_Gradient_Mutational_Effect(const unsigned int & ith_clone,  std::tuple<unsigned long long int, 
																							unsigned long long int, 
																							unsigned long long int, 
																							bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years);


	void Apply_Penalties_to_Mutants_V2(const unsigned int & ith_clone,  std::tuple<unsigned long long int, 
																									  unsigned long long int, 
																									  unsigned long long int, 
																									  bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years);

	void Update_Mutant_Mutational_Effects_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years, const std::tuple<unsigned long long int, 
																																					unsigned long long int, 
																																					unsigned long long int, 
																																					bool, bool> & Dying_and_Newborn);

	void Update_Mutant_Mutational_Effects_V2(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years,  std::tuple<unsigned long long int, 
																																								 unsigned long long int, 
																																								 unsigned long long int, 
																																								 bool, bool> & Dying_and_Newborn);

	void Update_Newborn_Parameters_V1(const std::vector<unsigned long long int> & NewBorn_Cells, std::tuple<unsigned long long int, 
																											unsigned long long int, 
																											unsigned long long int, 
																											bool, bool> & Dying_and_Newborn);

	void Basic_Clonal_Expansion_V1(const unsigned int & ith_clone, std::tuple<unsigned long long int, 
																			  unsigned long long int, 
																			  unsigned long long int, bool, bool> & Dying_and_Newborn);

	void Write_Tumour_Evolution(const unsigned long long int & elapsed_hours, const std::string & data_path);

	void Write_Final_Population_Sattus_V1( const std::string & data_path);
	void Write_Final_Population_Sattus_V2( const std::string & data_path);

	void Stop_Gowth_Condition_Counter(unsigned int & stop_growth_counter, const unsigned int & _STOP_AFTER_DIAGNOSIS_COUNTER, const unsigned long long int & _DETECTABLE_POPULATION_SIZE );

	void Reset_to_G1(const unsigned int & ith_clone);
	void Valid_Size_Heterogeneity_V1(const unsigned int & ith_clone);
	void Print_Clone_V2(void);
	void Print_Clones(void);
	void Update_Population_V2(const unsigned int & hours, const unsigned int & years);
	void Update_Population_V1(const unsigned int & hours, const unsigned int & years);
	void Compute_Tumour_Growth(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path);

	void Update_Hours(unsigned int & hours, unsigned int & years);
	void Get_Competitive_Clonal_Size(unsigned int & Effective_Tumour_Size );
	void Write_To_File(	 std::ofstream & te_file, unsigned long long int & elapsed_hours);
	void print_Status( const unsigned int & hours, const unsigned int & years,  const bool & each_100, int & myID);

	void set_PR_DR_Difference(void);
	void set_Penalty_Max_Pop_Size(const unsigned long long int _MAXIMUM_POPULATION_SIZE_GROWTH);

	void Setting_Population_Parameters(const std::map<std::string, std::string> &logic, 
										unsigned long long int & _MAXIMUM_POPULATION_SIZE_GROWTH, 
										unsigned long long int & _DETECTABLE_POPULATION_SIZE,
										unsigned int & _STOP_AFTER_DIAGNOSIS_COUNTER );

	void Logic_File_Kernel_Configurations_V1(const std::map<std::string, std::string> &logic, bool & print_time_step, bool & write_tumour_evolution, bool & each_100);
	void Logic_File_Kernel_Configurations_V2(const std::map<std::string, std::string> &logic, bool & print_time_step, bool & write_tumour_evolution, bool & each_100);
	void Configure_Sampling_Distribution_Structures_V2(const std::map<std::string, std::string> &logic);
	void Configure_Sampling_Distribution_Structures_V1(const std::map<std::string, std::string> &logic);
	void Compute_Tumour_Growth_V2(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path);

	void Create_Folder_TE(const std::string & data_path);
	void Check_Data_Path_Folder_Creation(const std::string & data_path);
	void Create_TimeStamp_Path( std::string & data_path);
	void setVersion_Folder( std::string & data_path, const std::map<std::string, std::string> &logic );
	void setTumour_Evolution_Folder(std::string & data_path, const std::map<std::string, std::string> &logic );
	void Set_Data_Storage_Folders(const std::string & data_path, std::string & growth_path, std::string & evolution_path, int & myID, unsigned int & replicates);
	void Path_Bcast_From_Master(const std::string & data_path);
	void Path_Bcast_From_Salves(std::string & data_path, int & myID);
	void Compute_Tumour_Growth_V1(const std::map<std::string, std::string> &logic, unsigned int & replicates, int & myID, std::string & data_path );

	
};

#endif