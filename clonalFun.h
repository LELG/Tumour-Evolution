#ifndef GLOBAL_CLONAL_FUN
#define GLOBAL_CLONAL_FUN

#include "Clone.h"
#include "ClonalExpansion.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <memory>           // std::unique_ptr
#include <string>
#include <sstream>          // std::istringstream
#include <tuple>
#include <algorithm>
#include <iterator>
#include <map>
#include <functional>
/**	FUNCTION get_Clone_DS()
		
	This function returns a 
	unique pointer of type Clone. Remmeber
	that unique pointers are memory efficient
	and guaratee a unique scope of the DS.

********************************************/

std::unique_ptr<Clone> get_Clone_DS();

std::unique_ptr<Clonal_Expansion> get_Clonal_Expasion_DS();

unsigned long long int Division_Model(unsigned long long int Clone_Size);
void Division_Model(unsigned long long int *Clone_Size);
bool check_Clonal_Extinction(unsigned long long int clone_size, unsigned int death_cells, unsigned int newborne_cells);
bool check_Clonal_Extinction(unsigned long long int *clone_size, unsigned int *death_cells, unsigned int *newborne_cells);
void Enviromental_Penalty(double *feedback, double *NewBorn_probability);
void Basic_Clonal_Expansion(std::unique_ptr<Clonal_Expansion>  & tmr, 
							int *clone_indx, 
							unsigned int *hours,
							unsigned long long int *dying_cells,
							unsigned long long int *newborn_cells,
							unsigned long long int *mutant_cells,
							bool *mutant_flag,
							bool *alive); // send data
std::tuple<unsigned long long int, unsigned long long int, unsigned long long int, bool, bool> Basic_Clonal_Expansion(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int hours);

void Generate_Clone_Generation_ID(std::string & newGeneration_ID,
								  int & Parent_Generation_ID_Counter, 
								  const std::string & Parent_ID, 
								  unsigned int & years, 
								  unsigned int & hours);

void Generate_Clone_Generation_ID( std::string *newGeneration_ID,
									int *Parent_Generation_ID_Counter, 
									const std::string &Parent_ID, 
									unsigned int *years, 
									unsigned int *hours);

void carcinogenesis(std::unique_ptr<Clonal_Expansion>  &tmr);
void carcinogenesis(std::unique_ptr<Clonal_Expansion>  &tmr, unsigned long long int &neoplastic_cells);

void carcinogensis_from_Driver(std::unique_ptr<Clonal_Expansion>  & tmr,
								int *clone_indx,
								unsigned int *years,
								unsigned int *hours);

double map_Model( double Current_Population_Size);
void map_Model( double *Current_Population_Size, double *Mapped);
void map_Feedback_Standard(double *Feedback, double *Population_Size );
void map_Feedback( std::unique_ptr<Clonal_Expansion>  & tmr, double scale  );
void Size_Dependent_Penalty(std::unique_ptr<Clonal_Expansion>  & tmr);
unsigned int add_multiple_mutations();
void mutations_in_cell_division(unsigned int *number_of_mutations);
unsigned long long int Mutational_Effects_Normal_Kernell(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int Number_of_Mutants, unsigned int years, unsigned int hours);
unsigned long long int Basic_Tester(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int Number_of_Mutants, unsigned int years, unsigned int hours);
unsigned long long int Size_Mapped(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned long long int Number_of_Mutants, unsigned int years, unsigned int hours);
unsigned long long int Laplace_Effect(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned long long int Number_of_Mutants, unsigned int years, unsigned int hours);
unsigned long long int Mutational_Effect_From_Mutants(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned long long int Number_of_Mutants, unsigned int years, unsigned int hours  );
void Check_for_size(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, signed long long int clone_size, signed long long int death_cells, signed long long int newborn_cells, signed long long int clonal_mutants, signed long long int total_mutants);
void compute_Mutations(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID ,std::tuple<unsigned long long  int, unsigned long long  int, unsigned long long int, bool, bool> Dying_and_Newborn, unsigned int years, unsigned int hours );
void Initial_Expasion_Mitosis(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int hours, unsigned int years);
void Mitosis(std::unique_ptr<Clonal_Expansion>  & tmr, int Generation_ID, unsigned int hours, unsigned int years);
void update_population(std::unique_ptr<Clonal_Expansion>  & tmr, unsigned int hours, unsigned int years);
void print_Status(std::unique_ptr<Clonal_Expansion>  & tmr, unsigned int hours, unsigned int years, int myID,  bool each_100 = false );
unsigned int Abort_Condition(unsigned long long int Population_Size, unsigned int times_to_wait);
void Compute_Tumour_Growth(std::unique_ptr<Clonal_Expansion>  & tmr, std::string BasePath, int myID);
void trim_path_to_vector(std::vector<std::string> & tokens, const std::string & s_path);
void build_core_path(std::vector<std::string> & tokens, std::string & path_to_core);
bool check_valid_folders(const std::string & path_to_core, const std::string & path_to_settings);
void build_settings_path(const std::string & path_to_core, std::string & path_to_settings);
void Load_Logic_File(std::map<std::string, std::string> & logic );
void print_logic_file(const std::map<std::string, std::string> & logic);

#endif