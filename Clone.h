/*
	Thisis a header file for the Clone Structure
*/

#ifndef CLONE_H
#define CLONE_H

#include <string>
#include <map>
#include <memory>           // std::unique_ptr
#include <vector>
#include <tuple>        // std::tuple, std::make_tuple, std::get
class Clone 
{

	/**
	    Variable Definitions
	**************************/

	/**
        NOTE: Defining as a const makes it read only
        consider this to model only non - Stem Cells.
    *************************************************/
public:
    unsigned int  Number_of_Memebers_to_Start_Heterogeneity;
	unsigned int Generation_ID_Counter;
	unsigned int Driver_10_fold_accumulation;
	bool clone_extinct;
	double  Mutation_Rate;

	/**
	    The Following variables can be modified at runtime.	
    ********************************************************/
	unsigned long long int  Number_of_Mutations;
	unsigned long long int Clone_Size; 
	bool Initiall_Expasion_Period;

	//Mitosis values
	bool In_G0_Phase;
	bool In_G1_Phase;
	bool In_S_Phase;
	bool In_G2_Phase;
	bool In_M_Phase;

	double P_Expansion[4]; 

	unsigned int Remaining_Time_in_G1_Phase;
	unsigned int Remaining_Time_in_S_Phase;
	unsigned int Remaining_Time_in_G2_Phase;
	unsigned int Remaining_Time_in_M_Phase;

	double init_PR;
	double max_PR;
	double final_PR;

	std::string Generation_ID;

	
	std::tuple<unsigned int, // G0{S} 
			   unsigned int, // G0{Dy}
			   unsigned int, // G0{->G1}
			   double, 		// px {Stay}
			   double> 		// py {Die}
			   			G0_status;
								
	std::tuple<unsigned int, // G1 {S}
			   unsigned int, // G1 {Dy}
			   unsigned int, // G1 {->G2}
			   double,		 // px
			   double>		 // py	 	
			   			G1_status;
	
	std::tuple<unsigned int, // G2{S} 
			   unsigned int, // G2{Dy}
			   unsigned int, // G2{->S}
			   double,		 // px
			   double> 		 // py
			   			G2_status;
	
	std::tuple<unsigned int, // S{S}
			   unsigned int, // S{Dy}
			   unsigned int, // S{->M} 
			   double,		 // px 
			   double> 		 // py 
			   			S_status;
	
	std::tuple<unsigned int, // M{S}
			   unsigned int, // M{Dy}
			   unsigned int, // M{Div}
			   double,       // px
			   double> 	     // py
			   			M_status;

	unsigned int G0_cells;
	unsigned int G1_cells;
	unsigned int G2_cells;
	unsigned int  S_cells;
	unsigned int  M_cells;

	unsigned int Available_cells;

public:
	Clone();
	~Clone();
	void setNumber_of_Memebers_to_Start_Heterogeneity(unsigned int _Number_of_Memebers_to_Start_Heterogeneity = 0);
	void printTest(void);
	void LoadElementsFromMap(std::map<std::string, std::string> *rawData);
	void FileReader(std::string path);
	void printValues (void);


};

#endif