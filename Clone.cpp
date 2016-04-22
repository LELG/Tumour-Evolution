#include "Clone.h"

#include "inout_funs.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <memory>           // std::unique_ptr
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>
#include <tuple>        // std::tuple, std::make_tuple, std::get

/***********************  CLONE DEFINITIONS ********************************/ 
/*
	Clone::Clone is the default constructor
	Clone::~Clone is the deconstructor

*/
/****************************************************************************/

// Deafult Parameters
Clone::Clone()
: 
	Number_of_Memebers_to_Start_Heterogeneity(500), 
	Generation_ID_Counter(0),
	Driver_10_fold_accumulation(10),
	clone_extinct(0),
	Mutation_Rate(MUT_RATE),
	Number_of_Mutations(0),
	Clone_Size(0),
	Initiall_Expasion_Period(0),
	In_G0_Phase(0),
	In_G1_Phase(0),
	In_S_Phase(0),
	In_G2_Phase(0),
	In_M_Phase(0),
	P_Expansion{DEATH_RATE, PROLIFERATION_RATE, MUT_RATE, 0.01}, //P_Expansion{DEATH_RATE, PROLIFERATION_RATE, MUT_RATE, 1.0-(DEATH_RATE + PROLIFERATION_RATE + MUT_RATE)},
	Remaining_Time_in_G1_Phase(0),
	Remaining_Time_in_S_Phase(0),
	Remaining_Time_in_G2_Phase(0),
	Remaining_Time_in_M_Phase(0),
	init_PR(0.0),
	max_PR(0.0),
	final_PR(0.0),
	Generation_ID(""),
	G0_status {0, 0, 0, 0.0, 0.0},
	G1_status {0, 0, 0, 0, 0.0, 0.0},			
	G2_status {0, 0, 0, 0.0, 0.0},
	S_status  {0, 0, 0, 0.0, 0.0},
	M_status  {0, 0, 0, 0.0, 0.0},
	G0_cells (0),
	G1_cells (0),
	G2_cells (0),
	S_cells  (0),
	M_cells  (0),
	Available_cells(0)
 { }

 Clone::~Clone()
 { }


//No need for this, not private
void Clone::setNumber_of_Memebers_to_Start_Heterogeneity(unsigned int _Number_of_Memebers_to_Start_Heterogeneity)
 {
 	Number_of_Memebers_to_Start_Heterogeneity = _Number_of_Memebers_to_Start_Heterogeneity;
 }

 void Clone::printTest(void)
 {
 	std::cout << " This should be zero " << Number_of_Memebers_to_Start_Heterogeneity << std::endl;
 }



void  Clone::LoadElementsFromMap(std::map<std::string, std::string> *rawData)
{

	double _Mutation_Rate = 0.0;
	double _Proliferation_Rate = 0.0;
	double _Death_Rate = 0.0;
	bool _Mutation_Rate_flag = false;
	bool _Proliferation_Rate_flag = false;
	bool _Death_Rate_flag = false;

	if(rawData -> count("Mutation_Rate") )
	{
	 	_Mutation_Rate = std::stod(rawData -> at ("Mutation_Rate"));
	 	Mutation_Rate = _Mutation_Rate;
	 	_Mutation_Rate_flag = true;
	 }

	 if(rawData -> count("Proliferation_Rate"))
	 {
	 	_Proliferation_Rate = std::stod(rawData -> at ("Proliferation_Rate"));
	 	_Proliferation_Rate_flag = true;
	 }

	 if(rawData -> count("Death_Rate"))
	 {
	 	 _Death_Rate = std::stod(rawData -> at ("Death_Rate"));
	 	_Death_Rate_flag = true;
	 }

	 if( rawData -> count("Number_of_Memebers_to_Start_Heterogeneity") )
	 	std::stringstream(rawData -> at("Number_of_Memebers_to_Start_Heterogeneity")) >> Number_of_Memebers_to_Start_Heterogeneity ;

	 if( rawData -> count("Generation_ID_Counter") )
	 	std::stringstream(rawData -> at("Generation_ID_Counter")) >> Generation_ID_Counter ;

	 if( rawData -> count("Driver_10_fold_accumulation") )
	 	std::stringstream(rawData -> at("Driver_10_fold_accumulation")) >> Driver_10_fold_accumulation ;

	 if( rawData -> count("clone_extinct") )
	 {
	 	if(rawData -> at("clone_extinct") == "0")
	 	{
	 		clone_extinct = false;
	 	}
	 	else
	 	{
	 		clone_extinct = true;
	 	}
	 }

	 if( rawData -> count("Number_of_Mutations") )
	 {
	 	Number_of_Mutations = std::stoull( rawData -> at("Number_of_Mutations") ) ;
	 }

	 if( rawData -> count("Clone_Size") )
	 {
	 	Clone_Size = std::stoull( rawData -> at("Clone_Size")  );
	 }

	 if( rawData -> count("Initiall_Expasion_Period") )
	 {
	 	if(rawData -> at("Initiall_Expasion_Period") == "0")
	 	{
	 		Initiall_Expasion_Period = false;
	 	}
	 	else
	 	{
	 		Initiall_Expasion_Period = true;
	 	}
	 }

	 if( rawData -> count("In_G0_Phase") )
	 {
	 	if(rawData -> at("In_G0_Phase") == "0")
	 	{
	 		In_G0_Phase = false;
	 	}
	 	else
	 	{
	 		In_G0_Phase = true;
	 	}
	 }

	  if( rawData -> count("In_G1_Phase") )
	  {
	 	if(rawData -> at("In_G1_Phase") == "0")
	 	{
	 		In_G1_Phase = false;
	 	}
	 	else
	 	{
	 		In_G1_Phase = true;
	 	}
	 }

	  if( rawData -> count("In_G2_Phase") )
	  {
	 	if(rawData -> at("In_G2_Phase") == "0")
	 	{
	 		In_G2_Phase = false;
	 	}
	 	else
	 	{
	 		In_G2_Phase = true;
	 	}
	 }

	 if( rawData -> count("In_S_Phase") )
	 {
	 	if(rawData -> at("In_S_Phase") == "0")
	 	{
	 		In_S_Phase = false;
	 	}
	 	else
	 	{
	 		In_S_Phase = true;
	 	}
	 }

	 if( rawData -> count("In_M_Phase") )
	 {
	 	if(rawData -> at("In_M_Phase") == "0")
	 	{
	 		In_M_Phase = false;
	 	}
	 	else
	 	{
	 		In_M_Phase = true;
	 	}
	 }


	 if( rawData -> count("Remaining_Time_in_G1_Phase") )
	 	std::stringstream(rawData -> at("Remaining_Time_in_G1_Phase")) >> Remaining_Time_in_G1_Phase ;

	 if( rawData -> count("Remaining_Time_in_S_Phase") )
	 	std::stringstream(rawData -> at("Remaining_Time_in_S_Phase")) >> Remaining_Time_in_S_Phase ;

	 if( rawData -> count("Remaining_Time_in_G2_Phase") )
	 	std::stringstream(rawData -> at("Remaining_Time_in_G2_Phase")) >> Remaining_Time_in_G2_Phase ;

	 if( rawData -> count("Remaining_Time_in_M_Phase") )
	 	std::stringstream(rawData -> at("Remaining_Time_in_M_Phase")) >> Remaining_Time_in_M_Phase ;

	 if(rawData -> count("init_PR"))
	 	init_PR = std::stod(rawData -> at ("init_PR"));

	 if(rawData -> count("max_PR"))
	 	max_PR = std::stod(rawData -> at ("max_PR"));

	 if(rawData -> count("final_PR"))
	 	final_PR = std::stod(rawData -> at ("final_PR"));

	 if( rawData -> count("Generation_ID") )
	 {
	 	if ( (rawData -> at("Generation_ID")) == "-" )
	 	{
	 		Generation_ID = "";
	 	}
	 	else
	 	{
	 		Generation_ID = rawData -> at("Generation_ID");
	 	}
	 }

	 if(_Mutation_Rate_flag)
	 {
	 	P_Expansion[2] = _Mutation_Rate;
	 	P_Expansion[3] = 1.0 - (_Mutation_Rate + DEATH_RATE + PROLIFERATION_RATE );
	 }

	 if(_Proliferation_Rate_flag)
	 {
	 	P_Expansion[1] = _Proliferation_Rate;
	 	P_Expansion[3] = 1.0 - (DEATH_RATE + _Proliferation_Rate + MUT_RATE);
	 }

	 if(_Death_Rate_flag)
	 {
	 	P_Expansion[0] = _Death_Rate;
	 	P_Expansion[3] = 1.0 - (_Death_Rate + PROLIFERATION_RATE + MUT_RATE);
	 }

	 if(_Mutation_Rate_flag && _Proliferation_Rate_flag && _Death_Rate_flag)
	 	P_Expansion[3] = 1.0 - (_Mutation_Rate + _Proliferation_Rate + _Death_Rate);

	 if(_Mutation_Rate_flag && _Proliferation_Rate_flag )
	 	P_Expansion[3] = 1.0 - (DEATH_RATE + _Proliferation_Rate + _Mutation_Rate);

	 if(_Mutation_Rate_flag && _Death_Rate_flag )
	 	P_Expansion[3] = 1.0 - (_Death_Rate + PROLIFERATION_RATE + _Mutation_Rate);

	 if(_Proliferation_Rate_flag && _Death_Rate_flag )
	 	P_Expansion[3] = 1.0 - (_Death_Rate + _Proliferation_Rate + MUT_RATE);

	 
}

 void Clone::FileReader(std::string path)
 {
 	if( path.empty() )
 	{
 		std::cout << "NO INPUT PROVIDED, reading default file section..." << std::endl;
 		//Read defulkat or call default parameters
 	}
 	else
 	{
 		std::cout << "[ CLONAL FILE ] " << path << std::endl;
 		if( FileExists(path) )
 		{
 			std::map<std::string, std::string> rawData  = ReadTextFile(path);
 			LoadElementsFromMap(&rawData);
 		}
 	}
 }

 void Clone::printValues (void)
 {
 	std::cout << "\n\n################# C L O N E   V A L U E S  #################" << "\n\n";
 	std::cout << "################# G E N E R A L  V A L U E S  #################" << "\n\n";
 	std::cout << "Number of Members to Start Seeding Heterogeneity: " << Number_of_Memebers_to_Start_Heterogeneity << "\n";
 	std::cout << "Generation Branch ID: " << Generation_ID_Counter << "\n";
 	std::cout << "Log(10) Penalty Accumulation: " << Driver_10_fold_accumulation << "\n";
 	std::cout << "Is the Clone extinct: " << clone_extinct << "\n";
 	std::cout << "Initial Proliferation Rate: " << init_PR << "\n";
 	std::cout << "Maximum Proliferation Rate Achieved: " << max_PR << "\n";
 	std::cout << "Final Proliferation Rate: " << final_PR << "\n";
 	std::cout << "\n";
 	std::cout << "################# M I T O S I S #################" << "\n\n";
 	std::cout << "Is In Initial Expansion Period: " << Initiall_Expasion_Period << "\n";
 	std::cout << "In G0: " << In_G0_Phase << "\n";
 	std::cout << "In G1: " << In_G1_Phase << "\n";
 	std::cout << "In S: " << In_S_Phase << "\n";
 	std::cout << "In G2: " << In_G2_Phase << "\n";
 	std::cout << "In M: " << In_M_Phase << "\n";
 	std::cout << "Remaining Time in G1: " << Remaining_Time_in_G1_Phase << "\n";
 	std::cout << "Remaining Time in G2: " << Remaining_Time_in_G2_Phase << "\n";
 	std::cout << "Remaining Time in S: " << Remaining_Time_in_S_Phase << "\n";
	std::cout << "Remaining Time in M: " << Remaining_Time_in_M_Phase << "\n";
 	std::cout << "\n";
 	std::cout << "################# E X P A N S I O N  #################" << "\n\n";
 	std::cout << "Clone Size: " << Clone_Size << "\n";
 	std::cout << "Death Rate : " << P_Expansion[0] << "\n";
 	std::cout << "Proliferation Rate : " << P_Expansion[1] << "\n";
 	std::cout << "Mutation Rate : " << P_Expansion[2] << "\n";
 	std::cout << "\n";
 	std::cout << "############################################" << "\n\n";

 }


 
 


