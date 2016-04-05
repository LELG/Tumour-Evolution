#include "Clone.h"
#include "config.h"
#include "inout_funs.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <memory>           // std::unique_ptr
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>


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
	P_Expansion{DEATH_RATE, PROLIFERATION_RATE, MUT_RATE, 1.0-(DEATH_RATE + PROLIFERATION_RATE + MUT_RATE)},
	Remaining_Time_in_G1_Phase(0),
	Remaining_Time_in_S_Phase(0),
	Remaining_Time_in_G2_Phase(0),
	Remaining_Time_in_M_Phase(0),
	init_PR(0.0),
	max_PR(0.0),
	final_PR(0.0),
	Generation_ID("")
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
	 	Mutation_Rate = Mutation_Rate;
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
	 	if(rawData -> at("clone_extinct") == "0")
	 		clone_extinct = false;
	 	else
	 		clone_extinct = true;

	 if( rawData -> count("Number_of_Mutations") )
	 	Number_of_Mutations = std::stoull( rawData -> at("Number_of_Mutations") ) ;

	 if( rawData -> count("Clone_Size") )
	 	Clone_Size = std::stoull( rawData -> at("Clone_Size")  );

	 if( rawData -> count("Initiall_Expasion_Period") )
	 	if(rawData -> at("Initiall_Expasion_Period") == "0")
	 		Initiall_Expasion_Period = false;
	 	else
	 		Initiall_Expasion_Period = true;

	 if( rawData -> count("In_G0_Phase") )
	 	if(rawData -> at("In_G0_Phase") == "0")
	 		In_G0_Phase = false;
	 	else
	 		In_G0_Phase = true;

	  if( rawData -> count("In_G1_Phase") )
	 	if(rawData -> at("In_G1_Phase") == "0")
	 		In_G1_Phase = false;
	 	else
	 		In_G1_Phase = true;

	  if( rawData -> count("In_G2_Phase") )
	 	if(rawData -> at("In_G2_Phase") == "0")
	 		In_G2_Phase = false;
	 	else
	 		In_G2_Phase = true;

	 if( rawData -> count("In_S_Phase") )
	 	if(rawData -> at("In_S_Phase") == "0")
	 		In_S_Phase = false;
	 	else
	 		In_S_Phase = true;

	 if( rawData -> count("In_M_Phase") )
	 	if(rawData -> at("In_M_Phase") == "0")
	 		In_M_Phase = false;
	 	else
	 		In_M_Phase = true;


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
	 		Generation_ID = "";
	 	else
	 		Generation_ID = rawData -> at("Generation_ID");
	 }
	 	

 // 	std::cout << Number_of_Memebers_to_Start_Heterogeneity << std::endl;
	// std::cout << "IN " << std::endl;
	// std::cout << _Mutation_Rate << std::endl;
	// std::cout << _Proliferation_Rate << std::endl;
	// std::cout << _Death_Rate << std::endl;
	// std::cout << Number_of_Mutations << std::endl;
	// std::cout << Remaining_Time_in_G1_Phase << std::endl;
	// std::cout << Generation_ID << std::endl;


}


// TODO: 
// Implement a read file to run withouth default conditions
 void Clone::FileReader(std::string path)
 {
 	if( path.empty() )
 	{
 		std::cout << "NO IMPUT PROVIDED, reading default file section..." << std::endl;
 		//Read defulkat or call default parameters
 	}
 	else
 	{
 		std::cout << "[ CLONAL FILE ] " << path << std::endl;
 		if( FileExists(path) )
 		{
 			//Validate content
 			std::map<std::string, std::string> rawData  = ReadTextFile(path);

 			std::cout << rawData["Number_of_Memebers_to_Start_Heterogeneity"] << std::endl;
 			LoadElementsFromMap(&rawData);
 		}
 		//std::cout << FileExists(path) 
 		// Validate file
 		// Check that everything opens and set data
 	}
 
 }


