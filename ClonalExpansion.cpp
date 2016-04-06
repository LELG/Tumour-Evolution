#include "ClonalExpansion.h"
#include "Clone.h"
#include "clonalFun.h"
#include <map>
#include <memory>           // std::unique_ptr
#include <string>
#include <vector>
#include <iostream>         // std::cout,std::endl


Clonal_Expansion::Clonal_Expansion()
 :
 	Population_Size(0),
 	feedback(0.0)
 {       }

Clonal_Expansion::~Clonal_Expansion()
 {       }

 void Clonal_Expansion::printParameters(void)
 {
 	std::cout << "\n\n##################################\n\n CLONAL EXPANSION DS \n{ " << std::endl;
 	std::cout << "\tSelective Pressure: " << feedback << std::endl;
 	std::cout << "\tPopulation Size: " << Population_Size << std::endl;
 	std::cout << "\tNumber of Clones: " << Tumour -> size() << std::endl;
 	std::cout << "}\n#################################" << std::endl; 
 }


 void Clonal_Expansion::add_Clone_to_Tumour(void)
 {
 	Tumour -> push_back ( get_Clone_DS() );
 }

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

	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[1];

	Population_Size++;
}