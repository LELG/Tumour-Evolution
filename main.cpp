#include "Clone.h"
#include "ClonalExpansion.h"
#include "clonalFun.h"
#include "config.h"
#include "Random.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <stdexcept>        // std::exception, std::runtime_error
#include <memory>           // std::unique_ptr

int main(int argc, char**argv)
{

	std::unique_ptr<Clone>  cln = get_Clone_DS();
  	
  	cln -> FileReader("settings.txt");
  	//cln -> printValues();

  	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
  
  	tmr -> printParameters();

  	//this will override input params, validate that
  	tmr -> carcinogenesis();

  	tmr -> printParameters();

  	tmr -> Tumour -> at (0) -> printValues();

  	r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  	long seed = time (NULL) * getpid();    
  	gsl_rng_set (r_global, seed);  

  	Random r;



  	std::cout << "TESTING pointer functions " << "\n";

  	unsigned long long int cs = 10;
  	std::cout << "STANDARD " << Division_Model(cs) << "\n";
  	unsigned long long int csp = 10;
  	Division_Model(&csp);
  	std::cout<<"POINTER " << csp << "\n";

	std::cout << "\nTESTING with bool " << "\n";

	unsigned int death_cells = 0;
	unsigned int newborne_cells = 100;

	std::cout << "NORMAL " << check_Clonal_Extinction(tmr -> Tumour -> at (0) -> Clone_Size, death_cells,   newborne_cells) << "\n";

	death_cells = 1000;
	std::cout << "PTR  " << check_Clonal_Extinction(&tmr -> Tumour -> at (0) -> Clone_Size, &death_cells,   &newborne_cells) << "\n";
  	

	std::cout << "Feedback " << tmr -> feedback << "\n";
	tmr -> feedback = 0.3;
	double NewBorn_probability = 0.2;
	std::cout << "BEFORE NewBorn_probability " << NewBorn_probability << "\n";

 	Enviromental_Penalty(&tmr -> feedback, &NewBorn_probability);

 	std::cout << "AFTER NewBorn_probability " << NewBorn_probability << "\n";

 	tmr -> feedback = 0.1;
 	NewBorn_probability = 0.2;

 	Enviromental_Penalty(&tmr -> feedback, &NewBorn_probability);

 	std::cout << "AFTER 2 NewBorn_probability " << NewBorn_probability << "\n";

 	std::cout << "\nTESTING BASIC CLONAL Clonal_Expansion " << "\n";


 	tmr -> feedback = 0.000000001;
 	NewBorn_probability = 0.03;

 	int clone_indx = 0;
 	unsigned int hours = 10;
 	unsigned long long int dying_cells = 10;
 	unsigned long long int newborn_cells = 500;
 	unsigned long long int mutant_cells = 40;
 	bool mutant_flag = false;
 	bool extinction = true;
 	tmr -> Tumour -> at (clone_indx) -> Clone_Size = 100000000;

 	Basic_Clonal_Expansion( tmr, 
							&clone_indx, 
							&hours,
							&dying_cells,
							&newborn_cells,
							&mutant_cells,
							&mutant_flag,
							&extinction);

std::string newGeneration_ID = "" ;
int Parent_Generation_ID_Counter = 0;
std::string Parent_ID = "P-0:0";
unsigned int years=0;

std::cout << "Parent ID " << Parent_ID << "\n";

 Generate_Clone_Generation_ID(&newGeneration_ID,
								  &Parent_Generation_ID_Counter, 
								  Parent_ID, 
								 &years, 
								  &hours);

std::cout << "NEW GEN ID " << newGeneration_ID << "\n";

Parent_ID = newGeneration_ID;
years = 10;
hours = 2;

Parent_Generation_ID_Counter = 5;

 Generate_Clone_Generation_ID(&newGeneration_ID,
								  &Parent_Generation_ID_Counter, 
								  Parent_ID, 
								 &years, 
								  &hours);


std::cout << "NEW GEN ID " << newGeneration_ID << "\n";

std::cout << "CARCINOGENSIS RROM DRIVER " << newGeneration_ID << "\n";

 carcinogensis_from_Driver(tmr,
						   &clone_indx,
						   &years,
						   &hours);

tmr -> Tumour -> at( tmr -> Tumour -> size() -1 ) -> printValues();

//delete tmr;
tmr.reset();
if(tmr)
std::cout<<"YES" << std::endl;
else
std::cout << "NO" << std::endl;

tmr.get();
tmr = get_Clonal_Expasion_DS();

if(tmr)
std::cout<<"YES" << std::endl;
else
std::cout << "NO" << std::endl;

int j = 0;
int gr = 0;
while(j < 100)
{
	Compute_Tumour_Growth( tmr, ".", 1);

	if(tmr -> Population_Size > 1)
	{
		std::cout <<  " Active Population Size: " << tmr -> Population_Size << std::endl;
		j++;
		tmr.reset();
		tmr.get();
		tmr = get_Clonal_Expasion_DS();
		std::cout <<  " Active Population Size Removed : " << tmr -> Population_Size << std::endl;
		j++;
		
	}
	else
	{
		std::cout <<  " \nTumour did not grow! " << std::endl;
		gr++;
		
	}
}
std::cout << "Successfull growth ratio: " << j << " : " << j+gr << std::endl;


  return 0;
}