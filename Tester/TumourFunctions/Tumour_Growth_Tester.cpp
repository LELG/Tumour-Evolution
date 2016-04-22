#include "Clone.h"
#include "ClonalExpansion.h"
#include "clonalFun.h"
#include "config.h"
#include "Random.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <stdexcept>        // std::exception, std::runtime_error
#include <memory>           // std::unique_ptr
#include <cassert>
#include <limits>


int main(int argc, char**argv)
{
	//***** SET UP Random ********//
	long seed;
	r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  	seed = time (NULL) * getpid();    
  	gsl_rng_set (r_global, seed);  
  	//***** SET UP Random ********//


	std::cout << "TESTING MITOSIS " << std::endl;

	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
	tmr -> carcinogenesis();

	tmr -> Tumour -> at(0) -> Available_cells = 5000;

	tmr -> Tumour -> at(0) -> G0_cells = 5;
	std::get<P_STAYING>( tmr -> Tumour -> at(0) -> G0_status ) = 0.40;
	std::get<P_DYING>( tmr -> Tumour -> at(0) -> G0_status ) = 0.10;

	tmr -> Tumour -> at(0) -> G1_cells = 5;
	std::get<P_STAYING_G1>( tmr -> Tumour -> at(0) -> G1_status ) = 0.40;
	std::get<P_DYING_G1>( tmr -> Tumour -> at(0) -> G1_status ) = 0.10;

	 tmr -> Tumour -> at(0) -> G2_cells = 5;
	 std::get<P_STAYING>( tmr -> Tumour -> at(0) -> G2_status ) = 0.40;
	 std::get<P_DYING>( tmr -> Tumour -> at(0) -> G2_status ) = 0.10;


	 tmr -> Tumour -> at(0) -> S_cells = 5;
	 std::get<P_STAYING>( tmr -> Tumour -> at(0) -> S_status ) = 0.40;
	 std::get<P_DYING>( tmr -> Tumour -> at(0) -> S_status ) = 0.10;

	 tmr -> Tumour -> at(0) -> M_cells = 5;
	 std::get<P_STAYING>( tmr -> Tumour -> at(0) -> M_status ) = 0.40;
	 std::get<P_DYING>( tmr -> Tumour -> at(0) -> M_status ) = 0.10;


	 tmr -> Tumour -> at(0) -> Clone_Size = static_cast<unsigned long long int>(tmr -> Tumour -> at(0) -> Available_cells + 
	 																			tmr -> Tumour -> at(0) -> G0_cells        + 
	 																			tmr -> Tumour -> at(0) -> G1_cells        + 
	 																			tmr -> Tumour -> at(0) -> G2_cells        + 
	 																			tmr -> Tumour -> at(0) -> S_cells         + 
	 																			tmr -> Tumour -> at(0) -> M_cells);


	
	std::cout << "G0 cells: " << tmr -> Tumour -> at(0) -> G0_cells << std::endl;
	std::cout << "G1 cells: " << tmr -> Tumour -> at(0) -> G1_cells << std::endl;
	std::cout << "G2 cells: " << tmr -> Tumour -> at(0) -> G2_cells << std::endl;
	std::cout << "S cells: " << tmr -> Tumour -> at(0) -> S_cells << std::endl;
	std::cout << "M cells: " << tmr -> Tumour -> at(0) -> M_cells << std::endl;
	std::cout << "Available_cells: " << tmr -> Tumour -> at(0) -> Available_cells << std::endl;
	std::cout << "Clone Size " << tmr -> Tumour -> at(0) -> Clone_Size << std::endl;
	std::cout << "P_G0(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> G0_status ) << std::endl;
	std::cout << "P_G0(Dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> G0_status ) << std::endl;
	std::cout << "P_G1(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> G1_status ) << std::endl;
	std::cout << "P_G1(dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> G1_status ) << std::endl;
	std::cout << "P_G2(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> G2_status ) << std::endl;
	std::cout << "P_G2(dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> G2_status ) << std::endl;
	std::cout << "P_S(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> S_status ) << std::endl;
	std::cout << "P_S(dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> S_status ) << std::endl;
	std::cout << "P_M(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> M_status ) << std::endl;
	std::cout << "P_M(staying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> M_status ) << std::endl;
	std::cout << "\n\n";

	//unsigned int years = 3;
	//unsigned int hours = 100;
	//unsigned int ith_clone = 0;
	//tmr -> Check_Mitosis_Network_Status( ith_clone, years, hours ) ;

	//std::cout << "\n\n";
	//unsigned int * Random::Update_S_Phase(unsigned int & S_cells, double & p_staying_S)

	unsigned int seconds = 0;
	unsigned int hours = 0;
	unsigned int years = 0;
	unsigned int elapsed_hours = 0;
	unsigned int times_to_wait = 0;
	unsigned int ith_clone = 0; 

	while( ((times_to_wait < STOP_GROWTH_AFTER_DIAGNOSIS) && (years < 60) ) && (tmr -> Population_Size > 0) )
	{

		seconds += dt;
		if(seconds == 3600)
		{
			seconds = 0; hours ++;
			for (ith_clone = 0; ith_clone < tmr -> Tumour -> size() ; ith_clone++) 			//(1) Is my clone not extinct?	
			{
				tmr -> Check_Mitosis_Network_Status( ith_clone, years, hours ) ;
				std::cout << "\n\n" << std::endl;
				getchar();
			}

		}

		if(hours == 8764)//8764
		{
			hours = 0; years ++;
			
		}

	}




	return 0;
}// End of main







