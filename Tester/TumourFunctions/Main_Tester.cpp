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
#include <limits>
#include <chrono>
#include <thread>
#include <mpi.h>	// Needs to be included in order to use MPI

int main(int argc, char**argv)
{
	//*************** SET UP Random ********//
	// long seed;
	// r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
 //  	seed = time (NULL) * getpid();    
 //  	gsl_rng_set (r_global, seed);  
  	//*************** SET UP Random ********//

  	int			myID;
	int			N_Procs;

	unsigned int replicates = 10;

  	
  	std::map<std::string, std::string> logic;
  	std::unique_ptr<Clonal_Expansion>  Primary ;
  	std::unique_ptr<Clonal_Expansion>  Metastasis ;

  	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,	&myID); 
	MPI_Comm_size(MPI_COMM_WORLD,	&N_Procs);
  	


  	Load_Logic_File(logic);

  	std::string data_path = "";

  	if(logic.size() >= 1)
	{
		print_logic_file(logic);
	}

	//std::cout << "VERSION TO RUN: " << logic["Version"] << std::endl;

	// Growth
for(unsigned int  i = 0; i < replicates; i++)
{
	if(myID == 0) {std::cout << "Type of Growth: " ;}
	
	if(logic["Growth"] == "Primary")
	{
		if(myID == 0) {std::cout << "Primary \n";}
		Primary = get_Clonal_Expasion_DS();
	}
	else
	{
		if(myID == 0) {std::cout << " Multiple \n";}
		//lop with metastasis or multiple primary tumours
	}

	// Input
	if(myID == 0){std::cout << "Input of Growth: " ;}
	if(logic["Input"] == "carcinogenesis")
	{
		if(myID == 0){std::cout << "carcinogenesis \n";}
		if(myID == 0){std::cout << "VERSION TO RUN: " << logic["Version"] << std::endl;
	}
		
		Primary -> setVersion_type(logic["Version"] );
		
		if(myID == 0){std::cout << Primary -> getVersion_type() << " associated with ID " <<  Primary -> Version << std::endl;}

		Primary -> Compute_Tumour_Growth(logic, i, myID, data_path );/// can take out of here
		
		//Primary -> Select_Carcinogenesis();
	}
	else
	{
		if(myID == 0){ std::cout << " File \n";}
	}

	Primary.reset();
	Primary.get();
}


	//std::unique_ptr<Clonal_Expansion>  Primary = get_Clonal_Expasion_DS(); // calling clonalfun 



	//Primary -> carcinogenesis(); 

	std::cout << "END DATA PROCESS: " << myID << std::endl;

	// tmr -> Tumour -> at(0) -> Available_cells = 0;

	// tmr -> Tumour -> at(0) -> G0_cells = 0;
	// std::get<P_STAYING>( tmr -> Tumour -> at(0) -> G0_status ) = 0.70;
	// std::get<P_DYING>( tmr -> Tumour -> at(0) -> G0_status ) = 0.10;

	// tmr -> Tumour -> at(0) -> G1_cells = 100;
	// std::get<P_STAYING_G1>( tmr -> Tumour -> at(0) -> G1_status ) = 0.70;
	// std::get<P_DYING_G1>( tmr -> Tumour -> at(0) -> G1_status ) = 0.016;

	//  tmr -> Tumour -> at(0) -> G2_cells = 0;
	//  std::get<P_STAYING>( tmr -> Tumour -> at(0) -> G2_status ) = 0.70;
	//  std::get<P_DYING>( tmr -> Tumour -> at(0) -> G2_status ) = 0.016;


	//  tmr -> Tumour -> at(0) -> S_cells = 0;
	//  std::get<P_STAYING>( tmr -> Tumour -> at(0) -> S_status ) = 0.70;
	//  std::get<P_DYING>( tmr -> Tumour -> at(0) -> S_status ) = 0.016;

	//  tmr -> Tumour -> at(0) -> M_cells = 0;
	//  std::get<P_STAYING>( tmr -> Tumour -> at(0) -> M_status ) = 0.70;
	//  std::get<P_DYING>( tmr -> Tumour -> at(0) -> M_status ) = 0.016;


	//  tmr -> Tumour -> at(0) -> Clone_Size = static_cast<unsigned long long int>(tmr -> Tumour -> at(0) -> Available_cells + 
	//  																			tmr -> Tumour -> at(0) -> G0_cells        + 
	//  																			tmr -> Tumour -> at(0) -> G1_cells        + 
	//  																			tmr -> Tumour -> at(0) -> G2_cells        + 
	//  																			tmr -> Tumour -> at(0) -> S_cells         + 
	//  																			tmr -> Tumour -> at(0) -> M_cells);



	//  std::cout << "NUMBER LIMITS\t"
 //              << std::numeric_limits<unsigned int>::lowest() << '\t'
 //              << std::numeric_limits<unsigned int>::max() << '\n';
	
	// std::cout << "G0 cells: " << tmr -> Tumour -> at(0) -> G0_cells << std::endl;
	// std::cout << "G1 cells: " << tmr -> Tumour -> at(0) -> G1_cells << std::endl;
	// std::cout << "G2 cells: " << tmr -> Tumour -> at(0) -> G2_cells << std::endl;
	// std::cout << "S cells: " << tmr -> Tumour -> at(0) -> S_cells << std::endl;
	// std::cout << "M cells: " << tmr -> Tumour -> at(0) -> M_cells << std::endl;
	// std::cout << "Available_cells: " << tmr -> Tumour -> at(0) -> Available_cells << std::endl;
	// std::cout << "Clone Size " << tmr -> Tumour -> at(0) -> Clone_Size << std::endl;
	// std::cout << "P_G0(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> G0_status ) << std::endl;
	// std::cout << "P_G0(Dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> G0_status ) << std::endl;
	// std::cout << "P_G1(staying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> G1_status ) << std::endl;
	// std::cout << "P_G1(dying) cells: " << std::get<5>( tmr -> Tumour -> at(0) -> G1_status ) << std::endl;
	// std::cout << "P_G2(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> G2_status ) << std::endl;
	// std::cout << "P_G2(dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> G2_status ) << std::endl;
	// std::cout << "P_S(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> S_status ) << std::endl;
	// std::cout << "P_S(dying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> S_status ) << std::endl;
	// std::cout << "P_M(staying) cells: " << std::get<3>( tmr -> Tumour -> at(0) -> M_status ) << std::endl;
	// std::cout << "P_M(staying) cells: " << std::get<4>( tmr -> Tumour -> at(0) -> M_status ) << std::endl;
	// std::cout << "\n\n";

	// //unsigned int years = 3;
	// //unsigned int hours = 100;
	// //unsigned int ith_clone = 0;
	// //tmr -> Check_Mitosis_Network_Status( ith_clone, years, hours ) ;

	// //std::cout << "\n\n";
	// //unsigned int * Random::Update_S_Phase(unsigned int & S_cells, double & p_staying_S)

	// unsigned int seconds = 0;
	// unsigned int hours = 0;
	// unsigned int years = 0;
	// //unsigned int elapsed_hours = 0;
	// //unsigned int times_to_wait = 0;
	// unsigned int ith_clone = 0; 

	// Random r;
	// std::ofstream te_file;
	// te_file.open ("./Mitosis_Model/te_file_31.txt");

	// unsigned long long int elapsed_hours = 0;

	// //getchar(); 
	// while( tmr -> Population_Size < 4000000000 && !(tmr -> Population_Size == 0))
	// {
	// 	std::cout << "WHILE IN G0 >  " <<  tmr -> Tumour -> at( 0 ) -> G0_cells << std::endl;

	// 	seconds += dt;
	// 	if(seconds == 3600)
	// 	{
	// 		seconds = 0; hours ++;
	// 		for (ith_clone = 0; ith_clone < tmr -> Tumour -> size() ; ith_clone++) 			//(1) Is my clone not extinct?	tmr -> Tumour -> size()
	// 		{	
	// 			std::cout << "TIME [Y:H] " << years << " : " << hours << std::endl;
	// 			tmr -> Check_Mitosis_Network_Status( ith_clone, years, hours, r) ;
	// 		}

	// 		tmr -> Update_Tumour_Size();
	// 		tmr -> map_Feedback();

	// 	}
	// 	//std::this_thread::sleep_for(std::chrono::nanoseconds(30000000));
	// 	if(hours == 8764)//8764
	// 	{
	// 		hours = 0; years ++;
			
	// 	}

	// 	//SDave to file 
	// 	//Pop size hours and heterogeneity
	// 	te_file << tmr -> Population_Size << "\t" << elapsed_hours << "\t"<< tmr -> Tumour -> size()   << "\n";
	// 	elapsed_hours++;

	// }

	// std::cout << "Clone Sizes: " << std::endl;
	// for (ith_clone = 0; ith_clone < tmr -> Tumour -> size() ; ith_clone++) 			//(1) Is my clone not extinct?	tmr -> Tumour -> size()
	// {	
	// 	std::cout << "C_S[ " << ith_clone << " ] = " << tmr -> Tumour -> at( ith_clone ) -> Clone_Size <<  "\t\t PR: " <<  tmr -> Tumour -> at (ith_clone) -> P_Expansion[1] << " MR: " << tmr -> Tumour -> at (ith_clone) -> P_Expansion[2]  <<std::endl;
	// }

	// te_file.close();


	MPI_Finalize();
	return 0;
}// End of main







