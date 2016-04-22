#include "ClonalExpansion.h"
#include "Clone.h"
#include "config.h"
#include "clonalFun.h"
#include "Random.h"
#include <map>
#include <memory>           // std::unique_ptr
#include <string>
#include <vector>
#include <iostream>         // std::cout,std::endl


Clonal_Expansion::Clonal_Expansion()
 :
 	Population_Size(0),
 	feedback(0.0),
 	mitosis(Clonal_Expansion::Standard)
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

 void Clonal_Expansion::printValues (std::unique_ptr<Clone> const & ith_clone)
 {
 	std::cout << "\n\n################# C L O N E   V A L U E S  #################" << "\n\n";
 	std::cout << "################# G E N E R A L  V A L U E S  #################" << "\n\n";
 	std::cout << "Number of Members to Start Seeding Heterogeneity: " << ith_clone -> Number_of_Memebers_to_Start_Heterogeneity << "\n";
 	std::cout << "Generation Branch ID: " << ith_clone ->  Generation_ID_Counter << "\n";
 	std::cout << "Log(10) Penalty Accumulation: " << ith_clone ->  Driver_10_fold_accumulation << "\n";
 	std::cout << "Is the Clone extinct: " << ith_clone ->  clone_extinct << "\n";
 	std::cout << "Initial Proliferation Rate: " << ith_clone ->  init_PR << "\n";
 	std::cout << "Maximum Proliferation Rate Achieved: " << ith_clone ->  max_PR << "\n";
 	std::cout << "Final Proliferation Rate: " << ith_clone -> final_PR << "\n";
 	std::cout << "\n";
 	std::cout << "################# M I T O S I S #################" << "\n\n";
 	std::cout << "Is In Initial Expansion Period: " <<ith_clone ->  Initiall_Expasion_Period << "\n";
 	std::cout << "In G0: " << ith_clone -> In_G0_Phase << "\n";
 	std::cout << "In G1: " << ith_clone -> In_G1_Phase << "\n";
 	std::cout << "In S: " << ith_clone -> In_S_Phase << "\n";
 	std::cout << "In G2: " << ith_clone -> In_G2_Phase << "\n";
 	std::cout << "In M: " << ith_clone -> In_M_Phase << "\n";
 	std::cout << "Remaining Time in G1: " << ith_clone -> Remaining_Time_in_G1_Phase << "\n";
 	std::cout << "Remaining Time in G2: " << ith_clone -> Remaining_Time_in_G2_Phase << "\n";
 	std::cout << "Remaining Time in S: " << ith_clone -> Remaining_Time_in_S_Phase << "\n";
	std::cout << "Remaining Time in M: " << ith_clone -> Remaining_Time_in_M_Phase << "\n";
 	std::cout << "\n";
 	std::cout << "################# E X P A N S I O N  #################" << "\n\n";
 	std::cout << "Clone Size: " << ith_clone -> Clone_Size << "\n";
 	std::cout << "Death Rate : " << ith_clone ->  P_Expansion[0] << "\n";
 	std::cout << "Proliferation Rate : " << ith_clone -> P_Expansion[1] << "\n";
 	std::cout << "Mutation Rate : " << ith_clone -> P_Expansion[2] << "\n";
 	std::cout << "\n";
 	std::cout << "############################################" << "\n\n";

 }


 void Clonal_Expansion::add_Clone_to_Tumour(void)
 {
 	Tumour -> push_back ( get_Clone_DS() );
 }

std::string Clonal_Expansion::getMitosis_type(void)
 {
 	std::string types[] = {"Standard", "Binomial", "Network"};
 	return types[mitosis]; 
 }

/**
	Adds a parent clone with size 1 in the tumour.
	Steps:
	1) Ask for a new Clone DS and append to vector Tumour
	2) From the recently added clone, initialise values

*******************************************/
void Clonal_Expansion::carcinogenesis(void)
{
	// 1) Add a Clone
	Random r;
	// Adding a new clone
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();

	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[P_Expansion_PR];
	Tumour -> back() -> init_PR  = Tumour -> back() -> P_Expansion[P_Expansion_PR];


	Population_Size++;
	mitosis = Clonal_Expansion::Standard; // Requires modification from input parameter
	
	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
}

void Clonal_Expansion::carcinogenesis_from_driver(unsigned int & ith_clone, unsigned int & years, unsigned int & hours)
{
	Random r;
	double mr = Tumour -> at( ith_clone ) -> Mutation_Rate;
	unsigned int AC = Tumour -> at( ith_clone ) -> Driver_10_fold_accumulation;
	unsigned int Accumulated_Drivers = Tumour -> at( ith_clone ) -> Driver_10_fold_accumulation;
	unsigned long long int NOM = Tumour -> at( ith_clone ) -> Number_of_Mutations;
	double pr = Tumour -> at( ith_clone ) -> P_Expansion[1];

	std::string cloneName = "";	
	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) ->Generation_ID_Counter;
	Generate_Clone_Generation_ID(cloneName, 
								 Parent_Generation_ID_Counter, 
								 Tumour -> at( ith_clone ) -> Generation_ID,
								 years, 
								 hours);

	Tumour -> push_back( get_Clone_DS() );

	Tumour -> back() -> Generation_ID = cloneName;
	Tumour -> back() -> Initiall_Expasion_Period = false;
	Tumour -> back() -> Clone_Size = 1;
	Tumour -> back() -> Number_of_Mutations = NOM + 1;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;
	
	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Generation_ID_Counter = 0;
	Tumour -> back() -> Mutation_Rate = r.Uniform_Mutation_Rate_2(mr);
	double PR = r.Update_Proliferation_Rate(pr);
	Tumour -> back() -> Driver_10_fold_accumulation = AC + 5;
	Tumour -> back() -> P_Expansion[1] =  PR;
	Tumour -> back() -> init_PR = PR;
	Tumour -> back() -> max_PR = PR;
	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
	

	Population_Size++;

	//Tumour -> back() -> printValues();


	//std::cout << "\n\n \t\t\t VS \n\n";

	//Tumour -> at( ith_clone ) -> printValues();
	
	//ith_clone = Tumour -> at(Tumour -> size() - 2 ) const;
	//printValues(Tumour -> at (Tumour -> size() - 2 ) );




}

/**
	Adds a parent clone with size neoplastic_cells in the tumour.
	Steps:
	1) Ask for a new Clone DS and append to vector Tumour
	2) From the recently added clone, initialise values

*******************************************/
void Clonal_Expansion::carcinogenesis(unsigned long long int &neoplastic_cells)
{
	Random r;

	Tumour -> push_back( get_Clone_DS() );
	Tumour -> back() -> Clone_Size = neoplastic_cells;

	Tumour -> back() -> In_G0_Phase = false;
	Tumour -> back() -> In_G1_Phase = true;
	Tumour -> back() -> In_S_Phase = false;
	Tumour -> back() -> In_G2_Phase = false;
	Tumour -> back() -> In_M_Phase = false;

	Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

	Tumour -> back() -> Generation_ID_Counter = 0; 
	Tumour -> back() -> Generation_ID = "P-0:0";
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[P_Expansion_PR];

	Population_Size = Tumour -> back() -> Clone_Size;
	mitosis = Clonal_Expansion::Standard;
		
	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
}

/**
	Check if the ith_clone is in G1 mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_G1_phase(std::unique_ptr<Clone> const & ith_clone)
{
	return (ith_clone -> In_G1_Phase && ith_clone -> Remaining_Time_in_G1_Phase !=0);
}


/**
	Check if the ith_clone is exiting the G1 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_G1_phase(std::unique_ptr<Clone> const & ith_clone)
{
	return (ith_clone -> In_G1_Phase && ith_clone -> Remaining_Time_in_G1_Phase == 0);
}

/**
	Check if the ith_clone is in S mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_S_phase(std::unique_ptr<Clone> const & ith_clone)
{
	return (ith_clone -> In_S_Phase && ith_clone -> Remaining_Time_in_S_Phase !=0);
}

/**
	Check if the ith_clone is exiting the G1 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_S_phase(std::unique_ptr<Clone> const & ith_clone)
{
	return (ith_clone -> In_S_Phase && ith_clone -> Remaining_Time_in_S_Phase == 0);
}

/**
	Check if the ith_clone is in G2 mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_G2_phase(std::unique_ptr<Clone> const & ith_clone)
{
	return (ith_clone -> In_G2_Phase && ith_clone -> Remaining_Time_in_G2_Phase !=0);
}

/**
	Check if the ith_clone is exiting the G2 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_G2_phase(std::unique_ptr<Clone> const & ith_clone)
{
	return (ith_clone -> In_G2_Phase && ith_clone -> Remaining_Time_in_G2_Phase == 0);
}

/**
	ith_Clone from G1 -> S phase
******************************************************/
void Clonal_Expansion::Transition_From_G1_S(std::unique_ptr<Clone> const & ith_clone)
{
	Random r;
	ith_clone -> In_G1_Phase  = false; 
	ith_clone -> Remaining_Time_in_G1_Phase = r.G1();
	ith_clone -> In_S_Phase = true;
}

/**
	ith_Clone from S -> G2 phase
******************************************************/
void Clonal_Expansion::Transition_From_S_G2(std::unique_ptr<Clone> const & ith_clone)
{
	Random r;
	ith_clone -> In_S_Phase  = false; 
	ith_clone -> Remaining_Time_in_S_Phase = r.S();
	ith_clone -> In_G2_Phase = true;
}

/**
	ith_Clone from G2 -> M phase
******************************************************/
void Clonal_Expansion::Transition_From_G2_M(std::unique_ptr<Clone> const & ith_clone)
{
	Random r;
	ith_clone -> In_G2_Phase  = false; 
	ith_clone -> Remaining_Time_in_G2_Phase = r.G2();
	ith_clone -> In_M_Phase = true;
}



// TODO Asymetric division
void Clonal_Expansion::Update_Population_After_Division(std::unique_ptr<Clone> const & ith_clone)
{

	unsigned long long int cells_after_division = ith_clone -> Clone_Size * 2;
	Population_Size =  (Population_Size -  ith_clone -> Clone_Size ) + cells_after_division;
	ith_clone -> Clone_Size = cells_after_division;
	std::cout << "Population_Size: " << Population_Size << std::endl;
}

void Clonal_Expansion::Estimate_Mutational_Effects(unsigned int & Mutant_Cells, unsigned long long int & Clone_Size) //years and hours
{

	Random r;

	unsigned int cs = static_cast<unsigned int>(Clone_Size);
	unsigned int * mutational_effect_buffer;
	double effect = 0.0;

	unsigned int kills = 0;
	unsigned int drivers = 0;
	double mutational_burden = 0.0;

	for(unsigned int i = 0; i < Mutant_Cells  ; i ++)
	{
		for(unsigned int j = 0; j < add_multiple_mutations() ; j++)
		{
			//std::cout << add_multiple_mutations() << " mutations to add for cell " << i << std::endl; 
			mutational_effect_buffer = r.Mutational_Effects();
			std::cout << "Killer Mutation " << mutational_effect_buffer[0]
					  << " Driver Mutation " <<  mutational_effect_buffer[1]
					  << " Gradient Effect "<<   mutational_effect_buffer[2]
					  << std::endl;


			if(mutational_effect_buffer[0] > 0)
			{
				std::cout << "Highly deleterious mutation... inducing cell death" << std::endl;
				kills++;
			}
			else if(mutational_effect_buffer[1] > 0)
			{
				std::cout << "Driver mutation... generating a new clone " << std::endl;
				drivers++;
			}
			else
			{
				r.Laplace(effect, cs );
				mutational_burden += effect;
				std::cout << " Mutation with effect of " << effect << std::endl;
			}

		}
	}

	std::cout<< " TOTAL KILLS " << kills << " DRIVERS " << drivers << " mutational burden " << mutational_burden << std::endl;  
	
}


void Clonal_Expansion::Check_for_Clonal_Extintion_at_min_size(unsigned int & ith_clone)
{

	if( Tumour -> at(ith_clone) -> Clone_Size == 1 )
	{
		Tumour -> at(ith_clone) -> Clone_Size = 0 ;
		Tumour -> at(ith_clone) -> clone_extinct = true; 
	}	

}



unsigned int * Clonal_Expansion::Update_Clonal_Mutational_Burden(unsigned int & ith_clone, unsigned int & Mutant_cells, unsigned int & years, unsigned int & hours )
{
	Random r;
	unsigned int * mutational_effect_buffer;
	unsigned int counter = 0;
	unsigned int death = 0;
	double effect_gain = static_cast<double>( Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] );

	unsigned int j = 0;


	while( (j < Mutant_cells  ) && !(Tumour -> at(ith_clone) -> clone_extinct) )
	{
		unsigned int i = 0;
		while( i < add_multiple_mutations() )
		{
			// check for mutations
			mutational_effect_buffer = r.Mutational_Effects();
			// buffer is size 3 --> need to check kills drivers and effects
			
			if ( mutational_effect_buffer[PASSENGER_CELLS] > 0 )// Passenger
			{
			 	//std::cout << " MUTATIONAL EFFECT " << std::endl; 
			 	unsigned int CS = static_cast<unsigned int>(Tumour -> at(ith_clone) -> Clone_Size);
			 	// std::cout << "Before Expansion " 
			 	// 		  << Tumour -> at(ith_clone) -> P_Expansion[1] << " " 
			 	// 		  << Tumour -> at(ith_clone) -> Clone_Size << std::endl;

			 	r.Laplace( Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR], CS );

			 	// std::cout << "After Expansion " 
			 	// 		  << Tumour -> at(ith_clone) -> P_Expansion[1] << " " 
			 	// 		  << Tumour -> at(ith_clone) -> Clone_Size << std::endl;

			 	counter++;
			}
			else if(mutational_effect_buffer[DRIVER_CELLS] > 0)// Driver
			{
				//std::cout << " DRIVER " << std::endl;
				carcinogenesis_from_driver( ith_clone,  years,  hours );

			}
			else 
			{
				//std::cout << " DEATH " << std::endl;
				Check_for_Clonal_Extintion_at_min_size( ith_clone );
				death++;
				if(Tumour -> at(ith_clone) -> clone_extinct)
					break;
			}
			i++;
		}//innner for
		j++;
	}// outer loop

	// VERY CAREFUL WITH THIS
	
	effect_gain = effect_gain -  Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] ;
	std::cout << "AFT Gain " << effect_gain << " counter " << counter << std::endl;
	std::cout << " Final PR " << Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] << std::endl;
	effect_gain = fabs(effect_gain);


	unsigned int MCC = Mutant_cells ;

	unsigned int * Mutations;

	
	//Mutations[0] = 0; Mutations[1] = 0; Mutations[2] = 0;

	if (!Tumour -> at(ith_clone) -> clone_extinct )
	{

		Mutations = r.Mutational_Proportions( effect_gain, Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] , MCC );
	
	}
	

	std::cout << "Values Muts: " << "Muts(G0) " << Mutations[MUTANTS_to_G0]
				  << "  Muts(G1) " << Mutations[MUTANTS_to_G1]
				  << " Muts(IDLE) " << Mutations[MUTANTS_to_IDLE]
				  << " Muts (Death) " << death 
				  << " Mutant cells "<< MCC 
				  << " effect gain " <<  effect_gain << std::endl;

	return Mutations;

}

void Clonal_Expansion::Probabilities_of_Cell_Division (unsigned int & ith_clone, double & p_idle, double & p_go_to_G0 )
{
	p_idle = 1.0 - Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR];
	p_go_to_G0 = fabs( Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR] - Tumour -> at( ith_clone ) -> max_PR);
}



//This can run in parallel 
// Mitosis
// https://www.youtube.com/watch?v=koudmJdil60
void Clonal_Expansion::Check_Mitosis_Network_Status(unsigned int & ith_clone, unsigned int & years, unsigned int & hours)
{
	// All of this runs
	Random r;
	

	// // for G0
	double comp_probability_G0 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) ); 
	double comp_probability_G1 = 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) ) - 0.1;
	double comp_probability_G2 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status ) );
	double comp_probability_S  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status ) );
	double comp_probability_M  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status ) );
	

	 unsigned int * G0_buffer = r.Update_G0_Phase( Tumour -> at( ith_clone ) -> G0_cells, 
	 							  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ), 
	 							  std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G0_status ),
	 							  comp_probability_G0);

	 unsigned int * G1_buffer = r.Update_G1_Phase( Tumour -> at( ith_clone ) -> G1_cells,
	 							   std::get<P_STAYING_G1> ( Tumour -> at( ith_clone ) -> G1_status ),
	 							   std::get<P_DYING_G1>   ( Tumour -> at( ith_clone ) -> G1_status),
	 							   comp_probability_G1);

	  unsigned int * G2_buffer = r.Update_G2_Phase( Tumour -> at( ith_clone ) -> G2_cells,
	  							   std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ),
	  							   std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G2_status ),
	  							   comp_probability_G2 );

	  unsigned int *  S_buffer = r.Update_S_Phase( Tumour -> at( ith_clone ) -> S_cells,
	  							  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ),
	  							  std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> S_status ),
	  							  comp_probability_S );

	 unsigned int *  M_buffer = r.Update_M_Phase( Tumour -> at( ith_clone ) -> M_cells,
	  							 std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status),
	  							 std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> M_status ),
	  							 comp_probability_M );


	 // TODO Change this
	 //double to_G0 = 0.001; // probability of expoansion 
	 unsigned int AC_at_t = Tumour -> at( ith_clone ) -> Available_cells;
	 unsigned int * Clonal_update_buffer = r.Clonal_Functions( Tumour -> at( ith_clone ) -> Available_cells,
	  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR],
	  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR], 
	  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR] ) ;

	

	// //for G0
	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[DIE];// Die 
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[NEXT_STAGE];// to G1
	//for G1
	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[DIE];// to Dying
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[NEXT_STAGE];// to G2
	std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[REVERT_TO_G0];// to G0
	//for G2
	 std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[STAY];// Stay
	 std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[DIE];// to Dying
	 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[NEXT_STAGE];
	// //for S
	 std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[STAY];// Stay
	 std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[DIE];// die
	 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[NEXT_STAGE];// to M
	// //for M
	 std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[STAY];// Stay
	 std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[DIE];// Dy
	 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[NEXT_STAGE];// to Division Model

	
	std::cout << "[ G0 =   " << Tumour -> at( ith_clone ) -> G0_cells << " ] "
							 << "P(st) = " << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; P(Dy) = " << std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; P(G1) = " << comp_probability_G0  
							 << " ; G0(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; G0(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) 
							 << " ; G0(Trans) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) << std::endl;


	std::cout << "[ G1  =  " << Tumour -> at ( ith_clone ) -> G1_cells << " ] "
							 << "P(st) = "  << std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) ->G1_status ) 
							 << " ; P(Dy) = " << std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status )
							 << " ; P(tra) = "<< comp_probability_G1
							 << " ; G1(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) 
							 << " ; G1(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) 
							 << " ; G1(Tr) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) 
							 << " ; G1(G0) = " << std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) << std::endl;

	std::cout << "[ G2 = "   << Tumour -> at ( ith_clone ) -> G2_cells << " ] "
							 << " P(St) = "  << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) 
							 << " ; P(Dy) = "  << std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status )
							 << " ; P(Tra) = " << comp_probability_G2
							 << " ; G2(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) 
							 << " ; G2(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) 
							 << " ; G2(Tr) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) << std::endl;

	std::cout << "[ S =  " << Tumour -> at ( ith_clone ) -> S_cells << " ] "
						   << " P(St) = "  << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ) 
						   << " ; P(Dy) = " << std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status )
						   << " ; P(Tra) = "<< comp_probability_S
						   << " ; S(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) 
						   << " ; S(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) 
						   << " ; S(Tr) = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) << std::endl;

	std::cout << "[ M = " << Tumour -> at( ith_clone ) -> M_cells << " ] " 
	 					  << " P(St) = "  << std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status ) 
						  << " ; P(Dy) = " << std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status )
						  << " ; P(Tra) = "<< comp_probability_M
						  << " ; M(St) = " << std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) ->M_status ) 
						  << " ; M(Dy) = " << std::get<DYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) 
						  << " ; M(Tr) =  " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) << std::endl;


	std::cout << "[ AC(t) = " << AC_at_t << " ] => "
	 						   << " ; Idling = " <<  Clonal_update_buffer[AV_CELLS_IDLING]
							   << " ; Dying = " << Clonal_update_buffer[AV_CELLS_DYING]
							   << " ; To G0 = " << Clonal_update_buffer[AV_CELLS_to_G0]
							   << " ; To G1 = " << Clonal_update_buffer[AV_CELLS_to_G1] << std::endl;

	

	double p_idle = 0.0;
	double p_go_to_G0 = 0.0;
	Probabilities_of_Cell_Division (ith_clone, p_idle, p_go_to_G0 );
	//Send this as part of the model and then soft  max 
	unsigned int * Division_Model = r.Newborn( std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ), 
											   p_idle, 
											   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR],
											   p_go_to_G0 
											  );
	std::cout << "[ DV = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) << " ] "
			  << "; PDiv(Idle) = "   << p_idle 
			  << " ; MR = "       <<  Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR]
			  << " ; Idling = "   <<  Division_Model[DIV_CELLS_IDLING]
			  << " ; Muts = "     << Division_Model[DIV_CELLS_MUTATING]
			  << " ; To G0 = "    << Division_Model[DIV_CELLS_to_G0]
			  << " ; To G1 = "    << Division_Model[DIV_CELLS_to_G1]
			  << std::endl;

	//std::cout << "[ DV = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status << " ] "
					    //    << "PDiv(Idle): "   << p_idle 
					    //    << " ; PDiv(G0) = " << p_go_to_G0 
					    //    << " ; MR = "       <<  Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR] 
						   // << " ; Idling = "   <<  Division_Model[DIV_CELLS_IDLING]
						   // << " ; Muts = "     << Division_Model[DIV_CELLS_MUTATING]
						   // << " ; To G0 = "    << Division_Model[DIV_CELLS_to_G0]
						   // << " ; To G1 = "    << Division_Model[DIV_CELLS_to_G1] 
						   //<< std::endl;


	// //update mutations
	// unsigned int Mutant_cells = 10; 
	// //Estimate_Mutational_Effects( Mutant_cells, ith_clone -> Clone_Size );
	unsigned int * Mutations = Update_Clonal_Mutational_Burden( ith_clone, Division_Model[DIV_CELLS_MUTATING], years, hours );

	if (!Tumour -> at(ith_clone) -> clone_extinct )
	{
		std::cout << "VALID BUFFER " << std::endl; 
	}

	//
	// UPDATE Values
	//

	std::cout << " G0(t) = " <<  Tumour -> at( ith_clone ) -> G0_cells << " to ";

	Tumour -> at( ith_clone ) -> G0_cells = std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) + 
											std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) + 
											Division_Model[DIV_CELLS_to_G0] + 
											Mutations[MUTANTS_to_G0] + 
											Clonal_update_buffer[AV_CELLS_to_G0];

	std::cout << " G0(t + dt) = " <<  Tumour -> at( ith_clone ) -> G0_cells  << std::endl;
											
	std::cout << " G1(t) = " <<  Tumour -> at( ith_clone ) -> G1_cells << " to "  ;

	Tumour -> at ( ith_clone ) -> G1_cells =  std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) + 
											  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G1_status ) + 
											  Division_Model[DIV_CELLS_to_G1] + 
											  Clonal_update_buffer[AV_CELLS_to_G1] + 
											  Mutations[MUTANTS_to_G1] ;
	
	std::cout << " G1(t + dt) = " <<  Tumour -> at( ith_clone ) -> G1_cells  << std::endl;

	std::cout << " S(t) = " <<  Tumour -> at( ith_clone ) -> S_cells  << " to ";

	Tumour -> at ( ith_clone ) -> S_cells = std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status )  +
											std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) ;

	std::cout << " S(t + dt) = " <<  Tumour -> at( ith_clone ) -> S_cells  << std::endl;

	std::cout << " G2(t) = " <<  Tumour -> at( ith_clone ) -> G2_cells  << " to " ;

	Tumour -> at ( ith_clone ) -> G2_cells = std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) +
											 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) ;

	std::cout << " G2(t + dt) = " <<  Tumour -> at( ith_clone ) -> G2_cells  << std::endl;

	std::cout << " M(t) = " <<  Tumour -> at( ith_clone ) -> M_cells  << " to " ;

	Tumour -> at ( ith_clone ) -> M_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) +
											std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) ;

	std::cout << " M(t + dt) = " <<  Tumour -> at( ith_clone ) -> M_cells  << std::endl;

	std::cout << " AC(t ) = " <<  Tumour -> at( ith_clone ) -> Available_cells  << " to ";

	Tumour -> at( ith_clone ) -> Available_cells = Clonal_update_buffer[AV_CELLS_IDLING] +
												   Division_Model[DIV_CELLS_IDLING] +
												   Mutations[MUTANTS_to_IDLE] ;

	std::cout << " AC(t + dt) = " <<  Tumour -> at( ith_clone ) -> Available_cells  << std::endl;





	

	

//Mutational model

//G0 = G0(S) +A(G0) + Div(G0)


}


void Clonal_Expansion::Non_Mutagenic_Mitosis_Standard(std::unique_ptr<Clone> const & ith_clone)
{
	if( check_G1_phase(ith_clone) )
		ith_clone -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone);
	else if( check_G1_phase(ith_clone) )
		ith_clone -> Remaining_Time_in_S_Phase--;
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone);
	else if( check_G2_phase(ith_clone) )
		ith_clone -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone);
	else if( ith_clone -> In_M_Phase)
	{
		std::cout << "SEND TO DIVISION MODEL " << std::endl;
		Update_Population_After_Division(ith_clone);
	}
}


