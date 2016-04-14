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
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[1];

	Population_Size++;
	mitosis = Clonal_Expansion::Standard; // Requires modification from input parameter
	
	std::cout << "\t\tC A R C I N O G E N E S I S" << std::endl;
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
	Tumour -> back() -> max_PR = Tumour -> back() -> P_Expansion[1];

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



//This can run in parallel 
void Clonal_Expansion::Check_Mitosis_Network_Status(std::unique_ptr<Clone> const & ith_clone)
{
	// All of this runs
	Random r;
	
	// // for G0
	double comp_probability_G0 = 1.0 - (std::get<P_STAYING>( ith_clone -> G0_status ) + std::get<P_DYING>( ith_clone -> G0_status ) ) ;
	double comp_probability_G1 = 1.0 - (std::get<P_STAYING>( ith_clone -> G1_status ) + std::get<P_DYING>( ith_clone -> G1_status ) );
	double comp_probability_G2 = 1.0 - (std::get<P_STAYING>( ith_clone -> G2_status ) + std::get<P_DYING>( ith_clone -> G2_status ) );
	double comp_probability_S  = 1.0 - (std::get<P_STAYING>( ith_clone -> S_status ) + std::get<P_DYING>( ith_clone -> S_status ) );
	double comp_probability_M  = 1.0 - (std::get<P_STAYING>( ith_clone -> M_status ) + std::get<P_DYING>( ith_clone -> M_status ) );
	

	unsigned int * G0_buffer = r.Update_G0_Phase(ith_clone -> G0_cells, 
								  std::get<P_STAYING>( ith_clone -> G0_status ), 
								  std::get<P_DYING>( ith_clone -> G0_status ),
								  comp_probability_G0);

	unsigned int * G1_buffer = r.Update_G1_Phase( ith_clone -> G1_cells,
								   std::get<P_STAYING>(ith_clone -> G1_status),
								   std::get<P_DYING>(ith_clone -> G1_status),
								   comp_probability_G1);

	 unsigned int * G2_buffer = r.Update_G2_Phase( ith_clone -> G2_cells,
	 							   std::get<P_STAYING>(ith_clone -> G2_status),
	 							   std::get<P_DYING>( ith_clone -> G2_status ),
	 							   comp_probability_G2 );

	 unsigned int *  S_buffer = r.Update_S_Phase( ith_clone -> S_cells,
	 							  std::get<P_STAYING>( ith_clone -> S_status),
	 							  std::get<P_DYING>( ith_clone -> S_status ),
	 							  comp_probability_S );

	unsigned int *  M_buffer = r.Update_M_Phase( ith_clone -> M_cells,
	 							 std::get<P_STAYING>( ith_clone -> M_status),
	 							 std::get<P_DYING> ( ith_clone -> M_status ),
	 							 comp_probability_M );

	double to_G0 = 0.001;

	unsigned int * Clonal_update_buffer = r.Clonal_Functions(ith_clone -> Available_cells,
	 									   to_G0,
	 									   ith_clone -> P_Expansion[1], 
	 									   ith_clone -> P_Expansion[0] ) ;

	

	//for G0
	std::get<STAYING_CELLS>( ith_clone -> G0_status ) = G0_buffer[0];// Stay
	std::get<DYING_CELLS>  ( ith_clone -> G0_status ) = G0_buffer[1];// Die 
	std::get<EXITING_CELLS>( ith_clone -> G0_status ) = G0_buffer[2];// to G1
	//for G1
	std::get<STAYING_CELLS>( ith_clone -> G1_status ) = G1_buffer[0];// Stay
	std::get<DYING_CELLS>  ( ith_clone -> G1_status ) = G1_buffer[1];// to Dying
	std::get<EXITING_CELLS>( ith_clone -> G1_status ) = G1_buffer[2];// to G2
	//for G2
	 std::get<STAYING_CELLS>( ith_clone -> G2_status ) = G2_buffer[0];// Stay
	 std::get<DYING_CELLS>  ( ith_clone -> G2_status ) = G2_buffer[1];// to Dying
	 std::get<EXITING_CELLS>( ith_clone -> G2_status ) = G2_buffer[2];
	// //for S
	 std::get<STAYING_CELLS>( ith_clone -> S_status ) = S_buffer[0];// Stay
	 std::get<DYING_CELLS>  ( ith_clone -> S_status ) = S_buffer[1];// to M
	 std::get<EXITING_CELLS>( ith_clone -> S_status ) = S_buffer[2];// to M
	// //for M
	 std::get<STAYING_CELLS>( ith_clone -> M_status ) = M_buffer[0];// Stay
	 std::get<DYING_CELLS>  ( ith_clone -> M_status ) = M_buffer[1];// Dy
	 std::get<EXITING_CELLS>( ith_clone -> M_status ) = M_buffer[2];// to Division Model

	
	std::cout << "Values G0: " << "P(st) " << std::get<P_STAYING>( ith_clone -> G0_status ) 
							   << " P(dy) " << std::get<P_DYING>( ith_clone -> G0_status ) 
							   << " p(trans) " << comp_probability_G0  
							   << " G0(St) " << std::get<STAYING_CELLS>( ith_clone -> G0_status ) 
							   << " G0(Dy) " << std::get<DYING_CELLS>( ith_clone -> G0_status ) 
							   << " G0(Trans) " << std::get<EXITING_CELLS>( ith_clone -> G0_status ) << std::endl;


	std::cout << "Values G1: " << "P(st) "  << std::get<P_STAYING>( ith_clone -> G1_status ) 
							   << " P(dy) " << std::get<P_DYING>( ith_clone -> G1_status )
							   << " P(tra) "<< comp_probability_G1
							   << " G1(St) " << std::get<STAYING_CELLS>( ith_clone -> G1_status ) 
							   << " G1(Dy) " << std::get<DYING_CELLS>( ith_clone -> G1_status ) 
							   << " G1(Tr) " << std::get<EXITING_CELLS>( ith_clone -> G1_status ) << std::endl;

	std::cout << "Values G2: " << "P(st) "  << std::get<P_STAYING>( ith_clone -> G2_status ) 
							   << " P(dy) " << std::get<P_DYING>( ith_clone -> G2_status )
							   << " P(tra) "<< comp_probability_G2
							   << " G2(St) " << std::get<STAYING_CELLS>( ith_clone -> G2_status ) 
							   << " G2(Dy) " << std::get<DYING_CELLS>( ith_clone -> G2_status ) 
							   << " G2(Tr) " << std::get<EXITING_CELLS>( ith_clone -> G2_status ) << std::endl;

	std::cout << "Values S: " << "P(st) "  << std::get<P_STAYING>( ith_clone -> S_status ) 
							   << " P(dy) " << std::get<P_DYING>( ith_clone -> S_status )
							   << " P(tra) "<< comp_probability_S
							   << " S(St) " << std::get<STAYING_CELLS>( ith_clone -> S_status ) 
							   << " S(Dy) " << std::get<DYING_CELLS>( ith_clone -> S_status ) 
							   << " S(Tr) " << std::get<EXITING_CELLS>( ith_clone -> S_status ) << std::endl;

	std::cout << "Values M: " << "P(st) "  << std::get<P_STAYING>( ith_clone -> M_status ) 
							   << " P(dy) " << std::get<P_DYING>( ith_clone -> M_status )
							   << " P(tra) "<< comp_probability_M
							   << " M(St) " << std::get<STAYING_CELLS>( ith_clone -> M_status ) 
							   << " M(Dy) " << std::get<DYING_CELLS>( ith_clone -> M_status ) 
							   << " M(Tr) " << std::get<EXITING_CELLS>( ith_clone -> M_status ) << std::endl;


	std::cout << "Clonal Values: " << "Cells idling " <<  Clonal_update_buffer[0]
								   << " Cells dying " << Clonal_update_buffer[1]
								   << " Cells going to G0 " << Clonal_update_buffer[2]
								   << " Cells entering Mitosis " << Clonal_update_buffer[3] << std::endl;

	double p_idle = 0.9;
	double p_go_to_G0 = 0.01;
	unsigned int * Division_Model = r.Newborn( std::get<EXITING_CELLS>( ith_clone -> M_status ), 
											   p_idle, 
											   ith_clone -> P_Expansion[2],
											    p_go_to_G0 
											  );

	std::cout << "Divided Cells: " << "Cells idling " <<  Division_Model[0]
								   << " Mutant Cells " << Division_Model[1]
								   << " Cells going to G0 " << Division_Model[2]
								   << " Cells entering Mitosis " << Division_Model[3] 
								   << " MR " << ith_clone -> P_Expansion[2] << std::endl;



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


