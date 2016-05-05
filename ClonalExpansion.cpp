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
#include <tuple>

Clonal_Expansion::Clonal_Expansion()
 :
 	Population_Size( 0 ),
 	feedback( 0.0 ),
 	mitosis( Clonal_Expansion::Standard ),
 	Version( Clonal_Expansion::Tester ),
 	SD_Penalty( Clonal_Expansion::standard ),
 	PR_Sampling( Clonal_Expansion::Uniform )
 {       }

Clonal_Expansion::~Clonal_Expansion()
 {       }

 unsigned int Clonal_Expansion::getSize()
 {
 	return Tumour -> size();
 }

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

std::string Clonal_Expansion::getVersion_type(void)
 {
 	std::string types[] = {"Tester", "V1", "V2", "V2R", "V3", "VM"};
 	return types[Version]; 
 }

 std::string Clonal_Expansion::getSD_Penalty(void)
 {
 	std::string types[] = {"standard", "mean", "quartile"};
 	return types[SD_Penalty]; 
 }

 // void Clonal_Expansion::getMitosis_type_str(void)
 // {
 // 	std::string types[] = {"Standard", "Binomial", "Network"};
 // 	return types[mitosis]; 
 // }

  void Clonal_Expansion::setVersion_type(const std::string & str_version)
 {

 	if(str_version == "Tester")
 	{
 		Version = Clonal_Expansion::Tester;
 	}
 	else if(str_version == "V1")
 	{
 		Version = Clonal_Expansion::V1;

 	}
 	else if(str_version == "V2")
 	{
 		Version = Clonal_Expansion::V2;

 	}
 	else if(str_version == "V2R")
 	{
 		Version = Clonal_Expansion::V2R;

 	}
 	else if(str_version == "V3")
 	{
 		Version = Clonal_Expansion::V3;

 	}
 	else if(str_version == "VM")
 	{
 		Version = Clonal_Expansion::VM;

 	}
 	else
 	{
 		Version = Clonal_Expansion::VM;
 	}

 }

void Clonal_Expansion::setSD_Penalty(const std::string & str_Penalty)
{
	if(str_Penalty == "standard")
	{
		SD_Penalty = Clonal_Expansion::standard;
	}
	else if(str_Penalty == "mean")
	{
		SD_Penalty = Clonal_Expansion::mean;
	}
	else if(str_Penalty == "quartile")
	{
		SD_Penalty = Clonal_Expansion::quartile;
	}
	else
	{
		SD_Penalty = Clonal_Expansion::standard;
	}

}


void Clonal_Expansion::Map_Feedback_Penalty()
{
	feedback =  0.0 + (DIFF - 0.0) * (( (double) Population_Size - 0.0) / ((double) PS - 0.0));
	//std::cout << "Feedback " << feedback << std::endl;
}

void Clonal_Expansion::Select_Size_Dependant_Penalty(void)
{
	switch(SD_Penalty)
	{
		case 0:
			Map_Feedback_Penalty();
			break;
		case 1:
			std::cout << " Mean Penalty " << std::endl;
			break;
		case 2:
			std::cout << " Quartile " << std::endl;
			break;
		default:
			std::cout << "STD default " << std::endl;
			break;

	}//switch
}


 void Clonal_Expansion::Select_Carcinogenesis(void)
 {
 	switch(Version)
 	{
 		case 0:
            std::cout << "0 Version: " << Version << std::endl;
            break;
        case 1:
            std::cout << "1 Version: " << Version << std::endl;
            carcinogenesis_V1();
            break;
        case 2:
            std::cout << "2 Version: " << Version << std::endl;
            break;
        case 3:
            std::cout << "3 Version: " << Version << std::endl;
            break;
        case 4:
            std::cout << "4 Version: " << Version << std::endl;
            break;
        case 5:
            std::cout << "5 Version: " << Version << std::endl;
            break;
        default:
        	 std::cout << "Def Version: " << Version << std::endl;
        	 break;
 	}//switch	
 
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
	Tumour -> back() -> clone_extinct = false;

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

/*
	Initialisation of Carcinogenesis of Version 1
*/
void Clonal_Expansion::carcinogenesis_V1(void)
{
	//std::cout << "Carcinogenesis Version 1 " << std::endl;
	// 1) Add a Clone
	Random r;
	// Adding a new clone
	Tumour -> push_back ( get_Clone_DS () );  
	Tumour -> back() -> Initiall_Expasion_Period = true;
	Tumour -> back() -> clone_extinct = false;
	Tumour -> back() -> Clone_Size = 1;
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

	Population_Size = Tumour -> back() -> Clone_Size;
}



void Clonal_Expansion::Generate_Clone_Generation_ID( std::string & newGeneration_ID,
								  const int & Parent_Generation_ID_Counter, 
								  const std::string & Parent_ID, 
								  const unsigned int & years, 
								  const unsigned int & hours)
{

	std::string heredity_pattern = Parent_ID.substr( 0, Parent_ID.find("-") );
	std::string Clone_ID ("");
		//cout << "I GET " << heredity_pattern << endl;

		if(heredity_pattern.compare("P") == 0)
		{
			//cout << "I GET 1 " << heredity_pattern << endl;
			std::string Clone_ID_tag ("PC");
			Clone_ID_tag.append(std::to_string(Parent_Generation_ID_Counter));Clone_ID_tag.append("-");
			Clone_ID_tag.append(std::to_string(years));Clone_ID_tag.append(":");Clone_ID_tag.append(std::to_string(hours));
			Clone_ID = Clone_ID_tag;
		}
		else
		{
			
			std::string Clone_ID_tag ("");
			Clone_ID_tag.append(heredity_pattern); Clone_ID_tag.append(","); Clone_ID_tag.append(std::to_string(Parent_Generation_ID_Counter));
			Clone_ID_tag.append("-");
			Clone_ID_tag.append(std::to_string(years));Clone_ID_tag.append(":");Clone_ID_tag.append(std::to_string(hours));
			Clone_ID = Clone_ID_tag;

		}

	//std::cout << "ID " << Parent_ID << std::endl;
//	std::cout << "HP " << heredity_pattern << std::endl;

	newGeneration_ID = Clone_ID;

}




void Clonal_Expansion::carcinogenesis_from_driver(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours, Random & r)
{
	

	double mr = Tumour -> at( ith_clone ) -> Mutation_Rate;
	unsigned int AC = Tumour -> at( ith_clone ) -> Driver_10_fold_accumulation;
	//unsigned int Accumulated_Drivers = Tumour -> at( ith_clone ) -> Driver_10_fold_accumulation;
	unsigned long long int NOM = Tumour -> at( ith_clone ) -> Number_of_Mutations;
	double pr = Tumour -> at( ith_clone ) -> P_Expansion[1];

	double p_st_g0 =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status );
	double p_dy_g0 =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status );

	double p_st_g1 =  std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status );
	double p_dy_g1 =  std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status );

	double p_st_g2 =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status );
	double p_dy_g2 =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status );

	double p_st_s =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status );
	double p_dy_s =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status );

	double p_st_m =  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status );
	double p_dy_m =  std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status );

	std::string cloneName = "";	
	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) -> Generation_ID_Counter;
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
	Tumour -> back() -> P_Expansion[2] = Tumour -> back() -> Mutation_Rate;
	double PR = r.Update_Proliferation_Rate(pr);
	Tumour -> back() -> Driver_10_fold_accumulation = AC + 5;
	Tumour -> back() -> P_Expansion[1] =  PR;
	Tumour -> back() -> init_PR = PR;
	Tumour -> back() -> max_PR = PR;
	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;

	Tumour -> back() -> G0_cells = 0;
	std::get<P_STAYING>( Tumour -> back() -> G0_status ) = p_st_g0;
	std::get<P_DYING>( Tumour -> back() -> G0_status ) = p_dy_g0;

	Tumour -> back() -> G1_cells = 1;
	std::get<P_STAYING_G1>( Tumour -> back() -> G1_status ) = p_st_g1;
	std::get<P_DYING_G1>( Tumour -> back() -> G1_status ) = p_dy_g1;

	 Tumour -> back() -> G2_cells = 0;
	 std::get<P_STAYING>( Tumour -> back() -> G2_status ) = p_st_g2;
	 std::get<P_DYING>( Tumour -> back() -> G2_status ) = p_dy_g2;


	 Tumour -> back() -> S_cells = 0;
	 std::get<P_STAYING>( Tumour -> back() -> S_status ) = p_st_s;
	 std::get<P_DYING>( Tumour -> back() -> S_status ) = p_dy_s;

	 Tumour -> back() -> M_cells = 0;
	 std::get<P_STAYING>( Tumour -> back() -> M_status ) = p_st_m;
	 std::get<P_DYING>( Tumour -> back() -> M_status ) = p_dy_m;
	

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
bool Clonal_Expansion::check_G1_phase( const unsigned int & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G1_Phase && Tumour -> at (ith_clone) -> Remaining_Time_in_G1_Phase !=0);
}


/**
	Check if the ith_clone is exiting the G1 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_G1_phase( const unsigned int  & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G1_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase == 0);
}

/**
	Check if the ith_clone is in S mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_S_phase( const unsigned int & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_S_Phase && Tumour -> at (ith_clone) -> Remaining_Time_in_S_Phase !=0);
}

/**
	Check if the ith_clone is exiting the G1 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_S_phase( const unsigned int  & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_S_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase == 0);
}

/**
	Check if the ith_clone is in G2 mitotic 
	phase
*******************************************/
bool Clonal_Expansion::check_G2_phase( const unsigned int  & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G2_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase !=0);
}

/**
	Check if the ith_clone is exiting the G2 mitotic 
	phase
******************************************************/
bool Clonal_Expansion::exiting_G2_phase( const unsigned int & ith_clone)
{
	return (Tumour -> at(ith_clone) -> In_G2_Phase && Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase == 0);
}

/**
	ith_Clone from G1 -> S phase
******************************************************/
void Clonal_Expansion::Transition_From_G1_S(const unsigned int & ith_clone, Random & r)
{
	Tumour -> at(ith_clone) -> In_G1_Phase  = false; 
	Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase = r.G1();
	Tumour -> at(ith_clone) -> In_S_Phase = true;
}

/**
	ith_Clone from S -> G2 phase
******************************************************/
void Clonal_Expansion::Transition_From_S_G2( const unsigned int  & ith_clone, Random & r)
{
	Tumour -> at (ith_clone) -> In_S_Phase  = false; 
	Tumour -> at (ith_clone) -> Remaining_Time_in_S_Phase = r.S();
	Tumour -> at (ith_clone) -> In_G2_Phase = true;
}

/**
	ith_Clone from G2 -> M phase
******************************************************/
void Clonal_Expansion::Transition_From_G2_M(const unsigned int  & ith_clone, Random & r)
{
	Tumour -> at(ith_clone) -> In_G2_Phase  = false; 
	Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase = r.G2();
	Tumour -> at(ith_clone) -> In_M_Phase = true;
}



// TODO Asymetric division
void Clonal_Expansion::Update_Population_After_Division_V1(const unsigned int  & ith_clone)
{

	unsigned long long int cells_after_division = Tumour -> at (ith_clone) -> Clone_Size * 2;
	Population_Size =  (Population_Size - Tumour -> at(ith_clone) -> Clone_Size ) + cells_after_division;
	Tumour -> at(ith_clone) -> Clone_Size = cells_after_division;
	//std::cout << "Population_Size: " << Population_Size << std::endl;
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


void Clonal_Expansion::Check_for_Clonal_Extintion_at_min_size(const unsigned int & ith_clone)
{

	if( Tumour -> at(ith_clone) -> Clone_Size == 1 )
	{
		Tumour -> at(ith_clone) -> Clone_Size = 0 ;
		Tumour -> at(ith_clone) -> clone_extinct = true; 
	}	

}



void Clonal_Expansion::Update_Clonal_Mutational_Burden(const unsigned int & ith_clone, const unsigned int & Mutant_cells, const unsigned int & years, const unsigned int & hours, std::vector<unsigned int> & Mutations, Random & r )
{
	
	unsigned int * mutational_effect_buffer;
	unsigned int counter = 0;
	unsigned int death = 0;
	double effect_gain = static_cast<double>( Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] );
	//bool mut_effect_flag = false;
	

	unsigned int j = 0;


	while( (j < Mutant_cells  ) && !(Tumour -> at(ith_clone) -> clone_extinct) )
	{
		unsigned int i = 0;
		while( i < add_multiple_mutations() )
		{
			// check for mutations
			mutational_effect_buffer = r.Mutational_Effects();
			// buffer is size 3 --> need to check kills drivers and effects
			std::cout << "P(Ki) = " << KILLER_PROBABILITY << "; P(Dr) = " << DRIVER_PROBABILITY << " ; P(P) =  " << 1.0 -(KILLER_PROBABILITY+DRIVER_PROBABILITY) << " [ " << i << " cells]:  "<< " KILL = " <<  mutational_effect_buffer[0] << " DR = " <<  mutational_effect_buffer[1] << " PASS= " << mutational_effect_buffer[2] << std::endl;
			if ( mutational_effect_buffer[PASSENGER_CELLS] > 0 )// Passenger
			{
			 	std::cout << " MUTATIONAL EFFECT " << std::endl; 
			 	unsigned int CS = static_cast<unsigned int>(Tumour -> at(ith_clone) -> Clone_Size);
			 	// std::cout << "Before Expansion " 
			 	// 		  << Tumour -> at(ith_clone) -> P_Expansion[1] << " " 
			 	// 		  << Tumour -> at(ith_clone) -> Clone_Size << std::endl;

			 	r.Laplace( Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR], CS );

			 	// std::cout << "After Expansion " 
			 	// 		  << Tumour -> at(ith_clone) -> P_Expansion[1] << " " 
			 	// 		  << Tumour -> at(ith_clone) -> Clone_Size << std::endl;

			 	counter++;
			 	//getchar();

			 	//mut_effect_flag = true;

			}
			else if(mutational_effect_buffer[DRIVER_CELLS] > 0)// Driver
			{
				std::cout << " DRIVER " << std::endl;
				carcinogenesis_from_driver( ith_clone,  years,  hours, r );
				
				//getchar();
			}
			else if(mutational_effect_buffer[0] > 0)
			{
				std::cout << " DEATH " << std::endl;
				Check_for_Clonal_Extintion_at_min_size( ith_clone );
				death++;
				//getchar();
				// if(Tumour -> at(ith_clone) -> clone_extinct)
				// 	break;
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


	//unsigned int MCC = Mutant_cells ;
	
	//Mutations[0] = 0; Mutations[1] = 0; Mutations[2] = 0;

	//if (!Tumour -> at(ith_clone) -> clone_extinct )
	//{

		//Mutations = 
		r.Mutational_Proportions( effect_gain, Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR] , Mutant_cells, Mutations );
	
	//}
	

	std::cout << "[ Muts  =  " << Mutant_cells<< " ] "
							 << "P(G0) = "  << effect_gain
							 << " ; P(G1) = " << Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR]
							 << " ; P(I) = " << 1.0 -(effect_gain + Tumour -> at(ith_clone) -> P_Expansion[P_Expansion_PR])
							 << " ; to G0 = " << Mutations[0]
							 << " ; to G1 = " << Mutations[1]
							 << " ; IDLE  = " << Mutations[2]
							 << " ; Deaths = " << death << std::endl;

	// std::cout << "Values Muts: " << "Muts(G0) " << Mutations[MUTANTS_to_G0]
	// 			  << "  Muts(G1) " << Mutations[MUTANTS_to_G1]
	// 			  << " Muts(IDLE) " << Mutations[MUTANTS_to_IDLE]
	// 			  << " Muts (Death) " << death 
	// 			  << " Mutant cells "<< Mutant_cells 
	// 			  << " effect gain " <<  effect_gain << std::endl;



	

}

void Clonal_Expansion::Probabilities_of_Cell_Division (const unsigned int & ith_clone, double & p_idle, double & p_go_to_G0 )
{
	p_go_to_G0 = fabs( Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR] - Tumour -> at( ith_clone ) -> max_PR);
	p_idle = 0.62 * (1.0 - ( Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR] + p_go_to_G0 ) );
	

}


void Clonal_Expansion::Complementary_Probability( double & comp_probability_G0, 
												  double & comp_probability_G1, 
												  double & comp_probability_G2,
												  double & comp_probability_S,
												  double & comp_probability_M,
												  const unsigned int & ith_clone )
{

	comp_probability_G0 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) ); 
	comp_probability_G1 = 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) ) - 0.1;
	comp_probability_G2 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status ) );
	comp_probability_S  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status ) );
	comp_probability_M  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status ) );

}


// Update G0
void Clonal_Expansion::Update_G0_Phase( const double & comp_probability_G0, const unsigned int & ith_clone, Random & r)
{
	unsigned int * G0_buffer = r.Update_G0_Phase( Tumour -> at( ith_clone ) -> G0_cells, 
	 							std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ), 
	 							std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G0_status ),
	 							comp_probability_G0);

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[DIE];// Die 
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[NEXT_STAGE];// to G1
}

//Update G1
void Clonal_Expansion::Update_G1_Phase(const double & comp_probability_G1, const unsigned int & ith_clone, Random & r)
{
	
	unsigned int * G1_buffer = r.Update_G1_Phase( Tumour -> at( ith_clone ) -> G1_cells,
	 							std::get<P_STAYING_G1> ( Tumour -> at( ith_clone ) -> G1_status ),
	 							std::get<P_DYING_G1>   ( Tumour -> at( ith_clone ) -> G1_status),
	 							comp_probability_G1);

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[DIE];// to Dying
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[NEXT_STAGE];// to G2
	std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[REVERT_TO_G0];// to G0
}

// Update G2
void Clonal_Expansion::Update_G2_Phase( const double & comp_probability_G2, const unsigned int & ith_clone, Random & r)
{
	unsigned int * G2_buffer = r.Update_G2_Phase( Tumour -> at( ith_clone ) -> G2_cells,
	  							   std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ),
	  							   std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G2_status ),
	  							   comp_probability_G2 );

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[DIE];// to Dying
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[NEXT_STAGE];
}

//Update S
void Clonal_Expansion::Update_S_Phase( const double & comp_probability_S, const unsigned int & ith_clone, Random & r)
{
	unsigned int *  S_buffer = r.Update_S_Phase( Tumour -> at( ith_clone ) -> S_cells,
	  							std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ),
	  							std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> S_status ),
	  							comp_probability_S );

	std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[STAY];// Stay
	std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[DIE];// die
	std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[NEXT_STAGE];// to M

}

//Update M

void Clonal_Expansion::Update_M_Phase( const double & comp_probability_M, const unsigned int & ith_clone, Random & r)
{
	unsigned int *  M_buffer = r.Update_M_Phase( Tumour -> at( ith_clone ) -> M_cells,
	  							 std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status),
	  							 std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> M_status ),
	  							 comp_probability_M );

	 std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[STAY];// Stay
	 std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[DIE];// Dy
	 std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[NEXT_STAGE];// to Division Model
}



//This can run in parallel 
// Mitosis
// https://www.youtube.com/watch?v=koudmJdil60
void Clonal_Expansion::Check_Mitosis_Network_Status(const unsigned int & ith_clone, const unsigned int & years, const unsigned int & hours, Random & r )
{
	// All of this runs
	//Random r;
	


	double comp_probability_G0 = 0;
	double comp_probability_G1 = 0;
	double comp_probability_G2 = 0;
	double comp_probability_S = 0;
	double comp_probability_M = 0;
	double p_idle = 0.0;
	double p_go_to_G0 = 0.0;

	std::vector<unsigned int> Division_Model;

	//std::vector<unsigned int> Division_Model = 

	//unsigned int * Mutations;

	std::vector<unsigned int> Mutations;

	Complementary_Probability( comp_probability_G0, 
							   comp_probability_G1, 
							   comp_probability_G2,
							   comp_probability_S,
							   comp_probability_M,
							   ith_clone );

	// // // for G0
	// // double comp_probability_G0 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G0_status ) ); 
	// // double comp_probability_G1 = 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) ) - 0.1;
	// // double comp_probability_G2 = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ) + std::get<P_DYING>( Tumour -> at( ith_clone ) -> G2_status ) );
	// // double comp_probability_S  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> S_status ) );
	// // double comp_probability_M  = 1.0 - ( std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status )  + std::get<P_DYING>( Tumour -> at( ith_clone ) -> M_status ) );
	
	// std::cout << "STATS G0 >  " <<  Tumour -> at( ith_clone ) -> G0_cells << std::endl;

	 Update_G0_Phase( comp_probability_G0, ith_clone, r );
	 Update_G1_Phase( comp_probability_G1, ith_clone, r );
	 Update_G2_Phase( comp_probability_G2, ith_clone, r );
	 Update_S_Phase(  comp_probability_S,  ith_clone, r );
	 Update_M_Phase(  comp_probability_M, ith_clone, r );

	 unsigned int AC_at_t = Tumour -> at( ith_clone ) -> Available_cells;
	 unsigned int * Clonal_update_buffer = r.Clonal_Functions( Tumour -> at( ith_clone ) -> Available_cells,
	   									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR],
	   									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR] - feedback, 
	   									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR] ) ;
	
	 Probabilities_of_Cell_Division ( ith_clone, p_idle, p_go_to_G0 );

	 r.Newborn( std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ), 
	 										   p_idle, 
	 										   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR],
	 										   p_go_to_G0,
	 										   Division_Model
	 										  );

	 

	


	//  // unsigned int * G0_buffer = r.Update_G0_Phase( Tumour -> at( ith_clone ) -> G0_cells, 
	//  // 							  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G0_status ), 
	//  // 							  std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G0_status ),
	//  // 							  comp_probability_G0);

	//  // unsigned int * G1_buffer = r.Update_G1_Phase( Tumour -> at( ith_clone ) -> G1_cells,
	//  // 							   std::get<P_STAYING_G1> ( Tumour -> at( ith_clone ) -> G1_status ),
	//  // 							   std::get<P_DYING_G1>   ( Tumour -> at( ith_clone ) -> G1_status),
	//  // 							   comp_probability_G1);

	//   // unsigned int * G2_buffer = r.Update_G2_Phase( Tumour -> at( ith_clone ) -> G2_cells,
	//   // 							   std::get<P_STAYING>( Tumour -> at( ith_clone ) -> G2_status ),
	//   // 							   std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> G2_status ),
	//   // 							   comp_probability_G2 );

	//   // unsigned int *  S_buffer = r.Update_S_Phase( Tumour -> at( ith_clone ) -> S_cells,
	//   // 							  std::get<P_STAYING>( Tumour -> at( ith_clone ) -> S_status ),
	//   // 							  std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> S_status ),
	//   // 							  comp_probability_S );

	//  // unsigned int *  M_buffer = r.Update_M_Phase( Tumour -> at( ith_clone ) -> M_cells,
	//  //  							 std::get<P_STAYING>( Tumour -> at( ith_clone ) -> M_status),
	//  //  							 std::get<P_DYING>  ( Tumour -> at( ith_clone ) -> M_status ),
	//  //  							 comp_probability_M );


	//  // TODO Change this
	//  //double to_G0 = 0.001; // probability of expoansion 
	//  // unsigned int AC_at_t = Tumour -> at( ith_clone ) -> Available_cells;
	//  // unsigned int * Clonal_update_buffer = r.Clonal_Functions( Tumour -> at( ith_clone ) -> Available_cells,
	//  //  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR],
	//  //  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR], 
	//  //  									   Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR] ) ;

	

	// // //for G0
	// // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[STAY];// Stay
	// // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[DIE];// Die 
	// // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) = G0_buffer[NEXT_STAGE];// to G1

	// std::cout << "ASSIGN G0 >  " <<  Tumour -> at( ith_clone ) -> G0_cells << std::endl;

	// //for G1
	// // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[STAY];// Stay
	// // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[DIE];// to Dying
	// // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[NEXT_STAGE];// to G2
	// // std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) = G1_buffer[REVERT_TO_G0];// to G0
	// //for G2
	//  // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[STAY];// Stay
	//  // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[DIE];// to Dying
	//  // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) = G2_buffer[NEXT_STAGE];
	// // //for S
	//  // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[STAY];// Stay
	//  // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[DIE];// die
	//  // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> S_status ) = S_buffer[NEXT_STAGE];// to M
	// // //for M
	//  // std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[STAY];// Stay
	//  // std::get<DYING_CELLS>  ( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[DIE];// Dy
	//  // std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) = M_buffer[NEXT_STAGE];// to Division Model

	
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
							 << " ; P(toG0) = " << 1.0 - ( std::get<P_STAYING_G1>( Tumour -> at( ith_clone ) ->G1_status ) + std::get<P_DYING_G1>( Tumour -> at( ith_clone ) -> G1_status ) + comp_probability_G1)
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


	std::cout << "[ DV = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status ) << " ] "
			  << "; PDiv(Idle) = "   << p_idle 
			  << "; P(G0) = " << p_go_to_G0
			  << " ; P(MR) = "       <<  Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR]
			  << " ; P(G1) = " 		<< 1.0 - ( p_idle + p_go_to_G0 + Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR])
			  << " : sum = " << ( 1.0 - ( p_idle + p_go_to_G0 + Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR])) + ( p_idle + p_go_to_G0 + Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR])
			  << " ; Idling = "   <<  Division_Model[DIV_CELLS_IDLING]
			  << " ; Muts = "     << Division_Model[DIV_CELLS_MUTATING]
			  << " ; To G0 = "    << Division_Model[DIV_CELLS_to_G0]
			  << " ; To G1 = "    << Division_Model[DIV_CELLS_to_G1]
			  << std::endl;


	std::cout << "[ AC(t) = " << AC_at_t << " ] => "
							   << " ; PR = " << Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_PR]
							   << " ; penalty = " << feedback
							   << " ; DR = " << Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_DR]
							   << " ; QR = " << Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_QR]
	 						   << " ; Idling = " <<  Clonal_update_buffer[AV_CELLS_IDLING]
							   << " ; Dying = " << Clonal_update_buffer[AV_CELLS_DYING]
							   << " ; To G0 = " << Clonal_update_buffer[AV_CELLS_to_G0]
							   << " ; To G1 = " << Clonal_update_buffer[AV_CELLS_to_G1] << std::endl;

	

	
	
	// //Send this as part of the model and then soft  max 
	

	

	//std::cout << "[ DV = " << std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> M_status << " ] "
					    //    << "PDiv(Idle): "   << p_idle 
					    //    << " ; PDiv(G0) = " << p_go_to_G0 
					    //    << " ; MR = "       <<  Tumour -> at( ith_clone ) -> P_Expansion[P_Expansion_MR] 
						   // << " ; Idling = "   <<  Division_Model[DIV_CELLS_IDLING]
						   // << " ; Muts = "     << Division_Model[DIV_CELLS_MUTATING]
						   // << " ; To G0 = "    << Division_Model[DIV_CELLS_to_G0]
						   // << " ; To G1 = "    << Division_Model[DIV_CELLS_to_G1] 
						   //<< std::endl;


	// // //update mutations
	// // unsigned int Mutant_cells = 10; 
	// // //Estimate_Mutational_Effects( Mutant_cells, ith_clone -> Clone_Size );



	
	//if( Division_Model[DIV_CELLS_MUTATING] >  0)
	//{
	 	Update_Clonal_Mutational_Burden( ith_clone, Division_Model[DIV_CELLS_MUTATING], years, hours, Mutations, r );
	//}

	// //
	// // UPDATE Values
	// //

	std::cout << " G0(t) = " <<  Tumour -> at( ith_clone ) -> G0_cells << " to ";

	Tumour -> at( ith_clone ) -> G0_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) + 
											std::get<FROM_G1_TO_G0>( Tumour -> at( ith_clone ) -> G1_status ) + 
											Division_Model[DIV_CELLS_to_G0] + 
											Mutations[MUTANTS_to_G0] + 
											Clonal_update_buffer[AV_CELLS_to_G0];

	std::cout << " G0(t + dt) = " <<  Tumour -> at( ith_clone ) -> G0_cells  << std::endl;
											
	std::cout << " G1(t) = " <<  Tumour -> at( ith_clone ) -> G1_cells << " to "  ;

	Tumour -> at ( ith_clone ) -> G1_cells =  std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G0_status ) + 
											  std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) + 
											  Division_Model[DIV_CELLS_to_G1] + 
											  Clonal_update_buffer[AV_CELLS_to_G1] + 
											  Mutations[MUTANTS_to_G1] ;
	
	std::cout << " G1(t + dt) = " <<  Tumour -> at( ith_clone ) -> G1_cells  << std::endl;

	std::cout << " S(t) = " <<  Tumour -> at( ith_clone ) -> S_cells  << " to ";

	Tumour -> at ( ith_clone ) -> S_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> S_status )  +
											std::get<EXITING_CELLS>( Tumour -> at( ith_clone ) -> G1_status ) ;

	std::cout << " S(t + dt) = " <<  Tumour -> at( ith_clone ) -> S_cells  << std::endl;

	std::cout << " G2(t) = " <<  Tumour -> at( ith_clone ) -> G2_cells  << " to " ;

	Tumour -> at ( ith_clone ) -> G2_cells = std::get<STAYING_CELLS>( Tumour -> at( ith_clone ) -> G2_status ) +
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

	Tumour -> at( ith_clone ) -> Clone_Size = static_cast<unsigned long long int>(Tumour -> at(ith_clone) -> Available_cells + 
	 																			  Tumour -> at(ith_clone) -> G0_cells        + 
	 																			  Tumour -> at(ith_clone) -> G1_cells        + 
	 																			  Tumour -> at(ith_clone) -> G2_cells        + 
	 																			  Tumour -> at(ith_clone) -> S_cells         + 
	 																			  Tumour -> at(ith_clone) -> M_cells);

	std::cout << " Clone_S( " <<  ith_clone << " ) = "<< Tumour -> at( ith_clone ) -> Clone_Size  << std::endl;
	
}




void Clonal_Expansion::map_Feedback(  )
{

	feedback =  0.0 + (DIFF - 0.0) * (( (double) Population_Size - 0.0) / ((double) PS - 0.0));	
	//std::cout << "Population Feedback " << feedback << " = " << DIFF << " * " << Population_Size << " / " << PS << std::endl;
}


void Clonal_Expansion::Update_Tumour_Size()
{
	unsigned int ith_clone = 0;
	unsigned long long int population = 0;
	for (ith_clone = 0; ith_clone < Tumour -> size() ; ith_clone++) 
	{
		population += Tumour -> at(ith_clone) -> Clone_Size;
	}
	Population_Size = population;

	//std::cout << "Tumour SIZE = " <<  Population_Size << std::endl;
}


// void Clonal_Expansion::Non_Mutagenic_Mitosis_Standard(std::unique_ptr<Clone> & ith_clone)
// {

// 	if( check_G1_phase(ith_clone) )
// 		ith_clone -> Remaining_Time_in_G1_Phase--;
// 	else if( exiting_G1_phase(ith_clone) )
// 		Transition_From_G1_S(ith_clone);
// 	else if( check_G1_phase(ith_clone) )
// 		ith_clone -> Remaining_Time_in_S_Phase--;
// 	else if( exiting_S_phase(ith_clone) )
// 		Transition_From_S_G2(ith_clone);
// 	else if( check_G2_phase(ith_clone) )
// 		ith_clone -> Remaining_Time_in_G2_Phase--;
// 	else if( exiting_G2_phase(ith_clone) )
// 		Transition_From_G2_M(ith_clone);
// 	else if( ith_clone -> In_M_Phase)
// 	{
// 		std::cout << "SEND TO DIVISION MODEL " << std::endl;
// 		Update_Population_After_Division(ith_clone);
// 	}
// }

void Clonal_Expansion::Reset_to_G1(const unsigned int & ith_clone)
{
	Tumour -> at (ith_clone) -> In_M_Phase = false;
	Tumour -> at (ith_clone) -> In_G1_Phase = true;
}

void Clonal_Expansion::Valid_Size_Heterogeneity_V1(const unsigned int & ith_clone)
{
	
	if( Tumour -> at(ith_clone) -> Clone_Size >= Tumour -> at(ith_clone) -> Number_of_Memebers_to_Start_Heterogeneity)
	{
		Tumour -> at(ith_clone) ->  Initiall_Expasion_Period = false;		
	}
			
	// Tumour -> at (ith_clone) -> In_M_Phase = false;
	// Tumour -> at (ith_clone) -> In_G1_Phase = true;
}


// V1 Final version

void Clonal_Expansion::Non_Mutagenic_Mitosis_Standard_V1(const unsigned int & ith_clone, Random & r)
{
	if( check_G1_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone, r);
	else if( check_S_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase--;	
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone, r);
	else if( check_G2_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone, r);
	else if(Tumour -> at (ith_clone) -> In_M_Phase)
	{
		//std::cout << "SEND TO DIVISION MODEL " << std::endl;
		Update_Population_After_Division_V1(ith_clone);
		Valid_Size_Heterogeneity_V1(ith_clone);
	}
	else
	{
		std::cout << "ERROR " << std::endl;
	}



}

void Clonal_Expansion::Checking_for_Mutants(std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn, Random & r, const unsigned int & ith_clone)
{
	unsigned int mutant_cells = 0;

	if(std::get<1>(Dying_and_Newborn) > 0)
	{	
		std::get<2>(Dying_and_Newborn) = r.Binomial_Mutants(
															std::get<1>(Dying_and_Newborn),
															Tumour -> at(ith_clone) -> Mutation_Rate
														   );
		// Announce that there are mutants
		mutant_cells = std::get<2>(Dying_and_Newborn); // Mutant cells

		//std::cout << "Mut Cells " << mutant_cells << std::endl;
		std::get<3>(Dying_and_Newborn) =  true;// flag that indicates that there are mutant cells

	
		if(mutant_cells > std::get<1>(Dying_and_Newborn))
		{
			std::cout << "More mutant cells than newborn? " << " MC: " << mutant_cells << " NB: " << std::get<1>(Dying_and_Newborn) << std::endl;
		}
		else
		{
			std::get<1>(Dying_and_Newborn) -= mutant_cells;
		}
	}

}

void Clonal_Expansion::Update_Newborn_Parameters_V1(const std::vector<unsigned int> & NewBorn_Cells, std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn)
{

	std::get<0>(Dying_and_Newborn) = NewBorn_Cells[0];// For Dying Clones 		
	std::get<1>(Dying_and_Newborn) = NewBorn_Cells[1]; // For Newborn Clones 

}

void Clonal_Expansion::Adjust_Feedback_PR(double & P_NB)
{
	if(feedback >= P_NB)
		P_NB = 0;
	else
		P_NB -= feedback;
}


void Clonal_Expansion::Check_Clonal_Extinction_BCE_V1(const unsigned int & ith_clone, std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn)
{
	int o = (int) ( Tumour -> at(ith_clone) -> Clone_Size + std::get<1>(Dying_and_Newborn) ) - (int) (std::get<0>(Dying_and_Newborn) + std::get<2>(Dying_and_Newborn));
	if(o <= 0)
	{
		Tumour -> at(ith_clone) -> clone_extinct = true;
		std::get<4>(Dying_and_Newborn) =  false;
		Population_Size -=  Tumour  -> at(ith_clone) -> Clone_Size ;
		Tumour -> at(ith_clone) -> Clone_Size = 0;
	}

}

void Clonal_Expansion::Print_Debug_Carcinogeneisis_From_Driver(const unsigned int & ith_clone)
{
	std::cout << " TUMOUR STATUS " << std::endl;
	printParameters();

	std::cout << "ID: " << Tumour -> at(ith_clone) -> Generation_ID << " to " <<Tumour -> back() -> Generation_ID << std::endl;
	std::cout << "CS: " << Tumour -> at(ith_clone) -> Clone_Size << " to " <<Tumour -> back() -> Clone_Size << std::endl;
	std::cout << "NOM: " << Tumour -> at(ith_clone) -> Number_of_Mutations << " to " <<Tumour -> back() -> Number_of_Mutations << std::endl;
	std::cout << "MR: " <<  Tumour -> at(ith_clone) -> Mutation_Rate << " to " <<Tumour -> back() -> Mutation_Rate << std::endl;
	std::cout << "PR: " << Tumour -> at(ith_clone) -> P_Expansion[1] << " to " <<Tumour -> back() -> P_Expansion[1] << std::endl;

	//getchar();


}


void Clonal_Expansion::Carcinogenesis_From_Driver_V1(const unsigned int & ith_clone, const unsigned int & hours, const unsigned int & years, Random & r)
{
	double mr = Tumour -> at(ith_clone) -> Mutation_Rate;
	unsigned long long int NOM = Tumour -> at (ith_clone) -> Number_of_Mutations;
	double pr = Tumour -> back() -> P_Expansion[1];
	
	std::string cloneName = "";	
	int Parent_Generation_ID_Counter = (int) Tumour -> at( ith_clone ) -> Generation_ID_Counter;
	Generate_Clone_Generation_ID(cloneName, 
								 Parent_Generation_ID_Counter, 
								 Tumour -> at( ith_clone ) -> Generation_ID,
								 years, 
								 hours);
	//std::cout << "NID: " << cloneName << std::endl;

	Tumour -> push_back( get_Clone_DS() );
	Tumour -> back() -> Generation_ID = cloneName;
	Tumour -> back() -> Initiall_Expasion_Period = false;//true
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
	Tumour -> back() -> P_Expansion[1] =  r.Update_Proliferation_Rate(pr);
	Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
	Population_Size++;
	
	//Print_Debug_Carcinogeneisis_From_Driver( ith_clone );

}

void Clonal_Expansion::Mutant_Effects_V1(const unsigned int & ith_clone, const std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn,  const unsigned int & hours, const unsigned int & years, Random & r)
{
	//std::cout << "Driver Mutant " << std::endl;
	unsigned int i =0;

	for( i = 0; i < std::get<2>(Dying_and_Newborn); i++)
		Carcinogenesis_From_Driver_V1( ith_clone, hours, years, r);
}

void Clonal_Expansion::Apply_Penalties_to_Mutants_V1(const unsigned int & ith_clone, const std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn, const unsigned int & hours, const unsigned int & years, Random & r)
{
	
	if( std::get<3>(Dying_and_Newborn) && std::get<2>(Dying_and_Newborn) > 0 )
	{
		Mutant_Effects_V1(ith_clone, Dying_and_Newborn, hours, years, r);
	}

	unsigned int New_poulation_Size = ( Tumour -> at(ith_clone) -> Clone_Size  - std::get<0>(Dying_and_Newborn) ) 
										+ 	std::get<1>(Dying_and_Newborn) ;

		// Update population size by
		// Pop size =Pop_Size + (CS[t] - CS[t-1]) 
		Population_Size +=  ( static_cast<unsigned long long int>(New_poulation_Size) - Tumour  -> at(ith_clone) -> Clone_Size) ;
		Tumour -> at(ith_clone) -> Clone_Size = static_cast<unsigned long long int>(New_poulation_Size) ;
}

void Clonal_Expansion::Update_Mutant_Mutational_Effects(const unsigned int & ith_clone, Random & r, const unsigned int & hours, const unsigned int & years, const std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn)
{

	if(std::get<4>(Dying_and_Newborn))
	{
		//compute_Mutations(CE, Generation_ID, Dying_and_Newborn, years, hours );
		
		Apply_Penalties_to_Mutants_V1( ith_clone, Dying_and_Newborn, hours,  years,  r);
		
	}
}

void Clonal_Expansion::Basic_Clonal_Expansion_V1(const unsigned int & ith_clone, std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> & Dying_and_Newborn, Random & r)
{
	std::vector<unsigned int> NewBorn_Cells;
	
	//1) Let's get basline parameters
	double P_DR =  Tumour -> at(ith_clone) -> P_Expansion[0];
	double P_NB =  Tumour -> at(ith_clone) -> P_Expansion[1];
	//2) Apply enviromental penalty to growth variables
	Adjust_Feedback_PR(P_NB);
	//3) Adjust Probability Mass
	double P_NT =  1.0 - (P_DR + P_NB );

	r.Basic_Clonal_Expansion_Sampling_V1( Tumour -> at (ith_clone) -> Clone_Size, P_DR, P_NB, P_NT, NewBorn_Cells);

	Update_Newborn_Parameters_V1(NewBorn_Cells, Dying_and_Newborn);

	Checking_for_Mutants( Dying_and_Newborn, r, ith_clone );

	Check_Clonal_Extinction_BCE_V1( ith_clone, Dying_and_Newborn );
	
}



void Clonal_Expansion::Dealyed_Mitosis_V1(const unsigned int & ith_clone, Random & r, const unsigned int & hours, const unsigned int & years)
{
	
	std::tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn (0, 0, 0, false, true);

	if( check_G1_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase--;
	else if( exiting_G1_phase(ith_clone) )
		Transition_From_G1_S(ith_clone, r);
	else if( check_S_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase--;	
	else if( exiting_S_phase(ith_clone) )
		Transition_From_S_G2(ith_clone, r);
	else if( check_G2_phase(ith_clone) )
		Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase--;
	else if( exiting_G2_phase(ith_clone) )
		Transition_From_G2_M(ith_clone, r);
	else if(Tumour -> at (ith_clone) -> In_M_Phase)
	{
		Basic_Clonal_Expansion_V1(ith_clone, Dying_and_Newborn, r);
		Update_Mutant_Mutational_Effects(ith_clone,  r,  hours, years, Dying_and_Newborn);
		Reset_to_G1(ith_clone);
	
	}

}



void Clonal_Expansion::Update_Population_V1(unsigned int & hours, unsigned int & years, Random & r)
{
	unsigned int ith_clone = 0;
	unsigned int number_of_clones = Tumour -> size();

	for(ith_clone = 0; ith_clone < number_of_clones ; ith_clone++)
	{
		if(!Tumour -> at (ith_clone) -> clone_extinct)// if not extinct
		{
			if( Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
			{
				Non_Mutagenic_Mitosis_Standard_V1(ith_clone, r);
			}
			else
			{
				Dealyed_Mitosis_V1(ith_clone, r, hours, years);	
			}
		}
	
	}//for
}



void Clonal_Expansion::print_Status( const unsigned int & hours, const unsigned int & years,  const bool & each_100)
{
	if(each_100 && (hours % 100 == 0))
	{	
		
		std::cout << " \n\n [ALIVE CELLS]: " <<  Population_Size  
			 << " [CLONES]: " << Tumour -> size() 
			 << "   H: " << hours  
			 << " Y: " << years 
			 << " FD: " << feedback  
			 << std::endl;
	}
	else
	{
		std::cout << " \n\n [ALIVE CELLS]: " <<  Population_Size  
			 << " [CLONES]: " << Tumour -> size() 
			 << "   H: " << hours  
			 << " Y: " << years 
			 << " FD: " << feedback  
			 << std::endl;
	}

}



// Load all properties function

//this will become the main
void Clonal_Expansion::Compute_Tumour_Growth_V1(const std::map<std::string, std::string> &logic)
{
	//finalise initilisasing parameters
	
	setSD_Penalty(logic.at("Penalty"));
	std::cout << "Penlaty type " << getSD_Penalty() << std::endl;
	std::cout << "PRINT TS: " << logic.at("Print_TS") << std::endl;
	bool print_time_step = false;
	if(logic.at("Print_TS") == "true")
		print_time_step = true;

	std::cout << "PRINT TS: " << print_time_step << std::endl;
	bool each_100 = true;
	getchar();


	Random r;
	unsigned int seconds = 0;
	unsigned int hours = 0;
	unsigned int years = 0;

	/* File stream */
	carcinogenesis_V1();
	printValues (Tumour -> at(0));
	while( Population_Size < 4000000000 && !( Population_Size == 0) )
	{
		seconds += dt;
		if(seconds == 3600)
		{
			seconds = 0; hours ++;
			Update_Population_V1(hours, years, r);
			Select_Size_Dependant_Penalty();
			
			if(print_time_step)
				print_Status( hours, years, each_100 );
		}
		if(hours == 8764)//8764
		{
			hours = 0; years ++;
		}
	}


}

void Clonal_Expansion::Compute_Tumour_Growth(const std::map<std::string, std::string> &logic)
{
	switch(Version)
 	{
 		case 0:
            std::cout << "0 Version: " << Version << std::endl;
            break;
        case 1:
            std::cout << "1 Version: " << Version << std::endl;
            Compute_Tumour_Growth_V1(logic);
            break;
        case 2:
            std::cout << "2 Version: " << Version << std::endl;
            break;
        case 3:
            std::cout << "3 Version: " << Version << std::endl;
            break;
        case 4:
            std::cout << "4 Version: " << Version << std::endl;
            break;
        case 5:
            std::cout << "5 Version: " << Version << std::endl;
            break;
        default:
        	 std::cout << "Def Version: " << Version << std::endl;
        	 break;
 	}//switch	
}


//void Clonal_Expansion::Grow_Tumour()
//{

//}


