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

void test_external_carcinogenesis();
void test_external_carcinogenesis_with_size();
void test_member_carcinogenesis();
void test_member_carcinogenesis_with_size();

int main(int argc, char**argv)
{
	//***** SET UP Random ********//
	long seed;
	r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  	seed = time (NULL) * getpid();    
  	gsl_rng_set (r_global, seed);  
  	//***** SET UP Random ********//

	test_external_carcinogenesis();
	test_external_carcinogenesis_with_size();
	test_member_carcinogenesis();
	test_member_carcinogenesis_with_size();

	std::cout << "ALL TESTS PASSED " << std::endl;

	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
	carcinogenesis(tmr);

	std::cout << "IN G1 Phase: " << tmr -> check_G1_phase( tmr -> Tumour-> at(0) ) << std::endl;
	std::cout << "EXITING G1 Phase: " << tmr -> exiting_G1_phase( tmr -> Tumour-> at(0) ) << std::endl;
	std::cout << "IN S Phase: " << tmr -> check_S_phase( tmr -> Tumour-> at(0) ) << std::endl;
	std::cout << "EXITING S Phase: " << tmr -> exiting_S_phase( tmr -> Tumour-> at(0) ) << std::endl;
	std::cout << "IN G2 Phase: " << tmr -> check_G2_phase( tmr -> Tumour-> at(0) ) << std::endl;
	std::cout << "EXITING G2 Phase: " << tmr -> exiting_G2_phase( tmr -> Tumour-> at(0) ) << std::endl;

	tmr -> Transition_From_G1_S(tmr -> Tumour-> at(0));
	std::cout << "IN S Phase: " << tmr -> check_S_phase( tmr -> Tumour-> at(0) ) << std::endl;
	tmr -> Transition_From_S_G2(tmr -> Tumour-> at(0));
	std::cout << "IN G2 Phase: " << tmr -> check_G2_phase( tmr -> Tumour-> at(0) ) << std::endl;
	tmr -> Transition_From_G2_M(tmr -> Tumour-> at(0));
	std::cout << "IN M Phase: " << tmr -> Tumour -> at(0) -> In_M_Phase << std::endl;

	tmr -> Non_Mutagenic_Mitosis_Standard( tmr -> Tumour -> at(0));

	for(int i = 0; i < 10 ; i++)
		tmr -> Update_Population_After_Division(tmr -> Tumour -> at(0));

	unsigned long long int test = 87932;

	unsigned int tc = static_cast<unsigned int>(test);
	unsigned int max_unsigned_int_size = std::numeric_limits<unsigned int>::max();

	std::cout << "test value: " << test << " \n" << "test casted: " << tc << "\n" << "max " << max_unsigned_int_size << std::endl;

	std::cout << "MS " << MAX_CLONE_SIZE << std::endl;

	Random r;

	unsigned int * buffer;
	double p = 0.002;
	buffer = r.Newborn_G0_and_G1(test, p);

	std::cout << "Values: " << buffer[0] << " " << buffer[1] << std::endl;

	unsigned int *buffer_K;
	double p_staying_G0 = 0.33; 
	double p_dying_in_G0 = 0.33;
	double p_exiting_G0 = 0.33;

	buffer_K = r.Update_G0_Phase(tc, p_staying_G0, p_dying_in_G0, p_exiting_G0 );
	std::cout << "Values G0: " << buffer_K[0] << " " << buffer_K[1] << " " << buffer_K[2] << std::endl;

	unsigned int *buffer_G1;
	double p_staying_G1 = 0.2;
	buffer_G1 = r.Update_G1_Phase(tc, p_staying_G1);
	std::cout << "Values G1: " << buffer_G1[0] << " " << buffer_G1[1] << std::endl;

	unsigned int *buffer_G2;
	double p_staying_G2 = 0.3;
	buffer_G2 = r.Update_G2_Phase(tc, p_staying_G2);
	std::cout << "Values G2: " << buffer_G2[0] << " " << buffer_G2[1] << std::endl;


	unsigned int *buffer_S;
	double p_staying_S = 0.5;
	buffer_S = r.Update_S_Phase(tc, p_staying_S);
	std::cout << "Values S: " << buffer_S[0] << " " << buffer_S[1] << std::endl;

	unsigned int *buffer_M;
	double p_staying_M = 0.1;
	buffer_M = r.Update_S_Phase(tc, p_staying_M);
	std::cout << "Values M: " << buffer_M[0] << " " << buffer_M[1] << std::endl;

	//unsigned int * Random::Update_S_Phase(unsigned int & S_cells, double & p_staying_S)

	return 0;
}// End of main




/**
	Test Functions for the multiple Types
	of carcinogenesis (member and free functions)
**************************************************/
void test_external_carcinogenesis()
{
	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
	carcinogenesis(tmr);
	std::cout << "\nPRINT DS EXTERANL POP SIZE 1" << std::endl; 
	tmr -> printParameters();
	 assert(tmr -> Population_Size == 1);
	std::cout << "\nTesting External Carcinogenesis with Population Size = 1: PASSED"  << std::endl;

}

void test_external_carcinogenesis_with_size()
{
	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
	unsigned long long int neoplastic_cells = 100;
	carcinogenesis(tmr, neoplastic_cells);
	std::cout << "\nPRINT DS EXTERANL POP SIZE 100" << std::endl; 
	tmr -> printParameters();
	assert(tmr -> Population_Size == 100);
	std::cout << "\nTesting External Carcinogenesis with Population Size = " << neoplastic_cells << " : PASSED"  << std::endl;
}

void test_member_carcinogenesis()
{
	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
	tmr -> carcinogenesis();
	std::cout << "\nPRINT DS CLASS POP SIZE 1" << std::endl; 
	tmr -> printParameters();
	assert(tmr -> Population_Size == 1);
	std::cout << "\nTesting Class Carcinogenesis with Population Size = 1 "  << " : PASSED"  << std::endl;
}

void test_member_carcinogenesis_with_size()
{
	std::unique_ptr<Clonal_Expansion>  tmr = get_Clonal_Expasion_DS();
	unsigned long long int neoplastic_cells = 700;
	tmr -> carcinogenesis(neoplastic_cells);
	std::cout << "\nPRINT DS CLASS POP SIZE 700" << std::endl; 
	tmr -> printParameters();
	assert(tmr -> Population_Size == 700);
	std::cout << "\nTesting Class Carcinogenesis with Population Size = 700 "  << " : PASSED"  << std::endl;
}





