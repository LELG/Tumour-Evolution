#include "config.h"
#include <map>
#include <string>

const double MUT_RATE           = 0.0000001; 
const double DEATH_RATE         = 0.02;
const double PROLIFERATION_RATE = 0.03;
// We can overide thios in a config file, this are default
const std::map<std::string, bool> CLONE_VARIABLES = {	{"Mutation_Rate", true}, \
														{"Proliferation_Rate", true}, \
														{"Death_Rate", true}, \
														{"Number_of_Memebers_to_Start_Heterogeneity", true},\
														{"Generation_ID_Counter", true}, \
														{"Driver_10_fold_accumulation", true}, \
														{"clone_extinct", true}, \
														{"Number_of_Mutations", true}, \
														{"Clone_Size", true}, \
														{"Initiall_Expasion_Period", true}, \
														{"In_G0_Phase", true}, \
														{"In_G1_Phase", true}, \
														{"In_S_Phase", true}, \
														{"In_G2_Phase", true}, \
														{"In_M_Phase", true}, \
														{"Remaining_Time_in_G1_Phase", true}, \
														{"Remaining_Time_in_S_Phase", true}, \
														{"Remaining_Time_in_G2_Phase", true}, \
														{"Remaining_Time_in_M_Phase", true}, \
														{"init_PR", true}, \
														{"max_PR", true}, \
														{"final_PR", true}, \
														{"Generation_ID", true}\
														 };														
														
