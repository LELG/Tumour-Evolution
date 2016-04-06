/*
	Thisis a header file for the Clone Structure
*/

#ifndef CLONAL_EXPANSION_H
#define CLONAL_EXPANSION_H

#include "Clone.h"
#include <string>
#include <map>
#include <memory>           // std::unique_ptr
#include <vector>

class Clonal_Expansion
{

	/**
	    Variable Definitions
	**************************/
public:
	unsigned long long int Population_Size;

	double feedback; // This is the output of a proportional control system
		
		/**
			This is the main Data Structure to store all clones, this will
			reflect the heterogeneity.

			NOTE: This is a map of of unorder maps, this will allow us to model
			micro enviroments. The layer can be computed in MPI
		*/
		//unordered_map<string, unique_ptr<Clone> > Tumour;
		//unordered_map<int, struct Clone*> Tumour;
		/* We can create as many constructures as we want */

		 std::vector<std::unique_ptr<Clone> > *Tumour = new std::vector<std::unique_ptr<Clone> >;



public:
	Clonal_Expansion();
	~Clonal_Expansion();
	void printParameters();
	void add_Clone_to_Tumour(void);
	void carcinogenesis(void);
};

#endif