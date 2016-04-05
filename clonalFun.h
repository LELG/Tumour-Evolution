#ifndef GLOBAL_CLONAL_FUN
#define GLOBAL_CLONAL_FUN

#include "Clone.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <memory>           // std::unique_ptr

/**	FUNCTION get_Clone_DS()
		
	This function returns a 
	unique pointer of type Clone. Remmeber
	that unique pointers are memory efficient
	and guaratee a unique scope of the DS.

********************************************/

std::unique_ptr<Clone> get_Clone_DS();



#endif