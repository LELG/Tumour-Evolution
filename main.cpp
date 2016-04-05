#include "Clone.h"
#include "clonalFun.h"
#include "config.h"
#include <iostream>         // std::cout,std::endl
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <stdexcept>        // std::exception, std::runtime_error
#include <memory>           // std::unique_ptr

int main(int argc, char**argv)
{

	std::unique_ptr<Clone> const c = get_Clone_DS();
	std::cout << "NO BREAK " << std::endl;
  	std::cout << MUT_RATE << std::endl;

  	c -> printTest();
  	c -> FileReader("settings.txt");

  return 0;
}