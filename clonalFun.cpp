#include "clonalFun.h"
#include "Clone.h"
#include <memory>           // std::unique_ptr
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <sys/stat.h>
#include <unistd.h>


std::unique_ptr<Clone> get_Clone_DS()
{
    return std::unique_ptr<Clone>( new Clone( ) );
} // end function


inline bool FileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}