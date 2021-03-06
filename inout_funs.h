#ifndef GLOBAL_IN_OUT_FUN
#define GLOBAL_IN_OUT_FUN

#include "config.h"
#include <iostream>         // std::cout,std::endl
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <map>

/**	FUNCTION get_Clone_DS()
		
	This function returns a 
	unique pointer of type Clone. Remmeber
	that unique pointers are memory efficient
	and guaratee a unique scope of the DS.

********************************************/


bool FileExists (const std::string& name);
bool FolderExists(const std::string& pathname);
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) ;
 void ltrim(std::string &s);
 void rtrim(std::string &s); 
 void trim(std::string &s);

 std::string ltrimmed(std::string s);
std::string rtrimmed(std::string s);
 std::string trimmed(std::string s) ;
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

std::map<std::string, std::string> ReadAndParse( std::string fileName);
std::vector<std::string> UniqueElements(std::vector<std::string> vec);
std::vector<std::string> getKeysFromMap( std::map<std::string, std::string> rawData);
void ValidateElements(std::vector<std::string> *keys, std::map<std::string, std::string> *rawData );

std::map<std::string, std::string> ReadTextFile( std::string fileName);



#endif