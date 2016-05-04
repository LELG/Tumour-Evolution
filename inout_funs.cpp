#include "config.h"
#include <iostream>         // std::cout,std::endl
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <map>

/*
	This function validates that the faile exists
*/
bool FileExists (const std::string& name) 
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

bool FolderExists(const std::string& path)
{
	struct stat info;
	bool exists = false;
	const char * pathname = path.c_str();

	if( stat( pathname, &info ) != 0 )
	{
    	std::cout << "cannot access" << pathname << std::endl;
	}
	else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
	{
   	 	std::cout << "is a directory: " <<  pathname << std::endl;
   	 	exists = true;
	}
	else
	{
    	std::cout << "is not a directory: " << pathname << std::endl;
	}

	return exists;
}


// trim from start (in place)
 void ltrim(std::string &s) 
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
 void rtrim(std::string &s) 
{
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

/*
	String trim both ends
*/

// trim from both ends (in place)
 void trim(std::string &s) 
{
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
 std::string ltrimmed(std::string s) 
{
    ltrim(s);
    return s;
}

// trim from end (copying)
 std::string rtrimmed(std::string s) 
{
    rtrim(s);
    return s;
}

// trim from both ends (copying)
 std::string trimmed(std::string s) 
{
    trim(s);
    return s;
}

/*
  Parsing options 
*/
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) 
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) 
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::map<std::string, std::string> ReadAndParse( std::string fileName)
{
	std::ifstream file(fileName);
	std::string line;
	std::string file_contents;
	std::map<std::string, std::string> rawData;

	while ( std::getline(file, line) )
	{
		std::vector<std::string> line_splited;
		if ( !line.empty() )
		{
			split(line, '=',line_splited) ;

			if( line_splited.size() == 2 )
				rawData[ trimmed(line_splited.at(0)) ] = trimmed(line_splited.at(1));
		}
	}

	return rawData;

}

std::vector<std::string> getKeysFromMap( std::map<std::string, std::string> rawData)
{
	std::vector<std::string> keys;
	for(std::map<std::string, std::string>::iterator it = rawData.begin(); it != rawData.end(); ++it)
		keys.push_back( it -> first );
		//std::cout << it -> first << std::endl;

	return keys;
}

std::vector<std::string> UniqueElements(std::vector<std::string> vec)
{
	sort( vec.begin(), vec.end() );
	vec.erase( unique( vec.begin(), vec.end() ), vec.end() );

	return vec;
}

void ValidateElements(std::vector<std::string> *keys, std::map<std::string, std::string> *rawData )
{
	for(std::vector<std::string>::iterator it = keys -> begin (); it != keys -> end (); ++it )
	{
		if(!CLONE_VARIABLES.count(*it))
		{
			std::cout << "This entry " << *it <<  " is not valid ... erasing" << std::endl;
			rawData -> erase(*it);
		}	
	}
}

/*
	This function reads a file and print its content
*/
std::map<std::string, std::string> ReadTextFile( std::string fileName)
{	
	std::map<std::string, std::string> rawData = ReadAndParse( fileName );
	std::vector<std::string> keys = UniqueElements( getKeysFromMap( rawData ) );
	ValidateElements(&keys, &rawData);

	return rawData;
}

