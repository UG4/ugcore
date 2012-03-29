//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#include <cstring>
#include "loader_util.h"

using namespace std;

namespace ug
{

void split_parameters(std::vector<std::string>& paramsOut, char* strParams, const char* delims)
{
	paramsOut.clear();
	char* param = strtok(strParams, delims);
	while(param != NULL)
	{
		if(*param != '\0')
			paramsOut.push_back(string(param));
		param = strtok(NULL, delims);
	}
}

}
