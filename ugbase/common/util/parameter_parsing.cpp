// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)
 
#include <cstring>
#include <cstdlib>
#include "parameter_parsing.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
int GetParamIndex(const char* param, int argc, const char * const * argv) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(param, argv[i]) == 0) {
			return i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
bool FindParam(const char* param, int argc, const char * const * argv) {
	return GetParamIndex(param, argc, argv) != -1;
}


////////////////////////////////////////////////////////////////////////
bool ParamToInt(int& iOut, const char* param, int argc, const char * const * argv) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	iOut = atoi(argv[i + 1]);
	return true;
}

////////////////////////////////////////////////////////////////////////
bool ParamToDouble(double &dOut, const char *param, int argc, const char * const * argv)
{
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	dOut = atof(argv[i + 1]);
	return true;
}


////////////////////////////////////////////////////////////////////////
bool ParamToString(const char ** strOut, const char* param, int argc, const char * const * argv) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	*strOut = argv[i + 1];
	return true;
}

int ParamToInt(const char* param, int argc, const char * const * argv, int iDefault)
{
	int i;
	if(ParamToInt(i, param, argc, argv) == true)
		return i;
	else return iDefault;
}
double ParamToDouble(const char *param, int argc, const char * const * argv, double dDefault)
{
	double d;
	if(ParamToDouble(d, param, argc, argv) == true)
		return d;
	else return dDefault;
}

}// end of namespace
