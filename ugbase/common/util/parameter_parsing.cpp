// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)
 
#include <cstring>
#include <cstdlib>
#include "parameter_parsing.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
int GetParamIndex(const char* param, int argc, char* argv[]) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(param, argv[i]) == 0) {
			return i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
bool FindParam(const char* param, int argc, char* argv[]) {
	return GetParamIndex(param, argc, argv) != -1;
}

////////////////////////////////////////////////////////////////////////
bool ParamToInt(int& iOut, const char* param, int argc, char* argv[]) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	iOut = atoi(argv[i + 1]);
	return true;
}

////////////////////////////////////////////////////////////////////////
bool ParamToString(char** strOut, const char* param, int argc, char* argv[]) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	*strOut = argv[i + 1];
	return true;
}

}// end of namespace
