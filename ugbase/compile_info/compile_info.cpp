#include "compile_info.h"

namespace ug
{
//	the following variables are defined through a script at build-time
//	and are stored in a file in your build folder (compile_info_vars.cpp)
extern const char *UG_SVN_REVISION;
extern const char *UG_BUILD_HOST;
extern const char *UG_COMPILE_DATE;


const char* UGSvnRevision()
{
	return UG_SVN_REVISION;
}

const char* UGBuildHost()
{
	return UG_BUILD_HOST;
}

const char* UGCompileDate()
{
	return UG_COMPILE_DATE;
}

}// end of namespace

