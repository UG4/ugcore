//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#ifndef __H__UTIL__LOADER_UTIL__
#define __H__UTIL__LOADER_UTIL__

#include <vector>
#include <string>

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

///	fills paramsOut with the parameters that are seperated by the characters given in delim.
/**
 * if you pass " ,.-" to delims, then paramsOut will contain all
 * parameters that are separated by ' ', ',', '.' or '-'.
 *
 * Please note that strParams is not const!
 */
void split_parameters(std::vector<std::string>& paramsOut, char* strParams, const char* delims = " ");

// end group ugbase_common_util
/// \}

}

#endif
