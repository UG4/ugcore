//	string_util.h
//	created by Andreas Vogel

#ifndef __H__COMMON_STRING_UTIL__
#define __H__COMMON_STRING_UTIL__

#include <string>
#include <vector>

namespace ug{

// help function to tokenize the parameter string
void TokenizeString(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);

// help function to remove white space from string
void RemoveWhitespaceFromString(std::string& string);

} // end namespace ug

#endif /*__H__COMMON_STRING_UTIL__*/
