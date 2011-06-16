//	string_util.h
//	created by Andreas Vogel

#ifndef __H__COMMON_STRING_UTIL__
#define __H__COMMON_STRING_UTIL__

#include <string>
#include <vector>
#include "hash.h"

namespace ug{

// help function to tokenize the parameter string
void TokenizeString(const std::string& str, std::vector<std::string>& tokens, const char delimiter);

// help function to remove white space from string
void RemoveWhitespaceFromString(std::string& string);

// help function to remove whitespace from front and end of string
std::string TrimString(const std::string& str);

/// returns the number of digits of an integer (expressed with base 10)
/**
 * This functions returns the number of digits for the passed number. A minus
 * sign is ignored.
 *
 * \param[in]	n 		number to count the number of digits
 */
int NumberOfDigits(int n);

///	appends a counter number to a string
/**
 * This functions appends to a string a counter preceded by some indicator. If
 * a maxCounter is passed, the field is adjusted to the maximum needed width
 * and additional space is filled by zeros.
 *
 * \param[in, out]	str			string to append the counter
 * \param[in]		indicator	some string preceding the counter
 * \param[in]		counter		counter added
 * \param[in]		maxCounter	maximum counter to be added
 */
void AppendCounterToString(std::string& str, std::string indicator,
                           int counter, int maxCounter = -1);

//sreiter
///	this template function creates a hash key for a string value.
template <> unsigned long hash_key(const std::string& key);

/// returns the part of the string after the last '/' character
std::string FilenameWithoutPath(const std::string& str);

} // end namespace ug

#endif /*__H__COMMON_STRING_UTIL__*/
