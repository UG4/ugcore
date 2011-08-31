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


/**
 * Replaces each substring of <code>target</code> string that is equal to
 * <code>oldstr</code> with <code>newstr</code>
 * @param target string to modify
 * @param oldstr string to raplace
 * @param newstr replacement string
 * @return a copy of the specified <code>target</code> string where
 *         all occurences of <code>oldstr</code> are replaced with
 *         <code>newstr</code>
 */
std::string ReplaceAll(
		std::string target,
		const std::string oldstr,
		const std::string newstr);


/**
 * Checks whether <code>str</code> starts with <code>search</code>.
 * @param str string
 * @param search string to search
 * @return <code>true</code> if <code>str</code> starts
 * with <code>search</code>; <code>false</code> otherwise
 */
bool StartsWith(std::string str, std::string search);

/**
 * Checks whether <code>str</code> contains <code>search</code>.
 * @param str string
 * @param search string to search
 * @return <code>true</code> if <code>str</code> contains
 * <code>search</code>; <code>false</code> otherwise
 */
bool Contains(std::string str, std::string search);

/**
 * Checks whether the specified name is a valid identifier name.
 * <p>
 * <b>Note:</b> identifiers starting with <code>F_</code>, <code>C_</code>,
 * <code>I_</code> or containing <code>__</code>
 *  (double underscore) are invalid.
 * </p>
 * @param name name to check
 * @return  <code>true</code> if the specified name is valid;
 *          <code>false</code> otherwise
 */
bool IdentifierIsValid(std::string name);

/**
 * Returns a message describing which identifiers are valid and which ar not.
 * @return message string
 */
std::string GetIdentifierMessage();

} // end namespace ug

#endif /*__H__COMMON_STRING_UTIL__*/
