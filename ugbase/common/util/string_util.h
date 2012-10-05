//	string_util.h
//	created by Andreas Vogel

#ifndef __H__COMMON_STRING_UTIL__
#define __H__COMMON_STRING_UTIL__

#include <string>
#include <vector>
#include <algorithm> 
#include <cctype>
#include "hash.h"
#include "common/ug_config.h"

namespace ug{

// help function to tokenize the parameter string
UG_API void TokenizeString(const std::string& str, std::vector<std::string>& tokens, const char delimiter);

// help function to remove white space from string
UG_API void RemoveWhitespaceFromString(std::string& string);

// help function to remove whitespace from front and end of string
UG_API std::string TrimString(const std::string& str);

/// returns the number of digits of an integer (expressed with base 10)
/**
 * This functions returns the number of digits for the passed number. A minus
 * sign is ignored.
 *
 * \param[in]	n 		number to count the number of digits
 */
UG_API int NumberOfDigits(int n);

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
UG_API void AppendCounterToString(std::string& str, std::string indicator,
                           int counter, int maxCounter = -1);

UG_API std::string AppendSpacesToString(std::string& str, int nSpaces);

//sreiter
///	this template function creates a hash key for a string value.
template <> UG_API unsigned long hash_key(const std::string& key);

/// returns the part of the string after the last '/' or '\' character (e.g. "/sw/bla.txt" -> "bla.txt")
UG_API std::string FilenameWithoutPath(const std::string& str);

/// returns the part of the string before the last '/' or '\' character (e.g. "/sw/bla.txt" -> "/sw/")
UG_API std::string PathFromFilename(const std::string &str);

/// returns the filename without path and with extension (e.g. "/sw/bla.txt" -> "bla")
UG_API std::string FilenameWithoutExtension(std::string str);

/// returns the extension of the filename (e.g. "/sw/bla.txt" -> "txt")
UG_API std::string GetFilenameExtension(const std::string &str);

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
UG_API std::string ReplaceAll(
		std::string target,
		const std::string& oldstr,
		const std::string& newstr);


/**
 * Checks whether <code>str</code> starts with <code>search</code>.
 * @param str string
 * @param search string to search
 * @return <code>true</code> if <code>str</code> starts
 * with <code>search</code>; <code>false</code> otherwise
 */
UG_API bool StartsWith(const std::string& str, const std::string& search);

/**
 * Checks whether <code>str</code> contains <code>search</code>.
 * @param str string
 * @param search string to search
 * @return <code>true</code> if <code>str</code> contains
 * <code>search</code>; <code>false</code> otherwise
 */
UG_API bool Contains(const std::string& str, const std::string& search);

/**
 * Returns a lower case version of the specified string.
 * <p><b>Note: </b>this function does not support custom locales. Thus,
 *                 only ascii strings shall be specified.</p> 
 * @param str string to convert
 * @return a lower case version of the specified string
 */
UG_API std::string ToLower(std::string str);

/**
 * Returns an upper case version of the specified string.
 * <p><b>Note: </b>this function does not support custom locales. Thus,
 *                 only ascii strings shall be specified.</p> 
 * @param str string to convert
 * @return an upper case version of the specified string
 */
UG_API std::string ToUpper(std::string str);

/**
 * Searches for duplicates in the specified vector and returns a vector 
 * containing all elements that occur multiple times.
 * @param vec vector to analyze
 * @return a vector containing all elements that occur multiple times
 */
UG_API std::vector<std::string> FindDuplicates(const std::vector<std::string>& vec);

/**
 * @param c  the character
 * @param nr number of times to repeat c
 * @return string with nr times c
 */
UG_API std::string repeat(char c, int nr);

/**
 * Levenshtein distance calculates the minimum number of edits to transform one string
 * into the other with allowable edit operations insertion, deletion, or substitution of a
 * single character.
 * @param s1 string 1
 * @param s2 string 2
 * @return minimum number of edits needed to transform one string into the other
 * see http://en.wikipedia.org/wiki/Levenshtein_distance
 */
UG_API size_t LevenshteinDistance( const std::string& s1, const std::string& s2 );


/**
 * \brief function to get some lines of a file
 * \param filename				file name
 * \param fromline
 * \param toline
 * \param includeLineNumbers	if true, add the line number in front of each line and a tab.
 *
 * \return lines fromline to toline of file filename.
 */
UG_API std::string GetFileLines(const char *filename, size_t fromline, size_t toline, bool includeLineNumbers);

/**
 * \brief function to get a line of a file
 * \param filename				file name
 * \param line
 * \return the line of the file.
 */
UG_API std::string GetFileLine(const char *filename, size_t line);

/**
 * IsLonger can be used to get the longest string in a vector of strings:
 * int maxLength = (*max_element(vecStr.begin(), vecStr.end(), IsLonger)).size();
 * @param a
 * @param b
 * @return true if b is longer then a
 */
UG_API bool IsLonger(const std::string &a, const std::string &b);

} // end namespace ug

#endif /*__H__COMMON_STRING_UTIL__*/
