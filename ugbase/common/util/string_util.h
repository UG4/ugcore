//	string_util.h
//	created by Andreas Vogel

#ifndef __H__COMMON_STRING_UTIL__
#define __H__COMMON_STRING_UTIL__

#include <string>
#include <vector>
#include <algorithm> 
#include <sstream>
#include <cctype>
#include "hash_function.h"
#include "common/ug_config.h"

namespace ug{

/**
 * \defgroup ugbase_common_util_strings String Utilities
 * \ingroup ugbase_common_util
 * \{
 */

/**
 * \brief splits the string into parts based on a separating char
 * \details the string parts using a separator char in order to indicate parts
 * \note any prior content of \c vToken will get deleted
 * \param[in]     str       original string
 * \param[in,out] vToken    tokenized parts
 * \param[in]     delimiter char used as separator
 */
UG_API void TokenizeString( const std::string& str, std::vector<std::string>& vToken, 
                            const char delimiter=',' );

/**
 * \brief splits the string into parts based on a separating character
 * \details returns the string parts using a separator char in order to indicate parts
 * \param[in] str       original string
 * \param[in] delimiter char used as separator
 * \return tokenized parts
 */
UG_API std::vector<std::string> TokenizeString( const std::string& str,
                                                const char delimiter=',' );

/**
 * \brief splits the string into parts based on a separating character
 * \details returns the string parts using a separator char in order to indicate parts
 * \param[in] str       original string
 * \param[in] delimiter char used as separator
 * \return tokenized parts
 */
UG_API std::vector<std::string> TokenizeString( const char* str,
                                                const char delimiter=',' );

/**
 * \brief splits the string into trimmed parts based on a separating char
 * \details returns the string parts separated by \c delimiter and trims all parts
 * \param[in] str       original string
 * \param[in] delimiter char used as separator
 * \return tokenized and trimmed parts
 */
UG_API std::vector<std::string> TokenizeTrimString( const std::string& str, 
                                                    const char delimiter=',' );

/**
 * \brief removes all white space from a string, also within the string
 * \param[in,out] string the string to modify
 */
UG_API void RemoveWhitespaceFromString(std::string& string);

/**
 * \brief removes all white space from the front and end of a string
 * \param[in] string the string to modify
 * \return the modified string
 */
UG_API std::string TrimString(const std::string& str);

/**
 * \brief creates a truncated string and may add truncation symbol at end
 * \param[in] string 		the string to modify
 * \param[in] totalSize		the total size of snippet
 * \param[in] replaceLast	the number of last chars to be replaced by symbol
 * \param[in] replace		the replace symbol
 * \return the modified string
 */
UG_API std::string SnipString(const std::string& str, size_t totalSize,
                              size_t replaceLast = 0, const char replace = '.');

/**
 * \brief creates a truncated string and may add truncation symbol at front
 * \param[in] string 		the string to modify
 * \param[in] totalSize		the total size of snippet
 * \param[in] replaceLast	the number of last chars to be replaced by symbol
 * \param[in] replace		the replace symbol
 * \return the modified string
 */
UG_API std::string SnipStringFront(const std::string& str, size_t totalSize,
                                   size_t replaceFront = 0, const char replace = '.');

/**
 * \brief returns the number of digits of an integer (expressed with base 10)
 * \details Determines the number of digits for the passed base-10 number.
 *   A minus sign is ignored.
 * \param[in] n number to count the number of digits
 * \returns number of digits
 */
UG_API int NumberOfDigits(int n);

/**
 * \brief appends a counter number to a string
 * \details This functions appends to a string a counter preceded by some 
 *   indicator.
 *   If a \c maxCounter is passed, the field is adjusted to the maximum needed 
 *   width and additional space is filled by zeros.
 * \param[in,out] str        string to append the counter
 * \param[in]     indicator  some string preceding the counter
 * \param[in]     counter    counter added
 * \param[in]     maxCounter maximum counter to be added
 */
UG_API void AppendCounterToString( std::string& str, std::string indicator,
                                   int counter, int maxCounter=-1 );

/**
 * \brief padding a string with spaces to predefined length
 * \details Appends spaces to the given string so that the resulting string has
 *   a predefined length of \c totalLength
 * \param[in] str         string to be padded
 * \param[in] totalLength desired total length of the string
 * \returns padded string
 */
UG_API std::string AppendSpacesToString(std::string& str, int totalLength);

/**
 * \brief creates a hash key from a string value
 * \details this template function creates a hash key for a string value
 * \param[in] str string to create hash for
 * \returns hash key for given \c key
 * \note Implementation is copied from some book or website. Can't remember... (sreiter)
 */
template <> UG_API size_t hash_key(const std::string& str);

/**
 * \brief determines last occurrence of '/' or '\'
 * \param[in] str string to lock in
 * \returns position of the last occurrence of '/' or '\' in \c str; 
 *   returns `std::string::npos` if none are found
 */
std::string::size_type GetDirectorySeperatorPos(const std::string &str);

/**
 * \brief returns best guess of a filename from a given string
 * \details returns the part of the string after the last '/' or '\' character 
 *   (e.g. `/sw/bla.txt` -> `bla.txt`)
 * \param[in] str to retrieve the filename from
 * \return best guess of the file name from given path; if no guess can be made
 *   the complete string is returned
 */
UG_API std::string FilenameWithoutPath(const std::string &str);

/**
 * \brief returns best guess of a path without a filename from a given string
 * \details returns the part of the string before the last '/' or '\' character 
 *   (e.g. `/sw/bla.txt` -> `/sw/`)
 * \param[in] str to retrieve the filename from
 * \return best guess of the file name from given path; if no guess can be made
 *   '.' is returned
 */
UG_API std::string PathFromFilename(const std::string &str);

/**
 * \brief returns the best guess of the filename from given string
 * \details returns the part of the string without path and extension
 *   (e.g. `/sw/bla.txt` -> `bla`)
 * \param[in] str to retrieve filename from
 * \returns best guess of the filename without path and extension; if no guess 
 *   can be made, the whole string is returned
 */
UG_API std::string FilenameWithoutExtension(std::string str);

/**
 * \brief returns the best guess of a file extensions from given string
 * \details returns the extension of the filename (e.g. `/sw/bla.txt` -> `txt`).
 *   Everything after the last dot ('.') of \c str is considered the file extension.
 * \param[in] str to retrieve file extension from
 * \returns best guess of the file extension; empty string if no guess can be made
 */
UG_API std::string GetFilenameExtension(const std::string &str);

/**
 * \brief Substitutes substrings of given string with other substrings
 * \details Replaces each substring of \c target string that is equal to \c oldstr 
 *   with \c newstr
 * \param[in] target string to modify
 * \param[in] oldstr string to raplace
 * \param[in] newstr replacement string
 * \return a copy of the specified \c target string where all occurences of 
 *   \c oldstr are replaced with \c newstr.
 */
UG_API std::string ReplaceAll( std::string target, const std::string& oldstr, 
                               const std::string& newstr );

/**
 * \brief checks whether a given string starts with a specified substring
 * \details Checks whether \c str starts with \c search.
 * \param[in] str    string
 * \param[in] search string to search for
 * \return \c true if \c str starts with \c search; \c false otherwise
 */
UG_API bool StartsWith(const std::string& str, const std::string& search);

/**
 * \brief Checks whether given string contains a specified substring
 * \details Checks whether \c str contains \c search.
 * \param[in] str    string
 * \param[in] search string to search for
 * \return \c true if \c str contains \c search; \c false otherwise
 */
UG_API bool Contains(const std::string& str, const std::string& search);

/**
 * \brief Returns a lower case version of the specified string.
 * \note this function does not support custom locales.
 *   Thus, only ascii strings shall be specified.
 * \param[in] str string to convert
 * \return a lower case version of the specified string
 */
UG_API std::string ToLower(std::string str);

/**
 * \brief Returns an upper case version of the specified string.
 * \note this function does not support custom locales.
 *   Thus, only ascii strings shall be specified.
 * \param[in] str string to convert
 * \return an upper case version of the specified string
 */
UG_API std::string ToUpper(std::string str);

/**
 * \brief Finds and returns all duplicate elements of given vector
 * \details Searches for duplicates in the specified vector and returns a vector 
 *   containing all elements that occur multiple times.
 * \param[in] vec vector to analyze
 * \return a vector containing all elements that occur multiple times
 */
UG_API std::vector<std::string> FindDuplicates(const std::vector<std::string>& vec);

/**
 * \brief Builds a string with specified repetitions of given character
 * \param[in] c  the character
 * \param[in] nr number of times to repeat \c c
 * \return string with \c nr times \c c
 */
UG_API std::string repeat(char c, int nr);

/**
 * \brief Calculate Levenshtein Distance of to strings
 * \details Levenshtein distance calculates the minimum number of edits to 
 *   transform one string into the other with allowable edit operations 
 *   insertion, deletion, or substitution of a single character.
 * \note taken from http://en.wikipedia.org/wiki/Levenshtein_distance
 *   (check copyright or recreate!)
 * \param[in] s1 string 1
 * \param[in] s2 string 2
 * \return minimum number of edits needed to transform one string into the other
 */
UG_API size_t LevenshteinDistance( const std::string& s1, const std::string& s2 );


/**
 * \brief get some specified lines of a file
 * \param[in] filename           file name
 * \param[in] fromline           line number to start from
 * \param[in] toline             line number to stop at
 * \param[in] includeLineNumbers if true, add the line number in front of each 
 *   line and a tab.
 * \return lines fromline to toline of file filename.
 */
UG_API std::string GetFileLines( const char *filename, size_t fromline, size_t toline, 
                                 bool includeLineNumbers );

/**
 * \brief get a specific line of a file
 * \param filename file name
 * \param line     line number to extract
 * \return the line of the file
 */
UG_API std::string GetFileLine(const char *filename, size_t line);

/**
 * \brief checks whether second string is longer than first string
 * \details This can be used to get the longest string in a vector of strings:
 * 
 *     int maxLength = (*max_element(vecStr.begin(), vecStr.end(), IsLonger)).size();
 * \param[in] a
 * \param[in] b
 * \return \c true if \c b is longer then \c a; \c false otherwise
 */
UG_API bool IsLonger(const std::string &a, const std::string &b);


/**
 * \brief Convert a object supporting '`std::cout << obj`' to a string
 * \tparam T type of the object; must support `std::ostream operator<<()`
 * \param[in] t object to convert to string
 * \return a string with the object as if you would use operator << (like `std::cout`)
 */
template<typename T>
inline std::string ToString(const T &t)
{
    std::stringstream out;
    out << t;
    return out.str();
}

/**
 * \brief returns a string suitable for XML files
 * this functions escapes the characters <, >, ', " and &
 * @sa http://www.hdfgroup.org/HDF5/XML/xml_escape_chars.htm
 * @param[in] s
 * @return escaped string
 */
UG_API std::string XMLStringEscape(std::string s);

/**
 * \brief wildcard matches like bla.* or *.bla or t?st
 * @param[in] str a string
 * @param[in] pattern a pattern with wildcards * or ?
 * @return true if match otherwise false
 */
UG_API bool WildcardMatch(const char *str, const char *pattern);

/**
 * this function replaces XML special characters with their escaped versions:
 * & -> &amp;
 * " -> &quot;
 * ' -> "&apos;
 * < -> &lt;
 * > -> &gt;
 * @param s a normal text
 * @return a text where special XML characters are escaped
 */
UG_API std::string XMLStringEscape(std::string s);

// end group ugbase_common_util_strings
/// \}

/**
 * \brief returns a "shifted" string
 * one-line strings are not shifted
 * two line strings are shifted like this: input:
 * "MyLine1\nMyLine2\n"
 * Output:
 * "\n | MyLine1\n | MyLine2"
 * note that they get an additional \n at the beginning, and
 * doubled \n and \n at the end are removed, so you can use
 * ConfigShift like this
 * \code
 * strstr << 	"MySubcomponent = " << ConfigShift(comp1.config_string()) << "\n"
 * 				"MySubcomponent2 = " << ConfigShift(comp2.config_string()) << "\n"
 * \endcode
 * Depending on comp1.config_string(), this results in
 * "MySubcomponent1 = sub1 ... " or "MySubcomponent =\n | sub1.1\n | sub1.2 ..."
 * @param[in] s
 * @return shifted string
 */
UG_API std::string ConfigShift(std::string s);

template<typename T>
inline std::string OstreamShift(const T &t)
{
	std::stringstream ss; ss << t;
	return ConfigShift(ss.str());
}

/**
 * Helper function to display byte sizes like 2411724 => 2,3 MB
 * @param[in] s size in bytes
 * @param[in] length if != 0, fixes the returned string length to this length (for tables etc.).
 * @return string describing the size s=1024 -> 1 kb
 */
std::string GetBytesSizeString(size_t s, int length=0);

inline const char *TrueFalseString(bool b)
{
	return b ? "TRUE" : "FALSE";
}

inline const char *OnOffString(bool b)
{
	return b ? "ON" : "OFF";
}

} // end namespace ug

#endif /*__H__COMMON_STRING_UTIL__*/
