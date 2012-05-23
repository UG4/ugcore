// author: Andreas Vogel

#include "string_util.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cctype>
#include "common/common.h"

namespace ug{

void RemoveWhitespaceFromString(std::string& string)
{
	string.erase(std::remove_if(string.begin(), string.end(), isspace), string.end());
}

// help function to tokenize the parameter string
void TokenizeString(const std::string& str, std::vector<std::string>& tokens, const char delimiter)
{
	tokens.clear();
	std::stringstream tokenstream;
	tokenstream << str;
	std::string token;

	while ( std::getline (tokenstream, token, delimiter ) )
	{
			tokens.push_back(token);
	}
}

std::string TrimString(const std::string& str)
{
	const size_t start = str.find_first_not_of(" \t");
	const size_t end = str.find_last_not_of(" \t");
	if(start == std::string::npos || end == std::string::npos) return "";
	return str.substr(start, end - start + 1);
}

/// returns the number of digits of an integer (expressed with base 10)
int NumberOfDigits(int n)
{
//	a 0 has 1 digit
	if (n == 0) return 1;

//	divide by 10 until
    int cnt = 0;
    while (n > 0)
    {
    	++cnt;
    	n /= 10;
    }

//	return number of digits
    return cnt;
}

///	appends a counter number to a string
void AppendCounterToString(std::string& str, std::string indicator,
                           int counter, int maxCounter)
{
//	check correct usage
	if(maxCounter >= 0)
		UG_ASSERT(counter <= maxCounter, "Wrong usage of function");

//	get number of digits needed for counter field (default is 4)
	int numDigits = 4;
	if(maxCounter >= 0) numDigits = NumberOfDigits(maxCounter);

//	create stringstream
	std::stringstream ss;

//	print stringstream
	ss << indicator << std::setfill ('0') << std::setw (numDigits) << counter;

//	append to string
	str.append(ss.str());
}

///	appends a number of spaces to a string, returns "padded" string
std::string AppendSpacesToString(std::string& str, int totalLength)
{

	int numSpaces = std::max((int)(totalLength-str.length()), 0);
	for(int i = 0; i < numSpaces; ++i) str.append(" ");

	return str;
}

std::string::size_type GetDirectorySeperatorPos(const std::string &str)
{
	std::string::size_type pos1 = str.rfind("/");
	std::string::size_type pos2 = str.rfind("\\");
	if(pos1 != std::string::npos)
	{
		// FIXME: pos2 > pos2
		if(pos2 != std::string::npos && pos2 > pos2)	return pos2;
		else return pos1;
	}
	else return pos2;
}

std::string FilenameWithoutPath(const std::string& str)
{
	std::string::size_type pos = GetDirectorySeperatorPos(str);
	if( pos != std::string::npos ) return str.substr( pos+1 );
	else return str;
}

std::string PathFromFilename(const std::string &str)
{
	std::string::size_type pos = GetDirectorySeperatorPos(str);
	if( pos != std::string::npos ) return str.substr(0, pos+1 );
	else return ".";
}

std::string FilenameWithoutExtension(std::string str)
{
	str = FilenameWithoutPath(str);
	std::string::size_type pos = str.rfind(".");
	if( pos != std::string::npos ) return str.substr(0, pos );
	else return str;
}

std::string GetFilenameExtension(const std::string &str)
{
	std::string::size_type pos = str.rfind(".");
	if( pos != std::string::npos ) return str.substr(pos+1);
	else return "";
}

std::string ReplaceAll(
		std::string target,
		const std::string& oldstr,
		const std::string& newstr) {

	// no substitution necessary
	if (oldstr == newstr) {
		return target;
	}

	for (size_t x = target.find(oldstr); x != std::string::npos; x = target.find(oldstr, x + newstr.size())) {
		target.erase(x, oldstr.length());
		target.insert(x, newstr);
	}

	return target;
}

bool StartsWith(const std::string& str, const std::string& begin) {
	return str.find(begin) == 0;
}

bool Contains(const std::string& str, const std::string& search) {
	return str.find(search) !=std::string::npos;
}

//sreiter - Implementation is copied from some book or website. Can't remember...
template <> unsigned long hash_key(const std::string& key)
{

	unsigned long hash = 5381;
	const char* str = key.c_str();
	int c;

	while((c = *str)){
		hash = hash * 33 + c;
		++str;
	}

	return hash;
}


std::string ToLower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

std::string ToUpper(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

std::vector<std::string> FindDuplicates(const std::vector<std::string>& vec) {
		
	std::vector<std::string> result;
	
	// search for duplicates
	for (size_t i = 0; i < vec.size();i++) {
		
		bool duplicateExists = false;
		
		for (size_t j = 0; j < vec.size();j++) {
			if (vec[i]==vec[j] && i!=j) {
				duplicateExists = true;
				break;
			}
		}
		
		// if not already added, add entry to result vec
		bool duplicateAddedToResult =
			std::find(result.begin(), result.end(), vec[i])!=result.end();
		
		if (duplicateExists && !duplicateAddedToResult) {
			result.push_back(vec[i]);
		}
	}
	
	return result;	
}


const unsigned int cost_del = 1;
const unsigned int cost_ins = 1;
const unsigned int cost_sub = 1;
/**
 * Levenshtein distance algorithm, taken from
 * http://www.freemedialibrary.com/index.php/Levenshtein_distance
 * (check copyright or recreate!)
 *
 */
unsigned int LevenshteinDistance( const std::string& s1, const std::string& s2 )
{
  unsigned int n1 = s1.length();
  unsigned int n2 = s2.length();

  unsigned int* p = new unsigned int[ n2+1 ];
  unsigned int* q = new unsigned int[ n2+1 ];
  unsigned int* r;

  p[ 0 ] = 0;
  for( unsigned int j = 1; j <= n2; ++j )
    p[ j ] = p[ j-1 ] + cost_ins;

  for( unsigned int i = 1; i <= n1; ++i )
    {
      q[ 0 ] = p[ 0 ] + cost_del;
      for( unsigned int j = 1; j <= n2; ++j )
        {
          unsigned int d_del = p[ j   ] + cost_del;
          unsigned int d_ins = q[ j-1 ] + cost_ins;
          unsigned int d_sub = p[ j-1 ] + ( s1[i-1] == s2[j-1] ? 0 : cost_sub );
          q[ j ] = std::min( std::min( d_del, d_ins ), d_sub );
      }
      r = p;
      p = q;
      q = r;
    }

  unsigned int tmp = p[ n2 ];
  delete[] p;
  delete[] q;

  return tmp;
}

std::string repeat(char c, int nr)
{
	if(nr > 0)
		return std::string(nr, c);
	else
		return std::string("");
}
}

