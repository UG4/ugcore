// author: Andreas Vogel

#include "string_util.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <fstream>
#include "common/common.h"
#include "common/assert.h"

namespace ug{

using namespace std;

void RemoveWhitespaceFromString(std::string& str)
{
	str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}


void TokenizeString(const string& str, vector<string>& vToken, const char delimiter)
{
	vToken.clear();
	stringstream tokenstream;
	tokenstream << str;
	string token;

	while ( getline (tokenstream, token, delimiter ) )
		vToken.push_back(token);
}

vector<string> TokenizeString(const string& str, const char delimiter)
{
	vector<string> vToken;
	TokenizeString(str, vToken, delimiter);
	return vToken;
}

vector<string> TokenizeString(const char* str, const char delimiter)
{
	vector<string> vToken;
	TokenizeString(string(str), vToken, delimiter);
	return vToken;
}

vector<string> TokenizeTrimString(const string& str, const char delimiter)
{
	vector<string> vToken;
	TokenizeString(str, vToken, delimiter);
	for(size_t i = 0; i < vToken.size(); ++i)
		vToken[i] = TrimString(vToken[i]);
	return vToken;
}

string TrimString(const string& str)
{
	const size_t start = str.find_first_not_of(" \t");
	const size_t end = str.find_last_not_of(" \t");
	if(start == string::npos || end == string::npos) return "";
	return str.substr(start, end - start + 1);
}

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

void AppendCounterToString( string& str, string indicator, int counter, 
                            int maxCounter )
{
//	check correct usage
	if(maxCounter >= 0)
		UG_ASSERT(counter <= maxCounter, "Wrong usage of function");

//	get number of digits needed for counter field (default is 4)
	int numDigits = 4;
	if(maxCounter >= 0) numDigits = NumberOfDigits(maxCounter);

//	create stringstream
	stringstream ss;

//	print stringstream
	ss << indicator << setfill ('0') << setw (numDigits) << counter;

//	append to string
	str.append(ss.str());
}

string AppendSpacesToString(string& str, int totalLength)
{
	int numSpaces = max((int)(totalLength-str.length()), 0);
	for(int i = 0; i < numSpaces; ++i) str.append(" ");

	return str;
}

string::size_type GetDirectorySeperatorPos(const string &str)
{
	string::size_type pos1 = str.rfind("/");
	string::size_type pos2 = str.rfind("\\");
	if(pos1 != string::npos)
	{
		if(pos2 != string::npos && pos2 > pos1) return pos2;
		else return pos1;
	}
	else return pos2;
}

string FilenameWithoutPath(const string &str)
{
	string::size_type pos = GetDirectorySeperatorPos(str);
	if( pos != string::npos ) return str.substr( pos+1 );
	else return str;
}

string PathFromFilename(const string &str)
{
	string::size_type pos = GetDirectorySeperatorPos(str);
	if( pos != string::npos ) return str.substr(0, pos+1 );
	else return ".";
}

string FilenameWithoutExtension(string str)
{
	str = FilenameWithoutPath(str);
	string::size_type pos = str.rfind(".");
	if( pos != string::npos ) return str.substr(0, pos );
	else return str;
}

string GetFilenameExtension(const string &str)
{
	string::size_type pos = str.rfind(".");
	if( pos != string::npos ) return str.substr(pos+1);
	else return "";
}

string ReplaceAll( string target, const string& oldstr, const string& newstr ) {
	// no substitution necessary
	if (oldstr == newstr) {
		return target;
	}

	for (size_t x = target.find(oldstr); x != string::npos; x = target.find(oldstr, x + newstr.size())) {
		target.erase(x, oldstr.length());
		target.insert(x, newstr);
	}

	return target;
}

bool StartsWith(const string& str, const string& begin) {
	return str.find(begin) == 0;
}

bool Contains(const string& str, const string& search) {
	return str.find(search) != string::npos;
}

//sreiter - Implementation is copied from some book or website. Can't remember...
template <> unsigned long hash_key(const string& key)
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

string ToLower(string str) {
    transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

string ToUpper(string str) {
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

vector<string> FindDuplicates(const vector<string>& vec) {
		
	vector<string> result;
	
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
			find(result.begin(), result.end(), vec[i])!=result.end();
		
		if (duplicateExists && !duplicateAddedToResult) {
			result.push_back(vec[i]);
		}
	}
	
	return result;	
}


const unsigned int cost_del = 1;
const unsigned int cost_ins = 1;
const unsigned int cost_sub = 1;
size_t LevenshteinDistance( const string& s1, const string& s2 )
{
  size_t n1 = s1.length();
  size_t n2 = s2.length();

  size_t* p = new size_t[ n2+1 ];
  size_t* q = new size_t[ n2+1 ];
  size_t* r;

  p[ 0 ] = 0;
  for( size_t j = 1; j <= n2; ++j )
    p[ j ] = p[ j-1 ] + cost_ins;

  for( size_t i = 1; i <= n1; ++i )
    {
      q[ 0 ] = p[ 0 ] + cost_del;
      for( size_t j = 1; j <= n2; ++j )
        {
          size_t d_del = p[ j   ] + cost_del;
          size_t d_ins = q[ j-1 ] + cost_ins;
          size_t d_sub = p[ j-1 ] + ( s1[i-1] == s2[j-1] ? 0 : cost_sub );
          q[ j ] = min( min( d_del, d_ins ), d_sub );
      }
      r = p;
      p = q;
      q = r;
    }

  size_t tmp = p[ n2 ];
  delete[] p;
  delete[] q;

  return tmp;
}

string repeat(char c, int nr)
{
	if(nr > 0)
		return string(nr, c);
	else
		return string("");
}

bool IsLonger(const string &a, const string &b)
{
	return b.size() > a.size();
}

string GetFileLines(const char *filename, size_t fromline, size_t toline, bool includeLineNumbers)
{
	if(filename[0] == '@') filename++;
	char buf[512];
	fstream file(filename, std::ios::in);
	if(file.is_open() == false) return string("");
	for(size_t i=0; i<fromline && !file.eof(); i++)
		file.getline(buf, 512);
	stringstream ss;
	if(includeLineNumbers)
		ss << fromline << "\t";
	ss << buf;
	for(; fromline < toline && !file.eof(); fromline++)
	{
		file.getline(buf, 512);
		ss << "\n";
		if(includeLineNumbers)
			ss << fromline+1 << "\t";
		ss << buf;
	}
	return ss.str();
}

string GetFileLine(const char *filename, size_t line)
{
	return GetFileLines(filename, line, line, false);
}

string XMLStringEscape(string s)
{
	s = ReplaceAll(s, "&", "&amp;");
	s = ReplaceAll(s, "\"", "&quot;");
	s = ReplaceAll(s, "\'", "&apos;");
	s = ReplaceAll(s, "<", "&lt;");
	s = ReplaceAll(s, ">", "&gt;");
	return s;
}

bool WildcardMatch(const char *str, const char *pattern)
{
	int strLen = strlen(str);
	int patternLen = strlen(pattern);
	int j=0;
	for(int i=0; i<patternLen; i++)
	{
		if(j >= strLen)
			return pattern[i] == '*';
		switch(pattern[i])
		{
		case '*':
			if(i == patternLen-1) return true;
			else
			{
				char nextSign = pattern[i+1];
				UG_ASSERT(nextSign != '*', "** is not a valid pattern");
				while(1)
				{
					for(;j<strLen; j++)
						if(str[j] == nextSign)
							break;
					if(j == strLen) return false;
					if(WildcardMatch(str+j, pattern+i+1)) return true;
					j++;
				}
			}
			break;

		case '?':
			j++;
			break;

		default:
			if(pattern[i] != str[j]) return false;
			j++;
			break;
		}
	}
	return true;
}


}

