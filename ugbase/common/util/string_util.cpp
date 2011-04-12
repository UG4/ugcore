// author: Andreas Vogel

#include "string_util.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <cctype>

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

//sreiter - Implementation is copied from some book or website. Can't remember...
template <> unsigned long hash_key(const char* key)
{

	unsigned long hash = 5381;
	int c;

	while((c = *key)){
		hash = hash * 33 + c;
		++key;
	}

	return hash;
}

}

