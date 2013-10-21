/*
 * demangle.cpp
 *
 *  Created on: 26.09.2013
 *      Author: mrupp
 */
#ifdef UG_POSIX
#include <cxxabi.h>
#endif

#include <sstream>
#include <string>
#include <stdlib.h>

using namespace std;
namespace ug{

#ifdef UG_POSIX
string demangle_block(const char *str)
{
	stringstream ss;
	char lastc=0x00;
    string s;
	int status;
	for(char c = *str; c != 0x00; c = *(++str))
	{
		// mangled names start with _ . consider only if last sign was space or tab or newline
		if(c == '_' && (lastc==' ' || lastc == '\n' || lastc == '\t'))
		{
			s = "_";
			// add all characters to the string until space, tab or newline
			for(c = *(++str); c != 0x00; c = *(++str))
			{
				if(c == ' ' || c == '\n' || c == '\t')
					break;
				s += c;
			}
			// some compilers add an additional _ in front. skip it.
			const char *p = s.c_str();
			if(s.length() > 2 && s[1] == '_') p = p+1;
			char *realname = abi::__cxa_demangle(p, 0, 0, &status);
			if(status==0)
				ss << realname; // demangle successfull
			else
				ss << s; // demangling failed, print normal string
			free(realname);
		}
		ss << c;
		lastc =c;
	}
	return ss.str();
}

string demangle(const char *str)
{
	int status;
	char *realname = abi::__cxa_demangle(str, 0, 0, &status);
	if(status==0)
		return realname; // demangle successfull
	else
		return str; // demangling failed, print normal string
}
#else
string demangle_block(const char *str)
{
	return str;
}
string demangle(const char *str)
{
	return str;
}
#endif
}
