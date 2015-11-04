/*
 * string_concat.h
 *
 *  Created on: 22.05.2013
 *      Author: mrupp
 */

#ifndef STRING_CONCAT_H_
#define STRING_CONCAT_H_
#include <string>
#include "string_util.h"

namespace ug{
template<typename T>
inline std::string operator << (std::string a, T t)
{
	return a + ug::ToString(t);
}

template<>
inline std::string operator << (std::string a, const char *s)
{
	return a + s;
}

template<>
inline std::string operator << (std::string a, std::string b)
{
	return a + b;
}

}

#endif /* STRING_CONCAT_H_ */
