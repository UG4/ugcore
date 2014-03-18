/*
 * stringify.h
 *
 *  Created on: 12.03.2014
 *      Author: mrupp
 */

#ifndef STRINGIFY_H_
#define STRINGIFY_H_

#include "smart_pointer.h"
#include <sstream>

class Stringify
{
public:
	SmartPtr<std::stringstream> ss;
	Stringify()
	{
		ss = make_sp(new std::stringstream);
	}
	template<typename T>
	Stringify &operator << (T t)
	{
		*ss << t;
		return *this;
	}

	std::string str() const
	{
		return ss->str();
	}

	operator std::string() const
	{
		return str();
	}
};

template<typename T>
inline T& operator << (T &t, const Stringify &s)
{
	t << s.str();
	return t;
}


#endif /* STRINGIFY_H_ */
