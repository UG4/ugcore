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

/**
 * use this class like this
 * std::string myTmpString = Stringify() << "My_" << iInteger << "_SuperString" << ".txt";
 * also in functions
 * CallStdStringFunction(Stringify() << "My_" << iInteger << "_SuperString" << ".txt");
 * Note that this can NOT work for const char * functions, since you need an
 * temporary object which is valid for the whol call of the function.
 * Then you'd do
 * std::string myTmpString = Stringify() << "My_" << iInteger << "_SuperString" << ".txt";
 * CallConstCharFunction(myTmpString.c_str());
 */
class Stringify
{
public:
	std::stringstream ss;
	Stringify()
	{

	}
	template<typename T>
	Stringify &operator << (T t)
	{
		ss << t;
		return *this;
	}

	std::string str() const
	{
		return ss.str();
	}

	operator std::string() const
	{
		return str();
	}
};


class ConstCharify
{
public:
	std::stringstream ss;
	ConstCharify()
	{

	}
	template<typename T>
	ConstCharify &operator << (T t)
	{
		ss << t;
		return *this;
	}

	std::string str() const
	{
		return ss.str();
	}

	operator const char*() const
	{
		return str().c_str();
	}
};

template<typename T>
inline T& operator << (T &t, const Stringify &s)
{
	t << s.str();
	return t;
}

template<typename T>
inline T& operator << (T &t, const ConstCharify &s)
{
	t << s.str();
	return t;
}


#endif /* STRINGIFY_H_ */
