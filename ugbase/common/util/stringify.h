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
 * \note be sure that the object c_str is referring to is valid.
 * According to the C++ standard,
 * "Temporary objects are destroyed as the last step in evaluating the full-expression (1.9) that (lexically) contains the point where they were created. [12.2/3]"
 * so you can also do
 * CallConstCharFunction((Stringify() << "bla").c_str());
 * however,
 * const char *str = (Stringify() << "bla").c_str());
 * will NOT work.
 *
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


template<typename T>
inline T& operator << (T &t, const Stringify &s)
{
	t << s.str();
	return t;
}

#endif /* STRINGIFY_H_ */
