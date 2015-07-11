/*
 * stringify.h
 *
 *  Created on: 12.03.2014
 *      Author: mrupp
 */

#ifndef STRINGIFY_H_
#define STRINGIFY_H_

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
 * You can use mkstr as a wrapper around Stringify like this:
 * \code
 * std::string s = mkstr("Hello " << 5 << " is a number");
 * \endcode
 *
 * \sa mkstr
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

///	Comfortable (but not necessarily efficient) string building
/** Based on Stringify, mkstr allows you to easily join strings and convert
 * strings to numbers like this:
 *
 * \code
 * std::string s = mkstr("Hello " << 5 << " is a number");
 * \endcode
 *
 * Note that mkstr is not highly efficient. It should only be used in code
 * which is not performance critical.*/
#define mkstr(s)	(Stringify() << s).str()


#endif /* STRINGIFY_H_ */
