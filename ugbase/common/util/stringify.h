/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
	Stringify()	= default;

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


#endif
