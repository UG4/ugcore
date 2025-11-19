/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

/*
 * simplified casting of overloaded functions/methods
 */

#ifndef UTIL_OVERLOADED_H_
#define UTIL_OVERLOADED_H_

/// for a function with signature int *MyFunction(char *name, OtherType &t) :
/// OVERLOADED_FUNCTION_PTR(int *, MyFunction, (char *name, OtherType &t ) )
#define OVERLOADED_FUNCTION_PTR(returnType, methodName, arguments) \
		static_cast<returnType  (*)arguments >(&methodName)

/// for a member function with signature bool MyClass::my_member_fu(double d) const :
/// OVERLOADED_CONST_METHOD_PTR(bool, MyClass, my_member_fu, (double d) )
#define OVERLOADED_CONST_METHOD_PTR(returnType, classType, methodName, arguments) \
		static_cast<returnType  (classType::*)arguments const>(&classType::methodName)

/// for a member function with signature void MyClass::my_member_fu() :
/// OVERLOADED_METHOD_PTR(void, MyClass, my_member_fu, () )
#define OVERLOADED_METHOD_PTR(returnType, classType, methodName, arguments) \
		static_cast<returnType  (classType::*)arguments>(&classType::methodName)

/// for adding of constructors with signatures.
/**
 * example:
 * \code
 * reg.add_class_<T, TBase>(name, grp)
 *			.add_constructor()
 *			. ADD_CONSTRUCTOR( (int parameter) )("parameter")
 *	\endcode
 *	note the double brackets ! ADD_CONSTRUCTOR ( (int a, int b) )
 */
#define ADD_CONSTRUCTOR(__theSignature) template add_constructor<void (*) __theSignature >

#endif
