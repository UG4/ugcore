/*
 * util_overloaded.h
 *
 *  Created on: 05.12.2013
 *      Author: mrupp
 *
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



#endif /* UTIL_OVERLOADED_H_ */
