/*
 * lib_algebra_bridge.h
 *
 *  Created on: 23.03.2011
 *      Author: mrupp
 */

#ifndef UG_BRIDGE_LIB_ALGEBRA_BRIDGE_H_
#define UG_BRIDGE_LIB_ALGEBRA_BRIDGE_H_


#include "bridge/bridge.h"
#include "lib_algebra/algebra_selector.h"
#include "lib_algebra/algebra_types.h"

namespace ug
{

namespace bridge
{


/**
 * RegisterAlgebraClass function
 * This template function is intended to simplify the registration of different algebras.
 * If you want to add an algebra type, you only need to add it here.
 *
 * If you want to register classes or functions, you can to that
 * - in lib_algebra_bridge.cpp's cRegisterAlgebraType::reg function,
 * - or you do it in another .cpp file.
 * There you need a class like CRegisterMyClass
 * with a function static bool reg(Registry &reg, const char *parentGroup)
 * Then you write a global function
 \code
  bool RegisterMyClass(Registry& reg, int algebra_type, const char* parentGroup)
  {
		return RegisterAlgebraClass<CRegisterMyClass>(reg, algebra_type, parentGroup);
  }
  \endcode
 * and add it to lib_algebra_bridge.cpp's RegisterDynamicLibAlgebraInterface.
 *
 **/
// add your RegisterMyClass function in lib_algebra_bridge.cpp.

template<template<typename> class TRegister>
bool RegisterAlgebraClass(Registry& reg, std::string parentGroup)
{
	try
	{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Algebra";
		std::string grp = groupString.str();

	//	register
		TRegister<CPUAlgebra >::reg(reg, grp);
//		TRegister<CPUBlockAlgebra<2> >::reg(reg, grp);
		TRegister<CPUBlockAlgebra<3> >::reg(reg, grp);
//		TRegister<CPUBlockAlgebra<4> >::reg(reg, grp);
//		TRegister<CPUVariableBlockAlgebra>::reg(reg, grp);
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibAlgebraInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}
	return true;
}

}

}
#endif /* UG_BRIDGE_LIB_ALGEBRA_BRIDGE_H_ */
