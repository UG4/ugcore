/*
 * lib_algebra_bridge.h
 *
 *  Created on: 23.03.2011
 *      Author: mrupp
 */

#ifndef UG_BRIDGE_LIB_ALGEBRA_BRIDGE_H_
#define UG_BRIDGE_LIB_ALGEBRA_BRIDGE_H_


#include "ug_bridge/ug_bridge.h"
#include "lib_algebra/algebra_chooser.h"

namespace ug
{

namespace bridge
{

template<template<typename> class TRegister>
bool RegisterAlgebraClass(Registry& reg, int algebra_type, const char* parentGroup)
{
	try
	{
		//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Algebra";
		std::string grp = groupString.str();

		// register algebra
		switch(algebra_type)
		{
		case eCPUAlgebra:		 		TRegister<CPUAlgebra >::reg(reg, grp.c_str()); break;
//		case eCPUBlockAlgebra2x2: 		TRegister<CPUBlockAlgebra<2> >::reg(reg, grp.c_str()); break;
//		case eCPUBlockAlgebra3x3: 		TRegister<CPUBlockAlgebra<3> >::reg(reg, grp.c_str()); break;
//		case eCPUBlockAlgebra4x4: 		TRegister<CPUBlockAlgebra<4> >::reg(reg, grp.c_str()); break;
//		case eCPUVariableBlockAlgebra: 	TRegister<CPUVariableBlockAlgebra>::reg(reg, grp.c_str()); break;
		default: UG_ASSERT(0, "In RegisterAlgebraClass: " << algebra_type << " is unsupported algebra type");
					UG_LOG("In RegisterAlgebraClass: " << algebra_type << " is unsupported algebra type");
					return false;
		}
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
