/*
 * lib_discretization_interface.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_interface.h"
#include "lib_discretization_interface.h"

namespace ug
{
namespace interface
{


void RegisterLibDiscretizationInterface(Registry& reg)
{
	reg.register_object<MGDomainObject<2> >();
	reg.register_object<MGDomainObject<3> >();
}

}//	end of namespace ug
}//	end of namespace interface
