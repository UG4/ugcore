/*
 * lib_discretization_interface.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_interface.h"
#include "lib_discretization_interface.h"

#include "../ugbridge/registry.h"

namespace ug
{
namespace interface
{


void f1()
{
	return;
}


void RegisterLibDiscretizationInterface(Registry& reg)
{
/*	InterfaceRegistry regist;
	regist.add_function("test", &f1);
*/

	reg.register_object<MGDomainObject<2> >();
	reg.register_object<MGDomainObject<3> >();

	reg.register_object<P1ConformFunctionPatternObject>();
	reg.register_object<ApproximationSpaceObject<2, MartinAlgebra> >();
}

}//	end of namespace ug
}//	end of namespace interface
