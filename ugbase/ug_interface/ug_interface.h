//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#ifndef __H__UG__INTERFACE__
#define __H__UG__INTERFACE__

#include "ugbridge/registry.h"

namespace ug{
namespace interface
{

///	registers all standard interfaces.
/**	This method is called by the constructor of Registry automatically.
 *	You don't have to call it yourself!
 */
void RegisterStandardInterfaces(InterfaceRegistry& reg);

///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterLibGridInterface(InterfaceRegistry& reg);

///	registers lib-discretization interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterLibDiscretizationInterface(InterfaceRegistry& reg);

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterTestInterface(InterfaceRegistry& reg);

}//	end of namespace 
}//	end of namespace 

#endif
