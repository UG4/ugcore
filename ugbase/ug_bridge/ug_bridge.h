//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#ifndef __H__UG_BRIDGE__UG_BRIDGE__
#define __H__UG_BRIDGE__UG_BRIDGE__

#include "registry.h"

namespace ug
{
namespace bridge
{

///	registers all standard interfaces.
/**	This method is called by the constructor of Registry automatically.
 *	You don't have to call it yourself!
 */
void RegisterStandardInterfaces(Registry& reg, const char* parentGroup = "/ug4");

///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterLibGridInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers lib-algebra interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterLibAlgebraInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers lib-discretization interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterLibDiscretizationInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterTestInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers info commands TypeInfo and ClassUsage
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterInfoCommands(Registry &reg, const char* parentGroup = "/ug4");

}//	end of namespace 
}//	end of namespace 

#endif /* __H__UG_BRIDGE__UG_BRIDGE__ */
