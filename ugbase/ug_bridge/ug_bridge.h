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

///	returns the default registry used in ug
Registry & GetUGRegistry();

///	registers all standard interfaces.
/**	This method is called by the constructor of Registry automatically.
 *	You don't have to call it yourself!
 */
bool RegisterStandardInterfaces(Registry& reg, const char* parentGroup = "/ug4");

///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterLibGridInterface(Registry& reg, const char* parentGroup = "/ug4");

///	Registers the domain object and related methods
bool RegisterDomainInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers lib-algebra interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterStaticLibAlgebraInterface(Registry& reg, const char* parentGroup = "/ug4");
bool RegisterDynamicLibAlgebraInterface(Registry& reg, int algebra_type, const char* parentGroup = "/ug4");

///	registers lib-discretization interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterStaticLibDiscretizationInterface(Registry& reg, const char* parentGroup = "/ug4");
bool RegisterDynamicLibDiscretizationInterface(Registry& reg, int algebra_type, const char* parentGroup = "/ug4");

/// registers user data
bool RegisterUserData(Registry& reg, const char* parentGroup = "/ug4");

///	registers user functions for the elder problem.
void RegisterElderUserFunctions(Registry& reg, const char* parentGroup);

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterTestInterface(Registry& reg, const char* parentGroup = "/ug4");

/// registers access to profiling functions at the registry
bool RegisterProfileFunctions(Registry &reg, const char* parentGroup = "/ug4");


}//	end of namespace 
}//	end of namespace 

#endif /* __H__UG_BRIDGE__UG_BRIDGE__ */
