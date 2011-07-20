//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#ifndef __H__UG_BRIDGE__UG_BRIDGE__
#define __H__UG_BRIDGE__UG_BRIDGE__

#include <string>
#include <sstream>
#include "registry/registry.h"

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

///	Registers types and functions for 1, 2, 3 and 4 dimensional vector math.
bool RegisterVecMathBridge(Registry& reg, const char* parentGroup = "/ug4");

///	registers some util methods like path-access and script-parsing.
bool RegisterUtilInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterLibGridInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers methods for a parallel environment
bool RegisterPCLInterface(Registry& reg, const char* parentGroup = "/ug4");

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterTestInterface(Registry& reg, const char* parentGroup = "/ug4");

/// registers access to profiling functions at the registry
bool RegisterProfileFunctions(Registry &reg, const char* parentGroup = "/ug4");

bool RegisterMiscFunctions(Registry &reg, const char* parentGroup = "/ug4");

///	Registers the domain object and related methods
bool RegisterDomainInterface(Registry& reg, const char* parentGroup = "/ug4");

///	Registers refiners and marking methods.
bool RegisterRefinementBridge(Registry& reg, const char* parentGroup = "/ug4");

/// Registers the element discretizations
bool RegisterLibDiscElemDisc(Registry& reg, const char* parentGroup = "/ug4");


#ifdef UG_ALGEBRA
///	registers lib-algebra interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterStaticLibAlgebraInterface(Registry& reg, const char* parentGroup = "/ug4");
bool RegisterDynamicLibAlgebraInterface(Registry& reg, int algebra_type, const char* parentGroup = "/ug4");

///	registers lib-discretization interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterStaticLibDiscInterface(Registry& reg, const char* parentGroup = "/ug4");
bool RegisterDynamicLibDiscDomain(Registry& reg, int algebra_type, const char* parentGroup = "/ug4");
bool RegisterDynamicLibDiscAlgebra(Registry& reg, int algebra_type, const char* parentGroup = "/ug4");

bool RegisterDynamicLibDiscretizationInterface(Registry& reg, int algebra_type, const char* parentGroup = "/ug4");

/// registers user data
bool RegisterUserData(Registry& reg, const char* parentGroup = "/ug4");

///	registers user functions for the elder problem.
void RegisterElderUserFunctions(Registry& reg, const char* parentGroup);
#endif

////////////////////////////////////////////////////////////////////////////////
//	Suffix and Tag - Section
////////////////////////////////////////////////////////////////////////////////

/// returns the dim-suffix for a domain (e.g. "3d")
template <typename TDomain>
std::string GetDomainSuffix()
{
//	extract dimension of domain
	static const int dim = TDomain::dim;

//	the dimension suffix
	std::stringstream ss; ss << dim << "d";

//	return the suffix
	return ss.str();
}

/// returns the dim-tag for a domain (e.g. "dim=3d")
template <typename TDomain>
std::string GetDomainTag()
{
//	the dimension suffix
	std::stringstream ss; ss << "dim="<< GetDomainSuffix<TDomain>();

//	return the suffix
	return ss.str();
}

}//	end of namespace 
}//	end of namespace 

#endif /* __H__UG_BRIDGE__UG_BRIDGE__ */
