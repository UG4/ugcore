//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#ifndef __H__UG_BRIDGE__UG_BRIDGE__
#define __H__UG_BRIDGE__UG_BRIDGE__

#include <string>
#include <sstream>
#include "registry/registry.h"
#include "lib_algebra/algebra_type.h"
#include "common/ug_config.h"

namespace ug{
namespace bridge{

///	returns the default registry used in ug
UG_API Registry & GetUGRegistry();

///	Sets the default classes of class-groups based on a tags using default DoFManager
UG_API void InitUG(int dim, const AlgebraType& algebraType);

/// calls RegisterStandardInterfaces and LoadPlugins if UG_PLUGINS is defined
UG_API bool InitBridge();

///	registers all standard interfaces.
/**	This method is called by the constructor of Registry automatically.
 *	You don't have to call it yourself!
 */
UG_API bool RegisterStandardInterfaces(Registry& reg, std::string parentGroup = "/ug4");

///	Registers types and functions for 1, 2, 3 and 4 dimensional vector math.
bool RegisterVecMathBridge(Registry& reg, std::string parentGroup = "/ug4");

///	registers some util methods like path-access and script-parsing.
bool RegisterUtilInterface(Registry& reg, std::string parentGroup = "/ug4");

///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterLibGridInterface(Registry& reg, std::string parentGroup = "/ug4");

///	registers methods for a parallel environment
bool RegisterPCLInterface(Registry& reg, std::string parentGroup = "/ug4");

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
bool RegisterTestInterface(Registry& reg, std::string parentGroup = "/ug4");

/// registers access to profiling functions at the registry
bool RegisterProfileFunctions(Registry &reg, std::string parentGroup = "/ug4");

bool RegisterMiscFunctions(Registry &reg, std::string parentGroup = "/ug4");

///	Registers the domain object and related methods
bool RegisterDomainInterface(Registry& reg, std::string parentGroup = "/ug4");

///	Registers refiners and marking methods.
bool RegisterRefinementBridge(Registry& reg, std::string parentGroup = "/ug4");

/// Registers the element discretizations
bool RegisterElemDiscs(Registry& reg, std::string parentGroup = "/ug4");

///	Registers the common part of lib_discretization
bool RegisterLibDisc_Common(Registry& reg, std::string parentGroup = "/ug4");

/// registers user data
bool RegisterUserData(Registry& reg, std::string parentGroup = "/ug4");

#ifdef UG_ALGEBRA
///	registers lib-algebra interface methods at the registry.
bool RegisterLibAlgebra(Registry& reg, std::string parentGroup = "/ug4");

///	registers lib-discretization interface methods at the registry.
bool RegisterLibDisc_Algebra(Registry& reg, std::string parentGroup = "/ug4");
bool RegisterLibDisc_Domain(Registry& reg, std::string parentGroup = "/ug4");
bool RegisterConstraints(Registry& reg, std::string parentGroup = "/ug4");
bool RegisterMultiGrid(Registry& reg, std::string parentGroup = "/ug4");
bool RegisterOutput(Registry& reg, std::string parentGroup = "/ug4");
bool RegisterAdaptiveTools(Registry& reg, std::string parentGroup = "/ug4");
bool RegisterIntegrate(Registry& reg, std::string parentGroup = "/ug4");
#endif

}//	end bridge
}//	end ug

#include "suffix_tag.h"

#endif /* __H__UG_BRIDGE__UG_BRIDGE__ */
