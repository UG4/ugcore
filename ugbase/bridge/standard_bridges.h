#ifndef __H__UG_BRIDGE__UG_STANDARD_BRIDGES__
#define __H__UG_BRIDGE__UG_STANDARD_BRIDGES__

#include "bridge.h"

namespace ug{
namespace bridge{

/**
 * \defgroup bridge Bridge
 * \ingroup ugbase
 * \{
 */

///	Registers types and functions for 1, 2, 3 and 4 dimensional vector math.
void RegisterBridge_VecMath(Registry& reg, std::string grp = UG4_GRP);

///	registers some util methods like path-access and script-parsing.
void RegisterBridge_Util(Registry& reg, std::string grp = UG4_GRP);

///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterBridge_Grid(Registry& reg, std::string grp = UG4_GRP);

///	registers methods for a parallel environment
void RegisterBridge_PCL(Registry& reg, std::string grp = UG4_GRP);

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
void RegisterBridge_Test(Registry& reg, std::string grp = UG4_GRP);

/// registers access to profiling functions at the registry
void RegisterBridge_Profiler(Registry &reg, std::string grp = UG4_GRP);

void RegisterBridge_Misc(Registry &reg, std::string grp = UG4_GRP);

///	Registers the domain object and related methods
void RegisterBridge_Domain(Registry& reg, std::string grp = UG4_GRP);

/// Registers periodic boundary identification
void RegisterBridge_PeriodicBoundary(Registry& reg, std::string grp = UG4_GRP);

///	Registers methods to perform ray-tracing on domains
void RegisterBridge_DomainRayTracing(Registry& reg, std::string grp = UG4_GRP);

///	Registers refiners and marking methods.
void RegisterBridge_Refinement(Registry& reg, std::string grp = UG4_GRP);

///	Registers methods to perform selections on the elements of a domain.
void RegisterBridge_Selection(Registry& reg, std::string grp = UG4_GRP);

///	Registers methods to transform the vertices of a domain.
void RegisterBridge_Transform(Registry& reg, std::string grp = UG4_GRP);

/// Registers the element discretizations
void RegisterBridge_ElemDiscs(Registry& reg, std::string grp = UG4_GRP);

///	Registers the common part of lib_discretization
void RegisterBridge_DiscCommon(Registry& reg, std::string grp = UG4_GRP);

/// registers user data
void RegisterBridge_UserData(Registry& reg, std::string grp = UG4_GRP);

///	registers LoadBalancer, partitioners, etc
void RegisterBridge_LoadBalancing(Registry& reg, std::string grp = UG4_GRP);


#ifdef UG_ALGEBRA
///	registers lib-algebra interface methods at the registry.
void RegisterBridge_AlgebraCommon(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Preconditioner(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Schur(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Obstacle(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_PILUT(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Solver(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Eigensolver(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_DomainDependentPreconditioner(Registry& reg, std::string grp = UG4_GRP);
///	registers restart functionality
void RegisterBridge_Restart(Registry& reg, std::string grp = UG4_GRP);

///	registers lib-discretization interface methods at the registry.
void RegisterBridge_DiscAlgebra(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Constraints(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_MultiGrid(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Output(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_AdaptiveTools(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_FiniteVolume(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Integrate(Registry& reg, std::string grp = UG4_GRP);

void RegisterBridge_DomainDisc(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_GridFunction(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Interpolate(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Evaluate(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_MaxError(Registry& reg, std::string grp = UG4_GRP);
void RegisterBridge_Ordering(Registry& reg, std::string grp = UG4_GRP);

void RegisterBridge_ManifoldUtil(Registry& reg, std::string grp = UG4_GRP);
#endif

// end group bridge
/// \}

}//	end bridge
}//	end ug

#endif /* __H__UG_BRIDGE__UG_STANDARD_BRIDGES__ */
