/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#ifndef __H__UG_BRIDGE__UG_STANDARD_BRIDGES__
#define __H__UG_BRIDGE__UG_STANDARD_BRIDGES__

#include "bridge.h"

namespace ug{


/**
 * \defgroup bridge Bridge
 * \ingroup ugbase
 * \{
 */
using bridge::UG4_GRP;

///	Registers types and functions for 1, 2, 3 and 4 dimensional vector math.
UG_REGISTRY_DECL(RegisterBridge_VecMath, Registry& reg, std::string grp = bridge::UG4_GRP);
//void RegisterBridge_VecMath(Registry& reg, std::string grp = UG4_GRP);

///	registers some util methods like path-access and script-parsing.
UG_REGISTRY_DECL(RegisterBridge_Util, Registry& reg, std::string grp = bridge::UG4_GRP);
// void RegisterBridge_Util(Registry& reg, std::string grp = UG4_GRP);


///	registers lib-grid interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
UG_REGISTRY_DECL(RegisterBridge_Grid, Registry& reg, std::string grp = bridge::UG4_GRP);

///	registers methods for a parallel environment
UG_REGISTRY_DECL(RegisterBridge_PCL, Registry& reg, std::string grp = bridge::UG4_GRP);

///	registers tests for the interface methods at the registry.
/**	This method is automatically invoked during the creation of the Registry.*/
UG_REGISTRY_DECL(RegisterBridge_Test, Registry& reg, std::string grp = bridge::UG4_GRP);

/// registers access to profiling functions at the registry
UG_REGISTRY_DECL(RegisterBridge_Profiler, Registry& reg, std::string grp = bridge::UG4_GRP);


UG_REGISTRY_DECL(RegisterBridge_Misc, Registry& reg, std::string grp = bridge::UG4_GRP);

///	Registers the domain object and related methods
UG_REGISTRY_DECL(RegisterBridge_Domain, Registry& reg, std::string grp = bridge::UG4_GRP);

/// Registers periodic boundary identification
UG_REGISTRY_DECL(RegisterBridge_PeriodicBoundary, Registry& reg, std::string grp = bridge::UG4_GRP);

///	Registers methods to perform ray-tracing on domains
UG_REGISTRY_DECL(RegisterBridge_DomainRayTracing, Registry& reg, std::string grp = bridge::UG4_GRP);

///	Registers refiners and marking methods.
UG_REGISTRY_DECL(RegisterBridge_Refinement, Registry& reg, std::string grp = bridge::UG4_GRP);

///	Registers methods to perform selections on the elements of a domain.
UG_REGISTRY_DECL(RegisterBridge_Selection, Registry& reg, std::string grp = bridge::UG4_GRP);

///	Registers methods to transform the vertices of a domain.
UG_REGISTRY_DECL(RegisterBridge_Transform, Registry& reg, std::string grp = bridge::UG4_GRP);

/// Registers the element discretizations
UG_REGISTRY_DECL(RegisterBridge_ElemDiscs, Registry& reg, std::string grp = bridge::UG4_GRP);

///	Registers the common part of lib_discretization
UG_REGISTRY_DECL(RegisterBridge_DiscCommon, Registry& reg, std::string grp = bridge::UG4_GRP);

/// registers user data
UG_REGISTRY_DECL(RegisterBridge_UserData, Registry& reg, std::string grp = bridge::UG4_GRP);

///	registers LoadBalancer, partitioners, etc
UG_REGISTRY_DECL(RegisterBridge_LoadBalancing, Registry& reg, std::string grp = bridge::UG4_GRP);

///	registers rasters, e.g. for 1,2,3 dimensional image data or density distributions
UG_REGISTRY_DECL(RegisterBridge_Raster, Registry& reg, std::string grp = bridge::UG4_GRP);

///	registers orthogonal polynomials
UG_REGISTRY_DECL(RegisterBridge_OrthoPoly, Registry& reg, std::string grp = bridge::UG4_GRP);

/// registers reference mapping test functionality (common)
UG_REGISTRY_DECL(RegisterBridge_ReferenceMappingTest, Registry& reg, std::string grp = bridge::UG4_GRP);

#ifdef UG_ALGEBRA
///	registers lib-algebra interface methods at the registry.
UG_REGISTRY_DECL(RegisterBridge_AlgebraCommon, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Preconditioner, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Schur, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Obstacle, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_PILUT, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_AlgebraOrdering, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Solver, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Eigensolver, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_DomainDependentPreconditioner, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_ConstrainedLinearIterator, Registry& reg, std::string grp = bridge::UG4_GRP);
///	registers restart functionality
UG_REGISTRY_DECL(RegisterBridge_Restart, Registry& reg, std::string grp = bridge::UG4_GRP);

///	registers lib-discretization interface methods at the registry.
UG_REGISTRY_DECL(RegisterBridge_DiscAlgebra, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Constraints, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_MultiGrid, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Output, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_AdaptiveTools, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_FiniteVolume, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Integrate, Registry& reg, std::string grp = bridge::UG4_GRP);

UG_REGISTRY_DECL(RegisterBridge_DomainDisc, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_GridFunction, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Interpolate, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Evaluate, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_MaxError, Registry& reg, std::string grp = bridge::UG4_GRP);
UG_REGISTRY_DECL(RegisterBridge_Ordering, Registry& reg, std::string grp = bridge::UG4_GRP);

UG_REGISTRY_DECL(RegisterBridge_ManifoldUtil, Registry& reg, std::string grp = bridge::UG4_GRP);
#endif

// end group bridge
/// \}

// }//	end bridge
}//	end ug

#endif /* __H__UG_BRIDGE__UG_STANDARD_BRIDGES__ */
