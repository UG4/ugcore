/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG_BRIDGE__UTIL__
#define __H__UG_BRIDGE__UTIL__

#include "registry/registry.h"
#include "suffix_tag.h"

namespace ug{
namespace bridge{

/// \addtogroup bridge
/// \{

template <typename Functionality>
void RegisterCommon(Registry& reg, std::string grp)
{
	Functionality::Common(reg,grp);
}

template <typename Functionality>
void RegisterDimensionDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_1
	Functionality::template Dimension<1>(reg,grp);
#endif
#ifdef UG_DIM_2
	Functionality::template Dimension<2>(reg,grp);
#endif
#ifdef UG_DIM_3
	Functionality::template Dimension<3>(reg,grp);
#endif
}

template <typename Functionality>
void RegisterDimension1dDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_1
	Functionality::template Dimension<1>(reg,grp);
#endif
}

template <typename Functionality>
void RegisterDimension2dDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_2
	Functionality::template Dimension<2>(reg,grp);
#endif
}

template <typename Functionality>
void RegisterDimension3dDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_3
	Functionality::template Dimension<3>(reg,grp);
#endif
}

template <typename Functionality>
void RegisterDimension2d3dDependent(Registry& reg, std::string grp)
{
	RegisterDimension2dDependent<Functionality>(reg, grp);
	RegisterDimension3dDependent<Functionality>(reg, grp);
}

// end group bridge
/// \}

} // end namespace bridge
} // end namespace ug

#define UG_REGISTRY_CATCH_THROW(grp)	\
		catch(ug::bridge::UGRegistryError& ex) {\
			UG_ERR_LOG("### ERROR while registering functionality at '"<<(grp)<<"'. "\
					"Registration failed (using name " << ex.name << ").\n");\
			throw(ex);}

#endif /* __H__UG_BRIDGE__UTIL__ */
