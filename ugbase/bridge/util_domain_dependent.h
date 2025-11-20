/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef UTIL_DOMAIN_DEPENDENT_H
#define	UTIL_DOMAIN_DEPENDENT_H



#include "registry/registry.h"
#include "lib_disc/domain.h"

#include "suffix_tag.h"

#include <boost/mpl/if.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>


namespace ug{
namespace bridge{

/// \addtogroup bridge
/// \{

////////////////////////////////////////////////////////////////////////////////
// 	Default Domain List
////////////////////////////////////////////////////////////////////////////////

using CompileDomainList = boost::mpl::list<
#ifdef UG_DIM_1
	Domain1d
#endif
#if defined UG_DIM_1 && (defined UG_DIM_2 || defined UG_DIM_3)
	,
#endif
#ifdef UG_DIM_2
	Domain2d
#endif
#if defined UG_DIM_2 && defined UG_DIM_3
	,
#endif
#ifdef UG_DIM_3
	Domain3d
#endif
>;

////////////////////////////////////////////////////////////////////////////////
//  Register invokers
////////////////////////////////////////////////////////////////////////////////

template <typename Functionality, typename List = CompileDomainList>
struct RegisterDomainDependent
{
	RegisterDomainDependent(Registry& reg, std::string grp)
	{
		static constexpr bool isEmpty = boost::mpl::empty<List>::value;
		typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (reg,grp);
	}
	struct RegEnd{ RegEnd(Registry& reg, std::string grp){} };
	struct RegNext
	{
		RegNext(Registry& reg, std::string grp)
		{
			using DomainType = typename boost::mpl::front<List>::type;
			using NextList = typename boost::mpl::pop_front<List>::type;
			Functionality::template Domain<DomainType>(reg,grp);
			RegisterDomainDependent<Functionality, NextList>(reg,grp);
		}
	};
};

template<typename Functionality>
void RegisterDomain1dDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_1
	RegisterDomainDependent<Functionality, boost::mpl::list<Domain1d> > (reg, grp);
#endif
}

template<typename Functionality>
void RegisterDomain2dDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_2
	RegisterDomainDependent<Functionality, boost::mpl::list<Domain2d> > (reg, grp);
#endif
}

template<typename Functionality>
void RegisterDomain3dDependent(Registry& reg, std::string grp)
{
#ifdef UG_DIM_3
	RegisterDomainDependent<Functionality, boost::mpl::list<Domain3d> > (reg, grp);
#endif
}

template<typename Functionality>
void RegisterDomain2d3dDependent(Registry& reg, std::string grp)
{
	RegisterDomain2dDependent<Functionality>(reg, grp);
	RegisterDomain3dDependent<Functionality>(reg, grp);
}

// end group bridge
/// \}

} // namespace bridge
} // namespace ug

#endif
