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

typedef boost::mpl::list<
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
> CompileDomainList;

////////////////////////////////////////////////////////////////////////////////
//  Register invokers
////////////////////////////////////////////////////////////////////////////////

template <typename Functionality, typename List = CompileDomainList, typename TRegistry=Registry>
struct RegisterDomainDependent
{
	RegisterDomainDependent(TRegistry& reg, std::string grp)
	{
		static const bool isEmpty = boost::mpl::empty<List>::value;
		typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (reg,grp);
	}
	struct RegEnd{ RegEnd(TRegistry& reg, std::string grp){} };
	struct RegNext
	{
		RegNext(TRegistry& reg, std::string grp)
		{
			typedef typename boost::mpl::front<List>::type DomainType;
			typedef typename boost::mpl::pop_front<List>::type NextList;
			Functionality::template Domain<DomainType>(reg,grp);
			RegisterDomainDependent<Functionality, NextList, TRegistry>(reg,grp);
		}
	};
};

template<typename Functionality, typename TRegistry=Registry>
void RegisterDomain1dDependent(TRegistry& reg, std::string grp)
{
#ifdef UG_DIM_1
	RegisterDomainDependent<Functionality, boost::mpl::list<Domain1d> > (reg, grp);
#endif
}

template<typename Functionality, typename TRegistry=Registry>
void RegisterDomain2dDependent(TRegistry& reg, std::string grp)
{
#ifdef UG_DIM_2
	RegisterDomainDependent<Functionality, boost::mpl::list<Domain2d> > (reg, grp);
#endif
}

template<typename Functionality, typename TRegistry=Registry>
void RegisterDomain3dDependent(TRegistry& reg, std::string grp)
{
#ifdef UG_DIM_3
	RegisterDomainDependent<Functionality, boost::mpl::list<Domain3d> > (reg, grp);
#endif
}

template<typename Functionality, typename TRegistry=Registry>
void RegisterDomain2d3dDependent(TRegistry& reg, std::string grp)
{
	RegisterDomain2dDependent<Functionality, TRegistry>(reg, grp);
	RegisterDomain3dDependent<Functionality, TRegistry>(reg, grp);
}

// end group bridge
/// \}

} // namespace bridge

#ifdef UG_USE_PYBIND11
namespace pybind{

template <typename TFunctionality>
void RegisterDomainDependent(RegistryAdapter& reg, std::string grp)
{
	typedef typename ug::pybind::RegistryAdapter TRegistry;
	typedef typename ug::bridge::CompileDomainList TDomainList;
	ug::bridge::RegisterDomainDependent<TFunctionality, TDomainList, TRegistry>(reg, grp);
}



}
#endif

} // namespace ug

#endif	/* UTIL_DOMAIN_DEPENDENT_H */

