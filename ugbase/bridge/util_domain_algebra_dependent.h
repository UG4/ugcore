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


#ifndef UTIL_DOMAIN_ALGEBRA_DEPENDENT_H
#define	UTIL_DOMAIN_ALGEBRA_DEPENDENT_H

#include "util_algebra_dependent.h"
#include "util_domain_dependent.h"

namespace ug{
namespace bridge{

/// \addtogroup bridge
/// \{

template <	typename Functionality,
			typename TRegistry = Registry,
			typename DomainList = CompileDomainList,
			typename AlgebraList = CompileAlgebraList>
struct RegisterDomainAlgebraDependent
{
	RegisterDomainAlgebraDependent(TRegistry& reg, std::string grp)
	{
		static const bool domainIsEmpty = boost::mpl::empty<DomainList>::value;
		typename boost::mpl::if_c<domainIsEmpty, RegEnd, RegNextDomain>::type (reg,grp);
	}
	struct RegEnd{ RegEnd(TRegistry& reg, std::string grp){} };

	template <typename CurrAlgebraList>
	struct RegNextDomainAlgebra
	{
		RegNextDomainAlgebra(TRegistry& reg, std::string grp)
		{
			typedef typename boost::mpl::front<DomainList>::type DomainType;
			typedef typename boost::mpl::front<CurrAlgebraList>::type AlgebraType;
			typedef typename boost::mpl::pop_front<CurrAlgebraList>::type NextAlgebraList;

			Functionality::template DomainAlgebra<DomainType, AlgebraType, TRegistry>(reg,grp);
			RegAlgebra<NextAlgebraList>(reg,grp);
		}
	};

	template <typename CurrAlgebraList>
	struct RegAlgebra
	{
		RegAlgebra(TRegistry& reg, std::string grp)
		{
			static const bool algebraIsEmpty = boost::mpl::empty<CurrAlgebraList>::value;
			typename boost::mpl::if_c<algebraIsEmpty, RegEnd, RegNextDomainAlgebra<CurrAlgebraList> >::type (reg,grp);
		}
	};

	struct RegNextDomain
	{
		RegNextDomain(TRegistry& reg, std::string grp)
		{
			typedef typename boost::mpl::pop_front<DomainList>::type NextDomainList;

			RegAlgebra<AlgebraList>(reg,grp);
			RegisterDomainAlgebraDependent<Functionality, TRegistry, NextDomainList, AlgebraList>(reg,grp);
		}
	};
};


template <	typename Functionality, typename AlgebraList = CompileAlgebraList, typename TRegistry=Registry>
struct RegisterDomain1dAlgebraDependent
{
	RegisterDomain1dAlgebraDependent(TRegistry& reg, std::string grp)
	{
#ifdef UG_DIM_1
		RegisterDomainAlgebraDependent<Functionality, boost::mpl::list<Domain1d>, AlgebraList>(reg, grp);
#endif
	}
};

template <	typename Functionality, typename AlgebraList = CompileAlgebraList, typename TRegistry=Registry>
struct RegisterDomain2dAlgebraDependent
{
	RegisterDomain2dAlgebraDependent(TRegistry& reg, std::string grp)
	{
#ifdef UG_DIM_2
		RegisterDomainAlgebraDependent<Functionality, boost::mpl::list<Domain2d>, AlgebraList>(reg, grp);
#endif
	}
};


template <	typename Functionality,	typename AlgebraList = CompileAlgebraList, typename TRegistry=Registry>
struct RegisterDomain3dAlgebraDependent
{
	RegisterDomain3dAlgebraDependent(TRegistry& reg, std::string grp)
	{
#ifdef UG_DIM_3
		RegisterDomainAlgebraDependent<Functionality, boost::mpl::list<Domain3d>, AlgebraList>(reg, grp);
#endif
	}
};

template <	typename Functionality,	typename AlgebraList = CompileAlgebraList, typename TRegistry=Registry>
struct RegisterDomain2d3dAlgebraDependent
{
	RegisterDomain2d3dAlgebraDependent(TRegistry& reg, std::string grp)
	{
		RegisterDomain2dAlgebraDependent<Functionality, AlgebraList, TRegistry>(reg, grp);
		RegisterDomain3dAlgebraDependent<Functionality, AlgebraList, TRegistry>(reg, grp);
	}
};

// end group bridge
/// \}

} // namespace bridge
} // namespace ug
#endif	/* UTIL_DOMAIN_ALGEBRA_DEPENDENT_H */

