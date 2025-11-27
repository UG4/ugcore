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


#include "registry/registry.h"

#include "lib_algebra/cpu_algebra_types.h"
#include "suffix_tag.h"

#include <boost/mpl/if.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>

#include "boost/mpl/size.hpp"

#include "common/util/end_boost_list.h"

#ifndef UTIL_ALGEBRA_DEPENDENT_H
#define	UTIL_ALGEBRA_DEPENDENT_H



namespace ug{
namespace bridge{

/// \addtogroup bridge
/// \{

////////////////////////////////////////////////////////////////////////////////
// 	Default Algebra List
////////////////////////////////////////////////////////////////////////////////

using CompileAlgebraList = boost::mpl::list<
#ifdef UG_GPU
		GPUAlgebra,
#endif

#ifdef UG_CPU_1
	CPUAlgebra,
#endif
		
#ifdef UG_CPU_2
	CPUBlockAlgebra<2>,
#endif
		
#ifdef UG_CPU_3
	CPUBlockAlgebra<3>,
#endif
		
#ifdef UG_CPU_4
		CPUBlockAlgebra<4>,
#endif

#ifdef UG_CPU_5
		CPUBlockAlgebra<5>,
#endif

#ifdef UG_CPU_6
		CPUBlockAlgebra<6>,
#endif

#ifdef UG_CPU_VAR
		CPUVariableBlockAlgebra,
#endif

	end_boost_list // see common/util/end_boost_list.h
>;


static constexpr size_t NUM_ALGEBRA_TYPES = boost::mpl::size<CompileAlgebraList>::type::value - 1;


struct AlgebraTypeIDProvider
{
	public:
		static AlgebraTypeIDProvider& instance()
		{
			static AlgebraTypeIDProvider instance;
			return instance;
		}

		// helper structs for storage of algebra type indices
		struct AlgebraIDBase
		{
			// needed - otherwise compiler will assume AlgebraIDBase is not polymorphic
			virtual ~AlgebraIDBase() = default;
		};

		template <typename TAlgebra>
		struct AlgebraID : AlgebraIDBase {};

		template <typename TAlgebra>
		void reg()
		{
			m_aid[n++] = (new AlgebraID<TAlgebra>());
		}

		template <typename TAlgebra>
		size_t id()
		{
			// maybe this can be done faster
			// however, unless number of algebra types grows significantly,
			// this should be fine
			for (size_t i = 0; i < NUM_ALGEBRA_TYPES; ++i)
			{
				AlgebraID<TAlgebra>* p_aid = dynamic_cast<AlgebraID<TAlgebra>*>(m_aid[i]);
				if (p_aid) return i;
			}
			UG_THROW("Cannot provide Algebra type index. Algebra type unknown.");
			return 0;
		}

	private:
		// helper struct for filling storage of indices via boost::mpl::list
		template <typename List = CompileAlgebraList>
		struct RegisterAlgebraIndices
		{
			RegisterAlgebraIndices(AlgebraTypeIDProvider& atidp)
			{
				static constexpr bool isEmpty = boost::mpl::empty<List>::value;
				(typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (atidp));
			}

			struct RegEnd
			{
				RegEnd(AlgebraTypeIDProvider&) {}
			};

			struct RegNext
			{
				RegNext(AlgebraTypeIDProvider& atidp)
				{
					using AlgebraType = typename boost::mpl::front<List>::type;
					using NextList = typename boost::mpl::pop_front<List>::type;
					atidp.reg<AlgebraType>();
					(RegisterAlgebraIndices<NextList> (atidp));
				}
			};
		};

		// constructor
		AlgebraTypeIDProvider()
		{
			n = 0;
			RegisterAlgebraIndices(*this);
		}

		// prevent copy constructor and assignment
		AlgebraTypeIDProvider(AlgebraTypeIDProvider const&);	// do not implement
		void operator = (AlgebraTypeIDProvider const&);			// do not implement

	private:
		// storage for algebra types (index in array corresponds to algebra type index)
		AlgebraIDBase* m_aid[NUM_ALGEBRA_TYPES];
		size_t n;
};


template <typename Functionality, typename List = CompileAlgebraList>
struct RegisterAlgebraDependent
{
	RegisterAlgebraDependent(Registry& reg, std::string grp)
	{
		static constexpr bool isEmpty = boost::mpl::empty<List>::value;
		typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (reg,grp);
	}

	struct RegEnd{
		RegEnd(Registry& reg, std::string grp)
		{
		}
	};
	struct RegNext
	{
		RegNext(Registry& reg, std::string grp)
		{
			using AlgebraType = typename boost::mpl::front<List>::type;
			using NextList = typename boost::mpl::pop_front<List>::type;
			Functionality::template Algebra<AlgebraType>(reg,grp);
			RegisterAlgebraDependent<Functionality, NextList>(reg,grp);
		}
	};
};

// end group bridge
/// \}

}
}
#endif
