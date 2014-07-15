/* 
 * File:   util_algebra_dependent.h
 *
 * Created on 2. November 2012, 10:48
 */

#include "registry/registry.h"

#include "lib_algebra/cpu_algebra_types.h"
#include "suffix_tag.h"

#include <boost/mpl/if.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>

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

typedef boost::mpl::list<
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

#ifdef UG_CPU_VAR
		CPUVariableBlockAlgebra,
#endif

	end_boost_list // see common/util/end_boost_list.h
> CompileAlgebraList;



template <typename Functionality, typename List = CompileAlgebraList>
struct RegisterAlgebraDependent
{
	RegisterAlgebraDependent(Registry& reg, std::string grp)
	{
		static const bool isEmpty = boost::mpl::empty<List>::value;
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
			typedef typename boost::mpl::front<List>::type AlgebraType;
			typedef typename boost::mpl::pop_front<List>::type NextList;
			Functionality::template Algebra<AlgebraType>(reg,grp);
			RegisterAlgebraDependent<Functionality, NextList>(reg,grp);
		}
	};
};

// end group bridge
/// \}

}
}
#endif	/* UTIL_ALGEBRA_DEPENDENT_H */

