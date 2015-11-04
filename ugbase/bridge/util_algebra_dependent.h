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


static const size_t NUM_ALGEBRA_TYPES = boost::mpl::size<CompileAlgebraList>::type::value - 1;


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
			virtual ~AlgebraIDBase() {};
		};

		template <typename TAlgebra>
		struct AlgebraID : public AlgebraIDBase {};

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
				static const bool isEmpty = boost::mpl::empty<List>::value;
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
					typedef typename boost::mpl::front<List>::type AlgebraType;
					typedef typename boost::mpl::pop_front<List>::type NextList;
					atidp.reg<AlgebraType>();
					(RegisterAlgebraIndices<NextList> (atidp));
				}
			};
		};

		// constructor
		AlgebraTypeIDProvider()
		{
			n = 0;
			RegisterAlgebraIndices<>(*this);
		}

		// prevent copy constructor and assignment
		AlgebraTypeIDProvider(AlgebraTypeIDProvider const&);	// do not implement
		void operator=(AlgebraTypeIDProvider const&);			// do not implement

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

