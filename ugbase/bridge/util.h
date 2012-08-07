/*
 * util.h
 *
 *  Created on: 24.05.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG_BRIDGE__UTIL__
#define __H__UG_BRIDGE__UTIL__

#include "registry/registry.h"
#include "lib_disc/domain.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "suffix_tag.h"

#include <boost/mpl/if.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>

namespace ug{
namespace bridge{

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
// 	Default Algebra List
////////////////////////////////////////////////////////////////////////////////

typedef boost::mpl::list<
#ifdef UG_CPU_1
		CPUAlgebra
#endif
#if defined UG_CPU_1 && (defined UG_CPU_2 || defined UG_CPU_3 || defined UG_CPU_4 || defined UG_CPU_VAR)
		,
#endif
#ifdef UG_CPU_2
		CPUBlockAlgebra<2>
#endif
#if defined UG_CPU_2 && (defined UG_CPU_3 || defined UG_CPU_4 || defined UG_CPU_VAR)
		,
#endif
#ifdef UG_CPU_3
		CPUBlockAlgebra<3>
#endif
#if defined UG_CPU_3 && (defined UG_CPU_4 || defined UG_CPU_VAR)
		,
#endif
#ifdef UG_CPU_4
		CPUBlockAlgebra<4>
#endif
#if defined UG_CPU_4 && defined UG_CPU_VAR
		,
#endif
#ifdef UG_CPU_VAR
		CPUVariableBlockAlgebra
#endif
> CompileAlgebraList;

////////////////////////////////////////////////////////////////////////////////
//  Register invokers
////////////////////////////////////////////////////////////////////////////////

template <typename Functionality, typename List = CompileDomainList>
struct RegisterDomainDependent
{
	RegisterDomainDependent(Registry& reg, std::string grp)
	{
		static const bool isEmpty = boost::mpl::empty<List>::value;
		typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (reg,grp);
	}
	struct RegEnd{ RegEnd(Registry& reg, std::string grp){} };
	struct RegNext
	{
		RegNext(Registry& reg, std::string grp)
		{
			typedef typename boost::mpl::front<List>::type DomainType;
			typedef typename boost::mpl::pop_front<List>::type NextList;
			Functionality::template Domain<DomainType>(reg,grp);
			RegisterDomainDependent<Functionality, NextList>(reg,grp);
		}
	};
};

template <typename Functionality, typename List = CompileAlgebraList>
struct RegisterAlgebraDependent
{
	RegisterAlgebraDependent(Registry& reg, std::string grp)
	{
		static const bool isEmpty = boost::mpl::empty<List>::value;
		typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (reg,grp);
	}

	struct RegEnd{ RegEnd(Registry& reg, std::string grp){} };
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

template <	typename Functionality,
			typename DomainList = CompileDomainList,
			typename AlgebraList = CompileAlgebraList>
struct RegisterDomainAlgebraDependent
{
	RegisterDomainAlgebraDependent(Registry& reg, std::string grp)
	{
		static const bool domainIsEmpty = boost::mpl::empty<DomainList>::value;
		typename boost::mpl::if_c<domainIsEmpty, RegEnd, RegNextDomain>::type (reg,grp);
	}
	struct RegEnd{ RegEnd(Registry& reg, std::string grp){} };

	template <typename CurrAlgebraList>
	struct RegNextDomainAlgebra
	{
		RegNextDomainAlgebra(Registry& reg, std::string grp)
		{
			typedef typename boost::mpl::front<DomainList>::type DomainType;
			typedef typename boost::mpl::front<CurrAlgebraList>::type AlgebraType;
			typedef typename boost::mpl::pop_front<CurrAlgebraList>::type NextAlgebraList;

			Functionality::template DomainAlgebra<DomainType, AlgebraType>(reg,grp);
			RegAlgebra<NextAlgebraList>(reg,grp);
		}
	};

	template <typename CurrAlgebraList>
	struct RegAlgebra
	{
		RegAlgebra(Registry& reg, std::string grp)
		{
			typedef typename boost::mpl::front<DomainList>::type DomainType;

			static const bool algebraIsEmpty = boost::mpl::empty<CurrAlgebraList>::value;
			typename boost::mpl::if_c<algebraIsEmpty, RegEnd, RegNextDomainAlgebra<CurrAlgebraList> >::type (reg,grp);
		}
	};

	struct RegNextDomain
	{
		RegNextDomain(Registry& reg, std::string grp)
		{
			typedef typename boost::mpl::front<DomainList>::type DomainType;
			typedef typename boost::mpl::pop_front<DomainList>::type NextDomainList;

			RegAlgebra<AlgebraList>(reg,grp);
			RegisterDomainAlgebraDependent<Functionality, NextDomainList, AlgebraList>(reg,grp);
		}
	};
};

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


} // end namespace bridge
} // end namespace ug

#define UG_REGISTRY_CATCH_THROW(grp)	\
		catch(UGRegistryError& ex) {\
			UG_ERR_LOG("### ERROR while registering functionality at '"<<(grp)<<"'. "\
					"Registration failed (using name " << ex.name << ").\n");\
			throw(ex);}

#endif /* __H__UG_BRIDGE__UTIL__ */
