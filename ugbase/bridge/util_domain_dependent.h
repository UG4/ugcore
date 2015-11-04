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

#endif	/* UTIL_DOMAIN_DEPENDENT_H */

