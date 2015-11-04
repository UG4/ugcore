#ifndef UTIL_DOMAIN_ALGEBRA_DEPENDENT_H
#define	UTIL_DOMAIN_ALGEBRA_DEPENDENT_H

#include "util_algebra_dependent.h"
#include "util_domain_dependent.h"

namespace ug{
namespace bridge{

/// \addtogroup bridge
/// \{

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
			static const bool algebraIsEmpty = boost::mpl::empty<CurrAlgebraList>::value;
			typename boost::mpl::if_c<algebraIsEmpty, RegEnd, RegNextDomainAlgebra<CurrAlgebraList> >::type (reg,grp);
		}
	};

	struct RegNextDomain
	{
		RegNextDomain(Registry& reg, std::string grp)
		{
			typedef typename boost::mpl::pop_front<DomainList>::type NextDomainList;

			RegAlgebra<AlgebraList>(reg,grp);
			RegisterDomainAlgebraDependent<Functionality, NextDomainList, AlgebraList>(reg,grp);
		}
	};
};


template <	typename Functionality, typename AlgebraList = CompileAlgebraList>
struct RegisterDomain1dAlgebraDependent
{
	RegisterDomain1dAlgebraDependent(Registry& reg, std::string grp)
	{
#ifdef UG_DIM_1
		RegisterDomainAlgebraDependent<Functionality, boost::mpl::list<Domain1d>, AlgebraList>(reg, grp);
#endif
	}
};

template <	typename Functionality, typename AlgebraList = CompileAlgebraList>
struct RegisterDomain2dAlgebraDependent
{
	RegisterDomain2dAlgebraDependent(Registry& reg, std::string grp)
	{
#ifdef UG_DIM_2
		RegisterDomainAlgebraDependent<Functionality, boost::mpl::list<Domain2d>, AlgebraList>(reg, grp);
#endif
	}
};


template <	typename Functionality,	typename AlgebraList = CompileAlgebraList>
struct RegisterDomain3dAlgebraDependent
{
	RegisterDomain3dAlgebraDependent(Registry& reg, std::string grp)
	{
#ifdef UG_DIM_3
		RegisterDomainAlgebraDependent<Functionality, boost::mpl::list<Domain3d>, AlgebraList>(reg, grp);
#endif
	}
};

template <	typename Functionality,	typename AlgebraList = CompileAlgebraList>
struct RegisterDomain2d3dAlgebraDependent
{
	RegisterDomain2d3dAlgebraDependent(Registry& reg, std::string grp)
	{
		RegisterDomain2dAlgebraDependent<Functionality, AlgebraList>(reg, grp);
		RegisterDomain3dAlgebraDependent<Functionality, AlgebraList>(reg, grp);
	}
};

// end group bridge
/// \}

}
}
#endif	/* UTIL_DOMAIN_ALGEBRA_DEPENDENT_H */

