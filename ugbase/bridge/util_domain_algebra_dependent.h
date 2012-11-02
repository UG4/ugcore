/* 
 * File:   util_domain_algebra_dependent.h
 *
 * Created on 2. November 2012, 10:52
 */

#ifndef UTIL_DOMAIN_ALGEBRA_DEPENDENT_H
#define	UTIL_DOMAIN_ALGEBRA_DEPENDENT_H

#include "util_algebra_dependent.h"
#include "util_domain_dependent.h"

namespace ug{
namespace bridge{
	
	
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


}
}
#endif	/* UTIL_DOMAIN_ALGEBRA_DEPENDENT_H */

