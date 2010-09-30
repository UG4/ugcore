/*
 * lib_discretization_interface.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "lib_discretization/lib_discretization.h"

//#include "../ugbridge/registry.h"

namespace ug
{
namespace interface
{

template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename)
{
	const char * p = strstr(filename, ".ugx");
	if(p == NULL)
	{
		UG_LOG("Currently only '.ugx' format supported for domains.\n");
		return false;
	}

	return LoadGridFromFile(domain.get_grid(), filename, domain.get_subset_handler());
}

bool AddP1Function(P1ConformFunctionPattern& pattern, std::string name, int dim)
{
	return pattern.add_discrete_function(name, LSFS_LAGRANGEP1, dim);
}

void RegisterLibDiscretizationInterface(InterfaceRegistry& reg)
{

//	Domain2d
	{
	typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
	reg.add_class_<domain_type>("Domain2d")
		.add_constructor()
		.add_method("get_subset_handler", (MGSubsetHandler& (domain_type::*)()) &domain_type::get_subset_handler)
		.add_method("get_grid", (MultiGrid& (domain_type::*)()) &domain_type::get_grid);

	reg.add_function("LoadDomain2d", &LoadDomain<domain_type>);
	}

//	Domain3d
	{
	typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
	reg.add_class_<domain_type>("Domain3d")
		.add_constructor()
		.add_method("get_subset_handler", (MGSubsetHandler& (domain_type::*)()) &domain_type::get_subset_handler)
		.add_method("get_grid", (MultiGrid& (domain_type::*)()) &domain_type::get_grid);

	reg.add_function("LoadDomain3d", &LoadDomain<domain_type>);
	}

//	FunctionPattern (Abstract Base Class)
	reg.add_class_<FunctionPattern>("FunctionPattern");

//	P1ConformFunctionPattern
	{
	typedef P1ConformFunctionPattern T;
	reg.add_class_<T, FunctionPattern>("P1ConformFunctionPattern")
		.add_constructor()
		.add_method("lock", &T::lock);
	}

//  Add discrete function to pattern
	reg.add_function("AddP1Function", &AddP1Function);

//  ApproximationSpace2d
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef ApproximationSpace<domain_type, P1ConformDoFDistribution, MartinAlgebra> T;
		reg.add_class_<T>("ApproximationSpace2d")
			.add_constructor()
			.add_method("assign_domain", &T::assign_domain)
			.add_method("assign_function_pattern", &T::assign_function_pattern);
	}

//  ApproximationSpace3d
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		typedef ApproximationSpace<domain_type, P1ConformDoFDistribution, MartinAlgebra> T;
		reg.add_class_<T>("ApproximationSpace3d")
			.add_constructor()
			.add_method("assign_domain", &T::assign_domain)
			.add_method("assign_function_pattern", &T::assign_function_pattern);
	}

}


}//	end of namespace ug
}//	end of namespace interface
