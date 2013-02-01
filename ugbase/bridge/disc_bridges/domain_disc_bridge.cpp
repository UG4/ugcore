/*
 * domain_discretization_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/parallelization/domain_distribution.h"

using namespace std;

namespace ug{
namespace bridge{
namespace DomainDisc{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedef
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	string approxGrp = grp; approxGrp.append("/ApproximationSpace");
	string domDiscGrp = grp; domDiscGrp.append("/SpatialDisc");

//	DomainDiscretization
	{
		typedef IDomainDiscretization<TAlgebra> TBase;
		typedef DomainDiscretization<TDomain, TAlgebra> T;
		string name = string("DomainDiscretization").append(suffix);
		reg.add_class_<T, TBase>(name, domDiscGrp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("ApproximationSpace")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IDomainConstraint<TDomain, TAlgebra> >)>(&T::add), "", "Post Process")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IDomainElemDisc<TDomain> >)>(&T::add), "", "Element Discretization")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IDiscretizationItem<TDomain, TAlgebra> >)>(&T::add), "", "DiscItem")
			.add_method("force_regular_grid", &T::force_regular_grid)
			.add_method("assemble_linear", static_cast<void (T::*)(typename TAlgebra::matrix_type&, GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra>&)>(&T::template assemble_linear<SurfaceDoFDistribution>))
			.add_method("assemble_rhs", static_cast<void (T::*)(typename TAlgebra::vector_type&, GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra>&)>(&T::template assemble_rhs<SurfaceDoFDistribution>))
			.add_method("assemble_rhs", static_cast<void (T::*)(GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra>&)>(&T::template assemble_rhs<SurfaceDoFDistribution>))
			.add_method("adjust_solution", static_cast<void (T::*)(GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra>&)>(&T::template adjust_solution<SurfaceDoFDistribution>))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DomainDiscretization", tag);
	}

//	IDiscretizationItem
	{
		typedef IDiscretizationItem<TDomain, TAlgebra> T;
		string name = string("IDiscretizationItem").append(suffix);
		reg.add_class_<T>(name, domDiscGrp);
		reg.add_class_to_group(name, "IDiscretizationItem", tag);
	}

//	DomainDistribution
	{
		typedef ug::GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;
		reg.add_function("DistributeDomain",
			static_cast<bool (*)(TDomain&, PartitionMap&, bool, std::vector<SmartPtr<TFct> >)>(
				&DistributeDomain<TDomain, TFct>),
			grp);
		//note that an overload of this method exists, regisdered in domain_bridges/domain_bridge
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

}

}; // end Functionality

}// namespace DomainDisc

void RegisterBridge_DomainDisc(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef DomainDisc::Functionality Functionality;

	try{
//		RegisterCommon<Functionality>(reg,grp);
//		RegisterDimensionDependent<Functionality>(reg,grp);
//		RegisterDomainDependent<Functionality>(reg,grp);
//		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
