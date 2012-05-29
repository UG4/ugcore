/*
 * domain_bridge.cpp
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

// lib_algebra includes
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/integrateDraft.h"
#include "lib_disc/dof_manager/cuthill_mckee.h"
#include "lib_disc/dof_manager/lexorder.h"

using namespace std;

namespace ug{
namespace bridge{
namespace DiscDomain{

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
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;

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

//	GridFunction
	{
		string name = string("GridFunction").append(suffix);
		reg.add_class_<TFct, vector_type>(name, approxGrp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("ApproximationSpace")
			.add_method("assign", static_cast<void (TFct::*)(const vector_type&)>(&TFct::assign),
						"Success", "Vector")
			.add_method("clone", &TFct::clone)
			.add_method("add_transfer", &TFct::add_transfer)
			.add_method("remove_transfer", &TFct::remove_transfer)
			.add_method("clear_transfers", &TFct::clear_transfers)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunction", tag);
	}

//	InterpolateFunction
	{
		reg.add_function("InterpolateFunction", static_cast<void (*)(number, TFct&, const char*, number)>(&InterpolateFunction<TFct>), grp);
		reg.add_function("InterpolateFunction", static_cast<void (*)(number, TFct&, const char*, number, const char*)>(&InterpolateFunction<TFct>), grp);

		reg.add_function("InterpolateFunction", static_cast<void (*)(IPData<number, dim>&, TFct&, const char*, number)>(&InterpolateFunction<TFct>), grp);
		reg.add_function("InterpolateFunction", static_cast<void (*)(IPData<number, dim>&, TFct&, const char*, number, const char*)>(&InterpolateFunction<TFct>),grp);

		#ifdef UG_FOR_LUA
		reg.add_function("InterpolateFunction", static_cast<void (*)(const char*, TFct&, const char*, number)>(&InterpolateFunction<TFct>), grp);
		reg.add_function("InterpolateFunction", static_cast<void (*)(const char*, TFct&, const char*, number, const char*)>(&InterpolateFunction<TFct>),grp);
		#endif
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
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	grp.append("/ApproximationSpace");

//  ApproximationSpace
	{
		typedef ApproximationSpace<TDomain> T;
		typedef IApproximationSpace TBase;
		string name = string("ApproximationSpace").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TDomain>)>("Domain")
			.add_method("domain", static_cast<SmartPtr<TDomain> (T::*)()>(&T::domain))
			.add_method("surface_view", static_cast<ConstSmartPtr<SurfaceView> (T::*)() const>(&T::surface_view))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ApproximationSpace", tag);
	}

//	Order Cuthill-McKee
	{
		reg.add_function("OrderCuthillMcKee", static_cast<void (*)(approximation_space_type&, bool)>(&OrderCuthillMcKee), grp);
	}

//	Order lexicographically
	{
		reg.add_function("OrderLex", static_cast<void (*)(approximation_space_type&, const char*)>(&OrderLex<TDomain>), grp);
	}

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
//	GridLevel
	reg.add_class_<GridLevel>("GridLevel", grp)
		.add_constructor()
		.add_constructor<void (*)(int)>("Level")
		.add_constructor<void (*)(int, std::string)>("Level, Type")
		.set_construct_as_smart_pointer(true);

//	IApproximationSpace
	{
	typedef IApproximationSpace T;
	reg.add_class_<T>("IApproximationSpace", grp)
		.add_method("print_statistic", static_cast<void (T::*)(int) const>(&T::print_statistic))
		.add_method("print_statistic", static_cast<void (T::*)() const>(&T::print_statistic))
		.add_method("print_layout_statistic", static_cast<void (T::*)(int) const>(&T::print_layout_statistic))
		.add_method("print_layout_statistic", static_cast<void (T::*)() const>(&T::print_layout_statistic))
		.add_method("print_local_dof_statistic", static_cast<void (T::*)(int) const>(&T::print_local_dof_statistic))
		.add_method("print_local_dof_statistic", static_cast<void (T::*)() const>(&T::print_local_dof_statistic))
		.add_method("num_levels", &T::num_levels)
		.add_method("init_levels", &T::init_levels)
		.add_method("init_surfaces", &T::init_surfaces)
		.add_method("init_top_surface", &T::init_top_surface)
		.add_method("defragment", &T::defragment)
		.add_method("clear", &T::clear)
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int, const char*)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order # Subsets", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order", "Adds a function to the Function Pattern",
					"currently no help available");

	}
}

}; // end Functionality

}// namespace DiscDomain

void RegisterBridge_DiscDomain(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef DiscDomain::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
