/*
 * grid_function_bridge.cpp
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
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"

using namespace std;

namespace ug{
namespace bridge{
namespace GridFunction{

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
	typedef ug::GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;

//	group string
	grp.append("/ApproximationSpace");

//	GridFunction
	{
		string name = string("GridFunction").append(suffix);
		reg.add_class_<TFct, vector_type>(name, grp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("ApproximationSpace")
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>, int)>("ApproximationSpace, Level")
			.add_method("assign", static_cast<void (TFct::*)(const vector_type&)>(&TFct::assign),
						"Success", "Vector")
			.add_method("clone", &TFct::clone)
			.add_method("add_transfer", &TFct::add_transfer)
			.add_method("remove_transfer", &TFct::remove_transfer)
			.add_method("clear_transfers", &TFct::clear_transfers)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunction", tag);
	}

//	GridFunctionNumberData
	{
		string name = string("GridFunctionNumberData").append(suffix);
		typedef GridFunctionNumberData<TFct> T;
		typedef UserData<number, dim> TBase;
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionNumberData", tag);
	}

//	GridFunctionVectorData
	{
		string name = string("GridFunctionVectorData").append(suffix);
		typedef GridFunctionVectorData<TFct> T;
		typedef UserData<MathVector<dim>, dim> TBase;
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Components")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionVectorData", tag);
	}

//	GridFunctionGradientData
	{
		string name = string("GridFunctionGradientData").append(suffix);
		typedef GridFunctionGradientData<TFct> T;
		typedef UserData<MathVector<dim>, dim> TBase;
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionGradientData", tag);
	}

//	GridFunctionGradientComponentData
	{
		string name = string("GridFunctionGradientComponentData").append(suffix);
		typedef GridFunctionGradientComponentData<TFct> T;
		typedef UserData<number, dim> TBase;
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*, size_t)>("GridFunction#Components")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionGradientComponentData", tag);
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
		.add_method("num_fct", &T::num_fct)
		.add_method("name", &T::name)
		.add_method("dim", &T::dim)
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int, const char*)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order # Subsets", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"crouzeix-raviart\",\"piecewise-constant\"] ", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const std::vector<std::string>&, const char*, int, const std::vector<std::string>&)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order # Subsets", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const std::vector<std::string>&, const char*, int)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const std::vector<std::string>&, const char*)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"crouzeix-raviart\",\"piecewise-constant\"] ", "Adds a function to the Function Pattern",
					"currently no help available");
	}
}

}; // end Functionality

}// namespace GridFunction

void RegisterBridge_GridFunction(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef GridFunction::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		//RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
		//RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
