/*
 * multigrid_bridge.cpp
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

#include "lib_disc/operator/linear_operator/projection_operator.h"
#include "lib_disc/operator/linear_operator/transfer_post_process.h"
#include "lib_disc/operator/linear_operator/prolongation_operator.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"
#include "lib_disc/operator/linear_operator/element_gauss_seidel/element_gauss_seidel.h"

using namespace std;

namespace ug{
namespace bridge{
namespace MultiGrid{

/**
 * \defgroup multigrid_bridge Multi Grid Bridge
 * \ingroup disc_bridge
 * \{
 */

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
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;

	grp.append("/MultiGrid");

//	Standard Transfer
	{
		typedef StdTransfer<TDomain, TAlgebra> T;
		typedef ITransferOperator<TAlgebra> TBase;
		string name = string("StdTransfer").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("Approximation Space")
			.add_method("set_restriction_damping", &T::set_restriction_damping)
			.add_method("add_constraint", &T::add_constraint)
			.add_method("set_debug", &T::set_debug)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StdTransfer", tag);
	}

//	Standard Injection
	{
		typedef InjectionTransfer<TDomain, TAlgebra> T;
		typedef ITransferOperator<TAlgebra> TBase;
		string name = string("InjectionTransfer").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("Approximation Space")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "InjectionTransfer", tag);
	}

//	Average Transfer Post Process
	{
		typedef AverageComponent<TDomain, TAlgebra> T;
		typedef ITransferPostProcess<TAlgebra> TBase;
		string name = string("AverageComponent").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>, const char*)>("Approximation Space#Components")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AverageComponent", tag);
	}

//	AssembledMultiGridCycle
	{
		typedef AssembledMultiGridCycle<TDomain, TAlgebra> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("GeometricMultiGrid").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("Approximation Space")
			.add_method("set_discretization", &T::set_discretization, "", "Discretization")
			.add_method("set_base_level", &T::set_base_level, "", "Base Level")
			.add_method("set_surface_level", &T::set_surface_level, "", "Surface Level")
			.add_method("set_parallel_base_solver", &T::set_parallel_base_solver,"", "Specifies if base solver works in parallel")
			.add_method("set_base_solver", &T::set_base_solver,"","Base Solver")
			.add_method("set_smoother", &T::set_smoother,"", "Smoother")
			.add_method("set_presmoother", &T::set_presmoother,"", "Smoother")
			.add_method("set_postsmoother", &T::set_postsmoother,"", "Smoother")
			.add_method("set_cycle_type", &T::set_cycle_type,"", "Cycle Type")
			.add_method("set_num_presmooth", &T::set_num_presmooth,"", "Number PreSmooth Steps")
			.add_method("set_num_postsmooth", &T::set_num_postsmooth,"", "Number PostSmooth Steps")
			.add_method("set_transfer", &T::set_transfer,"", "Transfer")
			.add_method("set_prolongation", &T::set_prolongation,"", "Prolongation")
			.add_method("set_restriction", &T::set_restriction,"", "Restriction")
			.add_method("set_projection", &T::set_projection,"", "Projection")
			.add_method("add_prolongation_post_process", &T::add_prolongation_post_process,"", "Prolongation Post Process")
			.add_method("add_restriction_post_process", &T::add_restriction_post_process,"", "Restriction Post Process")
			.add_method("set_debug", &T::set_debug)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GeometricMultiGrid", tag);
	}

	//	ElementGaussSeidel
	{
		typedef ElementGaussSeidel<TDomain, TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ElementGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Vanka Preconditioner")
		.add_constructor()
		.template add_constructor<void (*)(number)>("relax")
		.template add_constructor<void (*)(const std::string&)>("patch-type")
		.template add_constructor<void (*)(number, const std::string&)>("relax, patch-type")
		.add_method("set_relax", &T::set_relax, "", "relax")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ElementGaussSeidel", tag);
	}

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
	typedef typename TAlgebra::vector_type vector_type;
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	IProlongationOperator
	{
		typedef ITransferOperator<TAlgebra> T;
		string name = string("ITransferOperator").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ITransferOperator", tag);
	}

//	IProlongationOperator
	{
		typedef ITransferPostProcess<TAlgebra> T;
		string name = string("ITransferPostProcess").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ITransferPostProcess", tag);
	}
}
};

// end group multigrid_bridge
/// \}

}// namespace MultiGrid

/// \addtogroup multigrid_bridge
void RegisterBridge_MultiGrid(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef MultiGrid::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace ug
}//	end of namespace interface
