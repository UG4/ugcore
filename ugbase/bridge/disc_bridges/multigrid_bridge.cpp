/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

#include "lib_disc/operator/linear_operator/std_injection.h"
#include "lib_disc/operator/linear_operator/average_component.h"
#include "lib_disc/operator/linear_operator/std_transfer.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"
#include "lib_disc/operator/linear_operator/element_gauss_seidel/element_gauss_seidel.h"
#include "lib_disc/operator/linear_operator/element_gauss_seidel/component_gauss_seidel.h"

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
	typedef ApproximationSpace<TDomain> approximation_space_type;

	grp.append("/MultiGrid");


//	ITransferOperator
	{
		typedef ITransferOperator<TDomain, TAlgebra> T;
		string name = string("ITransferOperator").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ITransferOperator", tag);
	}

//	ITransferPostProcess
	{
		typedef ITransferPostProcess<TDomain, TAlgebra> T;
		string name = string("ITransferPostProcess").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ITransferPostProcess", tag);
	}

//	Standard Transfer
	{
		typedef StdTransfer<TDomain, TAlgebra> T;
		typedef ITransferOperator<TDomain, TAlgebra> TBase;
		string name = string("StdTransfer").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_restriction_damping", &T::set_restriction_damping)
			.add_method("add_constraint", &T::add_constraint)
			.add_method("set_debug", &T::set_debug)
			.add_method("set_use_transposed", &T::set_use_transposed)
			.add_method("enable_p1_lagrange_optimization", &T::enable_p1_lagrange_optimization)
			.add_method("p1_lagrange_optimization_enabled", &T::p1_lagrange_optimization_enabled)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StdTransfer", tag);
	}

//	Standard Injection
	{
		typedef StdInjection<TDomain, TAlgebra> T;
		typedef ITransferOperator<TDomain, TAlgebra> TBase;
		string name = string("StdInjection").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("Approximation Space")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StdInjection", tag);
	}

//	Average Transfer Post Process
	{
		typedef AverageComponent<TDomain, TAlgebra> T;
		typedef ITransferPostProcess<TDomain, TAlgebra> TBase;
		string name = string("AverageComponent").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const std::string&)>("Components")
			.template add_constructor<void (*)(const std::vector<std::string>&)>("Components")
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
					.template add_constructor<void (*)()>()
			.add_method("set_approximation_space", &T::set_approximation_space, "", "Approximation space")
			.add_method("set_discretization", &T::set_discretization, "", "Discretization")
			.add_method("set_base_level", &T::set_base_level, "", "Base Level")
			.add_method("set_surface_level", &T::set_surface_level, "", "Surface Level")
			.add_method("set_gathered_base_solver_if_ambiguous", &T::set_gathered_base_solver_if_ambiguous,"", "Specifies if gathered base solver used in case of Ambiguity")
			.add_method("set_base_solver", &T::set_base_solver,"","Base Solver")
			.add_method("set_smoother", &T::set_smoother,"", "Smoother")
			.add_method("set_presmoother", &T::set_presmoother,"", "Smoother")
			.add_method("set_postsmoother", &T::set_postsmoother,"", "Smoother")
			.add_method("set_cycle_type", static_cast<void (T::*)(int)>(&T::set_cycle_type),"", "Cycle Type")
			.add_method("set_cycle_type", static_cast<void (T::*)(const std::string&)>(&T::set_cycle_type),"", "Cycle Type")
			.add_method("set_num_presmooth", &T::set_num_presmooth,"", "Number PreSmooth Steps")
			.add_method("set_num_postsmooth", &T::set_num_postsmooth,"", "Number PostSmooth Steps")
			.add_method("set_transfer", &T::set_transfer,"", "Transfer")
			.add_method("set_prolongation", &T::set_prolongation,"", "Prolongation")
			.add_method("set_restriction", &T::set_restriction,"", "Restriction")
			.add_method("set_projection", &T::set_projection,"", "Projection")
			.add_method("add_prolongation_post_process", &T::add_prolongation_post_process,"", "Prolongation Post Process")
			.add_method("add_restriction_post_process", &T::add_restriction_post_process,"", "Restriction Post Process")
			.add_method("set_debug", &T::set_debug)
			.add_method("set_emulate_full_refined_grid", &T::set_emulate_full_refined_grid)
			.add_method("set_rap", &T::set_rap)
			.add_method("set_smooth_on_surface_rim", &T::set_smooth_on_surface_rim)
			.add_method("set_comm_comp_overlap", &T::set_comm_comp_overlap)
			.add_method("ignore_init_for_base_solver", static_cast<void (T::*)(bool)>(&T::ignore_init_for_base_solver), "", "ignore")
			.add_method("ignore_init_for_base_solver", static_cast<bool (T::*)() const>(&T::ignore_init_for_base_solver), "is ignored", "")
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
		.template add_constructor<void (*)(const std::string&)>("patch_type")
		.template add_constructor<void (*)(number, const std::string&)>("relax#patch_type")
		.add_method("set_relax", &T::set_relax, "", "relax")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ElementGaussSeidel", tag);
	}

	//	ComponentGaussSeidel
	{
		typedef ComponentGaussSeidel<TDomain, TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ComponentGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Vanka Preconditioner")
		.template add_constructor<void (*)(const std::vector<std::string>&)>("Cmps")
		.template add_constructor<void (*)(number, const std::vector<std::string>&)>("relax#Cmps")
		.template add_constructor<void (*)(number, const std::vector<std::string>&, const std::vector<int>&, const std::vector<number>&)>("relax#Cmps")
		.add_method("set_alpha", &T::set_alpha, "", "alpha")
		.add_method("set_beta", &T::set_beta, "", "beta")
		.add_method("set_weights", &T::set_weights, "", "weights")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ComponentGaussSeidel", tag);
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
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace ug
}//	end of namespace interface
