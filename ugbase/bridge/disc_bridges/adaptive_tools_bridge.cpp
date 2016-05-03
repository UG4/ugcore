/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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
#include "lib_disc/function_spaces/error_indicator.h"
#include "lib_disc/function_spaces/error_elem_marking_strategy.h"
#include "lib_disc/function_spaces/level_transfer.h"
#include "lib_disc/function_spaces/local_transfer.h"

// error-indicator implementation
#include "lib_grid/grid_objects/constraint_traits.h"
#include "lib_grid/algorithms/normal_calculation.h"
#include "lib_grid/algorithms/volume_calculation.h"

using namespace std;

namespace ug{

#ifdef UG_FOR_LUA
	template <typename TDomain, typename TAlgebra>
	static void MarkForAdaption_L2ErrorExactLUA(IRefiner& refiner,
	                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
	                                   const char* exactSolCallbackName,
	                                   const char* cmp,
	                                   number minL2Error,
	                                   number maxL2Error,
	                                   number refFrac,
	                                   int minLvl, int maxLvl,
	                                   number time, int quadOrder)
	{
		SmartPtr<UserData<number, TDomain::dim> > spCallback
			= make_sp(new LuaUserData<number, TDomain::dim>(exactSolCallbackName));
		MarkForAdaption_L2ErrorExact(refiner, u, spCallback, cmp, minL2Error,
									 maxL2Error, refFrac, minLvl, maxLvl, time, quadOrder);
	}

	template <typename TDomain, typename TAlgebra>
	static number MarkForAdaption_ResidualErrorP1AbsoluteLUA(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   const char* fCallbackName,
                                   const char* cmp,
                                   number time,
                                   number refTol,
                                   number coarsenTol,
                                   int maxLvl,
                                   int quadOrder, std::string quadType,
                                   bool refTopLvlOnly)
	{
		SmartPtr<UserData<number, TDomain::dim> > spCallback
			= make_sp(new LuaUserData<number, TDomain::dim>(fCallbackName));
		return MarkForAdaption_ResidualErrorP1Absolute(refiner, u, spCallback,
									cmp, time, refTol, coarsenTol, maxLvl, quadOrder,
									quadType, refTopLvlOnly);
	}

	template <typename TDomain, typename TAlgebra>
	static void MarkForAdaption_ResidualErrorP1RelativeLUA(IRefiner& refiner,
                                   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
                                   const char* fCallbackName,
                                   const char* cmp,
                                   number time,
                                   number refFrac,
                                   int minLvl, int maxLvl,
                                   int quadOrder, std::string quadType)
	{
		SmartPtr<UserData<number, TDomain::dim> > spCallback
			= make_sp(new LuaUserData<number, TDomain::dim>(fCallbackName));
		MarkForAdaption_ResidualErrorP1Relative(refiner, u, spCallback, cmp, time, refFrac,
												minLvl, maxLvl, quadOrder, quadType);
	}

#endif
	

namespace bridge{
namespace AdaptiveTools{

/**
 * \defgroup adaptivetools_bridge Adaptive Tools Bridge
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

//	Refinement
	{
	//	MarkForAdaption_GradientIndicator
		string grp("ug4/Refinement/");
		reg.add_function("MarkForAdaption_GradientIndicator",
						 &MarkForAdaption_GradientIndicator<TDomain, TAlgebra>, grp);
	
	//	MarkForAdaption_AbsoluteGradientIndicator
		reg.add_function("MarkForAdaption_AbsoluteGradientIndicator",
				 &MarkForAdaption_AbsoluteGradientIndicator<TDomain, TAlgebra>, grp);

	//	MarkForAdaption_GradientJumpIndicator
		reg.add_function("MarkForAdaption_GradientJumpIndicator",
						 &MarkForAdaption_GradientJumpIndicator<TDomain, TAlgebra>, grp);
						 
	//	MarkForAdaption_AbsoluteGradientJumpIndicator
		reg.add_function("MarkForAdaption_AbsoluteGradientJumpIndicator",
						 &MarkForAdaption_AbsoluteGradientJumpIndicator<TDomain, TAlgebra>, grp);

	//	MarkForAdaption_L2ErrorExact
		reg.add_function("MarkForAdaption_L2ErrorExact",
						 &MarkForAdaption_L2ErrorExact<TDomain, TAlgebra>, grp);
		
		#ifdef UG_FOR_LUA
			reg.add_function("MarkForAdaption_L2ErrorExact",
						 	 &MarkForAdaption_L2ErrorExactLUA<TDomain, TAlgebra>, grp);
		#endif

		reg.add_function("MarkForAdaption_GradientJump",
				 &MarkForAdaption_GradientJump<TDomain, TAlgebra>, grp);

		reg.add_function("MarkForAdaption_GradientAverage",
				&MarkForAdaption_GradientAverage<TDomain, TAlgebra>, grp);

	//	MarkForAdaption_ResidualErrorP1Absolute
		reg.add_function("MarkForAdaption_ResidualErrorP1Absolute",
						 &MarkForAdaption_ResidualErrorP1Absolute<TDomain, TAlgebra>, grp);
		
		#ifdef UG_FOR_LUA
			reg.add_function("MarkForAdaption_ResidualErrorP1Absolute",
						 	 &MarkForAdaption_ResidualErrorP1AbsoluteLUA<TDomain, TAlgebra>, grp);
		#endif

	//	MarkForAdaption_ResidualErrorP1Relative
		reg.add_function("MarkForAdaption_ResidualErrorP1Relative",
						 &MarkForAdaption_ResidualErrorP1Relative<TDomain, TAlgebra>, grp);
		
		#ifdef UG_FOR_LUA
			reg.add_function("MarkForAdaption_ResidualErrorP1Relative",
						 	 &MarkForAdaption_ResidualErrorP1RelativeLUA<TDomain, TAlgebra>, grp);
		#endif
	}

//	Prolongate
	{
		reg.add_function("Prolongate", &Prolongate<TDomain, TAlgebra>, grp);
	}

//	Restrict
	{
		reg.add_function("Restrict", &Restrict<TDomain, TAlgebra>, grp);
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

	//	group string
	grp.append("/Adaptive");

	//  IElementMarkingStrategy
	{
		typedef IElementMarkingStrategy<TDomain> T;
		string name = string("IElementMarkingStrategy").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("global_estimated_error", &T::global_estimated_error)
			.add_method("global_estimated_error_per_elem_max", &T::global_estimated_error_per_elem_max)
			.add_method("global_estimated_error_per_elem_min", &T::global_estimated_error_per_elem_min);
		reg.add_class_to_group(name, "IElementMarkingStrategy", tag);
	}


	//  StdRefinementStrategy
	{
		typedef StdRefinementMarkingStrategy<TDomain> T;
		typedef IElementMarkingStrategy<TDomain> TBase;
		string name = string("StdRefinementMarking").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(number, int)>("tolerated max error#max level")
			.add_method("set_tolerance", &T::set_tolerance)
			.add_method("set_max_level", &T::set_max_level)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StdRefinementMarking", tag);
	}


	//  GlobalMarking
	{
		typedef GlobalMarking<TDomain> T;
		typedef IElementMarkingStrategy<TDomain> TBase;
		string name = string("GlobalMarking").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(number, int)>("tolerated max error#max level")
			.add_method("set_tolerance", &T::set_tolerance)
			.add_method("set_max_level", &T::set_max_level)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GlobalMarking", tag);
	}


	//  StdCoarseningStrategy
	{
		typedef StdCoarseningMarkingStrategy<TDomain> T;
		typedef IElementMarkingStrategy<TDomain> TBase;
		string name = string("StdCoarseningMarking").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(number, number)>("tolerated max error#safety factor")
			.template add_constructor<void (*)(number)>("tolerated max error")
			.add_method("set_tolerance", &T::set_tolerance)
			.add_method("set_safety_factor", &T::set_safety_factor)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StdCoarseningMarking", tag);
	}

	//  MaximumMarking
	{
			typedef MaximumMarking<TDomain> T;
			typedef IElementMarkingStrategy<TDomain> TBase;
			string name = string("MaximumMarking").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
								   .template add_constructor<void (*)(number)>("theta")
								   .template add_constructor<void (*)(number, number)>("theta#eps")
								   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "MaximumMarking", tag);
	}

	//  MeanValueMarking
	{
			typedef MeanValueMarking<TDomain> T;
			typedef IElementMarkingStrategy<TDomain> TBase;
			string name = string("MeanValueMarking").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
								   .template add_constructor<void (*)(number, number)>("theta#factor")
								   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "MeanValueMarking", tag);
	}

	//  VarianceMarking
	{
			typedef VarianceMarking<TDomain> T;
			typedef IElementMarkingStrategy<TDomain> TBase;
			string name = string("VarianceMarking").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
									   .template add_constructor<void (*)(number)>("theta")
									   .template add_constructor<void (*)(number, number)>("theta#eps")
									   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VarianceMarking", tag);
	}


	//  VarianceMarking2
	{
			typedef VarianceMarkingEta<TDomain> T;
			typedef IElementMarkingStrategy<TDomain> TBase;
			string name = string("VarianceMarkingEta").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
						.template add_constructor<void (*)(number)>("theta")
						.template add_constructor<void (*)(number, number)>("theta#eps")
						.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VarianceMarkingEta", tag);
	}

	//  EquilibrationMarking
	{
		typedef EquilibrationMarkingStrategy<TDomain> T;
		typedef IElementMarkingStrategy<TDomain> TBase;
		string name = string("EquilibrationMarking").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
									.template add_constructor<void (*)(number)>("theta")
									.template add_constructor<void (*)(number, number)>("theta#eps")
									.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "EquilibrationMarking", tag);
	}

	//  AbsoluteMarking
	{
			typedef AbsoluteMarking<TDomain> T;
			typedef IElementMarkingStrategy<TDomain> TBase;
			string name = string("AbsoluteMarking").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
										.template add_constructor<void (*)(number)>("eta")
										.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "AbsoluteMarking", tag);
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
}

}; // end Functionality

// end group adaptivetools_bridge
/// \}

}// end AdaptiveTools

/// \addtogroup adaptivetools_bridge
void RegisterBridge_AdaptiveTools(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef AdaptiveTools::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace bridge
}// namespace ug
