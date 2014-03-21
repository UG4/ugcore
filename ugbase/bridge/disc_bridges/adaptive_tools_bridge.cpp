/*
 * adaptive_tools_bridge.cpp
 *
 *  Created on: 06.03.2012
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
#include "lib_disc/function_spaces/error_indicator.h"
#include "lib_disc/function_spaces/level_transfer.h"
#include "lib_disc/function_spaces/local_transfer.h"

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
		// reg.add_function("MarkForAdaption_L2ErrorExact",
		// 				 static_cast<void (*)(IRefiner&,
		// 										SmartPtr<GridFunction<TDomain, TAlgebra> >,
		// 										SmartPtr<UserData<number, TDomain::dim> >,
		// 										const char*, number, number,
		// 										int, int, number, int)>
		// 				 	(&MarkForAdaption_L2ErrorExact<TDomain, TAlgebra>),
		// 				 grp);

		// reg.add_function("MarkForAdaption_L2ErrorExact",
		// 				 static_cast<void (*)(IRefiner&,
		// 										SmartPtr<GridFunction<TDomain, TAlgebra> >,
		// 										const char*,
		// 										const char*, number, number,
		// 										int, int, number, int)>
		// 				 	(&MarkForAdaption_L2ErrorExact),
		// 			 	 grp);
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
