/*
 * rsamg_bridge.cpp
 *
 *  Created on: 23.03.2011
 *      Author: mrupp
 */


#include <iostream>
#include <sstream>
#include <string>

//#define UG_USE_AMG // temporary switch until AMG for systems works again

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

#include "lib_algebra/lib_algebra.h"

#include "lib_algebra/operator/preconditioner/amg/rsamg/rsamg.h"

using namespace std;

namespace ug{
namespace bridge{
namespace registerRSAMG{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

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
	string suffix = GetAlgebraSuffix<CPUAlgebra>();
	string tag = GetAlgebraTag<CPUAlgebra>();

//	typedefs for this algebra
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

#ifdef UG_PARALLEL
	reg.add_class_<	IParallelCoarsening > ("IParallelCoarsening", grp, "Parallel Coarsening Interface (RSAMG)");
	reg.add_function("GetFullSubdomainBlockingCoarsening", GetFullSubdomainBlockingCoarsening, grp);
	reg.add_function("GetColorCoarsening", GetColorCoarsening, grp);
	reg.add_function("GetRS3Coarsening", GetRS3Coarsening, grp);
	reg.add_function("GetCLJPCoarsening", GetCLJPCoarsening, grp);
	reg.add_function("GetFalgoutCoarsening", GetFalgoutCoarsening, grp);
	reg.add_function("GetMinimumSubdomainBlockingCoarsening", GetMinimumSubdomainBlockingCoarsening, grp);
	reg.add_function("GetCoarseGridClassificationCoarsening", GetCoarseGridClassificationCoarsening, grp);
	reg.add_function("GetSimpleParallelCoarsening", GetSimpleParallelCoarsening, grp);
#endif


	reg.add_class_<	RSAMG<algebra_type>, AMGBase<algebra_type> > (string("RSAMGPreconditioner").append(suffix), grp, "Ruge-Stueben Algebraic Multigrid Preconditioner")
		.add_constructor()
		.add_method("set_epsilon_strong", &RSAMG<algebra_type>::set_epsilon_strong, "", "epsilon_str", "sets epsilon_strong, a measure for strong connectivity")
		.add_method("get_epsilon_strong", &RSAMG<algebra_type>::get_epsilon_strong, "epsilon_strong", "")
		.add_method("set_prolongation_truncation", &RSAMG<algebra_type>::set_prolongation_truncation, "", "prolongation_tr", "sets the truncation of interpolation")
		.add_method("get_prolongation_truncation", &RSAMG<algebra_type>::get_prolongation_truncation, "prolongation_tr")

		.add_method("tostring", &RSAMG<algebra_type>::tostring)
		.add_method("enable_aggressive_coarsening_A", &RSAMG<algebra_type>::enable_aggressive_coarsening_A, "", "nrOfPaths", "enables aggressive coarsening (A1 or A2) on the first level.")
		.add_method("disable_aggressive_coarsening", &RSAMG<algebra_type>::disable_aggressive_coarsening, "", "", "disables aggressive coarsening")
		.add_method("is_aggressive_coarsening", &RSAMG<algebra_type>::is_aggressive_coarsening)
		.add_method("is_aggressive_coarsening_A", &RSAMG<algebra_type>::is_aggressive_coarsening_A)
#ifdef UG_PARALLEL
		.add_method("set_parallel_coarsening", &RSAMG<algebra_type>::set_parallel_coarsening)
#endif
		.set_construct_as_smart_pointer(true)
		;

	reg.add_class_to_group(string("RSAMGPreconditioner").append(suffix), "RSAMGPreconditioner", tag);

}
}; // end Functionality
}// end AMG

void RegisterBridge_RSAMG(Registry& reg, string grp)
{
	
	grp.append("/Algebra/Preconditioner");
	typedef registerRSAMG::Functionality Functionality;
	typedef boost::mpl::list<CPUAlgebra> AlgList;

	try{
		RegisterAlgebraDependent<Functionality, AlgList>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
