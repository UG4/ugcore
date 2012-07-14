/*
 * amg_bridge.cpp
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

#include "lib_algebra/operator/preconditioner/amg/famg/famg.h"

using namespace std;

namespace ug{
namespace bridge{
namespace registerFAMG{

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
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

	{
		typedef FAMG<algebra_type> T;
		typedef AMGBase<algebra_type> TBase;

		reg.add_class_<T, TBase> (string("FAMGPreconditioner").append(suffix), grp, "Filtering Algebraic Multigrid")
			.add_constructor()
			.add_method("tostring", &T::tostring)
			.add_method("set_aggressive_coarsening", &T::set_aggressive_coarsening)
			.add_method("set_delta", &T::set_delta, "", "delta", "\"Interpolation quality\" F may not be worse than this (F < m_delta). e.g. 0.5")
			.add_method("get_delta", &T::get_delta, "delta")
			.add_method("set_theta", &T::set_theta, "" , "theta", "with multiple parents paris, discard pairs with m_theta * F > min F. e.g. 0.9")
			.add_method("get_theta", &T::get_theta, "theta")


			.add_method("set_damping_for_smoother_in_interpolation_calculation",&T::set_damping_for_smoother_in_interpolation_calculation)

			.add_method("add_testvector", (void(T::*)(vector_type& c, double weight))&T::add_testvector, "testVector#weight",
					"adds a testvector with weight")
			.add_method("add_testvector", (void(T::*)(IVectorWriter<vector_type> *vw, double weight))&T::add_testvector, "testVector#weight",
					"adds a testvector with weight by using the IVectorWriter interface")
			.add_method("set_write_testvectors", &T::set_write_testvectors, "bWrite", "if true, write testvectors to path specified in set_matrix_write_path")
			.add_method("set_testvector_from_matrix_rows", &T::set_testvector_from_matrix_rows, "", "testvector is obtained by setting 1 for dirichlet nodes (nodes with only A(i,i) != 0) and 0 everywhere else")
			.add_method("set_testvector_smoother", &T::set_testvector_smoother, "smoother", "sets the smoother to smooth testvectors")
			.add_method("set_testvector_smooths", &T::set_testvector_smooths, "n", "number of smoothing steps to smooth testvectors")

			.add_method("reset_testvectors", &T::reset_testvectors, "", "removes all added testvectors")

			.add_method("set_prolongation_truncation", &T::set_prolongation_truncation, "", "prolongation_tr", "sets prolongation_truncation, a parameter used for truncation of interpolation. use with care! (like 1e-5)")
			.add_method("get_prolongation_truncation", &T::get_prolongation_truncation, "prolongation_tr")

			.add_method("set_galerkin_truncation", &T::set_galerkin_truncation, "", "galerkin_tr", "sets galerkin truncation, a parameter used to truncate the galerkin product. use with care! (like 1e-9)")
			.add_method("get_galerkin_truncation", &T::get_galerkin_truncation, "galerkin_tr")

			.add_method("set_H_reduce_interpolation_nodes_parameter", &T::set_H_reduce_interpolation_nodes_parameter, "", "HreduceParameter", "we can restrict the number of parent nodes by looking at the entries of H(i,j) to prevent high fill in rates (e.g. 1e-3)")
			.add_method("get_H_reduce_interpolation_nodes_parameter", &T::get_H_reduce_interpolation_nodes_parameter, "HreduceParameter")

			.add_method("set_prereduce_A_parameter", &T::set_prereduce_A_parameter, "", "prereduceA", "by setting this != 0.0, we reduce the matrix A before using it to its strong connections. (e.g. 1e-3)")
			.add_method("get_prereduce_A_parameter", &T::get_prereduce_A_parameter, "prereduceA")

			.add_method("set_external_coarsening", &T::set_external_coarsening, "bExternalCoarsening", "You need to set_parallel_coarsening in parallel.")
			.add_method("set_strong_connection_external", &T::set_strong_connection_external, "epsilon", "set strong_connection value for coarsening (like set_epsilon_strong in RSAMG)")
			.add_method("get_strong_connection_external", &T::get_strong_connection_external)

			.add_method("set_write_coarsening", &T::set_write_coarsening)
			.add_method("set_read_coarsening", &T::set_read_coarsening)

	#ifdef UG_PARALLEL
			.add_method("set_parallel_coarsening", &T::set_parallel_coarsening, "", "parallelCoarsening", "e.g. GetColorCoarsening()")
	#endif
			.add_method("set_use_precalculate", &T::set_use_precalculate, "", "bUsePrecalculate", "experimental way of coarsening. beta.")

			.add_method("set_debug_level_overlap", &T::set_debug_level_overlap)
			.add_method("set_debug_level_testvector_calc", &T::set_debug_level_testvector_calc)
			.add_method("set_debug_level_phase_3", &T::set_debug_level_phase_3)
			.add_method("set_debug_level_calculate_parent_pairs", &T::set_debug_level_calculate_parent_pairs)
			.add_method("set_debug_level_recv_coarsening", &T::set_debug_level_recv_coarsening)
			.add_method("set_debug_level_coloring", &T::set_debug_level_coloring)
			.add_method("set_debug_level_get_ratings", &T::set_debug_level_get_ratings)
			.add_method("set_debug_level_precalculate_coarsening", &T::set_debug_level_precalculate_coarsening)
			.add_method("set_debug_level_aggressive_coarsening", &T::set_debug_level_aggressive_coarsening)
			.add_method("set_debug_level_send_coarsening", &T::set_debug_level_send_coarsening)
			.add_method("set_debug_level_communicate_prolongation", &T::set_debug_level_communicate_prolongation)
			.add_method("set_debug_level_after_communciate_prolongation", &T::set_debug_level_after_communciate_prolongation)


			.add_method("set_write_f_values", &T::set_write_f_values)

			.add_method("check_testvector",	static_cast<bool(T::*)()>(&T::check_testvector))
			.add_method("check_testvector",	static_cast<bool(T::*)(size_t fromlevel)>(&T::check_testvector), "", "fromlevel")
			.add_method("check_testvector",	static_cast<bool(T::*)(size_t fromlevel, size_t tolevel)>(&T::check_testvector), "", "fromlevel#tolevel")

			.set_construct_as_smart_pointer(true)
			;
	}
	reg.add_class_to_group(string("FAMGPreconditioner").append(suffix), "FAMGPreconditioner", tag);

}
}; // end Functionality
}// end AMG

void RegisterBridge_FAMG(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef registerFAMG::Functionality Functionality;
	typedef boost::mpl::list<CPUAlgebra> AlgList;

	UG_LOG("register FAMG\n");
	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
