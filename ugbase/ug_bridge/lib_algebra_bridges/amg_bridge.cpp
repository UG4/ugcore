/*
 * amg_bridge.cpp
 *
 *  Created on: 23.03.2011
 *      Author: mrupp
 */


#include <iostream>
#include <sstream>

//#define UG_USE_AMG // temporary switch until AMG for systems works again

#include "lib_algebra_bridge.h"

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/operator_impl.h"

#include "lib_algebra/operator/preconditioner/amg/rsamg/rsamg.h"
#include "lib_algebra/operator/preconditioner/amg/famg/famg.h"

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
struct RegisterAMGClass
{
	static bool reg(Registry &reg, const char *parentGroup)
	{
		return true;
	}
};

template <>
struct RegisterAMGClass<CPUAlgebra>
{
	static bool reg(Registry &reg, const char *parentGroup)
	{
	//	typedefs for this algebra
		typedef CPUAlgebra algebra_type;
		typedef algebra_type::vector_type vector_type;
		typedef algebra_type::matrix_type matrix_type;

		//	get group string (use same as parent)
		std::string grp = std::string(parentGroup);
		std::stringstream grpSS2; grpSS2 << grp << "/Preconditioner";
		std::string grp2 = grpSS2.str();
	//	AMG
		reg.add_class_< amg_base<algebra_type>::LevelInformation > ("AMGLevelInformation", grp2.c_str())
			.add_method("get_creation_time_ms", &amg_base<algebra_type>::LevelInformation::get_creation_time_ms, "creation time of this level (in ms)")
			.add_method("get_nr_of_nodes", &amg_base<algebra_type>::LevelInformation::get_nr_of_nodes, "nr of nodes of this level, sum over all processors")
			.add_method("get_nr_of_nodes_min", &amg_base<algebra_type>::LevelInformation::get_nr_of_nodes_min, "nr of nodes of this level, minimum over all processors")
			.add_method("get_nr_of_nodes_max", &amg_base<algebra_type>::LevelInformation::get_nr_of_nodes_max, "nr of nodes of this level, maximum over all processors")
			.add_method("get_nnz", &amg_base<algebra_type>::LevelInformation::get_nnz, "nr of non-zeros, sum over all processors")
			.add_method("get_nnz_min", &amg_base<algebra_type>::LevelInformation::get_nnz_min, "nr of non-zeros, minimum over all processors")
			.add_method("get_nnz_max", &amg_base<algebra_type>::LevelInformation::get_nnz_max, "nr of non-zeros, maximum over all processors")
			.add_method("get_fill_in", &amg_base<algebra_type>::LevelInformation::get_fill_in, "nr of non-zeros / (nr of nodes)^2")
			.add_method("is_valid", &amg_base<algebra_type>::LevelInformation::is_valid, "true if this is a valid level information")
			.add_method("get_nr_of_interface_elements", &amg_base<algebra_type>::LevelInformation::get_nr_of_interface_elements, "nr of interface elements (including multiplicites)");

	//todo: existance of AMGPreconditioner class should not depend on defines.
		reg.add_class_<	amg_base<algebra_type>, IPreconditioner<algebra_type> > ("AMGBase", grp2.c_str())
			.add_method("set_num_presmooth", &amg_base<algebra_type>::set_num_presmooth, "", "nu1", "sets nr. of presmoothing steps (nu1)")
			.add_method("get_num_presmooth", &amg_base<algebra_type>::get_num_presmooth, "nr. of presmoothing steps (nu1)")

			.add_method("set_num_postsmooth", &amg_base<algebra_type>::set_num_postsmooth, "", "nu2", "sets nr. of postsmoothing steps (nu2)")
			.add_method("get_num_postsmooth", &amg_base<algebra_type>::get_num_postsmooth, "nr. of postsmoothing steps (nu2)")

			.add_method("set_cycle_type", &amg_base<algebra_type>::set_cycle_type, "", "gamma", "sets cycle type in multigrid cycle (gamma)")
			.add_method("get_cycle_type", &amg_base<algebra_type>::get_cycle_type, "cycle type in multigrid cycle (gamma)")

			.add_method("set_max_levels", &amg_base<algebra_type>::set_max_levels, "", "max_levels", "sets max nr of AMG levels")
			.add_method("get_max_levels", &amg_base<algebra_type>::get_max_levels,  "max nr of AMG levels")
			.add_method("get_used_levels", &amg_base<algebra_type>::get_used_levels, "used nr of AMG levels")

			.add_method("set_max_nodes_for_base", &amg_base<algebra_type>::set_max_nodes_for_base, "", "maxNrOfNodes", "sets the maximal nr of nodes for base solver")
			.add_method("get_max_nodes_for_base", &amg_base<algebra_type>::get_max_nodes_for_base, "maximal nr of nodes for base solver")

			.add_method("set_min_nodes_on_one_processor", &amg_base<algebra_type>::set_min_nodes_on_one_processor, "", "minNrOfNodes", "if the node number on one processor falls below this, agglomerate")
			.add_method("get_min_nodes_on_one_processor", &amg_base<algebra_type>::get_min_nodes_on_one_processor, "minNrOfNodes", "", "if the node number on one processor falls below this, agglomerate")

			.add_method("set_max_fill_before_base", &amg_base<algebra_type>::set_max_fill_before_base, "", "fillrate", "sets maximal fill rate before base solver is used")
			.add_method("get_max_fill_before_base", &amg_base<algebra_type>::get_max_fill_before_base, "maximal fill rate before base solver is used", "")

			.add_method("get_operator_complexity", &amg_base<algebra_type>::get_operator_complexity, "operator complexity c_A", "", "c_A = total nnz of all matrices divided by nnz of matrix A")
			.add_method("get_grid_complexity", &amg_base<algebra_type>::get_grid_complexity, "grid complexity c_G", "", "c_G = total number of nodes of all levels divided by number of nodes on level 0")
			.add_method("get_timing_whole_setup_ms", &amg_base<algebra_type>::get_timing_whole_setup_ms, "the time spent on the whole setup in ms")
			.add_method("get_timing_coarse_solver_setup_ms", &amg_base<algebra_type>::get_timing_coarse_solver_setup_ms, "the time spent in the coarse solver setup in ms")

			.add_method("get_level_information", &amg_base<algebra_type>::get_level_information, "information about the level L", "L")

			.add_method("set_presmoother", &amg_base<algebra_type>::set_presmoother, "", "presmoother")
			.add_method("set_postsmoother", &amg_base<algebra_type>::set_postsmoother, "", "postsmoother")
			.add_method("set_base_solver", &amg_base<algebra_type>::set_base_solver, "", "basesmoother")
			.add_method("check", &amg_base<algebra_type>::check, "", "x#b", "performs a check of convergence on all levels")
			.add_method("check2", &amg_base<algebra_type>::check2, "", "x#b", "performs a check of convergence on all levels")
			.add_method("set_matrix_write_path", &amg_base<algebra_type>::set_matrix_write_path, "", "matrixWritePath", "set the path where connectionviewer matrices of the levels are written")
			.add_method("set_fsmoothing", &amg_base<algebra_type>::set_fsmoothing, "", "enable", "")
			.add_method("get_fsmoothing", &amg_base<algebra_type>::get_fsmoothing, "f smoothing enabled", "")

			.add_method("set_position_provider2d", &amg_base<algebra_type>::set_position_provider2d, "", "prov", "needed for connectionviewer output")
			.add_method("set_position_provider3d", &amg_base<algebra_type>::set_position_provider3d, "", "prov", "needed for connectionviewer output")
			;

		reg.add_class_<	rsamg<algebra_type>, amg_base<algebra_type> > ("RSAMGPreconditioner", grp2.c_str(), "Ruge-Stueben Algebraic Multigrid Preconditioner")
			.add_constructor()
			.add_method("set_epsilon_strong", &rsamg<algebra_type>::set_epsilon_strong, "", "epsilon_str", "sets epsilon_strong, a measure for strong connectivity")
			.add_method("get_epsilon_strong", &rsamg<algebra_type>::get_epsilon_strong, "epsilon_strong", "")
			.add_method("set_epsilon_truncation", &rsamg<algebra_type>::set_epsilon_truncation, "", "epsilon_tr", "sets epsilon_truncation, a parameter used for truncation of interpolation")
			.add_method("get_epsilon_truncation", &rsamg<algebra_type>::get_epsilon_truncation, "epsilon_tr")

			.add_method("tostring", &rsamg<algebra_type>::tostring)
			.add_method("enable_aggressive_coarsening_A", &rsamg<algebra_type>::enable_aggressive_coarsening_A, "", "nrOfPaths", "enables aggressive coarsening (A1 or A2) on the first level.")
			.add_method("disable_aggressive_coarsening", &rsamg<algebra_type>::disable_aggressive_coarsening, "", "", "disables aggressive coarsening")
			.add_method("is_aggressive_coarsening", &rsamg<algebra_type>::is_aggressive_coarsening)
			.add_method("is_aggressive_coarsening_A", &rsamg<algebra_type>::is_aggressive_coarsening_A);

		reg.add_class_<	famg<algebra_type>, amg_base<algebra_type> > ("FAMGPreconditioner", grp2.c_str(), "Filtering Algebraic Multigrid")
			.add_constructor()
			.add_method("tostring", &famg<algebra_type>::tostring)
			.add_method("set_aggressive_coarsening", &famg<algebra_type>::set_aggressive_coarsening)
			.add_method("set_delta", &famg<algebra_type>::set_delta, "", "delta", "\"Interpolation quality\" F may not be worse than this (F < m_delta)")
			.add_method("get_delta", &famg<algebra_type>::get_delta, "delta")
			.add_method("set_theta", &famg<algebra_type>::set_theta, "" , "theta", "with multiple parents paris, discard pairs with m_theta * F > min F.")
			.add_method("get_theta", &famg<algebra_type>::get_theta, "theta")
			.add_method("set_damping_for_smoother_in_interpolation_calculation",
					&famg<algebra_type>::set_damping_for_smoother_in_interpolation_calculation)
			.add_method("set_testvector_damps", &famg<algebra_type>::set_testvector_damps)
			.add_method("reset_testvectors", &famg<algebra_type>::reset_testvectors)
			.add_method("add_testvector", &famg<algebra_type>::add_testvector)
			.add_method("add_vector_writer", &famg<algebra_type>::add_vector_writer)
			.add_method("set_epsilon_truncation", &famg<algebra_type>::set_epsilon_truncation, "", "epsilon_tr", "sets epsilon_truncation, a parameter used for truncation of interpolation")
			.add_method("get_epsilon_truncation", &famg<algebra_type>::get_epsilon_truncation, "epsilon_tr")

			.add_method("set_external_coarsening", &famg<algebra_type>::set_external_coarsening)
			.add_method("set_use_precalculate", &famg<algebra_type>::set_use_precalculate)

			.add_method("set_debug_level_overlap", &famg<algebra_type>::set_debug_level_overlap)
			.add_method("set_debug_level_testvector_calc", &famg<algebra_type>::set_debug_level_testvector_calc)
			.add_method("set_debug_level_phase_3", &famg<algebra_type>::set_debug_level_phase_3)
			.add_method("set_debug_level_calculate_parent_pairs", &famg<algebra_type>::set_debug_level_calculate_parent_pairs)
			.add_method("set_debug_level_recv_coarsening", &famg<algebra_type>::set_debug_level_recv_coarsening)
			.add_method("set_debug_level_coloring", &famg<algebra_type>::set_debug_level_coloring)
			.add_method("set_debug_level_get_ratings", &famg<algebra_type>::set_debug_level_get_ratings)
			.add_method("set_debug_level_precalculate_coarsening", &famg<algebra_type>::set_debug_level_precalculate_coarsening)
			.add_method("set_debug_level_aggressive_coarsening", &famg<algebra_type>::set_debug_level_aggressive_coarsening)
			.add_method("set_debug_level_send_coarsening", &famg<algebra_type>::set_debug_level_send_coarsening)
			.add_method("set_debug_level_communicate_prolongation", &famg<algebra_type>::set_debug_level_communicate_prolongation)
			.add_method("set_debug_level_after_communciate_prolongation", &famg<algebra_type>::set_debug_level_after_communciate_prolongation)
			;

		return true;
	}

};

bool RegisterAMG(Registry& reg, int algebra_type, const char* parentGroup)
{
	return RegisterAlgebraClass<RegisterAMGClass>(reg, algebra_type, parentGroup);
}


}
}
