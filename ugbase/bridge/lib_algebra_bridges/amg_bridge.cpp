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

#include "lib_algebra_bridge.h"

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/operator_impl.h"

#include "lib_algebra/operator/preconditioner/amg/rsamg/rsamg.h"
#include "lib_algebra/operator/preconditioner/amg/famg/famg.h"

using namespace std;

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
struct RegisterAMGClass
{
	static bool reg(Registry &reg, string parentGroup)
	{
		return true;
	}
};

template <>
struct RegisterAMGClass<CPUAlgebra>
{
	static bool reg(Registry &reg, string parentGroup)
	{
	//	typedefs for this algebra
		typedef CPUAlgebra algebra_type;
		typedef algebra_type::vector_type vector_type;
		typedef algebra_type::matrix_type matrix_type;

	//	get group string (use same as parent)
		std::string grp = std::string(parentGroup);
		std::stringstream grpSS2; grpSS2 << grp << "/Preconditioner";
		std::string grp2 = grpSS2.str();

	//	suffix and tag
		string algSuffix = GetAlgebraSuffix<CPUAlgebra>();
		string algTag = GetAlgebraTag<CPUAlgebra>();


#ifdef UG_PARALLEL
		reg.add_class_<	IParallelCoarsening > ("IParallelCoarsening", grp2.c_str(), "Parallel Coarsening Interface (RSAMG)");
		reg.add_function("GetFullSubdomainBlockingCoarsening", GetFullSubdomainBlockingCoarsening, grp2.c_str());
		reg.add_function("GetColorCoarsening", GetColorCoarsening, grp2.c_str());
		reg.add_function("GetRS3Coarsening", GetRS3Coarsening, grp2.c_str());
		reg.add_function("GetCLJPCoarsening", GetCLJPCoarsening, grp2.c_str());
		reg.add_function("GetFalgoutCoarsening", GetFalgoutCoarsening, grp2.c_str());
		reg.add_function("GetMinimumSubdomainBlockingCoarsening", GetMinimumSubdomainBlockingCoarsening, grp2.c_str());
		reg.add_function("GetCoarseGridClassificationCoarsening", GetCoarseGridClassificationCoarsening, grp2.c_str());
		reg.add_function("GetSimpleParallelCoarsening", GetSimpleParallelCoarsening, grp2.c_str());
#endif

	//	AMG

		reg.add_class_< AMGBase<algebra_type>::LevelInformation > (string("AMGLevelInformation").append(algSuffix), grp2.c_str())
			.add_method("get_creation_time_ms", &AMGBase<algebra_type>::LevelInformation::get_creation_time_ms, "creation time of this level (in ms)")
			.add_method("get_nr_of_nodes", &AMGBase<algebra_type>::LevelInformation::get_nr_of_nodes, "nr of nodes of this level, sum over all processors")
			.add_method("get_nr_of_nodes_min", &AMGBase<algebra_type>::LevelInformation::get_nr_of_nodes_min, "nr of nodes of this level, minimum over all processors")
			.add_method("get_nr_of_nodes_max", &AMGBase<algebra_type>::LevelInformation::get_nr_of_nodes_max, "nr of nodes of this level, maximum over all processors")
			.add_method("get_nnz", &AMGBase<algebra_type>::LevelInformation::get_nnz, "nr of non-zeros, sum over all processors")
			.add_method("get_nnz_min", &AMGBase<algebra_type>::LevelInformation::get_nnz_min, "nr of non-zeros, minimum over all processors")
			.add_method("get_nnz_max", &AMGBase<algebra_type>::LevelInformation::get_nnz_max, "nr of non-zeros, maximum over all processors")
			.add_method("get_fill_in", &AMGBase<algebra_type>::LevelInformation::get_fill_in, "nr of non-zeros / (nr of nodes)^2")
			.add_method("is_valid", &AMGBase<algebra_type>::LevelInformation::is_valid, "true if this is a valid level information")
			.add_method("get_nr_of_interface_elements", &AMGBase<algebra_type>::LevelInformation::get_nr_of_interface_elements, "nr of interface elements (including multiplicites)");

		reg.add_class_to_group(string("AMGLevelInformation").append(algSuffix), "AMGLevelInformation", algTag);

	//todo: existance of AMGPreconditioner class should not depend on defines.
		reg.add_class_<	AMGBase<algebra_type>, IPreconditioner<algebra_type> > (string("AMGBase").append(algSuffix), grp2.c_str())
			.add_method("set_num_presmooth", &AMGBase<algebra_type>::set_num_presmooth, "", "nu1", "sets nr. of presmoothing steps (nu1)")
			.add_method("get_num_presmooth", &AMGBase<algebra_type>::get_num_presmooth, "nr. of presmoothing steps (nu1)")

			.add_method("set_num_postsmooth", &AMGBase<algebra_type>::set_num_postsmooth, "", "nu2", "sets nr. of postsmoothing steps (nu2)")
			.add_method("get_num_postsmooth", &AMGBase<algebra_type>::get_num_postsmooth, "nr. of postsmoothing steps (nu2)")

			.add_method("set_cycle_type", &AMGBase<algebra_type>::set_cycle_type, "", "gamma", "sets cycle type in multigrid cycle (gamma)")
			.add_method("get_cycle_type", &AMGBase<algebra_type>::get_cycle_type, "cycle type in multigrid cycle (gamma)")

			.add_method("set_max_levels", &AMGBase<algebra_type>::set_max_levels, "", "max_levels", "sets max nr of AMG levels")
			.add_method("get_max_levels", &AMGBase<algebra_type>::get_max_levels,  "max nr of AMG levels")
			.add_method("get_used_levels", &AMGBase<algebra_type>::get_used_levels, "used nr of AMG levels")

			.add_method("set_max_nodes_for_base", &AMGBase<algebra_type>::set_max_nodes_for_base, "", "maxNrOfNodes", "sets the maximal nr of nodes for base solver")
			.add_method("get_max_nodes_for_base", &AMGBase<algebra_type>::get_max_nodes_for_base, "maximal nr of nodes for base solver")

			.add_method("set_min_nodes_on_one_processor", &AMGBase<algebra_type>::set_min_nodes_on_one_processor, "", "minNrOfNodes", "if the node number on one processor falls below this, agglomerate")
			.add_method("get_min_nodes_on_one_processor", &AMGBase<algebra_type>::get_min_nodes_on_one_processor, "minNrOfNodes", "", "if the node number on one processor falls below this, agglomerate")

			.add_method("set_preferred_nodes_on_one_processor", &AMGBase<algebra_type>::set_preferred_nodes_on_one_processor, "", "preferredNrOfNodes", "if we need to agglomerate, ensure all nodes have more than this number of unknowns")
			.add_method("get_preferred_nodes_on_one_processor", &AMGBase<algebra_type>::get_preferred_nodes_on_one_processor, "preferredNrOfNodes", "", "if we need to agglomerate, ensure all nodes have more than this number of unknowns")

			.add_method("set_max_fill_before_base", &AMGBase<algebra_type>::set_max_fill_before_base, "", "fillrate", "sets maximal fill rate before base solver is used")
			.add_method("get_max_fill_before_base", &AMGBase<algebra_type>::get_max_fill_before_base, "maximal fill rate before base solver is used", "")

			.add_method("get_operator_complexity", &AMGBase<algebra_type>::get_operator_complexity, "operator complexity c_A", "", "c_A = total nnz of all matrices divided by nnz of matrix A")
			.add_method("get_grid_complexity", &AMGBase<algebra_type>::get_grid_complexity, "grid complexity c_G", "", "c_G = total number of nodes of all levels divided by number of nodes on level 0")
			.add_method("get_timing_whole_setup_ms", &AMGBase<algebra_type>::get_timing_whole_setup_ms, "the time spent on the whole setup in ms")
			.add_method("get_timing_coarse_solver_setup_ms", &AMGBase<algebra_type>::get_timing_coarse_solver_setup_ms, "the time spent in the coarse solver setup in ms")

			.add_method("get_level_information", &AMGBase<algebra_type>::get_level_information, "information about the level L", "L")

			.add_method("set_presmoother", &AMGBase<algebra_type>::set_presmoother, "", "presmoother")
			.add_method("set_postsmoother", &AMGBase<algebra_type>::set_postsmoother, "", "postsmoother")
			.add_method("set_base_solver", &AMGBase<algebra_type>::set_base_solver, "", "basesmoother")
			.add_method("check", &AMGBase<algebra_type>::check, "", "x#b", "performs a check of convergence on all levels")
			.add_method("check2", &AMGBase<algebra_type>::check2, "", "x#b", "performs a check of convergence on all levels")
			.add_method("check_fsmoothing", &AMGBase<algebra_type>::check_fsmoothing, "", "", "")
			.add_method("set_matrix_write_path", &AMGBase<algebra_type>::set_matrix_write_path, "", "matrixWritePath", "set the path where connectionviewer matrices of the levels are written")
			.add_method("set_fsmoothing", &AMGBase<algebra_type>::set_fsmoothing, "", "enable", "")
			.add_method("get_fsmoothing", &AMGBase<algebra_type>::get_fsmoothing, "f smoothing enabled", "")

			.add_method("set_position_provider",
					(void(AMGBase<algebra_type>::*)(IPositionProvider<2> *))&AMGBase<algebra_type>::set_position_provider, "", "prov", "needed for connectionviewer output")
			.add_method("set_position_provider",
					(void(AMGBase<algebra_type>::*)(IPositionProvider<3> *))&AMGBase<algebra_type>::set_position_provider, "", "prov", "needed for connectionviewer output")
			.add_method("write_interfaces", &AMGBase<algebra_type>::write_interfaces)
			;
		reg.add_class_to_group(string("AMGBase").append(algSuffix), "AMGBase", algTag);

		reg.add_class_<	RSAMG<algebra_type>, AMGBase<algebra_type> > (string("RSAMGPreconditioner").append(algSuffix), grp2.c_str(), "Ruge-Stueben Algebraic Multigrid Preconditioner")
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
			;

		reg.add_class_to_group(string("RSAMGPreconditioner").append(algSuffix), "RSAMGPreconditioner", algTag);

		reg.add_class_<	FAMG<algebra_type>, AMGBase<algebra_type> > (string("FAMGPreconditioner").append(algSuffix), grp2.c_str(), "Filtering Algebraic Multigrid")
			.add_constructor()
			.add_method("tostring", &FAMG<algebra_type>::tostring)
			.add_method("set_aggressive_coarsening", &FAMG<algebra_type>::set_aggressive_coarsening)
			.add_method("set_delta", &FAMG<algebra_type>::set_delta, "", "delta", "\"Interpolation quality\" F may not be worse than this (F < m_delta)")
			.add_method("get_delta", &FAMG<algebra_type>::get_delta, "delta")
			.add_method("set_theta", &FAMG<algebra_type>::set_theta, "" , "theta", "with multiple parents paris, discard pairs with m_theta * F > min F.")
			.add_method("get_theta", &FAMG<algebra_type>::get_theta, "theta")
			.add_method("set_damping_for_smoother_in_interpolation_calculation",
					&FAMG<algebra_type>::set_damping_for_smoother_in_interpolation_calculation)

			.add_method("set_testvector_damps", &FAMG<algebra_type>::set_testvector_damps)
			.add_method("reset_testvectors", &FAMG<algebra_type>::reset_testvectors)
			.add_method("add_testvector", &FAMG<algebra_type>::add_testvector)
			.add_method("add_vector_writer", &FAMG<algebra_type>::add_vector_writer)
			.add_method("write_testvectors", &FAMG<algebra_type>::write_testvectors)
			.add_method("set_testvector_from_matrix_rows", &FAMG<algebra_type>::set_testvector_from_matrix_rows)

			.add_method("set_prolongation_truncation", &FAMG<algebra_type>::set_prolongation_truncation, "", "prolongation_tr", "sets prolongation_truncation, a parameter used for truncation of interpolation")
			.add_method("get_prolongation_truncation", &FAMG<algebra_type>::get_prolongation_truncation, "prolongation_tr")

			.add_method("set_galerkin_truncation", &FAMG<algebra_type>::set_galerkin_truncation, "", "galerkin_tr", "sets galerkin truncation, a parameter used to truncate the galerkin product. use with care! (like 1e-12)")
			.add_method("get_galerkin_truncation", &FAMG<algebra_type>::get_galerkin_truncation, "galerkin_tr")

			.add_method("set_H_reduce_interpolation_nodes_parameter", &FAMG<algebra_type>::set_H_reduce_interpolation_nodes_parameter, "", "HreduceParameter", "we can restrict the number of parent nodes by looking at the entries of H(i,j) to prevent high fill in rates")
			.add_method("get_H_reduce_interpolation_nodes_parameter", &FAMG<algebra_type>::get_H_reduce_interpolation_nodes_parameter, "HreduceParameter")

			.add_method("set_prereduce_A_parameter", &FAMG<algebra_type>::set_prereduce_A_parameter, "", "prereduceA", "by setting this != 0.0, we reduce the matrix A before using it to its strong connections. ex. parameter is 0.1")
			.add_method("get_prereduce_A_parameter", &FAMG<algebra_type>::get_prereduce_A_parameter, "prereduceA")

			.add_method("set_external_coarsening", &FAMG<algebra_type>::set_external_coarsening)
			.add_method("set_use_precalculate", &FAMG<algebra_type>::set_use_precalculate)

			.add_method("set_debug_level_overlap", &FAMG<algebra_type>::set_debug_level_overlap)
			.add_method("set_debug_level_testvector_calc", &FAMG<algebra_type>::set_debug_level_testvector_calc)
			.add_method("set_debug_level_phase_3", &FAMG<algebra_type>::set_debug_level_phase_3)
			.add_method("set_debug_level_calculate_parent_pairs", &FAMG<algebra_type>::set_debug_level_calculate_parent_pairs)
			.add_method("set_debug_level_recv_coarsening", &FAMG<algebra_type>::set_debug_level_recv_coarsening)
			.add_method("set_debug_level_coloring", &FAMG<algebra_type>::set_debug_level_coloring)
			.add_method("set_debug_level_get_ratings", &FAMG<algebra_type>::set_debug_level_get_ratings)
			.add_method("set_debug_level_precalculate_coarsening", &FAMG<algebra_type>::set_debug_level_precalculate_coarsening)
			.add_method("set_debug_level_aggressive_coarsening", &FAMG<algebra_type>::set_debug_level_aggressive_coarsening)
			.add_method("set_debug_level_send_coarsening", &FAMG<algebra_type>::set_debug_level_send_coarsening)
			.add_method("set_debug_level_communicate_prolongation", &FAMG<algebra_type>::set_debug_level_communicate_prolongation)
			.add_method("set_debug_level_after_communciate_prolongation", &FAMG<algebra_type>::set_debug_level_after_communciate_prolongation)
			.add_method("set_testvectorsmoother", &FAMG<algebra_type>::set_testvectorsmoother)
#ifdef UG_PARALLEL
			.add_method("set_parallel_coarsening", &FAMG<algebra_type>::set_parallel_coarsening)
#endif
			;
		reg.add_class_to_group(string("FAMGPreconditioner").append(algSuffix), "FAMGPreconditioner", algTag);




		return true;
	}

};

bool RegisterAMG(Registry& reg, string parentGroup)
{
	return RegisterAlgebraClass<RegisterAMGClass>(reg, parentGroup);
}


}
}
