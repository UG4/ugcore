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

#include "lib_algebra/operator/preconditioner/amg/amg_base.h"


using namespace std;

namespace ug{
namespace bridge{
namespace AMG{

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

//	AMG

	reg.add_class_< typename AMGBase<algebra_type>::LevelInformation > (string("AMGLevelInformation").append(suffix), grp)
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

	reg.add_class_to_group(string("AMGLevelInformation").append(suffix), "AMGLevelInformation", tag);


//todo: existance of AMGPreconditioner class should not depend on defines.
	reg.add_class_<	AMGBase<algebra_type>, IPreconditioner<algebra_type> > (string("AMGBase").append(suffix), grp)
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
		.add_method("set_nr_of_preiterations_at_check", &AMGBase<algebra_type>::set_nr_of_preiterations_at_check,
							"", "i", "nr of mg cycles performed before checking")
		.add_method("set_preiterations_mimum_defect_at_check", &AMGBase<algebra_type>::set_preiterations_mimum_defect_at_check, "",
						"d", "minimum defect for preiterations before checking")
		.add_method("set_matrix_write_path", &AMGBase<algebra_type>::set_matrix_write_path, "", "matrixWritePath", "set the path where connectionviewer matrices of the levels are written")
		.add_method("set_fsmoothing", &AMGBase<algebra_type>::set_fsmoothing, "", "enable", "")
		.add_method("get_fsmoothing", &AMGBase<algebra_type>::get_fsmoothing, "f smoothing enabled", "")
		.add_method("set_one_init", &AMGBase<algebra_type>::set_one_init, "", "b")

		.add_method("set_position_provider",
				(void(AMGBase<algebra_type>::*)(IPositionProvider<2> *))&AMGBase<algebra_type>::set_position_provider, "", "prov", "needed for connectionviewer output")
		.add_method("set_position_provider",
				(void(AMGBase<algebra_type>::*)(IPositionProvider<3> *))&AMGBase<algebra_type>::set_position_provider, "", "prov", "needed for connectionviewer output")
		.add_method("write_interfaces", &AMGBase<algebra_type>::write_interfaces)
		.add_method("set_checkLevel_post_iterations", &AMGBase<algebra_type>::set_checkLevel_post_iterations)
		.add_method("set_Y_cycle", &AMGBase<algebra_type>::set_Y_cycle)
		;
	reg.add_class_to_group(string("AMGBase").append(suffix), "AMGBase", tag);
}
}; // end Functionality
}// end AMG


void RegisterBridge_FAMG(Registry& reg, string grp);
void RegisterBridge_RSAMG(Registry& reg, string grp);

void RegisterBridge_AMG(Registry& reg, string grp)
{
	RegisterBridge_RSAMG(reg, grp);
	RegisterBridge_FAMG(reg, grp);

	grp.append("/Algebra/Preconditioner");

	try{
		RegisterAlgebraDependent<AMG::Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
