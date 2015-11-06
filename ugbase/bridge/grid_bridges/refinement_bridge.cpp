/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include "grid_bridges.h"
#include "lib_grid/algorithms/refinement/adaptive_regular_mg_refiner.h"
#include "lib_grid/algorithms/refinement/global_fractured_media_refiner.h"
#include "lib_grid/algorithms/refinement/global_multi_grid_refiner.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_grid.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_multi_grid.h"
#include "lib_grid/algorithms/refinement/refinement_projectors/loop_subdivision_projectors.h"
#include "lib_grid/algorithms/subdivision/subdivision_loop.h"
#include "lib_grid/algorithms/subdivision/subdivision_volumes.h"
#include "lib_grid/file_io/file_io.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/parallel_refinement/parallel_hanging_node_refiner_multi_grid.h"
#endif

using namespace std;

namespace ug{
namespace bridge{

void TestSubdivision(const char* fileIn, const char* fileOut, int numRefs)
{
	PROFILE_FUNC_GROUP("grid");
//todo: Callbacks have to make sure that their attachment is accessible in the grid.
//		even if they were initialized before the attachment was attached to the grid.
	MultiGrid mg;
	SubsetHandler sh(mg);
	SubdivisionLoopProjector<APosition> refCallback(mg, aPosition, aPosition);
	GlobalMultiGridRefiner ref(mg, &refCallback);
	
	if(LoadGridFromFile(mg, sh, fileIn)){
		for(int lvl = 0; lvl < numRefs; ++lvl){
			ref.refine();
		}

		ProjectToLimitPLoop(mg, aPosition, aPosition);
		SaveGridToFile(mg, mg.get_hierarchy_handler(), fileOut);

	}
	else{
		UG_LOG("Load failed. aborting...\n");
	}
}

bool CreateHierarchy(MultiGrid& mg, size_t numRefs)
{
	PROFILE_FUNC_GROUP("grid");

	GlobalMultiGridRefiner ref(mg);

	for(size_t lvl = 0; lvl < numRefs; ++lvl){
		ref.refine();
	}
	return true;
}

bool CreateSmoothHierarchy(MultiGrid& mg, size_t numRefs)
{
	PROFILE_FUNC_GROUP("grid");
	IRefinementCallback* refCallback = NULL;
//	we're only checking for the main attachments here.
//todo: improve this - add a domain-based hierarchy creator.
	if(mg.has_vertex_attachment(aPosition1))
		refCallback = new SubdivisionLoopProjector<APosition1>(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new SubdivisionLoopProjector<APosition2>(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new SubdivisionLoopProjector<APosition>(mg, aPosition, aPosition);
		
	if(!refCallback){
		UG_LOG("No standard position attachment found. Aborting.\n");
		return false;
	}
	
	GlobalMultiGridRefiner ref(mg, refCallback);

	for(size_t lvl = 0; lvl < numRefs; ++lvl){
		ref.refine();
	}

	if(mg.has_vertex_attachment(aPosition1))
		ProjectToLimitPLoop(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		ProjectToLimitPLoop(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		ProjectToLimitPLoop(mg, aPosition, aPosition);

	delete refCallback;
	return true;
}

bool CreateSmoothVolumeHierarchy(MultiGrid& mg, size_t numRefs, bool bPreserveBnd, bool bSubdivisionLoopBnd)
{
	PROFILE_FUNC_GROUP("grid");

	GlobalMultiGridRefiner ref(mg);

	for(size_t lvl = 0; lvl < numRefs; ++lvl){
		ref.refine();
		SubdivisionVolumes(mg, bPreserveBnd, bSubdivisionLoopBnd);
	}

	if(bSubdivisionLoopBnd)
		ProjectToLimitPLoop(mg, aPosition, aPosition);

	ProjectToLimitSubdivisionVolume(mg);

	return true;
}

bool CreateSemiSmoothHierarchy(MultiGrid& mg, size_t numRefs)
{
	PROFILE_FUNC_GROUP("grid");
	IRefinementCallback* refCallback = NULL;
//	we're only checking for the main attachments here.
//todo: improve this - add a domain-based hierarchy creator.
	if(mg.has_vertex_attachment(aPosition1))
		refCallback = new SubdivisionLoopBoundaryProjector<APosition1>(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new SubdivisionLoopBoundaryProjector<APosition2>(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new SubdivisionLoopBoundaryProjector<APosition>(mg, aPosition, aPosition);
		
	if(!refCallback){
		UG_LOG("No standard position attachment found. Aborting.\n");
		return false;
	}
	
	GlobalMultiGridRefiner ref(mg, refCallback);

	for(size_t lvl = 0; lvl < numRefs; ++lvl){
		ref.refine();
	}

	if(mg.has_vertex_attachment(aPosition1))
		ProjectToLimitSubdivBoundary(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		ProjectToLimitSubdivBoundary(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		ProjectToLimitSubdivBoundary(mg, aPosition, aPosition);

	delete refCallback;
	return true;
}


void RegisterGridBridge_Refinement(Registry& reg, string parentGroup)
{
	string grp = parentGroup;

	reg.add_class_<IRefinementCallback>("IRefinementCallback", grp);

//	IRefiner
	reg.add_class_<IRefiner>("IRefiner", grp)
		.add_method("refine", &IRefiner::refine)
		.add_method("coarsen", &IRefiner::coarsen)
		.add_method("save_marks_to_file", &IRefiner::save_marks_to_file, "", "filename")
		.add_method("set_adjusted_marks_debug_filename", &IRefiner::set_adjusted_marks_debug_filename, "", "filename")
		.add_method("mark_neighborhood",
					static_cast<void (IRefiner::*)(size_t)>(&IRefiner::mark_neighborhood),
					"", "numIterations")
		.add_method("clear_marks", &IRefiner::clear_marks)
		.add_method("set_refinement_callback", &IRefiner::set_refinement_callback)
		.add_method("enable_debugging", &IRefiner::enable_debugging)
		.add_method("num_marked_edges", static_cast<size_t (IRefiner::*)()>(&IRefiner::num_marked_edges))
		.add_method("num_marked_faces", static_cast<size_t (IRefiner::*)()>(&IRefiner::num_marked_faces))
		.add_method("num_marked_volumes", static_cast<size_t (IRefiner::*)()>(&IRefiner::num_marked_volumes))
		.add_method("num_marked_elements", static_cast<size_t (IRefiner::*)()>(&IRefiner::num_marked_elements));

//	HangingNodeRefiner
	reg.add_class_<HangingNodeRefiner_Grid, IRefiner>("HangingNodeRefiner_Grid", grp)
		.add_constructor()
		.add_method("assign_grid", &HangingNodeRefiner_Grid::assign_grid, "", "g")
		.set_construct_as_smart_pointer(true);

	reg.add_class_<HangingNodeRefiner_MultiGrid, IRefiner>("HangingNodeRefiner_MultiGrid", grp)
		.add_constructor()
		.add_method("assign_grid", &HangingNodeRefiner_MultiGrid::assign_grid, "", "mg")
		.set_construct_as_smart_pointer(true);

//	AdaptiveRegularMGRefiner
	reg.add_class_<AdaptiveRegularRefiner_MultiGrid, HangingNodeRefiner_MultiGrid>("AdaptiveRegularRefiner_MultiGrid", grp)
		.add_constructor()
		.add_method("assign_grid", &AdaptiveRegularRefiner_MultiGrid::assign_grid, "", "mg")
		.set_construct_as_smart_pointer(true);

//	GlobalMultiGridRefiner
	reg.add_class_<GlobalMultiGridRefiner, IRefiner>("GlobalMultiGridRefiner", grp)
		.add_constructor()
		.add_method("assign_grid", static_cast<void (GlobalMultiGridRefiner::*)(MultiGrid&)>(&GlobalMultiGridRefiner::assign_grid),
				"", "mg")
		.set_construct_as_smart_pointer(true);

//	FracturedMediaRefiner
	/*typedef FracturedMediaRefiner<typename TDomain::grid_type,
						  	  	  typename TDomain::position_attachment_type>	FracDomRef;
	reg.add_class_<FracDomRef, IRefiner>("FracturedMediumRefiner", grp)
		.add_constructor()
		.add_method("set_aspect_ratio_threshold", &FracDomRef::set_aspect_ratio_threshold);*/

//	GlobalFracturedDomainRefiner
	{
		typedef GlobalFracturedMediaRefiner cls;
		reg.add_class_<cls, IRefiner>("GlobalFracturedMediumRefiner", grp)
			.add_constructor()
			.add_method("assign_grid", static_cast<void (cls::*)(MultiGrid*)>(&cls::assign_grid), "", "g")
			.add_method("set_subset_handler", static_cast<void (cls::*)(ISubsetHandler*)>(&cls::set_subset_handler),
					"", "sh")
			.add_method("mark_as_fracture", &cls::mark_as_fracture, "", "subInd#bIsFracture")
			.add_method("is_fracture", &cls::is_fracture, "", "subInd")
			.set_construct_as_smart_pointer(true);
	}

//	parallel refinement
#ifdef UG_PARALLEL
	reg.add_class_<ParallelHangingNodeRefiner_MultiGrid, HangingNodeRefiner_MultiGrid>
		("ParallelHangingNodeRefiner_MultiGrid", grp)
		.add_constructor()
		.set_construct_as_smart_pointer(true);
/*Currently not directly usable. For domains, you may use the factory method
* GlobalFracturedDomainRefiner, which automatically creates a
* ParallelGlobalFracturedMediaRefiner, if required.
	reg.add_class_<ParallelGlobalFracturedMediaRefiner, GlobalFracturedMediaRefiner>
		("ParallelGlobalFracturedMediaRefiner", grp)
		.add_constructor()
		.set_construct_as_smart_pointer(true);
*/
#endif

//	refinement
	reg.add_function("TestSubdivision", &TestSubdivision, grp)
		.add_function("CreateHierarchy", &CreateHierarchy, grp)
		.add_function("CreateSmoothHierarchy", &CreateSmoothHierarchy, grp)
		.add_function("CreateSmoothVolumeHierarchy", &CreateSmoothVolumeHierarchy, grp)
		.add_function("CreateSemiSmoothHierarchy", &CreateSemiSmoothHierarchy, grp);

//	smooth volume subdivision
	reg.add_function("SubdivisionVolumes", &SubdivisionVolumes, grp);
}

}//	end of namespace
}//	end of namespace
