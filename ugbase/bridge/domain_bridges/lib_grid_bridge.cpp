//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include <iostream>
#include <sstream>
#include "registry/registry.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/tools/surface_view.h"
#include "common/profiler/profiler.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/tools/partition_map.h"
//todo: include this in algorithms.h
#include "lib_grid/algorithms/refinement/global_fractured_media_refiner.h"
#include "lib_grid/algorithms/refinement/adaptive_regular_mg_refiner.h"
#include "lib_grid/parallelization/util/partition_weighting_callbacks.h"

using namespace std;

namespace ug
{
namespace bridge
{
///	Wrapper object that simplifies script creation
class GridObject : public Grid
{
	public:
		GridObject() : Grid(GRIDOPT_STANDARD_INTERCONNECTION), m_sh(*this)	{}
		inline Grid& grid()	{return *this;}
		inline SubsetHandler& subset_handler()	{return m_sh;}
		
	protected:
		SubsetHandler m_sh;
};

bool LoadGridObject(GridObject& go, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return LoadGridFromFile(go.grid(), go.subset_handler(), filename);
}

bool SaveGridObject(GridObject& go, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(go.grid(), go.subset_handler(), filename);
}

GridObject* CreateGridObject(const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	GridObject* go = new GridObject;
	if(!LoadGridObject(*go, filename)){
		delete go;
		return NULL;
	}
	return go;
}

bool CreateFractal(Grid& grid, HangingNodeRefiner_Grid& href,
					number scaleFac, size_t numIterations)
{
	PROFILE_FUNC_GROUP("grid");
//	HangingNodeRefiner_IR1 href(grid);
	return CreateFractal_NormalScale(grid, href, scaleFac, numIterations);
//	return true;
}

bool LoadGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return LoadGridFromFile(grid, sh, filename);
}

bool LoadGrid(Grid& grid, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return LoadGridFromFile(grid, filename);
}

bool SaveGrid(Grid& grid, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(grid, filename);
}

bool SaveGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(grid, sh, filename);
}

bool SaveGrid(Grid& grid, const ISubsetHandler& sh, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(grid, *const_cast<ISubsetHandler*>(&sh), filename);
}

bool SaveGridHierarchy(MultiGrid& mg, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(mg, mg.get_hierarchy_handler(), filename);
}


void TestSubdivision(const char* fileIn, const char* fileOut, int numRefs)
{
	PROFILE_FUNC_GROUP("grid");
//todo: Callbacks have to make sure that their attachment is accessible in the grid.
//		even if they were initialized before the attachment was attached to the grid.
	MultiGrid mg;
	SubsetHandler sh(mg);
	RefinementCallbackSubdivisionLoop<APosition> refCallback(mg, aPosition, aPosition);
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
		refCallback = new RefinementCallbackSubdivisionLoop<APosition1>(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition2>(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition>(mg, aPosition, aPosition);
		
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

bool CreateSemiSmoothHierarchy(MultiGrid& mg, size_t numRefs)
{
	PROFILE_FUNC_GROUP("grid");
	IRefinementCallback* refCallback = NULL;
//	we're only checking for the main attachments here.
//todo: improve this - add a domain-based hierarchy creator.
	if(mg.has_vertex_attachment(aPosition1))
		refCallback = new RefinementCallbackSubdivBoundary<APosition1>(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new RefinementCallbackSubdivBoundary<APosition2>(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new RefinementCallbackSubdivBoundary<APosition>(mg, aPosition, aPosition);
		
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

template <class TElem>
void MarkForRefinement(MultiGrid& mg,
					  IRefiner& refiner,
					  float percentage)
{
	PROFILE_FUNC_GROUP("grid");
	typedef typename geometry_traits<TElem>::iterator iterator;
	for(iterator iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter)
	{
		if(urand<float>(0, 99) < percentage){
			refiner.mark(*iter);
		}
	}

/*
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	TElem* elem = FindByCoordinate<TElem>(vector3(-0.00001, -0.00001, -0.00001),
									mg.begin<TElem>(mg.num_levels()-1),
									mg.end<TElem>(mg.num_levels()-1),
									aaPos);

	if(elem)
		refiner.mark(elem);
	else{
		UG_LOG("No element found for refinement.\n");
	}*/


}

bool TestHangingNodeRefiner_MultiGrid(const char* filename,
									  const char* outFilename,
									  int numIterations,
									  float percentage)
{
	PROFILE_FUNC_GROUP("grid");
	MultiGrid mg;
	MGSubsetHandler sh(mg);
	HangingNodeRefiner_MultiGrid refiner(mg);

	PROFILE_BEGIN(PROFTEST_loading);
	if(!LoadGridFromFile(mg, sh, filename)){
		UG_LOG("  could not load " << filename << endl);
		return false;
	}
	PROFILE_END();

	{
		PROFILE_BEGIN(PROFTEST_refining);
		for(int i = 0; i < numIterations; ++i){
			UG_LOG("refinement step " << i+1 << endl);

			if(mg.num<Volume>() > 0)
				MarkForRefinement<Volume>(mg, refiner, percentage);
			else if(mg.num<Face>() > 0)
				MarkForRefinement<Face>(mg, refiner, percentage);
			else
				MarkForRefinement<EdgeBase>(mg, refiner, percentage);

			refiner.refine();
		}
	}

	UG_LOG("saving to " << outFilename << endl);
	SaveGridHierarchy(mg, outFilename);

	UG_LOG("grid element numbers:\n");
	PrintGridElementNumbers(mg);

	return true;
}


////////////////////////////////////////////////////////////////////////
///	A helper class for ExpandLayers.
/**	This class should never be publicly available, especially since
 * deriving from std::vector is a bad idea (compare 'Effective C++').
 * However, it is very useful in this situation.
 *
 * The class simply extends std::vector<FractureInfo> by an add_layer method.
 */
class ExpandLayersDesc : public std::vector<FractureInfo>
{
	public:
		ExpandLayersDesc() {}

		void add_layer(int subsetInd, int newSubsetInd, number width)
		{
			push_back(FractureInfo(subsetInd, newSubsetInd, width));
		}
};


////////////////////////////////////////////////////////////////////////
void RegisterBridge_Grid(Registry& reg, string parentGroup)
{
//	get group string
	stringstream groupString; groupString << parentGroup << "/Grid";
	string grp = groupString.str();
	try{
	//	Grid
		reg.add_class_<Grid>("Grid", grp)
			.add_constructor()
			.add_method("clear", static_cast<void (Grid::*)()>(&Grid::clear))
			.add_method("clear_geometry", &Grid::clear_geometry)
			.add_method("num_vertices", &Grid::num_vertices)
			.add_method("num_edges", &Grid::num_edges)
			.add_method("num_faces", &Grid::num_faces)
			.add_method("num_triangles", &Grid::num<Triangle>)
			.add_method("num_quadrilaterals", &Grid::num<Quadrilateral>)
			.add_method("num_volumes", &Grid::num_volumes)
			.add_method("num_tetrahedrons", &Grid::num<Tetrahedron>)
			.add_method("num_pyramids", &Grid::num<Pyramid>)
			.add_method("num_prisms", &Grid::num<Prism>)
			.add_method("num_hexahedrons", &Grid::num<Hexahedron>)
			.add_method("reserve_vertices", &Grid::reserve<VertexBase>)
			.add_method("reserve_edges", &Grid::reserve<EdgeBase>)
			.add_method("reserve_faces", &Grid::reserve<Face>)
			.add_method("reserve_volumes", &Grid::reserve<Volume>)
			.set_construct_as_smart_pointer(true);

	//	MultiGrid
		reg.add_class_<MultiGrid, Grid>("MultiGrid", grp)
			.add_constructor()
			.add_method("num_levels", &MultiGrid::num_levels)

			.add_method("num_vertices", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<VertexBase>)
			.add_method("num_edges", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<EdgeBase>)
			.add_method("num_faces", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Face>)
			.add_method("num_triangles", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Triangle>)
			.add_method("num_quadrilaterals", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Quadrilateral>)
			.add_method("num_volumes", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Volume>)
			.add_method("num_tetrahedrons", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Tetrahedron>)
			.add_method("num_pyramids", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Pyramid>)
			.add_method("num_prisms", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Prism>)
			.add_method("num_hexahedrons", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Hexahedron>)
			.set_construct_as_smart_pointer(true);



	////////////////////////
	//	SUBSET HANDLERS

	//  ISubsetHandler
		reg.add_class_<ISubsetHandler>("ISubsetHandler", grp)
			.add_method("num_subsets", &ISubsetHandler::num_subsets)
			.add_method("get_subset_name", &ISubsetHandler::get_subset_name)
			.add_method("set_subset_name", &ISubsetHandler::set_subset_name)
			.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(const char*) const>(
											&ISubsetHandler::get_subset_index));
		
	//	SubsetHandler
		reg.add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler", grp)
			.add_constructor()
			.add_method("assign_grid", static_cast<void (SubsetHandler::*)(Grid&)>(&SubsetHandler::assign_grid))
			.set_construct_as_smart_pointer(true);


	//	MGSubsetHandler
		reg.add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler", grp)
			.add_constructor()
			.add_method("assign_grid", &MGSubsetHandler::assign_grid)
			.set_construct_as_smart_pointer(true);

	//	SurfaceView
		reg.add_class_<SurfaceView>("SurfaceView", grp)
			.add_method("subset_handler", static_cast<ConstSmartPtr<MGSubsetHandler> (SurfaceView::*)() const>(
											&SurfaceView::subset_handler));

	//	Selector
		reg.add_class_<ISelector>("ISelector", grp);

		reg.add_class_<Selector, ISelector>("Selector", grp)
			.add_constructor<void (*)(Grid&)>()

			.add_method("num_vertices", static_cast<size_t (Selector::*)() const>(&Selector::num<VertexBase>))
			.add_method("num_edges", static_cast<size_t (Selector::*)() const>(&Selector::num<EdgeBase>))
			.add_method("num_faces", static_cast<size_t (Selector::*)() const>(&Selector::num<Face>))
			.add_method("num_triangles", static_cast<size_t (Selector::*)() const>(&Selector::num<Triangle>))
			.add_method("num_quadrilaterals", static_cast<size_t (Selector::*)() const>(&Selector::num<Quadrilateral>))
			.add_method("num_volumes", static_cast<size_t (Selector::*)() const>(&Selector::num<Volume>))
			.add_method("num_tetrahedrons", static_cast<size_t (Selector::*)() const>(&Selector::num<Tetrahedron>))
			.add_method("num_pyramids", static_cast<size_t (Selector::*)() const>(&Selector::num<Pyramid>))
			.add_method("num_prisms", static_cast<size_t (Selector::*)() const>(&Selector::num<Prism>))
			.add_method("num_hexahedrons", static_cast<size_t (Selector::*)() const>(&Selector::num<Hexahedron>))
			.set_construct_as_smart_pointer(true);

	////////////////////////
	//	REFINEMENT

	//	IRefiner
		reg.add_class_<IRefiner>("IRefiner", grp)
			.add_method("refine", &IRefiner::refine)
			.add_method("coarsen", &IRefiner::coarsen)
			.add_method("save_marks_to_file", &IRefiner::save_marks_to_file)
			.add_method("set_adjusted_marks_debug_filename", &IRefiner::set_adjusted_marks_debug_filename)
			.add_method("clear_marks", &IRefiner::clear_marks);

	//	HangingNodeRefiner
		reg.add_class_<HangingNodeRefiner_Grid, IRefiner>("HangingNodeRefiner_Grid", grp)
			.add_constructor()
			.add_method("assign_grid", &HangingNodeRefiner_Grid::assign_grid)
			.set_construct_as_smart_pointer(true);

		reg.add_class_<HangingNodeRefiner_MultiGrid, IRefiner>("HangingNodeRefiner_MultiGrid", grp)
			.add_constructor()
			.add_method("assign_grid", &HangingNodeRefiner_MultiGrid::assign_grid)
			.set_construct_as_smart_pointer(true);

	//	AdaptiveRegularMGRefiner
		reg.add_class_<AdaptiveRegularRefiner_MultiGrid, HangingNodeRefiner_MultiGrid>("AdaptiveRegularRefiner_MultiGrid", grp)
			.add_constructor()
			.add_method("assign_grid", &AdaptiveRegularRefiner_MultiGrid::assign_grid)
			.set_construct_as_smart_pointer(true);

	//	GlobalMultiGridRefiner
		reg.add_class_<GlobalMultiGridRefiner, IRefiner>("GlobalMultiGridRefiner", grp)
			.add_constructor()
			.add_method("assign_grid", static_cast<void (GlobalMultiGridRefiner::*)(MultiGrid&)>(&GlobalMultiGridRefiner::assign_grid))
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
				.add_method("assign_grid", static_cast<void (cls::*)(MultiGrid*)>(&cls::assign_grid))
				.add_method("set_subset_handler", static_cast<void (cls::*)(ISubsetHandler*)>(&cls::set_subset_handler))
				.add_method("mark_as_fracture", &cls::mark_as_fracture)
				.add_method("is_fracture", &cls::is_fracture)
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

	// partition weighting in metis partitioning
	reg.add_class_<PartitionWeighting>("PartitionWeighting", grp)
			.add_constructor()
			.add_method("set_default_weights", &PartitionWeighting::set_default_weights)
			.set_construct_as_smart_pointer(true);
	reg.add_class_<InterSubsetPartitionWeighting, PartitionWeighting>("InterSubsetPartitionWeighting", grp)
			.add_constructor()
			.add_method("set_inter_subset_weight", &InterSubsetPartitionWeighting::set_inter_subset_weight)
			.set_construct_as_smart_pointer(true);

	//	GridObject
		reg.add_class_<GridObject, Grid>("GridObject", grp)
			.add_constructor()
			.add_method("grid", &GridObject::grid)
			.add_method("subset_handler", &GridObject::subset_handler)
			.set_construct_as_smart_pointer(true);

	//	Grid functions
		reg.add_function("CreateFractal", &CreateFractal, grp)
			.add_function("PrintAttachmentInfo", &PrintAttachmentInfo, grp);

	//	UGXFileInfo
		reg.add_class_<UGXFileInfo>("UGXFileInfo", grp)
			.add_constructor()
			.add_method("parse_file", &UGXFileInfo::parse_file)
			.add_method("num_grids", &UGXFileInfo::num_grids)
			.add_method("num_subset_handlers", &UGXFileInfo::num_subset_handlers)
			.add_method("num_subsets", &UGXFileInfo::num_subsets)
			.add_method("grid_name", &UGXFileInfo::grid_name)
			.add_method("subset_handler_name", &UGXFileInfo::subset_handler_name)
			.add_method("subset_name", &UGXFileInfo::subset_name)
			.add_method("grid_world_dimension", &UGXFileInfo::grid_world_dimension)
			.set_construct_as_smart_pointer(true);

	//  GridObject functions
		reg.add_function("LoadGrid", static_cast<bool (*)(Grid&, ISubsetHandler&, const char*)>(&LoadGrid), grp)
			.add_function("LoadGrid", static_cast<bool (*)(Grid&, const char*)>(&LoadGrid), grp)
			.add_function("SaveGrid", static_cast<bool (*)(Grid&, const ISubsetHandler&, const char*)>(&SaveGrid), grp)
			.add_function("SaveGrid", static_cast<bool (*)(Grid&, ISubsetHandler&, const char*)>(&SaveGrid), grp)
			.add_function("SaveGrid", static_cast<bool (*)(Grid&, const char*)>(&SaveGrid), grp)
			.add_function("LoadGridObject", &LoadGridObject, grp)
			.add_function("SaveGridObject", &SaveGridObject, grp)
			.add_function("SaveGridHierarchy", &SaveGridHierarchy, grp)
			.add_function("SaveGridHierarchyTransformed",
						  static_cast<bool (*)(MultiGrid&, ISubsetHandler&, const char*, number)>(
								  &SaveGridHierarchyTransformed),
						  grp)
			.add_function("SaveGridHierarchyTransformed",
						  static_cast<bool (*)(MultiGrid&, const char*, number)>(
								  &SaveGridHierarchyTransformed),
						  grp)
			.add_function("SaveParallelGridLayout", &SaveParallelGridLayout)
			.add_function("SaveSurfaceViewTransformed", &SaveSurfaceViewTransformed)
			.add_function("CreateGridObject", &CreateGridObject, grp)
			.add_function("PrintGridElementNumbers", static_cast<void (*)(MultiGrid&)>(&PrintGridElementNumbers), grp)
			.add_function("PrintGridElementNumbers", static_cast<void (*)(Grid&)>(&PrintGridElementNumbers), grp);

	//	refinement
		reg.add_function("TestSubdivision", &TestSubdivision, grp)
			.add_function("TestHangingNodeRefiner_MultiGrid", &TestHangingNodeRefiner_MultiGrid, grp)
			.add_function("CreateHierarchy", &CreateHierarchy, grp)
			.add_function("CreateSmoothHierarchy", &CreateSmoothHierarchy, grp)
			.add_function("CreateSemiSmoothHierarchy", &CreateSemiSmoothHierarchy, grp);
		
	//	subset util
		reg.add_function("AdjustSubsetsForSimulation",
						static_cast<void (*)(SubsetHandler&, bool)>(
						&AdjustSubsetsForSimulation<SubsetHandler>), grp)
			.add_function("AdjustSubsetsForSimulation",
						static_cast<void (*)(MGSubsetHandler&, bool)>(
						&AdjustSubsetsForSimulation<MGSubsetHandler>), grp)
			.add_function("AssignSubsetsByElementType", &AssignSubsetsByElementType)
			.add_function("AssignSubsetColors", &AssignSubsetColors);

	//	PartitionMap
		reg.add_class_<PartitionMap>("PartitionMap", grp)
			.add_constructor()
			.add_method("clear", &PartitionMap::clear)
			.add_method("get_partition_handler", &PartitionMap::get_partition_handler)
			.add_method("add_target_proc", &PartitionMap::add_target_proc)
			.add_method("add_target_procs", &PartitionMap::add_target_procs)
			.add_method("num_target_procs", &PartitionMap::num_target_procs)
			.add_method("get_target_proc", &PartitionMap::get_target_proc)
			.add_method("shift_target_procs", &PartitionMap::shift_target_procs)
			.set_construct_as_smart_pointer(true);

	//	ExpandLayers
		typedef std::vector<FractureInfo> FracInfoVec;
		reg.add_class_<FracInfoVec>("FractureInfoVec", grp);

		reg.add_class_<ExpandLayersDesc, FracInfoVec>("ExpandLayersDesc", grp)
			.add_constructor()
			.add_method("add_layer", &ExpandLayersDesc::add_layer)
			.set_construct_as_smart_pointer(true);

		// \todo: this is uncommented, since in conflict with new std::vector
		//	handling of registry. Should be adapted.
//		reg.add_function("ExpandLayers2d", &ExpandFractures2d, grp)
//			.add_function("ExpandLayers3d", &ExpandFractures3d, grp);
			
	//	Debugging
		reg.add_function("CheckHangingNodeConsistency", static_cast<bool (*)(MultiGrid&)>(&CheckHangingNodeConsistency), grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace 
}//	end of namespace 
