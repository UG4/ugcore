//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../registry.h"
#include "../ug_bridge.h"
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"
#include <iostream>
#include <sstream>

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
		inline Grid& get_grid()	{return *this;}
		inline SubsetHandler& get_subset_handler()	{return m_sh;}
		
	protected:
		SubsetHandler m_sh;
};

bool LoadGridObject(GridObject& go, const char* filename)
{
	return LoadGridFromFile(go.get_grid(), filename, go.get_subset_handler());
}

bool SaveGridObject(GridObject& go, const char* filename)
{
	return SaveGridToFile(go.get_grid(), filename, go.get_subset_handler());
}

GridObject* CreateGridObject(const char* filename)
{
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
	PROFILE_FUNC();
//	HangingNodeRefiner_IR1 href(grid);
	return CreateFractal_NormalScale(grid, href, scaleFac, numIterations);
//	return true;
}







bool LoadGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	return LoadGridFromFile(grid, filename, sh);
}

bool SaveGrid(Grid& grid, SubsetHandler& sh, const char* filename)
{
	return SaveGridToFile(grid, filename, sh);
}

bool SaveGridHierarchy(MultiGrid& mg, const char* filename)
{
	return SaveGridToFile(mg, filename, mg.get_hierarchy_handler());
}


void TestSubdivision(const char* fileIn, const char* fileOut, int numRefs)
{
//todo: Callbacks have to make sure that their attachment is accessible in the grid.
//		even if they were initialized before the attachment was attached to the grid.
	MultiGrid mg;
	SubsetHandler sh(mg);
	RefinementCallbackSubdivisionLoop<APosition> refCallback(mg, aPosition, aPosition);
	GlobalMultiGridRefiner ref(mg, &refCallback);
	
	if(LoadGridFromFile(mg, fileIn, sh)){
		for(int lvl = 0; lvl < numRefs; ++lvl){
			ref.refine();
		}

		ProjectToLimitPLoop(mg, aPosition, aPosition);
		SaveGridToFile(mg, fileOut, mg.get_hierarchy_handler());

	}
	else{
		UG_LOG("Load failed. aborting...\n");
	}
}

bool CreateSmoothHierarchy(MultiGrid& mg, size_t numRefs)
{
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
					  HangingNodeRefiner_MultiGrid& refiner,
					  float percentage)
{
/*
	typedef typename geometry_traits<TElem>::iterator iterator;
	for(iterator iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter)
	{
		if(urand<float>(0, 99) < percentage){
			refiner.mark_for_refinement(*iter);
		}
	}
*/

	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	TElem* elem = FindByCoordinate<TElem>(vector3(-0.00001, -0.00001, -0.00001),
									mg.begin<TElem>(mg.num_levels()-1),
									mg.end<TElem>(mg.num_levels()-1),
									aaPos);

	if(elem)
		refiner.mark_for_refinement(elem);
	else{
		UG_LOG("No element found for refinement.\n");
	}


}

bool TestHangingNodeRefiner_MultiGrid(const char* filename,
									  const char* outFilename,
									  int numIterations,
									  float percentage)
{
	MultiGrid mg;
	MGSubsetHandler sh(mg);
	HangingNodeRefiner_MultiGrid refiner(mg);

	if(!LoadGridFromFile(mg, filename, sh)){
		UG_LOG("  could not load " << filename << endl);
		return false;
	}

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

	UG_LOG("saving to " << outFilename << endl;)
	SaveGridHierarchy(mg, outFilename);

	UG_LOG("grid element numbers:\n");
	PrintGridElementNumbers(mg);

//	create a surface view
	SurfaceView surfView(mg);

	if(mg.num<Volume>() > 0)
		CreateSurfaceView<Volume>(surfView, mg, sh);
	else if(mg.num<Face>() > 0)
		CreateSurfaceView<Face>(surfView, mg, sh);
	else
		CreateSurfaceView<EdgeBase>(surfView, mg, sh);

	SaveGridToFile(mg, "surface_view.ugx", surfView);

	UG_LOG("surface view element numbers:\n");
	PrintGridElementNumbers(surfView);

	return true;
}


template <class TAPos, class vector_t>
void MarkForRefinement_VerticesInSphere(Grid& grid, IRefiner& refiner,
										const vector_t& center,
										number radius, TAPos& aPos)
{
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("WARNING in MarkForRefinement_VerticesInSphere: position attachment missing.\n");
		return;
	}

	Grid::VertexAttachmentAccessor<TAPos> aaPos(grid, aPos);

//	we'll compare against the square radius.
	number radiusSq = radius * radius;

//	we'll store associated edges, faces and volumes in those containers
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all vertices of the grid. If a vertex is inside the given sphere,
//	then we'll mark all associated elements.
	for(VertexBaseIterator iter = grid.begin<VertexBase>();
		iter != grid.end<VertexBase>(); ++iter)
	{
		if(VecDistanceSq(center, aaPos[*iter]) <= radiusSq){
			CollectAssociated(vEdges, grid, *iter);
			CollectAssociated(vFaces, grid, *iter);
			CollectAssociated(vVols, grid, *iter);

			refiner.mark_for_refinement(vEdges.begin(), vEdges.end());
			refiner.mark_for_refinement(vFaces.begin(), vFaces.end());
			refiner.mark_for_refinement(vVols.begin(), vVols.end());
		}
	}
}

////////////////////////////////////////////////////////////////////////
void MarkForRefinement_VerticesInSphere(IRefiner& refiner, number centerX,
							number centerY, number centerZ, number radius)
{
//	get the associated grid of the refiner.
	Grid* pGrid = refiner.get_associated_grid();
	if(!pGrid){
		UG_LOG("WARNING: Circular refinement failed: no grid assigned.\n");
		return;
	}

//	depending on the position attachment, we'll call different versions of
//	MarkForRefinement_VerticesInSphere.
	if(pGrid->has_vertex_attachment(aPosition1)){
		MarkForRefinement_VerticesInSphere(*pGrid, refiner, vector1(centerX),
											radius, aPosition1);
	}
	else if(pGrid->has_vertex_attachment(aPosition2)){
		MarkForRefinement_VerticesInSphere(*pGrid, refiner,
											vector2(centerX, centerY),
											radius, aPosition2);
	}
	else if(pGrid->has_vertex_attachment(aPosition)){
		MarkForRefinement_VerticesInSphere(*pGrid, refiner,
											vector3(centerX, centerY, centerZ),
											radius, aPosition);
	}
	else{
		UG_LOG("WARNING in MarkForRefinement_VerticesInSphere: No Position attachment found. Aborting.\n");
	}
}


////////////////////////////////////////////////////////////////////////
bool RegisterLibGridInterface(Registry& reg, const char* parentGroup)
{
	try
	{
//	get group string
	std::stringstream groupString; groupString << parentGroup << "/Grid";
	std::string grp = groupString.str();

//	Grid
	reg.add_class_<Grid>("Grid", grp.c_str())
		.add_constructor()
		.add_method("clear", &Grid::clear)
		.add_method("num_vertices", &Grid::num_vertices)
		.add_method("num_edges", &Grid::num_edges)
		.add_method("num_faces", &Grid::num_faces)
		.add_method("num_volumes", &Grid::num_volumes);
		
//	MultiGrid
	reg.add_class_<MultiGrid, Grid>("MultiGrid", grp.c_str())
		.add_constructor()
		.add_method("num_vertices_on_level", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<VertexBase>)
		.add_method("num_edges_on_level", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<EdgeBase>)
		.add_method("num_faces_on_level", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Face>)
		.add_method("num_volumes_on_level", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Volume>);


////////////////////////
//	SUBSET HANDLERS

//  ISubsetHandler
	reg.add_class_<ISubsetHandler>("ISubsetHandler", grp.c_str())
		.add_method("num_subsets", &ISubsetHandler::num_subsets)
		.add_method("get_subset_name", &ISubsetHandler::get_subset_name)
		.add_method("set_subset_name", &ISubsetHandler::set_subset_name);
	
//	SubsetHandler
	reg.add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler", grp.c_str())
		.add_constructor()
		.add_method("assign_grid", &SubsetHandler::assign_grid);

//	MGSubsetHandler
	reg.add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler", grp.c_str())
		.add_constructor()
		.add_method("assign_grid", &MGSubsetHandler::assign_grid);

//	SurfaceView

	#ifndef FOR_VRL

	reg.add_class_<SurfaceView, SubsetHandler>("SurfaceView", grp.c_str())
		.add_constructor()
		.add_method("assign_grid", (void (SurfaceView::*)(MultiGrid&)) &SurfaceView::assign_grid);

	#endif

////////////////////////
//	REFINEMENT

//	IRefiner
	reg.add_class_<IRefiner>("IRefiner", grp.c_str())
		.add_method("refine", &IRefiner::refine);

//	HangingNodeRefiner
	reg.add_class_<HangingNodeRefiner_Grid, IRefiner>("HangingNodeRefiner_Grid", grp.c_str())
		.add_constructor()
		.add_method("assign_grid", &HangingNodeRefiner_Grid::assign_grid);

	reg.add_class_<HangingNodeRefiner_MultiGrid, IRefiner>("HangingNodeRefiner_MultiGrid", grp.c_str())
		.add_constructor()
		.add_method("assign_grid", &HangingNodeRefiner_MultiGrid::assign_grid);

//	GlobalMultiGridRefiner
	reg.add_class_<GlobalMultiGridRefiner, IRefiner>("GlobalMultiGridRefiner", grp.c_str())
		.add_constructor()
		.add_method("assign_grid", (void (GlobalMultiGridRefiner::*)(MultiGrid&)) &GlobalMultiGridRefiner::assign_grid);

//	parallel refinement
#ifdef UG_PARALLEL
	reg.add_class_<ParallelHangingNodeRefiner_MultiGrid, HangingNodeRefiner_MultiGrid>
		("ParallelHangingNodeRefiner_MultiGrid", grp.c_str())
		.add_constructor();
#endif

//	GridObject
	reg.add_class_<GridObject, Grid>("GridObject", grp.c_str())
		.add_constructor()
		.add_method("get_grid", &GridObject::get_grid)
		.add_method("get_subset_handler", &GridObject::get_subset_handler);

//	Grid functions
	reg.add_function("CreateFractal", &CreateFractal, grp.c_str());
	
//  GridObject functions
	reg.add_function("LoadGrid", &LoadGrid, grp.c_str())
		.add_function("SaveGrid", &SaveGrid, grp.c_str())
		.add_function("LoadGridObject", &LoadGridObject, grp.c_str())
		.add_function("SaveGridObject", &SaveGridObject, grp.c_str())
		.add_function("CreateGridObject", &CreateGridObject, grp.c_str());
		
//	refinement
	reg.add_function("TestSubdivision", &TestSubdivision)
		.add_function("TestHangingNodeRefiner_MultiGrid", &TestHangingNodeRefiner_MultiGrid)
		.add_function("CreateSmoothHierarchy", &CreateSmoothHierarchy, grp.c_str())
		.add_function("CreateSemiSmoothHierarchy", &CreateSemiSmoothHierarchy, grp.c_str())
		.add_function("SaveGridHierarchy", &SaveGridHierarchy, grp.c_str())
		.add_function("MarkForRefinement_VerticesInSphere", (void (*)(IRefiner&, number, number, number, number))
															&MarkForRefinement_VerticesInSphere, grp.c_str());
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibGridInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

}//	end of namespace 
}//	end of namespace 
