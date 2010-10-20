//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../registry.h"
#include "../ug_bridge.h"
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"
#include <iostream>
#include <sstream>

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

bool CreateFractal(Grid& grid, HangingNodeRefiner_IR1& href,
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
	RefinementCallbackSubdivisionLoop<APosition> refCallback(mg, aPosition);
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
		refCallback = new RefinementCallbackSubdivisionLoop<APosition1>(mg, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition2>(mg, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition>(mg, aPosition);
		
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
		refCallback = new RefinementCallbackSubdivBoundary<APosition1>(mg, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new RefinementCallbackSubdivBoundary<APosition2>(mg, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new RefinementCallbackSubdivBoundary<APosition>(mg, aPosition);
		
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


////////////////////////////////////////////////////////////////////////
void RegisterLibGridInterface(Registry& reg, const char* parentGroup)
{
//	get group string
	std::stringstream groupString; groupString << parentGroup << "/Grid";
	const char* grp = groupString.str().c_str();

//	Grid
	reg.add_class_<Grid>("Grid", grp)
		.add_constructor()
		.add_method("clear", &Grid::clear)
		.add_method("num_vertices", &Grid::num_vertices)
		.add_method("num_edges", &Grid::num_edges)
		.add_method("num_faces", &Grid::num_faces)
		.add_method("num_volumes", &Grid::num_volumes);
		
//	MultiGrid
	reg.add_class_<MultiGrid, Grid>("MultiGrid", grp)
		.add_constructor();

//  ISubsetHandler
	reg.add_class_<ISubsetHandler>("ISubsetHandler", grp)
		.add_method("num_subsets", &ISubsetHandler::num_subsets)
		.add_method("get_subset_name", &ISubsetHandler::get_subset_name)
		.add_method("set_subset_name", &ISubsetHandler::set_subset_name);
	
//	SubsetHandler
	reg.add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler", grp)
		.add_constructor()
		.add_method("assign_grid", &SubsetHandler::assign_grid);

//	MGSubsetHandler
	reg.add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler", grp)
		.add_constructor()
		.add_method("assign_grid", &MGSubsetHandler::assign_grid);

//	HangingNodeRefiner
	reg.add_class_<HangingNodeRefiner_IR1>("HangingNodeRefiner", grp)
		.add_constructor()
		.add_method("assign_grid", &HangingNodeRefiner_IR1::assign_grid);

//	GlobalMultiGridRefiner
	reg.add_class_<GlobalMultiGridRefiner>("GlobalMultiGridRefiner", grp)
		.add_constructor()
		.add_method("refine", &GlobalMultiGridRefiner::refine)
		.add_method("assign_grid", (void (GlobalMultiGridRefiner::*)(MultiGrid&)) &GlobalMultiGridRefiner::assign_grid);

//	GridObject
	reg.add_class_<GridObject, Grid>("GridObject", grp)
		.add_constructor()
		.add_method("get_grid", &GridObject::get_grid)
		.add_method("get_subset_handler", &GridObject::get_subset_handler);

//	Grid functions
	reg.add_function("CreateFractal", &CreateFractal, grp);
	
//  GridObject functions
	reg.add_function("LoadGrid", &LoadGrid, grp)
		.add_function("SaveGrid", &SaveGrid, grp)
		.add_function("LoadGridObject", &LoadGridObject, grp)
		.add_function("SaveGridObject", &SaveGridObject, grp)
		.add_function("CreateGridObject", &CreateGridObject, grp);
		
//	refinement
	reg.add_function("TestSubdivision", &TestSubdivision)
		.add_function("CreateSmoothHierarchy", &CreateSmoothHierarchy, grp)
		.add_function("CreateSemiSmoothHierarchy", &CreateSemiSmoothHierarchy, grp)
		.add_function("SaveGridHierarchy", &SaveGridHierarchy, grp);
}

}//	end of namespace 
}//	end of namespace 
