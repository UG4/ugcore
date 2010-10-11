//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../registry.h"
#include "../ug_bridge.h"
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"

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




void TestSubdivision(const char* fileIn, const char* fileOut, int numRefs)
{
//todo: Callbacks have to make sure that their attachment is accessible in the grid.
//		even if they were initialized before the attachment was attached to the grid.
	MultiGrid mg;
	SubsetHandler sh(mg);
	
	if(LoadGridFromFile(mg, fileIn, sh)){
		RefinementCallbackSubdivisionLoop refCallback(mg);
		GlobalMultiGridRefiner ref(mg, &refCallback);
		for(int lvl = 0; lvl < numRefs; ++lvl){
			ref.refine();
		}

		//SaveGridToFile(mg, fileOut, mg.get_hierarchy_handler());

		APosition aProjected;
		mg.attach_to_vertices(aProjected);
		ProjectToLimitPLoop(mg, aPosition, aProjected);
		SaveGridToFile(mg, fileOut, mg.get_hierarchy_handler(), aProjected);

	}
	else{
		UG_LOG("Load failed. aborting...\n");
	}
}


////////////////////////////////////////////////////////////////////////
void RegisterLibGridInterface(Registry& reg)
{

//	Grid
	reg.add_class_<Grid>("Grid")
		.add_constructor()
		.add_method("clear", &Grid::clear)
		.add_method("num_vertices", &Grid::num_vertices)
		.add_method("num_edges", &Grid::num_edges)
		.add_method("num_faces", &Grid::num_faces)
		.add_method("num_volumes", &Grid::num_volumes);
		
//	MultiGrid
	reg.add_class_<MultiGrid, Grid>("MultiGrid")
		.add_constructor();

//  ISubsetHandler
	reg.add_class_<ISubsetHandler>("ISubsetHandler")
		.add_method("num_subsets", &ISubsetHandler::num_subsets);
	
//	SubsetHandler
	reg.add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler")
		.add_constructor()
		.add_method("assign_grid", &SubsetHandler::assign_grid);

//	MGSubsetHandler
	reg.add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler")
		.add_constructor()
		.add_method("assign_grid", &MGSubsetHandler::assign_grid);

//	HangingNodeRefiner
	reg.add_class_<HangingNodeRefiner_IR1>("HangingNodeRefiner")
		.add_constructor()
		.add_method("assign_grid", &HangingNodeRefiner_IR1::assign_grid);

//	GridObject
	reg.add_class_<GridObject, Grid>("GridObject")
		.add_constructor()
		.add_method("get_grid", &GridObject::get_grid)
		.add_method("get_subset_handler", &GridObject::get_subset_handler);

//	Grid functions
	reg.add_function("CreateFractal", &CreateFractal);
	
//  GridObject functions
	reg.add_function("LoadGrid", &LoadGrid)
		.add_function("SaveGrid", &SaveGrid)
		.add_function("LoadGridObject", &LoadGridObject)
		.add_function("SaveGridObject", &SaveGridObject)
		.add_function("CreateGridObject", &CreateGridObject);
		
	reg.add_function("TestSubdivision", &TestSubdivision);
}

}//	end of namespace 
}//	end of namespace 
