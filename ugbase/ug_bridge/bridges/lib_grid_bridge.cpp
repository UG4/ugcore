//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../registry.h"
#include "../ug_bridge.h"
#include "lib_grid/lib_grid.h"

namespace ug
{
namespace bridge
{

bool LoadGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	return LoadGridFromFile(grid, filename, sh);
}

bool SaveGrid(Grid& grid, SubsetHandler& sh, const char* filename)
{
	return SaveGridToFile(grid, filename, sh);
}

void RegisterLibGridInterface(Registry& reg)
{

//	Grid
	reg.add_class_<Grid>("Grid")
		.add_constructor()
		.add_method("clear", &Grid::clear);
		
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

//  LoadGrid
	reg.add_function("LoadGrid", &LoadGrid);

//  SaveGrid
	reg.add_function("SaveGrid", &SaveGrid);
}

}//	end of namespace 
}//	end of namespace 
