// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "grid_bridges.h"
#include "common/space_partitioning/ntree_traverser.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/fractals.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_grid.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_grid/file_io/file_io.h"

using namespace std;

namespace ug{
namespace bridge{

bool CreateFractal(Grid& grid, HangingNodeRefiner_Grid& href,
					number scaleFac, size_t numIterations)
{
	PROFILE_FUNC_GROUP("grid");
//	HangingNodeRefiner_IR1 href(grid);
	return CreateFractal_NormalScale(grid, href, scaleFac, numIterations);
//	return true;
}


bool TestNTree(const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	Grid g;
	SubsetHandler sh(g);
	APosition aPos = aPosition;

	PROFILE_BEGIN(ntree_loading);
	if(!LoadGridFromFile(g, sh, filename, aPos)){
		UG_LOG("  could not load " << filename << endl);
		return false;
	}
	PROFILE_END();


	typedef lg_ntree<3, 3, Volume>	tree_t;
	tree_t	tree(g, aPos);

	PROFILE_BEGIN(ntree_creating_tree);
	tree.create_tree(g.volumes_begin(), g.volumes_end());
	PROFILE_END();

	size_t lvl = FindLowestLeafNodeLevel(tree);

	UG_LOG("Lowest leaf-node-level: " << lvl << endl);

	pair<size_t, size_t> minMax = GetMinMaxNumElements(tree, lvl);
	UG_LOG("Min num elements: " << minMax.first << endl);
	UG_LOG("Max num elements: " << minMax.second << endl);


	const size_t numPicks = 1000000;

	UG_LOG("Picking elements for " << numPicks << " random points:\n");
	Volume* e = NULL;
	size_t numSuccesses = 0;

	PROFILE_BEGIN(ntree_picking_elements);
	for(size_t i = 0; i < numPicks; ++i){
		//vector2 p(urand<number>(-1, 1), urand<number>(-1, 1));
		vector3 p(urand<number>(-1, 1), urand<number>(-1, 1), urand<number>(-1, 1));
		if(FindContainingElement(e, tree, p)){
			++numSuccesses;
		}
	}
	PROFILE_END();

	UG_LOG("  successes: " << numSuccesses << "\n");
	UG_LOG("  failures: " << numPicks - numSuccesses << "\n");
	return true;
}


void RegisterGridBridge_Misc(Registry& reg, string parentGroup)
{
	string grp = parentGroup;
	reg.add_function("CreateFractal", &CreateFractal, grp);

	reg.add_function("PrintGridElementNumbers", static_cast<void (*)(MultiGrid&)>(&PrintGridElementNumbers), grp)
		.add_function("PrintGridElementNumbers", static_cast<void (*)(Grid&)>(&PrintGridElementNumbers), grp)
		.add_function("PrintAttachmentInfo", &PrintAttachmentInfo, grp);

	reg.add_function("TestNTree", &TestNTree, grp);
}

}//	end of namespace
}//	end of namespace
