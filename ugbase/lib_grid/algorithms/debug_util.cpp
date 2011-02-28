// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 13.01.2011 (m,d,y)
 
#include "debug_util.h"
using namespace std;

namespace ug{

void PrintElementNumbers(const GeometricObjectCollection& goc)
{
	UG_LOG("grid element numbers:\n");
	for(size_t i = 0; i < goc.num_levels(); ++i)
	{
		if(goc.num_levels() > 1){
			UG_LOG("level " << i << endl);
		}
		UG_LOG("  vertices total:\t" << goc.num<VertexBase>(i) << endl);
		UG_LOG("    normal vrts:\t" << goc.num<Vertex>(i) << endl);
		UG_LOG("    hanging vrts:\t" << goc.num<HangingVertex>(i) << endl);

		UG_LOG("  edges total:\t\t" << goc.num<EdgeBase>(i) << endl);
		UG_LOG("    normal edges:\t" << goc.num<Edge>(i) << endl);
		UG_LOG("    constraining edges:\t" << goc.num<ConstrainingEdge>(i) << endl);
		UG_LOG("    constrained edges:\t" << goc.num<ConstrainedEdge>(i) << endl);

		UG_LOG("  faces total:\t\t" << goc.num<Face>(i) << endl);
		UG_LOG("    normal triangles:\t" << goc.num<Triangle>(i) << endl);
		UG_LOG("    constraining tris:\t" << goc.num<ConstrainingTriangle>(i) << endl);
		UG_LOG("    constrained tris:\t" << goc.num<ConstrainedTriangle>(i) << endl);

		UG_LOG("    normal quads:\t" << goc.num<Quadrilateral>(i) << endl);
		UG_LOG("    constraining quads:\t" << goc.num<ConstrainingQuadrilateral>(i) << endl);
		UG_LOG("    constrained quads:\t" << goc.num<ConstrainedQuadrilateral>(i) << endl);

		UG_LOG("  volumes total:\t" << goc.num<Volume>(i) << endl);

		UG_LOG(endl);
	}
}

void PrintGridElementNumbers(Grid& grid)
{
	PrintElementNumbers(grid.get_geometric_object_collection());
}

void PrintGridElementNumbers(MultiGrid& mg)
{
	PrintElementNumbers(mg.get_geometric_object_collection());
}

void PrintGridElementNumbers(GridSubsetHandler& sh)
{
	PrintElementNumbers(sh.get_geometric_object_collection());
}

}// end of namespace
