#include "common/util/vec_for_each.h"
#include "misc_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	EraseConnectingElements
///	erases all elements that connect v1 and v2
void EraseConnectingElements(Grid& grid, Vertex* v1, Vertex* v2)
{
//	check volumes for direct connection
	if(grid.num_volumes() > 0)
	{
	//	iterate through associated volumes
		std::vector<Volume*> vols;
		CollectVolumes(vols, grid, v1);
		for_each_in_vec(Volume* v, vols){
			uint numVrts = v->num_vertices();
			bool gotOne = false;
			for(uint i = 0; i < numVrts; ++i)
			{
				if(v->vertex(i) == v2)
				{
					gotOne = true;
					break;
				}
			}

			if(gotOne == true)
				grid.erase(v);
		}end_for;
	}
	
//	check faces for direct connection.
	if(grid.num_faces() > 0)
	{
	//	iterate through associated faces
		std::vector<Face*> faces;
		CollectFaces(faces, grid, v1);
		for_each_in_vec(Face* f, faces){
			uint numVrts = f->num_vertices();
			bool gotOne = false;
			for(uint i = 0; i < numVrts; ++i)
			{
				if(f->vertex(i) == v2)
				{
					gotOne = true;
					break;
				}
			}

			if(gotOne == true)
				grid.erase(f);
		}end_for;
	}

//	check edges
	if(grid.num_edges() > 0)
	{
	//	iterate through associated edges
		std::vector<Edge*> edges;
		CollectEdges(edges, grid, v1);
		for_each_in_vec(Edge* e, edges){
		//	if e contains v2 we have to remove it.
			if((e->vertex(0) == v2) || (e->vertex(1) == v2))
			{
				grid.erase(e);
			}
		}end_for;
	}
}

}//	end of namespace
