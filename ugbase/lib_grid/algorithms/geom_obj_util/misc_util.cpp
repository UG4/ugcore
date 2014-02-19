//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d03

#include "misc_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	EraseConnectingElements
///	erases all elements that connect v1 and v2
void EraseConnectingElements(Grid& grid, Vertex* v1, Vertex* v2)
{
//	check edges
	if(grid.num_edges() > 0)
	{
	//	iterate through associated edges
		Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(v1);
		Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v1);
		while(iter != iterEnd)
		{
			EdgeBase* e = *iter;
			++iter;
		//	if e contains v2 we have to remove it.
			if((e->vertex(0) == v2) || (e->vertex(1) == v2))
			{
				grid.erase(e);
			}
		}
	}

//	check faces for direct connection.
	if(grid.num_faces() > 0)
	{
	//	if FACEOPT_AUTOGENERATE_EDGES is enabled, all faces that connect
	//	v1 and v2 have been already removed on edge-erase.
		if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		{
			//UG_LOG(" ++++ this piece of code shouldn't be executed +++ EraseConnectingElements-Faces \n");
		//	iterate through associated faces
			Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v1);
			Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v1);
			while(iter != iterEnd)
			{
				Face* f = *iter;
				++iter;
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
			}
		}
	}

//	check volumes for direct connection
	if(grid.num_volumes() > 0)
	{
	//	connecting volumes can at this point only exist if neither
	//	VOLOPT_AUTOGENERATE_EDGES nor VOLOPT_AUTOGENERATE_FACES is enabled.
		if(!(grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES) ||
			grid.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)))
		{
			//UG_LOG(" ++++ this piece of code shouldn't be executed +++ EraseConnectingElements-VOLUMES \n");
		//	iterate through associated volumes
			Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(v1);
			Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(v2);
			while(iter != iterEnd)
			{
				Volume* v = *iter;
				++iter;
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
			}
		}
	}
}

}//	end of namespace
