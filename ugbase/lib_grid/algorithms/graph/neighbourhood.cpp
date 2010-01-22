//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d18

#include "neighbourhood.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
void CollectNeighbours(std::vector<EdgeBase*>& vNeighboursOut, EdgeBase* e,
					   Grid& grid, NeighbourhoodType nbhType)
{
//	clear the container
	vNeighboursOut.clear();
	
//	default neighbourhood:
	if(nbhType == NHT_DEFAULT)
		nbhType = NHT_VERTEX_NEIGHBOURS;

//	edges can only be vertex-neighbours
	if(nbhType != NHT_VERTEX_NEIGHBOURS)
		return;

//	begin marking
	grid.begin_marking();
	
//	mark the edge
	grid.mark(e);
	
//	mark the vertices of the edge
	grid.mark(e->vertex(0));
	grid.mark(e->vertex(1));	
	
//	iterate over all edges that are connected to the vertices.
//	if the edge is not yet marked, we have to push it to vNeighboursOut.
	for(uint i = 0; i < 2; ++i)
	{
		EdgeBaseIterator iterEnd = grid.associated_edges_end(e->vertex(i));
		for(EdgeBaseIterator iter = grid.associated_edges_begin(e->vertex(i));
			iter != iterEnd; ++iter)
		{
			if(!grid.is_marked(*iter))
			{
				vNeighboursOut.push_back(*iter);
				grid.mark(*iter);
			}
		}
	}

//	end marking
	grid.end_marking();
}


////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
void CollectNeighbours(std::vector<Face*>& vNeighboursOut, Face* f,
					   Grid& grid, NeighbourhoodType nbhType)
{
//	clear the container
	vNeighboursOut.clear();
	
//	default neighbourhood:
	if(nbhType == NHT_DEFAULT)
		nbhType = NHT_EDGE_NEIGHBOURS;

//	faces can't be face-neighbours
	if(nbhType == NHT_FACE_NEIGHBOURS)
		return;

//	begin marking
	grid.begin_marking();
	
//	mark the face
	grid.mark(f);
	
//	mark the vertices of the face
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
		grid.mark(f->vertex(i));
	
//	in order to get the maximum speed-up, we'll try to use
//	associated elements in grid.
	if((nbhType == NHT_EDGE_NEIGHBOURS)
		&& grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES
								  | EDGEOPT_STORE_ASSOCIATED_FACES
								  | FACEOPT_AUTOGENERATE_EDGES))
	{
	//	iterate through associated edges
		EdgeBaseIterator eEnd = grid.associated_edges_end(f);
		for(EdgeBaseIterator eIter = grid.associated_edges_begin(f);
			eIter != eEnd; ++eIter)
		{
		//	iterate through associated folumes of the eace
			FaceIterator fEnd = grid.associated_faces_end(*eIter);
			for(FaceIterator iter = grid.associated_faces_begin(*eIter);
				iter != fEnd; ++iter)
			{
			//	if the face is not yet marked, then add it to the neighbours
				if(!grid.is_marked(*iter))
				{
					grid.mark(*iter);
					vNeighboursOut.push_back(*iter);
				}
			}
		}
	//	we're done in here. end-marking and return.
		grid.end_marking();
		return;
	}


//	iterate over all faces that are connected to the vertices.
//	if the face shares the elements as required by nbhType and
//	it is not yet marked, we have to push it to vNeighboursOut.
//	to optimize speed we'll check both valid nbhTypes separately.
//	the first case indeed is a subcase of the second
//	(compare counted marked vertices against nbhType)
	switch(nbhType)
	{
	case NHT_VERTEX_NEIGHBOURS:
		for(uint i = 0; i < numVrts; ++i)
		{
			FaceIterator iterEnd = grid.associated_faces_end(f->vertex(i));
			for(FaceIterator iter = grid.associated_faces_begin(f->vertex(i));
				iter != iterEnd; ++iter)
			{
				if(!grid.is_marked(*iter))
				{
					vNeighboursOut.push_back(*iter);
					grid.mark(*iter);
				}
			}
		}
		break;

	case NHT_EDGE_NEIGHBOURS:
		for(uint i = 0; i < numVrts; ++i)
		{
			FaceIterator iterEnd = grid.associated_faces_end(f->vertex(i));
			for(FaceIterator iter = grid.associated_faces_begin(f->vertex(i));
				iter != iterEnd; ++iter)
			{
				Face* nf = *iter;
				if(!grid.is_marked(nf))
				{
				//	count the marked vertices that are contained by *iter
				//	if there are more than 1, the faces share an edge
				//	(at least in a regular grid)
					int counter = 0;
					
					for(uint j = 0; j < nf->num_vertices(); ++j)
					{
						if(grid.is_marked(nf->vertex(j)))
						{
							++counter;
							if(counter > 1)
							{
								vNeighboursOut.push_back(nf);
								grid.mark(nf);
								break;
							}
						}
					}
				}
			}
		}
		break;
	}

//	end marking
	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
void CollectNeighbours(std::vector<Volume*>& vNeighboursOut, Volume* v,
					   Grid& grid, NeighbourhoodType nbhType)
{
//	clear the container
	vNeighboursOut.clear();
	
//	default neighbourhood:
	if(nbhType == NHT_DEFAULT)
		nbhType = NHT_FACE_NEIGHBOURS;

//	begin marking
	grid.begin_marking();
	
//	mark the volume
	grid.mark(v);
	
//	mark the vertices of the volume
	uint numVrts = v->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
		grid.mark(v->vertex(i));

//	in order to get the maximum speed-up, we'll try to use
//	associated elements in grid.
	if((nbhType == NHT_FACE_NEIGHBOURS)
		&& grid.option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES
								  | FACEOPT_STORE_ASSOCIATED_VOLUMES
								  | VOLOPT_AUTOGENERATE_FACES))
	{
	//	iterate through associated faces
		FaceIterator fEnd = grid.associated_faces_end(v);
		for(FaceIterator fIter = grid.associated_faces_begin(v);
			fIter != fEnd; ++fIter)
		{
		//	iterate through associated volumes of the face
			VolumeIterator vEnd = grid.associated_volumes_end(*fIter);
			for(VolumeIterator iter = grid.associated_volumes_begin(*fIter);
				iter != vEnd; ++iter)
			{
			//	if the volume is not yet marked, then add it to the neighbours
				if(!grid.is_marked(*iter))
				{
					grid.mark(*iter);
					vNeighboursOut.push_back(*iter);
				}
			}
		}
	//	we're done in here. end-marking and return.
		grid.end_marking();
		return;
	}

//	iterate over all volumes that are connected to the vertices.
//	if the volume shares the elements as required by nbhType and
//	it is not yet marked, we have to push it to vNeighboursOut.
//	to optimize speed we'll check both valid nbhTypes separately.
//	the first case indeed is a subcase of the second
	if(nbhType == NHT_VERTEX_NEIGHBOURS)
		for(uint i = 0; i < numVrts; ++i)
		{
			VolumeIterator iterEnd = grid.associated_volumes_end(v->vertex(i));
			for(VolumeIterator iter = grid.associated_volumes_begin(v->vertex(i));
				iter != iterEnd; ++iter)
			{
				if(!grid.is_marked(*iter))
				{
					vNeighboursOut.push_back(*iter);
					grid.mark(*iter);
				}
			}
		}
	else
	{
		for(uint i = 0; i < numVrts; ++i)
		{
			VolumeIterator iterEnd = grid.associated_volumes_end(v->vertex(i));
			for(VolumeIterator iter = grid.associated_volumes_begin(v->vertex(i));
				iter != iterEnd; ++iter)
			{
				Volume* nv = *iter;
				if(!grid.is_marked(nv))
				{
				//	count the marked vertices that are contained by *iter
				//	if there as many as in nbhTypes, the volume is a neighbour.
				//	(at least in a regular grid)
					int counter = 0;
					
					for(uint j = 0; j < nv->num_vertices(); ++j)
					{
						if(grid.is_marked(nv->vertex(j)))
						{
							++counter;
							if(counter >= nbhType)
							{
								vNeighboursOut.push_back(nv);
								grid.mark(nv);
								break;
							}
						}
					}
				}
			}
		}
	}

//	end marking
	grid.end_marking();
}

}//	end of namespace
