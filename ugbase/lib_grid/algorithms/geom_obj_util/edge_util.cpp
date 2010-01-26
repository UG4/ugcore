//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#include <vector>
#include "edge_util.h"
#include "lib_grid/grid/grid_util.h"
#include "vertex_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
int GetEdgeIndex(Face* f, EdgeBase* e)
{
	uint numEdges = f->num_edges();
	EdgeDescriptor ed;
	for(uint i = 0; i < numEdges; ++i)
	{
		f->edge(i, ed);
		if(CompareVertices(e, &ed))
			return (int)i;
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
int GetEdgeIndex(Volume* vol, EdgeBase* e)
{
	uint numEdges = vol->num_edges();
	EdgeDescriptor ed;
	for(uint i = 0; i < numEdges; ++i)
	{
		vol->edge(i, ed);
		if(CompareVertices(e, &ed))
			return (int)i;
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
bool IsBoundaryEdge2D(Grid& grid, EdgeBase* e)
{
//	get the number of connected faces. if only one face is connected then
//	the edge is considered to be a boundary edge.
	int counter = 0;
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		for(FaceIterator iter = grid.associated_faces_begin(e);
			iter != grid.associated_faces_end(e); ++iter)
			++counter;

		if(counter == 1)
			return true;
	}
	else
	{
	//	fill a vector using a helper function
		vector<Face*> vFaces;
		CollectFaces(vFaces, grid, e, false);
		if(vFaces.size() == 1)
			return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
bool CollapseEdge(Grid& grid, EdgeBase* e, VertexBase* newVrt)
{
//	prepare the grid, so that we may perform Grid::replace_vertex.
//	create collapse geometries first and delete old geometries.
	{
	//	collect adjacent faces
		vector<Face*> vFaces;
		CollectFaces(vFaces, grid, e);

	//	a vector thats is used to store the collapse-geometries.
		vector<Face*> vNewFaces;

		for(uint i = 0; i < vFaces.size(); ++i)
		{
			Face* f = vFaces[i];
		//	create the collapse-geometry
			int eInd = GetEdgeIndex(f, e);
			f->collapse_edge(vNewFaces, eInd, newVrt);

		//	register the new faces
			for(uint j = 0; j < vNewFaces.size(); ++j)
				grid.register_element(vNewFaces[j], f);

		//	erase f
			grid.erase(f);
		}
	}

	{
	//	collect adjacent volumes
		vector<Volume*> vVolumes;
		CollectVolumes(vVolumes, grid, e);

	//	a vector thats used to store the collapse-geometries.
		vector<Volume*> vNewVolumes;

		for(uint i = 0; i < vVolumes.size(); ++i)
		{
			Volume* v = vVolumes[i];
		//	create the collapse-geometry
			int eInd = GetEdgeIndex(v, e);
			v->collapse_edge(vNewVolumes, eInd, newVrt);

		//	register the new volumes
			for(uint j = 0; j < vNewVolumes.size(); ++j)
				grid.register_element(vNewVolumes[i], v);

		//	erase v
			grid.erase(v);
		}
	}

//	store the end-points of e
	VertexBase* v[2];
	v[0] = e->vertex(0);
	v[1] = e->vertex(1);

//	erase e
	grid.erase(e);

//	perform replace_vertex for both endpoints
	for(uint i = 0; i < 2; ++i)
	{
		if(v[i] != newVrt)
			grid.replace_vertex(v[i], newVrt);
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
bool EdgeCollapseIsValid(Grid& grid, EdgeBase* e)
{
//TODO: Make sure that this is sufficient for geometries including Quads.
//TODO: Check validity for volumes.
//	in order to not destroy the topology of the grid, a collapse is
//	only valid if no three edges build a triangle that does not exist
//	in the grid.

//	first we need all vertices that are connected with end-points of e.
	vector<VertexBase*> vVertices1;
	vector<VertexBase*> vVertices2;

	CollectNeighbours(vVertices1, grid, e->vertex(0)); // e->vertex(0) is not contained in vVertices1!
	CollectNeighbours(vVertices2, grid, e->vertex(1)); // e->vertex(1) is not contained in vVertices2!

//	we need access to the triangles connected with e.
	vector<Face*> vFaces;
	CollectFaces(vFaces, grid, e);

//	this face descriptor will be needed during the algorithm.
	FaceDescriptor fd(3);
	fd.set_vertex(0, e->vertex(0));
	fd.set_vertex(1, e->vertex(1));

//	now check for each vertex in vVertices1 if it also exists in vVertices2.
	for(uint i = 0; i < vVertices1.size(); ++i)
	{
		VertexBase* v1 = vVertices1[i];

		if(v1 != e->vertex(1))
		{
			for(uint j = 0; j < vVertices2.size(); ++j)
			{
				VertexBase* v2 = vVertices2[j];

				if(v1 == v2)
				{
				//	v1 and v2 exist in both arrays.
				//	check if a triangle exists that connects e with v1.
				//	it is sufficient to search in vFaces.
					bool bGotOne = false;
					fd.set_vertex(2, v1);
					for(uint k = 0; k < vFaces.size(); ++k)
					{
						if(CompareVertices(vFaces[k], &fd))
						{
							bGotOne = true;
							break;
						}
					}

				//	if there was none, the collapse would be illegal.
					if(!bGotOne)
						return false;
				}
			}
		}
	}

//	everything is fine.
	return true;
}

////////////////////////////////////////////////////////////////////////
bool CreateEdgeSplitGeometry(Grid& destGrid, Grid& srcGrid, EdgeBase* e,
							VertexBase* newVertex, AVertexBase* paAssociatedVertices)
{

	if((paAssociatedVertices == NULL) && (&destGrid != &srcGrid))
		return false;

	vector<VertexBase*> vVrts;

//	If paAssociatedVertices is specified, we'll have to find the vertices
//	in destGrid for each element that is adjacent to e. We then have
//	to store these vertices in vVrts and pass them to the
//	split-routine of the element.

//	the attachment-accessor
	Grid::VertexAttachmentAccessor<AVertexBase> aaAssociatedVertices;

	if(paAssociatedVertices != NULL)
	{
		AVertexBase& aAssociatedVertices = *paAssociatedVertices;

	//	check if aVertexBase is properly attached.
		if(!srcGrid.has_vertex_attachment(aAssociatedVertices))
		//	attach it and initialize its values.
			srcGrid.attach_to_vertices_dv(aAssociatedVertices, NULL, false);

	//	initialize the attachment-accessor
		aaAssociatedVertices.access(srcGrid, aAssociatedVertices);
	}

//	we will store the substitute-vertices in this vector - if they are needed.
	vector<VertexBase*> vSubstituteVertices;

//	split the edge
	{
	//	simply create two new edges
		EdgeBase* parent = e;
	//	the grids do not match then we can't pass e as a parent
		if(&srcGrid != &destGrid)
			parent = NULL;
			
		if(paAssociatedVertices){
			destGrid.create<Edge>(EdgeDescriptor(aaAssociatedVertices[e->vertex(0)], newVertex), parent);
			destGrid.create<Edge>(EdgeDescriptor(newVertex, aaAssociatedVertices[e->vertex(1)]), parent);
		}
		else{
			destGrid.create<Edge>(EdgeDescriptor(e->vertex(0), newVertex), parent);
			destGrid.create<Edge>(EdgeDescriptor(newVertex, e->vertex(1)), parent);
		}
	}

//	split faces
	{
	//	collect all faces which are associated with e
		vector<Face*> vOldFaces;
		CollectFaces(vOldFaces, srcGrid, e, false);

	//	we will collect new faces in this vector
		vector<Face*> vNewFaces;

	//	iterate through those faces and split each.
	//	If vertices are missing in destGrid, create them first.
		for(vector<Face*>::iterator oldFaceIter = vOldFaces.begin();
			oldFaceIter != vOldFaces.end(); ++oldFaceIter)
		{
			Face* oldFace = *oldFaceIter;
			uint numVrts = oldFace->num_vertices();

		//	get the index of e in oldFace
			int edgeIndex = GetEdgeIndex(oldFace, e);

		//	clear the new faces-vec
			vNewFaces.clear();

		//	get the substitute-vertices if they are required
			if(paAssociatedVertices != NULL)
			{
				vSubstituteVertices.resize(numVrts);

			//	check for each vertex in oldFace, if an associated vertex exists in destGrid.
			//	if not create a new one by cloning the associated one in srcGrid.
			//	store the vertices in vSubstituteVertices
				{
					for(uint i = 0; i < numVrts; ++i)
					{
						if(aaAssociatedVertices[oldFace->vertex(i)] == NULL)
							aaAssociatedVertices[oldFace->vertex(i)] = *destGrid.create_by_cloning(oldFace->vertex(i));
						vSubstituteVertices[i] = aaAssociatedVertices[oldFace->vertex(i)];
					}
				}

			//	create the new faces by splitting the old face. use substitutes.
				oldFace->create_faces_by_edge_split(edgeIndex, newVertex, vNewFaces, &vSubstituteVertices);
			}
			else
			{
			//	create the new faces by splitting the old face.
			//	no substitutes required
				oldFace->create_faces_by_edge_split(edgeIndex, newVertex, vNewFaces);
			}

		//	register all new faces at destGrid
			{
				Face* pParent = NULL;
				if(&srcGrid == &destGrid)
					pParent = oldFace;

				for(vector<Face*>::iterator iter = vNewFaces.begin();
					iter != vNewFaces.end(); ++iter)
				{
					destGrid.register_element(*iter, pParent);
					if(pParent)
						destGrid.pass_on_values(pParent, *iter);
				}
			}
		}
	}

	return true;
}

}//	end of namespace
