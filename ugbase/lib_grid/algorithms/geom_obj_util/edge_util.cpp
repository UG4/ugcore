//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#include <vector>
#include "edge_util.h"
#include "lib_grid/grid/grid_util.h"
#include "vertex_util.h"
#include "lib_grid/algorithms/refinement/regular_refinement.h"

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
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(e);
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
//	GetAssociatedFaces
int GetAssociatedFaces(Face** facesOut, Grid& grid,
						EdgeBase* e, int maxNumFaces)
{
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		int counter = 0;
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(e);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(e);
			iter != grid.associated_faces_end(e); ++iter)
		{
			Face* tf = *iter;

			if(counter < maxNumFaces)
				facesOut[counter] = tf;

			counter++;
		}
		return counter;
	}
	else
	{
	//	we're using grid::mark for maximal speed.
		grid.begin_marking();
	//	mark the end-points of the edge
		grid.mark(e->vertex(0));
		grid.mark(e->vertex(1));

	//	we have to find the triangles 'by hand'
	//	iterate over all associated faces of vertex 0
		int counter = 0;
		VertexBase* v = e->vertex(0);
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v);
			iter != iterEnd; ++iter)
		{
			Face* tf = *iter;
			uint numVrts = tf->num_vertices();
			int numMarked = 0;
			for(uint i = 0; i < numVrts; ++i)
				if(grid.is_marked(tf->vertex(i)))
					numMarked++;
			if(numMarked > 1)
			{
			//	the face is connected with the edge
				if(counter < maxNumFaces)
					facesOut[counter] = tf;
				counter++;
			}
		}

	//	done with marking
		grid.end_marking();

		return counter;
	}
}

////////////////////////////////////////////////////////////////////////
int CalculateNormal(vector3& vNormOut, Grid& grid, EdgeBase* e,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Grid::FaceAttachmentAccessor<ANormal>* paaNormFACE)
{
	Face* f[2];
	
	int numFaces = GetAssociatedFaces(f, grid, e, 2);
	
	switch(numFaces){
	
	case 0:{ //	if there are no associated faces.
		//	we'll assume that the edge lies in the xy plane and return its normal
		vector3 dir;
		VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		vNormOut.x = dir.y;
		vNormOut.y = -dir.y;
		vNormOut.z = 0;
		VecNormalize(vNormOut, vNormOut);
		}break;
	
	case 1: //	if there is one face, the normal will be set to the faces normal
		if(paaNormFACE)
			vNormOut = (*paaNormFACE)[f[0]];
		else{
			CalculateNormal(vNormOut, f[0], aaPos);
		}
		break;
	
	default: //	there are at least 2 associated faces
		if(paaNormFACE)
			VecAdd(vNormOut, (*paaNormFACE)[f[0]], (*paaNormFACE)[f[1]]);
		else{
			vector3 fn0, fn1;
			CalculateNormalNoNormalize(fn0, f[0], aaPos);
			CalculateNormalNoNormalize(fn1, f[1], aaPos);
			VecAdd(vNormOut, fn0, fn1);
		}
		VecNormalize(vNormOut, vNormOut);
		break;
	}
	
	return numFaces;
}

int CalculateNormalNoNormalize(vector3& vNormOut, Grid& grid, EdgeBase* e,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Grid::FaceAttachmentAccessor<ANormal>* paaNormFACE)
{
	Face* f[2];
	
	int numFaces = GetAssociatedFaces(f, grid, e, 2);
	
	switch(numFaces){
	
	case 0:{ //	if there are no associated faces.
		//	we'll assume that the edge lies in the xy plane and return its normal
		vector3 dir;
		VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		vNormOut.x = dir.y;
		vNormOut.y = -dir.y;
		vNormOut.z = 0;
		}break;
	
	case 1: //	if there is one face, the normal will be set to the faces normal
		if(paaNormFACE)
			vNormOut = (*paaNormFACE)[f[0]];
		else{
			CalculateNormalNoNormalize(vNormOut, f[0], aaPos);
		}
		break;
	
	default: //	there are at least 2 associated faces
		if(paaNormFACE)
			VecAdd(vNormOut, (*paaNormFACE)[f[0]], (*paaNormFACE)[f[1]]);
		else{
			vector3 fn0, fn1;
			CalculateNormalNoNormalize(fn0, f[0], aaPos);
			CalculateNormalNoNormalize(fn1, f[1], aaPos);
			VecAdd(vNormOut, fn0, fn1);
			VecScale(vNormOut, vNormOut, 0.5);
		}
		VecNormalize(vNormOut, vNormOut);
		break;
	}
	
	return numFaces;
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
				if(numVrts > vSubstituteVertices.size())
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
				oldFace->create_faces_by_edge_split(edgeIndex, newVertex, vNewFaces, &vSubstituteVertices.front());
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

EdgeBase* SwapEdge(Grid& grid, EdgeBase* e)
{
//	get the associated faces.
//	Only 2 associated faces are allowed.
	Face* f[2];
	if(GetAssociatedFaces(f, grid, e, 2) != 2)
		return NULL;

//	make sure that both faces are triangles
	if((f[0]->num_vertices() != 3) || (f[1]->num_vertices() != 3))
		return NULL;

//	begin marking
	grid.begin_marking();

//	mark the end-points of the edge
	grid.mark(e->vertex(0));
	grid.mark(e->vertex(1));	

//	get the two vertices that will be connected by the new edge
	VertexBase* v[2];
	int vrtInd[2];
	for(int j = 0; j < 2; ++j){
		v[j] = NULL;
		for(int i = 0; i < 3; ++i){
			VertexBase* vrt = f[j]->vertex(i);
			if(!grid.is_marked(vrt)){
				vrtInd[j] = i;
				v[j] = vrt;
				break;
			}
		}
	}

//	we're done with marking
	grid.end_marking();

//	make sure that both vertices have been found.
	if(!(v[0] && v[1]))
		return NULL;

//	make sure that no edge exists between v[0] and v[1]
	if(grid.get_edge(v[0], v[1]))
		return NULL;

//	the indices of the marked points in the first triangle
	int ind1 = (vrtInd[0] + 1) % 3;
	int ind2 = (vrtInd[0] + 2) % 3;

//	create the new edge
	EdgeBase* nEdge = *grid.create_by_cloning(e, EdgeDescriptor(v[0], v[1]), e);

//	create the new faces
	grid.create<Triangle>(TriangleDescriptor(v[0], f[0]->vertex(ind1), v[1]), f[0]);
	grid.create<Triangle>(TriangleDescriptor(f[0]->vertex(ind2), v[0], v[1]), f[1]);

//	erase the old faces
	grid.erase(f[0]);
	grid.erase(f[1]);

//	erase the old edge
	grid.erase(e);
	
	return nEdge;
}

////////////////////////////////////////////////////////////////////////
bool CutEdgesWithPlane(Selector& sel, const vector3& p, const vector3& n,
						APosition& aPos)
{
	if(!sel.get_assigned_grid()){
		UG_LOG("ERROR in CutEdgesWithPlane: sel has to be assigned to a grid.\n");
		return false;
	}
	
	Grid& grid = *sel.get_assigned_grid();
	
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("ERROR in CutEdgesWithPlane: aPos has to be attached to the vertices of the grid.\n");
		return false;
	}
	
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	
//	used for plane-intersection
	number t;
	vector3 v;
	
//	iterate through all edges and deselect all that do not intersect the plane
//	deselect all vertices, faces and volumes, too.
	sel.clear<VertexBase>();
	sel.clear<Face>();
	sel.clear<Volume>();
		
	EdgeBaseIterator iter = sel.begin<EdgeBase>();
	while(iter != sel.end<EdgeBase>()){
		EdgeBase* e = *iter;
		iter++;
	
	//	check whether the edge intersects the plane
		vector3 rayDir;
		VecSubtract(rayDir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		
		bool bIntersect = RayPlaneIntersection(v, t, aaPos[e->vertex(0)],
												rayDir, p, n);
		if(!bIntersect || t < SMALL || t > 1. - SMALL)
		{
		//	it doesn't. Remove it from the selector
			sel.deselect(e);
		}
	}
	
//	refine all selected edges. RefinementCallbackEdgePlaneCut will insert
//	new vertices on the plane.
	RefinementCallbackEdgePlaneCut refCallbackEdgePlaneCut(grid, p, n, aPos);
	return Refine(grid, sel, &refCallbackEdgePlaneCut);
}

}//	end of namespace
