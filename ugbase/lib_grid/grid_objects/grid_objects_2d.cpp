//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d04

#include <vector>
#include <algorithm>
#include "grid_objects_2d.h"
//#include "../algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	TOOLS
///	helpful if a local vertex-order is required
/**
 * cornersOut and cornersIn both have to be of size numCorners.
 * After termination cornersOut will contain the vertices of
 * cornersIn, starting from firstCorner, taking vertices modulo numCorners.
 * If cornersOut == cornersIn, the method will fail! This is ok since
 * the method is used locally and has been created for a special case.
 */
static inline
bool ReorderCornersCCW(Vertex** cornersOut, Vertex** const cornersIn,
					   int numCorners, int firstCorner)
{
	cornersOut[0] = cornersIn[firstCorner];
	for(int i = 1; i < numCorners; ++i)
		cornersOut[i] = cornersIn[(firstCorner + i) % numCorners];
	return true;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	FACES

////////////////////////////////////////////////////////////////////////
//	TriangleDescriptor
TriangleDescriptor::TriangleDescriptor(const TriangleDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
}

TriangleDescriptor::TriangleDescriptor(Vertex* v1, Vertex* v2, Vertex* v3)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
}

////////////////////////////////////////////////////////////////////////
//	CustomTriangle
template <class ConcreteTriangleType, class BaseClass>
CustomTriangle<ConcreteTriangleType, BaseClass>::
CustomTriangle(const TriangleDescriptor& td)
{
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
}

template <class ConcreteTriangleType, class BaseClass>
CustomTriangle<ConcreteTriangleType, BaseClass>::
CustomTriangle(Vertex* v1, Vertex* v2, Vertex* v3)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
}

template <class ConcreteTriangleType, class BaseClass>
std::pair<GridBaseObjectId, int>
CustomTriangle<ConcreteTriangleType, BaseClass>::
get_opposing_object(Vertex* vrt) const
{
	for(int i = 0; i < 3; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(EDGE, (i + 1) % 3);
		}
	}

	UG_THROW("The given vertex is not contained in the given face.");
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
refine(std::vector<Face*>& vNewFacesOut,
		Vertex** newFaceVertexOut,
		Vertex** newEdgeVertices,
		Vertex* newFaceVertex,
		Vertex** pSubstituteVertices)
{
//TODO: complete triangle refine

	*newFaceVertexOut = newFaceVertex;
	vNewFacesOut.clear();

//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	get the number of new vertices.
	uint numNewVrts = 0;
	for(uint i = 0; i < 3; ++i)
	{
		if(newEdgeVertices[i] != NULL)
			++numNewVrts;
	}

//	if newFaceVertex is specified, then create three sub-triangles and
//	refine each. If not then refine the triangle depending
	if(newFaceVertex)
	{
		if(numNewVrts > 0){
			assert(!"Problem in CustomTriangle::refine: refine with newFaceVertex and newEdgeVertices is not yet supported.");
			return false;
		}

	//	create three new triangles
		vNewFacesOut.push_back(new ConcreteTriangleType(vrts[0], vrts[1], newFaceVertex));
		vNewFacesOut.push_back(new ConcreteTriangleType(vrts[1], vrts[2], newFaceVertex));
		vNewFacesOut.push_back(new ConcreteTriangleType(vrts[2], vrts[0], newFaceVertex));
		return true;
	}
	else
	{
		switch(numNewVrts)
		{
			case 1:
			{
			//	get the index of the edge that will be refined
				int iNew = -1;;
				for(int i = 0; i < 3; ++i){
					if(newEdgeVertices[i]){
						iNew = i;
						break;
					}
				}
				
			//	the corners. The first corner is the corner on the opposite side of iNew.
			//	Other follow in ccw order
				int iCorner[3];
				iCorner[0] = (iNew + 2) % 3;
				iCorner[1] = (iCorner[0] + 1) % 3;
				iCorner[2] = (iCorner[1] + 1) % 3;
					
			//	create the new triangles.
				vNewFacesOut.push_back(new ConcreteTriangleType(vrts[iCorner[0]], vrts[iCorner[1]],
																newEdgeVertices[iNew]));
				vNewFacesOut.push_back(new ConcreteTriangleType(vrts[iCorner[0]], newEdgeVertices[iNew],
																vrts[iCorner[2]]));
																
				return true;
			}

			case 2:
			{
			//	get the index of the edge that won't be refined
				int iFree = -1;
				for(int i = 0; i < 3; ++i){
					if(!newEdgeVertices[i]){
						iFree = i;
						break;
					}
				}
				
			//	the refined edges
				int iNew[2];
				iNew[0] = (iFree + 1) % 3;
				iNew[1] = (iFree + 2) % 3;
				
			//	the corners
				int iCorner[3];
				iCorner[0] = iFree;
				iCorner[1] = (iFree + 1) % 3;
				iCorner[2] = (iFree + 2) % 3;
				
			//	create the faces
				vNewFacesOut.push_back(new ConcreteTriangleType(newEdgeVertices[iNew[0]],
																vrts[iCorner[2]],
																newEdgeVertices[iNew[1]]));
				vNewFacesOut.push_back(new Quadrilateral(vrts[iCorner[0]], vrts[iCorner[1]],
														newEdgeVertices[iNew[0]], newEdgeVertices[iNew[1]]));
				return true;
			}

			case 3:
			{
			//	perform regular refine.
				vNewFacesOut.push_back(new ConcreteTriangleType(vrts[0], newEdgeVertices[0], newEdgeVertices[2]));
				vNewFacesOut.push_back(new ConcreteTriangleType(vrts[1], newEdgeVertices[1], newEdgeVertices[0]));
				vNewFacesOut.push_back(new ConcreteTriangleType(vrts[2], newEdgeVertices[2], newEdgeVertices[1]));
				vNewFacesOut.push_back(new ConcreteTriangleType(newEdgeVertices[0], newEdgeVertices[1], newEdgeVertices[2]));
				return true;
			}

		}
	}

	return false;
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
collapse_edge(std::vector<Face*>& vNewFacesOut,
				int edgeIndex, Vertex* newVertex,
				Vertex** pSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();
	return true;
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<Vertex*>& vNewEdgeVertices,
				Vertex** pSubstituteVertices)
{
	if(vNewEdgeVertices.size() > BaseClass::num_edges())
	{
		assert(!"WARNING in Triangle::collapse_edges(...): bad number of newEdgeVertices.");
		LOG("WARNING in Triangle::collapse_edges(...): bad number of newEdgeVertices.");
		return false;
	}

//	check if there is a valid entry in vNewEdgeVertices
	bool bGotOne = false;
	for(uint i = 0; i < vNewEdgeVertices.size(); ++i)
	{
		if(vNewEdgeVertices[i] != NULL)
		{
			bGotOne = true;
			break;
		}
	}

	if(!bGotOne)
	{
		assert(!"WARNING in Triangle::collapse:edges(...): no new vertex was specified.");
		LOG("WARNING in Triangle::collapse:edges(...): no new vertex was specified.");
		return false;
	}

//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();
	return true;
}

//	BEGIN Depreciated
template <class ConcreteTriangleType, class BaseClass>
void
CustomTriangle<ConcreteTriangleType, BaseClass>::
create_faces_by_edge_split(int splitEdgeIndex,
								Vertex* newVertex,
								std::vector<Face*>& vNewFacesOut,
								Vertex** pSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 3)) && "ERROR in Triangle::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	we have to find the indices ind0, ind1, ind2, where
//	ind0 is the index of the vertex on e before newVertex,
//	ind1 is the index of the vertex on e after newVertex
//	and ind2 is the index of the vertex not located on e.

	int ind0 = splitEdgeIndex;
	int ind1 = (ind0 + 1) % 3;
	int ind2 = (ind1 + 1) % 3;

	vNewFacesOut.push_back(new Triangle(vrts[ind0], newVertex, vrts[ind2]));
	vNewFacesOut.push_back(new Triangle(newVertex, vrts[ind1], vrts[ind2]));
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	QUADRILATERALS

////////////////////////////////////////////////////////////////////////
//	QuadrilateralDescriptor
QuadrilateralDescriptor::QuadrilateralDescriptor(const QuadrilateralDescriptor& qd)
{
	m_vertex[0] = qd.vertex(0);
	m_vertex[1] = qd.vertex(1);
	m_vertex[2] = qd.vertex(2);
	m_vertex[3] = qd.vertex(3);
}

QuadrilateralDescriptor::QuadrilateralDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
}

////////////////////////////////////////////////////////////////////////
//	Quad

template <class ConcreteQuadrilateralType, class BaseClass>
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
CustomQuadrilateral(const QuadrilateralDescriptor& qd)
{
	m_vertices[0] = qd.vertex(0);
	m_vertices[1] = qd.vertex(1);
	m_vertices[2] = qd.vertex(2);
	m_vertices[3] = qd.vertex(3);
}

template <class ConcreteQuadrilateralType, class BaseClass>
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
CustomQuadrilateral(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
}

template <class ConcreteQuadrilateralType, class BaseClass>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
get_opposing_side(EdgeVertices* e, EdgeDescriptor& edOut) const
{
	int localInd = Face::get_local_side_index(e);
	if(localInd == -1){
		return false;
	}

	edge_desc((localInd + 2) % 4, edOut);
	return true;
}

template <class ConcreteQuadrilateralType, class BaseClass>
std::pair<GridBaseObjectId, int>
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
get_opposing_object(Vertex* vrt) const
{
	for(int i = 0; i < 4; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(VERTEX, (i + 2) % 4);
		}
	}

	UG_THROW("The given vertex is not contained in the given face.");
}

template <class ConcreteQuadrilateralType, class BaseClass>
void
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
create_faces_by_edge_split(int splitEdgeIndex,
							Vertex* newVertex,
							std::vector<Face*>& vNewFacesOut,
							Vertex** pSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 4)) && "ERROR in Quadrilateral::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	we have to find the indices ind0, ind1, ind2, where
//	ind0 is the index of the vertex on e before newVertex,
//	ind1 is the index of the vertex on e after newVertex
//	and ind2 and ind3 are the indices of the vertices not located on e.
	int ind0 = splitEdgeIndex; //edge-index equals the index of its first vertex.
	int ind1 = (ind0 + 1) % 4;
	int ind2 = (ind0 + 2) % 4;
	int ind3 = (ind0 + 3) % 4;

	TriangleDescriptor td;

//	edge-split generates 3 triangles
	vNewFacesOut.push_back(new Triangle(vrts[ind0], newVertex, vrts[ind3]));
	vNewFacesOut.push_back(new Triangle(newVertex, vrts[ind1], vrts[ind2]));
	vNewFacesOut.push_back(new Triangle(vrts[ind3], newVertex, vrts[ind2]));

//	we're done.
}

////////////////////////////////////////////////////////////////////////
//	Quadrilateral::refine
template <class ConcreteQuadrilateralType, class BaseClass>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
refine(std::vector<Face*>& vNewFacesOut,
		Vertex** newFaceVertexOut,
		Vertex** edgeVrts,
		Vertex* newFaceVertex,
		Vertex** pSubstituteVertices)
{
//TODO: complete quad refine
	*newFaceVertexOut = newFaceVertex;
	vNewFacesOut.clear();
	
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	check which edges have to be refined and perform the required operations.
//	get the number of new vertices.
	uint numNewVrts = 0;
	for(uint i = 0; i < 4; ++i)
	{
		if(edgeVrts[i] != NULL)
			++numNewVrts;
	}

	switch(numNewVrts)
	{
		case 0:
		//	create four new triangles
			vNewFacesOut.push_back(new Triangle(vrts[0], vrts[1], newFaceVertex));
			vNewFacesOut.push_back(new Triangle(vrts[1], vrts[2], newFaceVertex));
			vNewFacesOut.push_back(new Triangle(vrts[2], vrts[3], newFaceVertex));
			vNewFacesOut.push_back(new Triangle(vrts[3], vrts[0], newFaceVertex));
			return true;
			
		case 1:
		{
			if(newFaceVertex){
				assert(!"Problem in CustomQuadrilateral::refine: refine with newFaceVertex and one newEdgeVertex is not yet supported.");
				return false;
			}
			
			int iNew = -1;
			for(int i = 0; i < 4; ++i){
				if(edgeVrts[i]){
					iNew = i;
					break;
				}
			}
			
		//	the corners in a local ordering relative to iNew. Keeping ccw order.
			Vertex* corner[4];
			ReorderCornersCCW(corner, vrts, 4, (iNew + 3) % 4);

		//	create the new elements
			vNewFacesOut.push_back(new Triangle(corner[0], corner[1], edgeVrts[iNew]));
			vNewFacesOut.push_back(new Triangle(corner[0], edgeVrts[iNew], corner[3]));
			vNewFacesOut.push_back(new Triangle(corner[3], edgeVrts[iNew], corner[2]));

			return true;
		}

		case 2:
		{
			if(newFaceVertex){
				assert(!"Problem in CustomQuadrilateral::refine: refine with newFaceVertex and two newEdgeVertices is not yet supported.");
				return false;
			}

		//	there are two cases (refined edges are adjacent or opposite to each other)
			int iNew[2];
			int counter = 0;
			for(int i = 0; i < 4; ++i){
				if(edgeVrts[i]){
					iNew[counter] = i;
					++counter;
				}
			}

		//	corners will be filled later on
			Vertex* corner[4];

		//	check which case applies
			if(iNew[1] - iNew[0] == 2){
			//	edges are opposite to each other
			//	the corners in a local ordering relative to iNew[0]. Keeping ccw order.
				ReorderCornersCCW(corner, vrts, 4, (iNew[0] + 3) % 4);
				
			//	create new faces
				vNewFacesOut.push_back(new ConcreteQuadrilateralType(corner[0], corner[1],
													edgeVrts[iNew[0]], edgeVrts[iNew[1]]));

				vNewFacesOut.push_back(new ConcreteQuadrilateralType(edgeVrts[iNew[1]], edgeVrts[iNew[0]],
														corner[2], corner[3]));
			}
			else{
			//	edges are adjacent
			//	we want that (iNew[0] + 1) % 4 == iNew[1]
				if((iNew[0] + 1) % 4 != iNew[1])
					swap(iNew[0], iNew[1]);
			
			//	the corners in a local ordering relative to iNew[0]. Keeping ccw order.
				ReorderCornersCCW(corner, vrts, 4, (iNew[0] + 3) % 4);

			//	create new faces
				vNewFacesOut.push_back(new Triangle(corner[0], corner[1], edgeVrts[iNew[0]]));
				vNewFacesOut.push_back(new Triangle(edgeVrts[iNew[0]], corner[2], edgeVrts[iNew[1]]));
				vNewFacesOut.push_back(new Triangle(corner[3], corner[0], edgeVrts[iNew[1]]));
				vNewFacesOut.push_back(new Triangle(corner[0], edgeVrts[iNew[0]], edgeVrts[iNew[1]]));
			}

			return true;
		}

		case 3:
		{
			if(newFaceVertex){
				assert(!"Problem in CustomQuadrilateral::refine: refine with newFaceVertex and three newEdgeVertices is not yet supported.");
				return false;
			}

			int iFree = -1;
			for(int i = 0; i < 4; ++i){
				if(!edgeVrts[i]){
					iFree = i;
					break;
				}
			}
			
		//	the vertices on the edges:
			Vertex* nvrts[3];
			nvrts[0] = edgeVrts[(iFree + 1) % 4];
			nvrts[1] = edgeVrts[(iFree + 2) % 4];
			nvrts[2] = edgeVrts[(iFree + 3) % 4];

		//	the corners in a local ordering relative to iNew. Keeping ccw order.
			Vertex* corner[4];
			ReorderCornersCCW(corner, vrts, 4, (iFree + 1) % 4);

		//	create the faces
			vNewFacesOut.push_back(new ConcreteQuadrilateralType(corner[0], nvrts[0], nvrts[2], corner[3]));
			vNewFacesOut.push_back(new Triangle(corner[1], nvrts[1], nvrts[0]));
			vNewFacesOut.push_back(new Triangle(corner[2], nvrts[2], nvrts[1]));
			vNewFacesOut.push_back(new Triangle(nvrts[0], nvrts[1], nvrts[2]));

			return true;
		}

		case 4:
		{
		//	we'll create 4 new quads. create a new center if required.
			if(!newFaceVertex)
				newFaceVertex = new RegularVertex;

			*newFaceVertexOut = newFaceVertex;
		
			vNewFacesOut.push_back(new ConcreteQuadrilateralType(vrts[0], edgeVrts[0], newFaceVertex, edgeVrts[3]));
			vNewFacesOut.push_back(new ConcreteQuadrilateralType(vrts[1], edgeVrts[1], newFaceVertex, edgeVrts[0]));
			vNewFacesOut.push_back(new ConcreteQuadrilateralType(vrts[2], edgeVrts[2], newFaceVertex, edgeVrts[1]));
			vNewFacesOut.push_back(new ConcreteQuadrilateralType(vrts[3], edgeVrts[3], newFaceVertex, edgeVrts[2]));
			return true;
		}
	}

	return false;
}

template <class ConcreteQuadrilateralType, class BaseClass>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
collapse_edge(std::vector<Face*>& vNewFacesOut,
				int edgeIndex, Vertex* newVertex,
				Vertex** pSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();

//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	the collapsed edge connects vertices at ind0 and ind1.
	int ind0 = edgeIndex; //edge-index equals the index of its first vertex.
	int ind1 = (ind0 + 1) % 4;
	int ind2 = (ind1 + 1) % 4;
	int ind3 = (ind2 + 1) % 4;

//	ind0 and ind1 will be replaced by newVertex.
	vNewFacesOut.push_back(new Triangle(newVertex, vrts[ind2], vrts[ind3]));
	return true;
}

template <class ConcreteQuadrilateralType, class BaseClass>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass>::
collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<Vertex*>& vNewEdgeVertices,
				Vertex** pSubstituteVertices)
{
	if(vNewEdgeVertices.size() > BaseClass::num_edges())
	{
		assert(!"WARNING in Quadrilateral::collapse_edges(...): bad number of newEdgeVertices.");
		LOG("WARNING in Quadrilateral::collapse_edges(...): bad number of newEdgeVertices.");
		return false;
	}

	vNewFacesOut.clear();

//	check if there is a valid entry in vNewEdgeVertices
	int collapseIndex = -1;
	int numCollapses = 0;
	for(uint i = 0; i < 4; ++i)
	{
		if(i < vNewEdgeVertices.size())
		{
			if(vNewEdgeVertices[i] != NULL)
			{
				++numCollapses;
				collapseIndex = i;
			}
		}
	}

//	if more than 1 edge is collapsed, nothing remains.
	if(numCollapses == 0)
	{
		assert(!"WARNING in Quadrilateral::collapse:edges(...): no new vertex was specified.");
		LOG("WARNING in Quadrilateral::collapse:edges(...): no new vertex was specified.");
		return false;
	}
	else if(numCollapses == 1)
	{
	//	call collapse_edge with the edges index.
		collapse_edge(vNewFacesOut, collapseIndex, vNewEdgeVertices[collapseIndex], pSubstituteVertices);
	}

	return true;
}

//	explicit instantiation
template class CustomTriangle<Triangle, Face>;
template class CustomTriangle<ConstrainedTriangle, ConstrainedFace>;
template class CustomTriangle<ConstrainingTriangle, ConstrainingFace>;

template class CustomQuadrilateral<Quadrilateral, Face>;
template class CustomQuadrilateral<ConstrainedQuadrilateral, ConstrainedFace>;
template class CustomQuadrilateral<ConstrainingQuadrilateral, ConstrainingFace>;


////////////////////////////////////////////////////////////////////////////////
//	CONSTRAINING FACE
template <> size_t
ConstrainingFace::
num_constrained<Vertex>() const
{
	return num_constrained_vertices();
}

template <> size_t
ConstrainingFace::
num_constrained<Edge>() const
{
	return num_constrained_edges();
}

template <> size_t
ConstrainingFace::
num_constrained<Face>() const
{
	return num_constrained_faces();
}


template <> Vertex*
ConstrainingFace::
constrained<Vertex>(size_t ind) const
{
	return constrained_vertex(ind);
}

template <> Edge*
ConstrainingFace::
constrained<Edge>(size_t ind) const
{
	return constrained_edge(ind);
}

template <> Face*
ConstrainingFace::
constrained<Face>(size_t ind) const
{
	return constrained_face(ind);
}

}//	end of namespace
