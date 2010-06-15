//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d04

#include <vector>
#include <algorithm>
#include "geometric_objects.h"
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
bool ReorderCornersCCW(VertexBase** cornersOut, VertexBase** const cornersIn,
					   int numCorners, int firstCorner)
{
	cornersOut[0] = cornersIn[firstCorner];
	for(int i = 1; i < numCorners; ++i)
		cornersOut[i] = cornersIn[(firstCorner + i) % numCorners];
	return true;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	EDGES

////////////////////////////////////////////////////////////////////////
//	Edge::refine
bool Edge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex,
				VertexBase** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new Edge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new Edge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new Edge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new Edge(newVertex, vertex(1)));
	}
	return true;
}

bool Edge::refine(std::vector<Edge*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge::refine
bool ConstrainedEdge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new ConstrainedEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new ConstrainedEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, vertex(1)));
	}
	return true;
}

bool ConstrainedEdge::refine(std::vector<ConstrainedEdge*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut),
				newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge::refine
bool ConstrainingEdge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new ConstrainingEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new ConstrainingEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, vertex(1)));
	}
	return true;
}

bool ConstrainingEdge::refine(std::vector<ConstrainingEdge*>& vNewEdgesOut, VertexBase* newVertex,
				  VertexBase** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
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

TriangleDescriptor::TriangleDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3)
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
	BaseClass::set_num_vertices(3);
	BaseClass::m_vertices[0] = td.vertex(0);
	BaseClass::m_vertices[1] = td.vertex(1);
	BaseClass::m_vertices[2] = td.vertex(2);
}

template <class ConcreteTriangleType, class BaseClass>
CustomTriangle<ConcreteTriangleType, BaseClass>::
CustomTriangle(VertexBase* v1, VertexBase* v2, VertexBase* v3)
{
	BaseClass::set_num_vertices(3);
	BaseClass::m_vertices[0] = v1;
	BaseClass::m_vertices[1] = v2;
	BaseClass::m_vertices[2] = v3;
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
refine(std::vector<Face*>& vNewFacesOut,
		VertexBase** newFaceVertexOut,
		VertexBase** newEdgeVertices,
		VertexBase* newFaceVertex,
		VertexBase** pSubstituteVertices)
{
//TODO: complete triangle refine

	*newFaceVertexOut = NULL;
	vNewFacesOut.clear();

//	handle substitute vertices.
	VertexBase** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = BaseClass::m_vertices;

//	if newFaceVertex is specified, then create three sub-triangles and
//	refine each. If not then refine the triangle depending
	if(newFaceVertex)
	{
		assert(!"PROBLEM in Triangle::refine(...): refine with newFaceVertex not yet implemented.");
		return false;
	}
	else
	{
	//	get the number of new vertices.
		uint numNewVrts = 0;
		for(uint i = 0; i < 3; ++i)
		{
			if(newEdgeVertices[i] != NULL)
				++numNewVrts;
		}

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
refine_regular(std::vector<Face*>& vNewFacesOut,
				VertexBase** newFaceVertexOut,
				std::vector<VertexBase*>& vNewEdgeVertices,
				VertexBase* newFaceVertex,
				const VertexBase& prototypeVertex,
				VertexBase** pSubstituteVertices)
{
	assert(vNewEdgeVertices.size() == 3 && "wrong number of newEdgeVertices.");
	assert(newFaceVertex == NULL && "regular triangle refine doesn't use newFaceVertex.");
	*newFaceVertexOut = NULL;
	
	return refine(vNewFacesOut, newFaceVertexOut, &vNewEdgeVertices.front(), *newFaceVertexOut, pSubstituteVertices);		
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
collapse_edge(std::vector<Face*>& vNewFacesOut,
				int edgeIndex, VertexBase* newVertex,
				VertexBase** pSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();
	return true;
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<VertexBase*>& vNewEdgeVertices,
				VertexBase** pSubstituteVertices)
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
								VertexBase* newVertex,
								std::vector<Face*>& vNewFacesOut,
								VertexBase** pSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 3)) && "ERROR in Triangle::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	VertexBase** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = BaseClass::m_vertices;

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

//	explicit instantiation
template class CustomTriangle<Triangle, Face>;
template class CustomTriangle<ConstrainedTriangle, ConstrainedFace>;
template class CustomTriangle<ConstrainingTriangle, ConstrainingFace>;


/*
////////////////////////////////////////////////////////////////////////
//	Triangle
Triangle::Triangle(const TriangleDescriptor& td) : CustomFace<3, 0>()
{
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
}

Triangle::Triangle(VertexBase* v1, VertexBase* v2, VertexBase* v3) : CustomFace<3, 0>()
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	Triangle
void Triangle::create_faces_by_edge_split(int splitEdgeIndex,
							VertexBase* newVertex,
							std::vector<Face*>& vNewFacesOut,
							std::vector<VertexBase*>* pvSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 3)) && "ERROR in Triangle::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	vector<VertexBase*> localVertices;
	if(!pvSubstituteVertices)
	{
		pvSubstituteVertices = &localVertices;
		localVertices.push_back(vertex(0));
		localVertices.push_back(vertex(1));
		localVertices.push_back(vertex(2));
	}
	else
		assert(pvSubstituteVertices->size() == 3 && "ERROR in Triangle::create_faces_by_edge_split(...): wrong size of pvSubstituteVertices.");

	vector<VertexBase*>& vVrts = *pvSubstituteVertices;

//	we have to find the indices ind0, ind1, ind2, where
//	ind0 is the index of the vertex on e before newVertex,
//	ind1 is the index of the vertex on e after newVertex
//	and ind2 is the index of the vertex not located on e.

	int ind0 = splitEdgeIndex;
	int ind1 = (ind0 + 1) % 3;
	int ind2 = (ind1 + 1) % 3;

	vNewFacesOut.push_back(new Triangle(vVrts[ind0], newVertex, vVrts[ind2]));
	vNewFacesOut.push_back(new Triangle(newVertex, vVrts[ind1], vVrts[ind2]));

//	we're done.
}

////////////////////////////////////////////////////////////////////////
//	Triangle::refine
bool Triangle::refine(std::vector<Face*>& vNewFacesOut,
						std::vector<VertexBase*>& vNewEdgeVertices,
						VertexBase* newFaceVertex,
						std::vector<VertexBase*>* pvSubstituteVertices)
{
//TODO: complete triangle refine

	vNewFacesOut.clear();

//	handle substitute vertices.
	vector<VertexBase*>	vTmpVrts;
	if(!pvSubstituteVertices)
	{
		vTmpVrts.resize(3);
		vTmpVrts[0] = vertex(0);
		vTmpVrts[1] = vertex(1);
		vTmpVrts[2] = vertex(2);
		pvSubstituteVertices = &vTmpVrts;
	}

//	use this vertex-vector during this algorithm.
	vector<VertexBase*>& vVrts = *pvSubstituteVertices;

//	if newFaceVertex is specified, then create three sub-triangles and
//	refine each. If not then refine the triangle depending
	if(newFaceVertex)
	{
		assert(!"PROBLEM in Triangle::refine(...): refine with newFaceVertex not yet implemented.");
		return false;
	}
	else
	{
	//	get the number of new vertices.
		uint numNewVrts = 0;
		for(uint i = 0; (i < 3) && (i < vNewEdgeVertices.size()); ++i)
		{
			if(vNewEdgeVertices[i] != NULL)
				++numNewVrts;
		}

		switch(numNewVrts)
		{
			case 1:
			{
				assert(!"PROBLEM in Triangle::refine(...): refine with 1 new edge vertex not yet implemented.");
				return false;
			}

			case 2:
			{
				assert(!"PROBLEM in Triangle::refine(...): refine with 2 new edge vertices not yet implemented.");
				return false;
			}

			case 3:
			{
			//	perform regular refine.
				vNewFacesOut.push_back(new Triangle(vVrts[0], vNewEdgeVertices[0], vNewEdgeVertices[2]));
				vNewFacesOut.push_back(new Triangle(vVrts[1], vNewEdgeVertices[1], vNewEdgeVertices[0]));
				vNewFacesOut.push_back(new Triangle(vVrts[2], vNewEdgeVertices[2], vNewEdgeVertices[1]));
				vNewFacesOut.push_back(new Triangle(vNewEdgeVertices[0], vNewEdgeVertices[1], vNewEdgeVertices[2]));
				return true;
			}

		}
	}

	return false;
}

bool Triangle::collapse_edge(std::vector<Face*>& vNewFacesOut,
					int edgeIndex, VertexBase* newVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();
	return true;
}

bool Triangle::collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<VertexBase*>& vNewEdgeVertices,
				std::vector<VertexBase*>* pvSubstituteVertices)
{
	if(vNewEdgeVertices.size() > num_edges())
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
*/
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

QuadrilateralDescriptor::QuadrilateralDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
}

////////////////////////////////////////////////////////////////////////
//	Quad
Quadrilateral::Quadrilateral(const QuadrilateralDescriptor& qd) : CustomFace<4, 1>()
{
	m_vertices[0] = qd.vertex(0);
	m_vertices[1] = qd.vertex(1);
	m_vertices[2] = qd.vertex(2);
	m_vertices[3] = qd.vertex(3);
}

Quadrilateral::Quadrilateral(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4) : CustomFace<4, 1>()
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
}

void Quadrilateral::create_faces_by_edge_split(int splitEdgeIndex,
							VertexBase* newVertex,
							std::vector<Face*>& vNewFacesOut,
							VertexBase** pSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 4)) && "ERROR in Quadrilateral::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	VertexBase** vrts;
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
bool Quadrilateral::refine(std::vector<Face*>& vNewFacesOut,
							VertexBase** newFaceVertexOut,
							VertexBase** edgeVrts,
							VertexBase* newFaceVertex,
							VertexBase** pSubstituteVertices)
{
//TODO: complete quad refine
	*newFaceVertexOut = NULL;
	vNewFacesOut.clear();
	
//	handle substitute vertices.
	VertexBase** vrts;
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
		case 1:
		{
			int iNew = -1;
			for(int i = 0; i < 4; ++i){
				if(edgeVrts[i]){
					iNew = i;
					break;
				}
			}
			
		//	the corners in a local ordering relative to iNew. Keeping ccw order.
			VertexBase* corner[4];
			ReorderCornersCCW(corner, vrts, 4, (iNew + 3) % 4);

		//	create the new elements
			vNewFacesOut.push_back(new Triangle(corner[0], corner[1], edgeVrts[iNew]));
			vNewFacesOut.push_back(new Triangle(corner[0], edgeVrts[iNew], corner[3]));
			vNewFacesOut.push_back(new Triangle(corner[3], edgeVrts[iNew], corner[2]));

			return true;
		}

		case 2:
		{
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
			VertexBase* corner[4];

		//	check which case applies
			if(iNew[1] - iNew[0] == 2){
			//	edges are opposite to each other
			//	the corners in a local ordering relative to iNew[0]. Keeping ccw order.
				ReorderCornersCCW(corner, vrts, 4, (iNew[0] + 3) % 4);
				
			//	create new faces
				vNewFacesOut.push_back(new Quadrilateral(corner[0], corner[1],
													edgeVrts[iNew[0]], edgeVrts[iNew[1]]));

				vNewFacesOut.push_back(new Quadrilateral(edgeVrts[iNew[1]], edgeVrts[iNew[0]],
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
			int iFree = -1;
			for(int i = 0; i < 4; ++i){
				if(!edgeVrts[i]){
					iFree = i;
					break;
				}
			}
			
		//	the vertices on the edges:
			VertexBase* nvrts[3];
			nvrts[0] = edgeVrts[(iFree + 1) % 4];
			nvrts[1] = edgeVrts[(iFree + 2) % 4];
			nvrts[2] = edgeVrts[(iFree + 3) % 4];

		//	the corners in a local ordering relative to iNew. Keeping ccw order.
			VertexBase* corner[4];
			ReorderCornersCCW(corner, vrts, 4, (iFree + 1) % 4);

		//	create the faces
			vNewFacesOut.push_back(new Quadrilateral(corner[0], nvrts[0], nvrts[2], corner[3]));
			vNewFacesOut.push_back(new Triangle(corner[1], nvrts[1], nvrts[0]));
			vNewFacesOut.push_back(new Triangle(corner[2], nvrts[2], nvrts[1]));
			vNewFacesOut.push_back(new Triangle(nvrts[0], nvrts[1], nvrts[2]));

			return true;
		}

		case 4:
		{
		//	we'll create 4 new quads. create a new center if required.
			if(!newFaceVertex)
				newFaceVertex = new Vertex;

			*newFaceVertexOut = newFaceVertex;
		
			vNewFacesOut.push_back(new Quadrilateral(vrts[0], edgeVrts[0], newFaceVertex, edgeVrts[3]));
			vNewFacesOut.push_back(new Quadrilateral(vrts[1], edgeVrts[1], newFaceVertex, edgeVrts[0]));
			vNewFacesOut.push_back(new Quadrilateral(vrts[2], edgeVrts[2], newFaceVertex, edgeVrts[1]));
			vNewFacesOut.push_back(new Quadrilateral(vrts[3], edgeVrts[3], newFaceVertex, edgeVrts[2]));
			return true;
		}
	}

	return false;
}

bool Quadrilateral::refine_regular(std::vector<Face*>& vNewFacesOut,
									VertexBase** newFaceVertexOut,
									std::vector<VertexBase*>& vNewEdgeVertices,
									VertexBase* newFaceVertex,
									const VertexBase& prototypeVertex,
									VertexBase** pSubstituteVertices)
{
	assert(vNewEdgeVertices.size() == 4 && "wrong number of newEdgeVertices.");

//	call refine with the correct parameters.
	if(newFaceVertex == NULL)
	{
		*newFaceVertexOut = newFaceVertex = dynamic_cast<VertexBase*>(prototypeVertex.create_empty_instance());
		assert(newFaceVertex && "create_empty_instance of a derivate of VertexBase did not create a derivate of VertexBase");
	}
	else
		*newFaceVertexOut = NULL;


	return refine(vNewFacesOut, newFaceVertexOut, &vNewEdgeVertices.front(), *newFaceVertexOut, pSubstituteVertices);
}

bool Quadrilateral::collapse_edge(std::vector<Face*>& vNewFacesOut,
					int edgeIndex, VertexBase* newVertex,
					VertexBase** pSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();

//	handle substitute vertices.
	VertexBase** vrts;
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

bool Quadrilateral::collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<VertexBase*>& vNewEdgeVertices,
				VertexBase** pSubstituteVertices)
{
	if(vNewEdgeVertices.size() > num_edges())
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

}//	end of namespace
