//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d04

#include <vector>
#include "geometric_objects.h"
#include "algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	EDGES

////////////////////////////////////////////////////////////////////////
//	Edge::refine
bool Edge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex)
{
	vNewEdgesOut.clear();
	vNewEdgesOut.push_back(new Edge(vertex(0), newVertex));
	vNewEdgesOut.push_back(new Edge(newVertex, vertex(1)));
	return true;
}

bool Edge::refine(std::vector<Edge*>& vNewEdgesOut, VertexBase* newVertex)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut), newVertex);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge::refine
bool ConstrainedEdge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex)
{
	vNewEdgesOut.clear();
	vNewEdgesOut.push_back(new ConstrainedEdge(vertex(0), newVertex));
	vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, vertex(1)));
	return true;
}

bool ConstrainedEdge::refine(std::vector<ConstrainedEdge*>& vNewEdgesOut, VertexBase* newVertex)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut), newVertex);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge::refine
bool ConstrainingEdge::refine(std::vector<EdgeBase*>& vNewEdgesOut, VertexBase* newVertex)
{
	vNewEdgesOut.clear();
	vNewEdgesOut.push_back(new ConstrainingEdge(vertex(0), newVertex));
	vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, vertex(1)));
	return true;
}

bool ConstrainingEdge::refine(std::vector<ConstrainingEdge*>& vNewEdgesOut, VertexBase* newVertex)
{
	return refine(reinterpret_cast<std::vector<EdgeBase*>&>(vNewEdgesOut), newVertex);
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
	BaseClass::m_vertices.resize(3);
	BaseClass::m_vertices[0] = td.vertex(0);
	BaseClass::m_vertices[1] = td.vertex(1);
	BaseClass::m_vertices[2] = td.vertex(2);
}

template <class ConcreteTriangleType, class BaseClass>
CustomTriangle<ConcreteTriangleType, BaseClass>::
CustomTriangle(VertexBase* v1, VertexBase* v2, VertexBase* v3)
{
	BaseClass::m_vertices.resize(3);
	BaseClass::m_vertices[0] = v1;
	BaseClass::m_vertices[1] = v2;
	BaseClass::m_vertices[2] = v3;
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
refine(std::vector<Face*>& vNewFacesOut,
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
		vTmpVrts[0] = BaseClass::vertex(0);
		vTmpVrts[1] = BaseClass::vertex(1);
		vTmpVrts[2] = BaseClass::vertex(2);
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
				vNewFacesOut.push_back(new ConcreteTriangleType(vVrts[0], vNewEdgeVertices[0], vNewEdgeVertices[2]));
				vNewFacesOut.push_back(new ConcreteTriangleType(vVrts[1], vNewEdgeVertices[1], vNewEdgeVertices[0]));
				vNewFacesOut.push_back(new ConcreteTriangleType(vVrts[2], vNewEdgeVertices[2], vNewEdgeVertices[1]));
				vNewFacesOut.push_back(new ConcreteTriangleType(vNewEdgeVertices[0], vNewEdgeVertices[1], vNewEdgeVertices[2]));
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
				std::vector<VertexBase*>* pvSubstituteVertices)
{
	assert(vNewEdgeVertices.size() == 3 && "wrong number of newEdgeVertices.");
	assert(newFaceVertex == NULL && "regular triangle refine doesn't use newFaceVertex.");
	*newFaceVertexOut = NULL;
	return refine(vNewFacesOut, vNewEdgeVertices, newFaceVertex, pvSubstituteVertices);
}

template <class ConcreteTriangleType, class BaseClass>
bool
CustomTriangle<ConcreteTriangleType, BaseClass>::
collapse_edge(std::vector<Face*>& vNewFacesOut,
				int edgeIndex, VertexBase* newVertex,
				std::vector<VertexBase*>* pvSubstituteVertices)
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
				std::vector<VertexBase*>* pvSubstituteVertices)
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
		localVertices.push_back(BaseClass::vertex(0));
		localVertices.push_back(BaseClass::vertex(1));
		localVertices.push_back(BaseClass::vertex(2));
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
							std::vector<VertexBase*>* pvSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 4)) && "ERROR in Quadrilateral::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	vector<VertexBase*> localVertices;
	if(!pvSubstituteVertices)
	{
		localVertices.reserve(4);
		pvSubstituteVertices = &localVertices;
		localVertices.push_back(vertex(0));
		localVertices.push_back(vertex(1));
		localVertices.push_back(vertex(2));
		localVertices.push_back(vertex(3));
	}
	else
		assert(pvSubstituteVertices->size() == 4 && "ERROR in Quadrilateral::create_faces_by_edge_split(...): wrong size of pvSubstituteVertices.");

	vector<VertexBase*>& vVrts = *pvSubstituteVertices;

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
	vNewFacesOut.push_back(new Triangle(vVrts[ind0], newVertex, vVrts[ind3]));
	vNewFacesOut.push_back(new Triangle(newVertex, vVrts[ind1], vVrts[ind2]));
	vNewFacesOut.push_back(new Triangle(vVrts[ind3], newVertex, vVrts[ind2]));

//	we're done.
}

////////////////////////////////////////////////////////////////////////
//	Quadrilateral::refine
bool Quadrilateral::refine(std::vector<Face*>& vNewFacesOut,
						std::vector<VertexBase*>& vNewEdgeVertices,
						VertexBase* newFaceVertex,
						std::vector<VertexBase*>* pvSubstituteVertices)
{
//TODO: complete quad refine

	vNewFacesOut.clear();

//	handle substitute vertices.
	vector<VertexBase*>	vTmpVrts;
	if(!pvSubstituteVertices)
	{
		vTmpVrts.resize(4);
		vTmpVrts[0] = vertex(0);
		vTmpVrts[1] = vertex(1);
		vTmpVrts[2] = vertex(2);
		vTmpVrts[3] = vertex(3);
		pvSubstituteVertices = &vTmpVrts;
	}

//	use this vertex-vector during this algorithm.
	vector<VertexBase*>& vVrts = *pvSubstituteVertices;

//	check which edges have to be refined and perform the required operations.
	{
	//	get the number of new vertices.
		uint numNewVrts = 0;
		for(uint i = 0; (i < 4) && (i < vNewEdgeVertices.size()); ++i)
		{
			if(vNewEdgeVertices[i] != NULL)
				++numNewVrts;
		}

		switch(numNewVrts)
		{
			case 1:
			{
				assert(!"PROBLEM in Quadrilateral::refine(...): refine with 1 new edge vertex not yet implemented.");
				return false;
			}

			case 2:
			{
				assert(!"PROBLEM in Quadrilateral::refine(...): refine with 2 new edge vertices not yet implemented.");
				return false;
			}

			case 3:
			{
				assert(!"PROBLEM in Quadrilateral::refine(...): refine with 3 new edge vertices not yet implemented.");
				return false;
			}

			case 4:
			{
			//	if a center-vertex has been specified we'll create 4 new quads.
			//	if not, one quad and 4 triangles will be created.
				if(newFaceVertex)
				{
					vNewFacesOut.push_back(new Quadrilateral(vVrts[0], vNewEdgeVertices[0], newFaceVertex, vNewEdgeVertices[3]));
					vNewFacesOut.push_back(new Quadrilateral(vVrts[1], vNewEdgeVertices[1], newFaceVertex, vNewEdgeVertices[0]));
					vNewFacesOut.push_back(new Quadrilateral(vVrts[2], vNewEdgeVertices[2], newFaceVertex, vNewEdgeVertices[1]));
					vNewFacesOut.push_back(new Quadrilateral(vVrts[3], vNewEdgeVertices[3], newFaceVertex, vNewEdgeVertices[2]));
					return true;
				}
				else
				{
					assert(!"PROBLEM in Quadrilateral::refine(...): refine with 4 new edge vertices but no new face vertex not yet implemented.");
					return false;
				}
			}
		}
	}

	return false;
}

bool Quadrilateral::refine_regular(std::vector<Face*>& vNewFacesOut,
									VertexBase** newFaceVertexOut,
									std::vector<VertexBase*>& vNewEdgeVertices,
									VertexBase* newFaceVertex,
									const VertexBase& prototypeVertex,
									std::vector<VertexBase*>* pvSubstituteVertices)
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

	return refine(vNewFacesOut, vNewEdgeVertices, newFaceVertex, pvSubstituteVertices);
}

bool Quadrilateral::collapse_edge(std::vector<Face*>& vNewFacesOut,
					int edgeIndex, VertexBase* newVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();

//	handle substitute vertices.
	vector<VertexBase*>	vTmpVrts;
	if(!pvSubstituteVertices)
	{
		vTmpVrts.resize(4);
		vTmpVrts[0] = vertex(0);
		vTmpVrts[1] = vertex(1);
		vTmpVrts[2] = vertex(2);
		vTmpVrts[3] = vertex(3);
		pvSubstituteVertices = &vTmpVrts;
	}

//	use this vertex-vector during this algorithm.
	vector<VertexBase*>& vVrts = *pvSubstituteVertices;

//	the collapsed edge connects vertices at ind0 and ind1.
	int ind0 = edgeIndex; //edge-index equals the index of its first vertex.
	int ind1 = (ind0 + 1) % 4;
	int ind2 = (ind1 + 1) % 4;
	int ind3 = (ind2 + 1) % 4;

//	ind0 and ind1 will be replaced by newVertex.
	vNewFacesOut.push_back(new Triangle(newVertex, vVrts[ind2], vVrts[ind3]));
	return true;
}

bool Quadrilateral::collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<VertexBase*>& vNewEdgeVertices,
				std::vector<VertexBase*>* pvSubstituteVertices)
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
		collapse_edge(vNewFacesOut, collapseIndex, vNewEdgeVertices[collapseIndex], pvSubstituteVertices);
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	VOLUMES

////////////////////////////////////////////////////////////////////////
//	TetrahedronDescriptor
TetrahedronDescriptor::TetrahedronDescriptor(const TetrahedronDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
	m_vertex[3] = td.vertex(3);
}

TetrahedronDescriptor::TetrahedronDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
}

////////////////////////////////////////////////////////////////////////
//	Tetrahedron
Tetrahedron::Tetrahedron(const TetrahedronDescriptor& td)
{
	m_vertices.resize(4);
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
}

Tetrahedron::Tetrahedron(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4)
{
	m_vertices.resize(4);
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
}

EdgeDescriptor Tetrahedron::edge(int index) const
{
	EdgeDescriptor ed;
	edge(index, ed);
	return ed;
}

void Tetrahedron::edge(int index, EdgeDescriptor& edOut) const
{
	switch(index)
	{
		case 0: edOut.set_vertices(m_vertices[0], m_vertices[1]);
				break;
		case 1: edOut.set_vertices(m_vertices[1], m_vertices[2]);
				break;
		case 2: edOut.set_vertices(m_vertices[2], m_vertices[0]);
				break;
		case 3: edOut.set_vertices(m_vertices[3], m_vertices[0]);
				break;
		case 4: edOut.set_vertices(m_vertices[3], m_vertices[1]);
				break;
		case 5: edOut.set_vertices(m_vertices[3], m_vertices[2]);
				break;
	}
}

uint Tetrahedron::num_edges() const
{
	return 6;
}

FaceDescriptor Tetrahedron::face(int index) const
{
	FaceDescriptor fd;
	face(index, fd);
	return fd;
}

void Tetrahedron::face(int index, FaceDescriptor& fdOut) const
{
	fdOut.set_num_vertices(3);
	switch(index)
		{
			case 0:
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[2]);
				fdOut.set_vertex(2, m_vertices[1]);
				break;
			case 1:
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[1]);
				fdOut.set_vertex(2, m_vertices[3]);
				break;
			case 2:
				fdOut.set_vertex(0, m_vertices[1]);
				fdOut.set_vertex(1, m_vertices[2]);
				fdOut.set_vertex(2, m_vertices[3]);
				break;
			case 3:
				fdOut.set_vertex(0, m_vertices[2]);
				fdOut.set_vertex(1, m_vertices[0]);
				fdOut.set_vertex(2, m_vertices[3]);
				break;
		}
}

uint Tetrahedron::num_faces() const
{
	return 4;
}

EdgeBase* Tetrahedron::create_edge(int index)
{
	EdgeDescriptor ed;
	edge(index, ed);
	return new Edge(ed);
}

Face* Tetrahedron::create_face(int index)
{
	FaceDescriptor fd;
	face(index, fd);
	return new Triangle(fd.vertex(0), fd.vertex(1), fd.vertex(2));
}

bool Tetrahedron::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, VertexBase* newVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//	if an edge of a tetrahedron is collapsed, nothing remains.
	vNewVolumesOut.clear();
	return true;
}

bool Tetrahedron::refine(std::vector<Volume*>& vNewVolumesOut,
					VertexBase** ppNewVertexOut,
					std::vector<VertexBase*>& vNewEdgeVertices,
					std::vector<VertexBase*>& vNewFaceVertices,
					VertexBase* newVolumeVertex,
					const VertexBase& prototypeVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//TODO: complete this refine method.
	vNewVolumesOut.clear();
	*ppNewVertexOut = NULL;

//	handle substitute vertices.
	vector<VertexBase*>	vTmpVrts;
	if(!pvSubstituteVertices)
	{
		vTmpVrts.resize(4);
		vTmpVrts[0] = vertex(0);
		vTmpVrts[1] = vertex(1);
		vTmpVrts[2] = vertex(2);
		vTmpVrts[3] = vertex(3);
		pvSubstituteVertices = &vTmpVrts;
	}

//	use this vertex-vector during this algorithm.
	vector<VertexBase*>& vVrts = *pvSubstituteVertices;

//	check which edges have to be refined and perform the required operations.
	{
	//	get the number of new vertices.
		uint numNewVrts = 0;
		for(uint i = 0; (i < 6) && (i < vNewEdgeVertices.size()); ++i)
		{
			if(vNewEdgeVertices[i] != NULL)
				++numNewVrts;
		}

		uint numNewFaceVrts = 0;
		for(uint i = 0; (i < 4) && (i < vNewFaceVertices.size()); ++i)
		{
			if(vNewFaceVertices[i] != NULL)
				++numNewFaceVrts;
		}

		assert(numNewFaceVrts == 0 && "PROBLEM in Tetrahedron::refine(): new face vertices are not yet supported!");

		switch(numNewVrts)
		{
			case 1:
			{
				assert(!"PROBLEM in Tetrahedron::refine(...): refine with 1 new edge vertex not yet implemented.");
				return false;
			}

			case 2:
			{
				assert(!"PROBLEM in Tetrahedron::refine(...): refine with 2 new edge vertices not yet implemented.");
				return false;
			}

			case 3:
			{
				assert(!"PROBLEM in Tetrahedron::refine(...): refine with 3 new edge vertices not yet implemented.");
				return false;
			}
			case 4:
			{
				assert(!"PROBLEM in Tetrahedron::refine(...): refine with 4 new edge vertices not yet implemented.");
				return false;
			}
			case 5:
			{
				assert(!"PROBLEM in Tetrahedron::refine(...): refine with 5 new edge vertices not yet implemented.");
				return false;
			}
			case 6:
			{
				if(!newVolumeVertex)
				{
					vNewVolumesOut.reserve(8);
					vNewVolumesOut.push_back(new Tetrahedron(vVrts[0], vNewEdgeVertices[0], vNewEdgeVertices[2], vNewEdgeVertices[3]));
					vNewVolumesOut.push_back(new Tetrahedron(vVrts[1], vNewEdgeVertices[1], vNewEdgeVertices[0], vNewEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(vVrts[2], vNewEdgeVertices[2], vNewEdgeVertices[1], vNewEdgeVertices[5]));
					vNewVolumesOut.push_back(new Tetrahedron(vNewEdgeVertices[0], vNewEdgeVertices[1], vNewEdgeVertices[2], vNewEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(vNewEdgeVertices[2], vNewEdgeVertices[0], vNewEdgeVertices[4], vNewEdgeVertices[3]));
					vNewVolumesOut.push_back(new Tetrahedron(vNewEdgeVertices[3], vNewEdgeVertices[2], vNewEdgeVertices[5], vNewEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(vNewEdgeVertices[2], vNewEdgeVertices[1], vNewEdgeVertices[5], vNewEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(vNewEdgeVertices[3], vNewEdgeVertices[4], vNewEdgeVertices[5], vVrts[3]));
					return true;
				}
				else
				{
					assert(!"PROBLEM in Tetrahedron::refine(...): refine with 4 new edge vertices and new new volume vertex not yet implemented.");
					return false;
				}
			}
		}
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
//	HexahedronDescriptor
HexahedronDescriptor::HexahedronDescriptor(const HexahedronDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
	m_vertex[3] = td.vertex(3);
	m_vertex[4] = td.vertex(4);
	m_vertex[5] = td.vertex(5);
	m_vertex[6] = td.vertex(6);
	m_vertex[7] = td.vertex(7);
}

HexahedronDescriptor::HexahedronDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4,
											VertexBase* v5, VertexBase* v6, VertexBase* v7, VertexBase* v8)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
	m_vertex[4] = v5;
	m_vertex[5] = v6;
	m_vertex[6] = v7;
	m_vertex[7] = v8;
}

////////////////////////////////////////////////////////////////////////
//	Hexahedron
Hexahedron::Hexahedron(const HexahedronDescriptor& td)
{
	m_vertices.resize(8);
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
	m_vertices[4] = td.vertex(4);
	m_vertices[5] = td.vertex(5);
	m_vertices[6] = td.vertex(6);
	m_vertices[7] = td.vertex(7);
}

Hexahedron::Hexahedron(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4,
						VertexBase* v5, VertexBase* v6, VertexBase* v7, VertexBase* v8)
{
	m_vertices.resize(8);
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
	m_vertices[4] = v5;
	m_vertices[5] = v6;
	m_vertices[6] = v7;
	m_vertices[7] = v8;
}

EdgeDescriptor Hexahedron::edge(int index) const
{
	EdgeDescriptor ed;
	edge(index, ed);
	return ed;
}

void Hexahedron::edge(int index, EdgeDescriptor& edOut) const
{
	if(index < 4)//base edges
		edOut.set_vertices(m_vertices[index], m_vertices[(index + 1) % 4]);
	else if(index < 8)//side edges
		edOut.set_vertices(m_vertices[index - 4], m_vertices[index]);
	else//top edges
		edOut.set_vertices(m_vertices[index - 4], m_vertices[(index - 7) % 4 + 4]);
}

uint Hexahedron::num_edges() const
{
	return 12;
}

FaceDescriptor Hexahedron::face(int index) const
{
	FaceDescriptor fd;
	face(index, fd);
	return fd;
}

void Hexahedron::face(int index, FaceDescriptor& fdOut) const
{
	fdOut.set_num_vertices(4);
	switch(index)
		{
			case 0://bottom
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[3]);
				fdOut.set_vertex(2, m_vertices[2]);
				fdOut.set_vertex(3, m_vertices[1]);
				break;
			case 1:
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[1]);
				fdOut.set_vertex(2, m_vertices[5]);
				fdOut.set_vertex(3, m_vertices[4]);
				break;
			case 2:
				fdOut.set_vertex(0, m_vertices[1]);
				fdOut.set_vertex(1, m_vertices[2]);
				fdOut.set_vertex(2, m_vertices[6]);
				fdOut.set_vertex(3, m_vertices[5]);
				break;
			case 3:
				fdOut.set_vertex(0, m_vertices[2]);
				fdOut.set_vertex(1, m_vertices[3]);
				fdOut.set_vertex(2, m_vertices[7]);
				fdOut.set_vertex(3, m_vertices[6]);
				break;
			case 4:
				fdOut.set_vertex(0, m_vertices[3]);
				fdOut.set_vertex(1, m_vertices[0]);
				fdOut.set_vertex(2, m_vertices[4]);
				fdOut.set_vertex(3, m_vertices[7]);
				break;
			case 5:
				fdOut.set_vertex(0, m_vertices[4]);
				fdOut.set_vertex(1, m_vertices[5]);
				fdOut.set_vertex(2, m_vertices[6]);
				fdOut.set_vertex(3, m_vertices[7]);
				break;
		}
}

uint Hexahedron::num_faces() const
{
	return 6;
}

EdgeBase* Hexahedron::create_edge(int index)
{
	EdgeDescriptor ed;
	edge(index, ed);
	return new Edge(ed);
}

Face* Hexahedron::create_face(int index)
{
	FaceDescriptor fd;
	face(index, fd);
	return new Quadrilateral(fd.vertex(0), fd.vertex(1), fd.vertex(2), fd.vertex(3));
}

bool Hexahedron::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, VertexBase* newVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//	NOT YET SUPPORTED!
//TODO: implement Hexahedron::collapse_edge
	assert(!"edge-collapse for hexahedrons not yet implemented... sorry");
	return false;
}

bool Hexahedron::refine(std::vector<Volume*>& vNewVolumesOut,
					VertexBase** ppNewVertexOut,
					std::vector<VertexBase*>& vNewEdgeVertices,
					std::vector<VertexBase*>& vNewFaceVertices,
					VertexBase* newVolumeVertex,
					const VertexBase& prototypeVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//TODO: refine for hexahedrons not yet implemented.
	assert(!"refine for hexahedrons not yet implemented... sorry");
	return false;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	PrismDescriptor
PrismDescriptor::PrismDescriptor(const PrismDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
	m_vertex[3] = td.vertex(3);
	m_vertex[4] = td.vertex(4);
	m_vertex[5] = td.vertex(5);
}

PrismDescriptor::PrismDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3,
									VertexBase* v4, VertexBase* v5, VertexBase* v6)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
	m_vertex[4] = v5;
	m_vertex[5] = v6;
}

////////////////////////////////////////////////////////////////////////
//	Prism
Prism::Prism(const PrismDescriptor& td)
{
	m_vertices.resize(6);
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
	m_vertices[4] = td.vertex(4);
	m_vertices[5] = td.vertex(5);
}

Prism::Prism(VertexBase* v1, VertexBase* v2, VertexBase* v3,
						VertexBase* v4, VertexBase* v5, VertexBase* v6)
{
	m_vertices.resize(6);
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
	m_vertices[4] = v5;
	m_vertices[5] = v6;
}

EdgeDescriptor Prism::edge(int index) const
{
	EdgeDescriptor ed;
	edge(index, ed);
	return ed;
}

void Prism::edge(int index, EdgeDescriptor& edOut) const
{
	if(index < 3)//base edges
		edOut.set_vertices(m_vertices[index], m_vertices[(index + 1) % 3]);
	else if(index < 6)//side edges
		edOut.set_vertices(m_vertices[index - 3], m_vertices[index]);
	else//top edges
		edOut.set_vertices(m_vertices[index - 3], m_vertices[(index - 5) % 3 + 3]);
}

uint Prism::num_edges() const
{
	return 9;
}

FaceDescriptor Prism::face(int index) const
{
	FaceDescriptor fd;
	face(index, fd);
	return fd;
}

void Prism::face(int index, FaceDescriptor& fdOut) const
{
	switch(index)
		{
			case 0://bottom
				fdOut.set_num_vertices(3);
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[2]);
				fdOut.set_vertex(2, m_vertices[1]);
				break;
			case 1:
				fdOut.set_num_vertices(4);
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[1]);
				fdOut.set_vertex(2, m_vertices[4]);
				fdOut.set_vertex(3, m_vertices[3]);
				break;
			case 2:
				fdOut.set_num_vertices(4);
				fdOut.set_vertex(0, m_vertices[1]);
				fdOut.set_vertex(1, m_vertices[2]);
				fdOut.set_vertex(2, m_vertices[5]);
				fdOut.set_vertex(3, m_vertices[4]);
				break;
			case 3:
				fdOut.set_num_vertices(4);
				fdOut.set_vertex(0, m_vertices[0]);
				fdOut.set_vertex(1, m_vertices[3]);
				fdOut.set_vertex(2, m_vertices[5]);
				fdOut.set_vertex(3, m_vertices[2]);
				break;
			case 4:
				fdOut.set_num_vertices(3);
				fdOut.set_vertex(0, m_vertices[3]);
				fdOut.set_vertex(1, m_vertices[4]);
				fdOut.set_vertex(2, m_vertices[5]);
				break;
		}
}

uint Prism::num_faces() const
{
	return 5;
}

EdgeBase* Prism::create_edge(int index)
{
	EdgeDescriptor ed;
	edge(index, ed);
	return new Edge(ed);
}

Face* Prism::create_face(int index)
{
	FaceDescriptor fd;
	face(index, fd);
	if(fd.num_vertices() == 3)
		return new Triangle(fd.vertex(0), fd.vertex(1), fd.vertex(2));
	else
		return new Quadrilateral(fd.vertex(0), fd.vertex(1), fd.vertex(2), fd.vertex(3));
}

bool Prism::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, VertexBase* newVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//	NOT YET SUPPORTED!
//TODO: implement prism::collapse_edge
	assert(!"edge-collapse for prism not yet implemented... sorry");
	return false;
}

bool Prism::refine(std::vector<Volume*>& vNewVolumesOut,
					VertexBase** ppNewVertexOut,
					std::vector<VertexBase*>& vNewEdgeVertices,
					std::vector<VertexBase*>& vNewFaceVertices,
					VertexBase* newVolumeVertex,
					const VertexBase& prototypeVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//TODO: refine for prisms not yet implemented.
	assert(!"refine for prisms not yet implemented... sorry");
	return false;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	PyramidDescriptor
PyramidDescriptor::PyramidDescriptor(const PyramidDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
	m_vertex[3] = td.vertex(3);
	m_vertex[4] = td.vertex(4);
}

PyramidDescriptor::PyramidDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3,
									VertexBase* v4, VertexBase* v5)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
	m_vertex[4] = v5;
}

////////////////////////////////////////////////////////////////////////
//	Pyramid
Pyramid::Pyramid(const PyramidDescriptor& td)
{
	m_vertices.resize(5);
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
	m_vertices[4] = td.vertex(4);
}

Pyramid::Pyramid(VertexBase* v1, VertexBase* v2, VertexBase* v3,
				VertexBase* v4, VertexBase* v5)
{
	m_vertices.resize(5);
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
	m_vertices[4] = v5;
}

EdgeDescriptor Pyramid::edge(int index) const
{
	EdgeDescriptor ed;
	edge(index, ed);
	return ed;
}

void Pyramid::edge(int index, EdgeDescriptor& edOut) const
{
	if(index < 4)//base edges
		edOut.set_vertices(m_vertices[index], m_vertices[(index + 1) % 4]);
	else
		edOut.set_vertices(m_vertices[index - 4], m_vertices[4]);
}

uint Pyramid::num_edges() const
{
	return 8;
}

FaceDescriptor Pyramid::face(int index) const
{
	FaceDescriptor fd;
	face(index, fd);
	return fd;
}

void Pyramid::face(int index, FaceDescriptor& fdOut) const
{
	switch(index)
	{
		case 0://bottom
			fdOut.set_num_vertices(4);
			fdOut.set_vertex(0, m_vertices[0]);
			fdOut.set_vertex(1, m_vertices[3]);
			fdOut.set_vertex(2, m_vertices[2]);
			fdOut.set_vertex(3, m_vertices[1]);
			break;
		case 1:
			fdOut.set_num_vertices(3);
			fdOut.set_vertex(0, m_vertices[0]);
			fdOut.set_vertex(1, m_vertices[1]);
			fdOut.set_vertex(2, m_vertices[4]);
			break;
		case 2:
			fdOut.set_num_vertices(3);
			fdOut.set_vertex(0, m_vertices[1]);
			fdOut.set_vertex(1, m_vertices[2]);
			fdOut.set_vertex(2, m_vertices[4]);
			break;
		case 3:
			fdOut.set_num_vertices(3);
			fdOut.set_vertex(0, m_vertices[2]);
			fdOut.set_vertex(1, m_vertices[3]);
			fdOut.set_vertex(2, m_vertices[4]);
			break;
		case 4:
			fdOut.set_num_vertices(3);
			fdOut.set_vertex(0, m_vertices[3]);
			fdOut.set_vertex(1, m_vertices[0]);
			fdOut.set_vertex(2, m_vertices[4]);
			break;
	}
}

uint Pyramid::num_faces() const
{
	return 5;
}

EdgeBase* Pyramid::create_edge(int index)
{
	EdgeDescriptor ed;
	edge(index, ed);
	return new Edge(ed);
}

Face* Pyramid::create_face(int index)
{
	FaceDescriptor fd;
	face(index, fd);
	if(fd.num_vertices() == 3)
		return new Triangle(fd.vertex(0), fd.vertex(1), fd.vertex(2));
	else
		return new Quadrilateral(fd.vertex(0), fd.vertex(1), fd.vertex(2), fd.vertex(3));
}

bool Pyramid::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, VertexBase* newVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//	NOT YET SUPPORTED!
//TODO: implement pyramids::collapse_edge
	assert(!"edge-collapse for pyramids not yet implemented... sorry");
	return false;
}

bool Pyramid::refine(std::vector<Volume*>& vNewVolumesOut,
					VertexBase** ppNewVertexOut,
					std::vector<VertexBase*>& vNewEdgeVertices,
					std::vector<VertexBase*>& vNewFaceVertices,
					VertexBase* newVolumeVertex,
					const VertexBase& prototypeVertex,
					std::vector<VertexBase*>* pvSubstituteVertices)
{
//TODO: refine for pyramids not yet implemented.
	assert(!"refine for pyramids not yet implemented... sorry");
	return false;
}

}//	end of namespace
