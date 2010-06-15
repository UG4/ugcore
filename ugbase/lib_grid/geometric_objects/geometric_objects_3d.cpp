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
							VertexBase** newEdgeVertices,
							VertexBase** newFaceVertices,
							VertexBase* newVolumeVertex,
							const VertexBase& prototypeVertex,
							VertexBase** pSubstituteVertices)
{
//TODO: complete this refine method.
	vNewVolumesOut.clear();
	*ppNewVertexOut = NULL;

//	handle substitute vertices.
	VertexBase** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = &BaseClass::m_vertices.front();

//	check which edges have to be refined and perform the required operations.
	{
	//	get the number of new vertices.
		uint numNewVrts = 0;
		for(uint i = 0; i < 6; ++i)
		{
			if(newEdgeVertices[i] != NULL)
				++numNewVrts;
		}
/*
		uint numNewFaceVrts = 0;
		for(uint i = 0; i < 4; ++i)
		{
			if(vNewFaceVertices[i] != NULL)
				++numNewFaceVrts;
		}

		assert(numNewFaceVrts == 0 && "PROBLEM in Tetrahedron::refine(): new face vertices are not yet supported!");
*/
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
					vNewVolumesOut.push_back(new Tetrahedron(vrts[0], newEdgeVertices[0], newEdgeVertices[2], newEdgeVertices[3]));
					vNewVolumesOut.push_back(new Tetrahedron(vrts[1], newEdgeVertices[1], newEdgeVertices[0], newEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(vrts[2], newEdgeVertices[2], newEdgeVertices[1], newEdgeVertices[5]));
					vNewVolumesOut.push_back(new Tetrahedron(newEdgeVertices[0], newEdgeVertices[1], newEdgeVertices[2], newEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(newEdgeVertices[2], newEdgeVertices[0], newEdgeVertices[4], newEdgeVertices[3]));
					vNewVolumesOut.push_back(new Tetrahedron(newEdgeVertices[3], newEdgeVertices[2], newEdgeVertices[5], newEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(newEdgeVertices[2], newEdgeVertices[1], newEdgeVertices[5], newEdgeVertices[4]));
					vNewVolumesOut.push_back(new Tetrahedron(newEdgeVertices[3], newEdgeVertices[4], newEdgeVertices[5], vrts[3]));
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
	UG_LOG("edge-collapse for hexahedrons not yet implemented... sorry\n");
	return false;
}

bool Hexahedron::refine(std::vector<Volume*>& vNewVolumesOut,
						VertexBase** ppNewVertexOut,
						VertexBase** newEdgeVertices,
						VertexBase** newFaceVertices,
						VertexBase* newVolumeVertex,
						const VertexBase& prototypeVertex,
						VertexBase** pSubstituteVertices)
{
//TODO: complete this refine method.
	vNewVolumesOut.clear();
	*ppNewVertexOut = NULL;

//	handle substitute vertices.
	VertexBase** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = &BaseClass::m_vertices.front();

//	check which edges have to be refined and perform the required operations.
	{
		VertexBase** evrts = newEdgeVertices;
		VertexBase** fvrts = newFaceVertices;
		
	//	get the number of new vertices.
		uint numNewEdgeVrts = 0;
		for(uint i = 0; i < 12; ++i)
		{
			if(newEdgeVertices[i] != NULL)
				++numNewEdgeVrts;
		}

		uint numNewFaceVrts = 0;
		for(uint i = 0; i < 6; ++i)
		{
			if(newFaceVertices[i] != NULL)
				++numNewFaceVrts;
		}

		switch(numNewEdgeVrts)
		{
			case 1:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 1 new edge vertex not yet implemented.");
				return false;
			}

			case 2:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 2 new edge vertices not yet implemented.");
				return false;
			}

			case 3:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 3 new edge vertices not yet implemented.");
				return false;
			}
			case 4:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 4 new edge vertices not yet implemented.");
				return false;
			}
			case 5:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 5 new edge vertices not yet implemented.");
				return false;
			}
			case 6:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 6 new edge vertices not yet implemented.");
				return false;
			}
			case 7:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 7 new edge vertices not yet implemented.");
				return false;
			}
			case 8:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 8 new edge vertices not yet implemented.");
				return false;
			}
			case 9:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 9 new edge vertices not yet implemented.");
				return false;
			}
			case 10:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 10 new edge vertices not yet implemented.");
				return false;
			}
			case 11:
			{
				UG_LOG("PROBLEM in Hexahedron::refine(...): refine with 11 new edge vertices not yet implemented.");
				return false;
			}
			case 12:
			{
				if(numNewFaceVrts != 6){
					assert(!"PROBLEM in Hexahedron::refine(...): All faces have to contain an inner vertex if all edges are marked.");
					return false;
				}
				
				if(!newVolumeVertex)
					newVolumeVertex = reinterpret_cast<VertexBase*>(prototypeVertex.create_empty_instance());
					
				VertexBase* vvrt = newVolumeVertex;
				*ppNewVertexOut = vvrt;
			
				vNewVolumesOut.reserve(8);
			//	left bottom front
				vNewVolumesOut.push_back(new Hexahedron(vrts[0], evrts[0], fvrts[0], evrts[3],
														evrts[4], fvrts[1], vvrt, fvrts[4]));
			//	right bottom front
				vNewVolumesOut.push_back(new Hexahedron(evrts[0], vrts[1], evrts[1], fvrts[0],
														fvrts[1], evrts[5], fvrts[2], vvrt));
			//	right bottom back
				vNewVolumesOut.push_back(new Hexahedron(fvrts[0], evrts[1], vrts[2], evrts[2],
														vvrt, fvrts[2], evrts[6], fvrts[3]));
			//	left bottom back
				vNewVolumesOut.push_back(new Hexahedron(evrts[3], fvrts[0], evrts[2], vrts[3],
														fvrts[4], vvrt, fvrts[3], evrts[7]));
			//	left top front
				vNewVolumesOut.push_back(new Hexahedron(evrts[4], fvrts[1], vvrt, fvrts[4],
														vrts[4], evrts[8], fvrts[5], evrts[11]));
			//	right top front
				vNewVolumesOut.push_back(new Hexahedron(fvrts[1], evrts[5], fvrts[2], vvrt,
														evrts[8], vrts[5], evrts[9], fvrts[5]));
			//	right top back
				vNewVolumesOut.push_back(new Hexahedron(vvrt, fvrts[2], evrts[6], fvrts[3],
														fvrts[5], evrts[9], vrts[6], evrts[10]));
			//	left top back
				vNewVolumesOut.push_back(new Hexahedron(fvrts[4], vvrt, fvrts[3], evrts[7],
														evrts[11], fvrts[5], evrts[10], vrts[7]));
				return true;
			}
		}
	}

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
	UG_LOG("edge-collapse for prism not yet implemented... sorry\n");
	return false;
}

bool Prism::refine(std::vector<Volume*>& vNewVolumesOut,
					VertexBase** ppNewVertexOut,
					VertexBase** newEdgeVertices,
					VertexBase** newFaceVertices,
					VertexBase* newVolumeVertex,
					const VertexBase& prototypeVertex,
					VertexBase** pSubstituteVertices)
{
//TODO: refine for prisms not yet implemented.
	UG_LOG("refine for prisms not yet implemented... sorry\n");
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
	UG_LOG("edge-collapse for pyramids not yet implemented... sorry\n");
	return false;
}

bool Pyramid::refine(std::vector<Volume*>& vNewVolumesOut,
						VertexBase** ppNewVertexOut,
						VertexBase** newEdgeVertices,
						VertexBase** newFaceVertices,
						VertexBase* newVolumeVertex,
						const VertexBase& prototypeVertex,
						VertexBase** pSubstituteVertices)
{
//TODO: refine for pyramids not yet implemented.
	UG_LOG("refine for pyramids not yet implemented... sorry\n");
	return false;
}

}//	end of namespace
