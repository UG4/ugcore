/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <vector>
#include <algorithm>
#include <cassert>
#include "grid_objects.h"
#include "common/common.h"
//ø #include "common/assert.h"

#include "tetrahedron_rules.h"
#include "octahedron_rules.h"
#include "pyramid_rules.h"
#include "prism_rules.h"
#include "hexahedron_rules.h"
#include "grid_object_ids.h"

//#include "../algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{

class TetrahedronClass{
	public:
		enum{
			NUM_VERTICES = tet_rules::NUM_VERTICES,
			NUM_EDGES = tet_rules::NUM_EDGES,
			NUM_FACES = tet_rules::NUM_FACES,
			MAX_NUM_INDS_OUT = tet_rules::MAX_NUM_INDS_OUT
		};
};

class OctahedronClass{
	public:
		enum{
			NUM_VERTICES = oct_rules::NUM_VERTICES,
			NUM_EDGES = oct_rules::NUM_EDGES,
			NUM_FACES = oct_rules::NUM_FACES,
			MAX_NUM_INDS_OUT = oct_rules::MAX_NUM_INDS_OUT
		};
};

class PyramidClass{
	public:
		enum{
			NUM_VERTICES = pyra_rules::NUM_VERTICES,
			NUM_EDGES = pyra_rules::NUM_EDGES,
			NUM_FACES = pyra_rules::NUM_FACES,
			MAX_NUM_INDS_OUT = pyra_rules::MAX_NUM_INDS_OUT
		};
};

class PrismClass{
	public:
		enum{
			NUM_VERTICES = prism_rules::NUM_VERTICES,
			NUM_EDGES = prism_rules::NUM_EDGES,
			NUM_FACES = prism_rules::NUM_FACES,
			MAX_NUM_INDS_OUT = prism_rules::MAX_NUM_INDS_OUT
		};
};

class HexahedronClass{
	public:
		enum{
			NUM_VERTICES = hex_rules::NUM_VERTICES,
			NUM_EDGES = hex_rules::NUM_EDGES,
			NUM_FACES = hex_rules::NUM_FACES,
			MAX_NUM_INDS_OUT = hex_rules::MAX_NUM_INDS_OUT
		};
};

struct GridObjectInfo{
	public:
		static size_t num_vertices(int gridObjectID)	{return inst().m_numVertices[gridObjectID];}

	private:
		static GridObjectInfo& inst(){
			static GridObjectInfo goi;
			return goi;
		}

		GridObjectInfo(){
			for(size_t i = 0; i <GridObjectID:: GOID_NUM_GRID_OBJECT_IDS; ++i)
				m_numVertices[i] = 0;

			m_numVertices[GridObjectID::GOID_TETRAHEDRON] = 4;
			m_numVertices[GridObjectID::GOID_PYRAMID] = 5;
			m_numVertices[GridObjectID::GOID_PRISM] = 6;
			m_numVertices[GridObjectID::GOID_OCTAHEDRON] = 6;
			m_numVertices[GridObjectID::GOID_HEXAHEDRON] = 8;
		}

		size_t	m_numVertices[GridObjectID::GOID_NUM_GRID_OBJECT_IDS];
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	TOOLS
static
void CreateVolumesFromElementIndexList (
		vector<Volume*>& volsOut,
		int* elemIndexList,
		int elemIndexListSize,
		Vertex** vrts)
{
	VolumeDescriptor vd;
	volsOut.clear();

	for(int i = 0; i < elemIndexListSize;){
		int gridObjectID = elemIndexList[i++];
		size_t num = GridObjectInfo::num_vertices(gridObjectID);
		vd.set_num_vertices(num);
		for(size_t j = 0; j < num; ++j){
			assert(vrts[elemIndexList[i]]);
			vd.set_vertex(j, vrts[elemIndexList[i++]]);
		}

		switch(gridObjectID){
			case GridObjectID::GOID_TETRAHEDRON: volsOut.push_back(new Tetrahedron(TetrahedronDescriptor(vd)));	break;
			case GridObjectID::GOID_PYRAMID: volsOut.push_back(new Pyramid(PyramidDescriptor(vd)));		break;
			case GridObjectID::GOID_PRISM: volsOut.push_back(new Prism(PrismDescriptor(vd))); 		break;
			case GridObjectID::GOID_HEXAHEDRON: volsOut.push_back(new Hexahedron(HexahedronDescriptor(vd)));	break;
			case GridObjectID::GOID_OCTAHEDRON: volsOut.push_back(new Octahedron(OctahedronDescriptor(vd)));	break;
		}
	}
}

/**	This refinement helper is called by the different refine implementations.
 * The last parameter is the actual refinement procedure as defined in
 * ug::tet_rules, ug::pyra_rules, ug::hex_rules or ug::prism_rules.
 */
template <typename TElemClass>
static bool Refine(std::vector<Volume*>& vNewVolumesOut,
					Vertex** ppNewVertexOut,
					Vertex** newEdgeVertices,
					Vertex** newFaceVertices,
					Vertex* newVolumeVertex,
					const Vertex& prototypeVertex,
					Vertex** vrts,
					int (*funcRefine)(int*, int*, bool&, vector3*, bool*),
					vector3* corners = nullptr,
					bool* isSnapPoint = nullptr)
{
	vNewVolumesOut.clear();
	*ppNewVertexOut = nullptr;

//	allVrts is an array holding both, the vertices and the new edge-vertices.
//	we will index it with the results later on to create the new elements.
	const int allVrtsSize = TElemClass::NUM_VERTICES + TElemClass::NUM_EDGES
							+ TElemClass::NUM_FACES + 1;
	Vertex* allVrts[allVrtsSize];
	for(int i = 0; i < TElemClass::NUM_VERTICES; ++i)
		allVrts[i] = vrts[i];

//	check which edge has to be refined, and which not
	int newEdgeVrts[TElemClass::NUM_EDGES];
	for(int i = 0; i < TElemClass::NUM_EDGES; ++i){
		allVrts[TElemClass::NUM_VERTICES + i] = newEdgeVertices[i];
		if(newEdgeVertices[i])
			newEdgeVrts[i] = 1;
		else
			newEdgeVrts[i] = 0;
	}

//	copy new face vertices to the allVrts array
	if(newFaceVertices){
		for(int i = 0; i < TElemClass::NUM_FACES; ++i){
			allVrts[TElemClass::NUM_VERTICES + TElemClass::NUM_EDGES + i] =
					newFaceVertices[i];
		}
	}

//	in this array we'll receive the new indices
	int newElemInds[TElemClass::MAX_NUM_INDS_OUT];

//	perform refine
	bool centerVrtRequired = false;
	int numElemInds = funcRefine(newElemInds, newEdgeVrts, centerVrtRequired, corners, isSnapPoint);

	assert(numElemInds != 0 && "PROBLEM in Refine(...): "
								"refine with 1 new edge vertex failed.");

	if(numElemInds == 0)
		return false;

//	if a new center vertex is required, then we'll create one now.
	if(centerVrtRequired){
		if(!newVolumeVertex)
			newVolumeVertex = static_cast<Vertex*>(
								prototypeVertex.create_empty_instance());
		*ppNewVertexOut = newVolumeVertex;
		allVrts[allVrtsSize - 1] = *ppNewVertexOut;
	}

//	debug log of inds
/*
	UG_LOG("newElemInds:");
	for(int i = 0; i < numElemInds; ++i){
		UG_LOG(" " << newElemInds[i]);
	}
	UG_LOG(endl);

	UG_LOG("allVrts, (allVrtsSize: " << allVrtsSize << ") ; ");
	for(int i = 0; i < allVrtsSize; ++i){
		UG_LOG(" " << allVrts[i]);
	}
	UG_LOG(endl);
*/

	CreateVolumesFromElementIndexList(vNewVolumesOut, newElemInds, numElemInds, allVrts);
	// for(int i = 0; i < numElemInds;){
	// 	int gridObjectID = newElemInds[i++];
	// 	size_t num = GridObjectInfo::num_vertices(gridObjectID);
	// 	vd.set_num_vertices(num);
	// 	for(size_t j = 0; j < num; ++j){
	// 		assert(allVrts[newElemInds[i]]);
	// 		vd.set_vertex(j, allVrts[newElemInds[i++]]);
	// 	}

	// 	switch(gridObjectID){
	// 		case GOID_TETRAHEDRON:	vNewVolumesOut.push_back(new Tetrahedron(vd));	break;
	// 		case GOID_PYRAMID:		vNewVolumesOut.push_back(new Pyramid(vd));		break;
	// 		case GOID_PRISM:		vNewVolumesOut.push_back(new Prism(vd)); 		break;
	// 		case GOID_HEXAHEDRON:	vNewVolumesOut.push_back(new Hexahedron(vd));	break;
	// 		case GOID_OCTAHEDRON:	vNewVolumesOut.push_back(new Octahedron(vd));	break;
	// 	}
	// }

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

TetrahedronDescriptor::TetrahedronDescriptor(const VolumeVertices& vv)
{
	assert((vv.num_vertices() == 4) &&	"Bad number of vertices in volume-descriptor. Should be 4.");
	m_vertex[0] = vv.vertex(0);
	m_vertex[1] = vv.vertex(1);
	m_vertex[2] = vv.vertex(2);
	m_vertex[3] = vv.vertex(3);
}

TetrahedronDescriptor::TetrahedronDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
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
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
}

Tetrahedron::Tetrahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
}

EdgeDescriptor Tetrahedron::edge_desc(int index) const
{
	EdgeDescriptor ed;
	edge_desc(index, ed);
	return ed;
}

void Tetrahedron::edge_desc(int index, EdgeDescriptor& edOut) const
{
	using namespace tet_rules;
	assert(index >= 0 && index <= NUM_EDGES);
	edOut.set_vertices(m_vertices[EDGE_VRT_INDS[index][0]],
					   m_vertices[EDGE_VRT_INDS[index][1]]);
}

uint Tetrahedron::num_edges() const
{
	return 6;
}

FaceDescriptor Tetrahedron::face_desc(int index) const
{
	FaceDescriptor fd;
	face_desc(index, fd);
	return fd;
}

void Tetrahedron::face_desc(int index, FaceDescriptor& fdOut) const
{
	using namespace tet_rules;
	assert(index >= 0 && index < NUM_FACES);

	fdOut.set_num_vertices(3);
	fdOut.set_vertex(0, m_vertices[FACE_VRT_INDS[index][0]]);
	fdOut.set_vertex(1, m_vertices[FACE_VRT_INDS[index][2]]);
	fdOut.set_vertex(2, m_vertices[FACE_VRT_INDS[index][1]]);
}

uint Tetrahedron::num_faces() const
{
	return 4;
}

Edge* Tetrahedron::create_edge(int index)
{
	using namespace tet_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	return new RegularEdge(m_vertices[e[0]], m_vertices[e[1]]);
}

Face* Tetrahedron::create_face(int index)
{
	using namespace tet_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
	return new Triangle(m_vertices[f[0]], m_vertices[f[2]], m_vertices[f[1]]);
}

void Tetrahedron::
get_vertex_indices_of_edge (size_t& ind1Out,
							size_t& ind2Out,
							size_t edgeInd) const
{
	assert(edgeInd >= 0 && edgeInd < 6);
	ind1Out = tet_rules::EDGE_VRT_INDS[edgeInd][0];
	ind2Out = tet_rules::EDGE_VRT_INDS[edgeInd][1];
}
											  
void Tetrahedron::
get_vertex_indices_of_face (std::vector<size_t>& indsOut,
							size_t side) const
{
	assert(side >= 0 && side < 4);
	indsOut.resize(3);
	indsOut[0] = tet_rules::FACE_VRT_INDS[side][0];
	indsOut[1] = tet_rules::FACE_VRT_INDS[side][2];
	indsOut[2] = tet_rules::FACE_VRT_INDS[side][1];
}

int Tetrahedron::
get_edge_index_from_vertices(	const size_t vi0,
								const size_t vi1) const
{
	return tet_rules::EDGE_FROM_VRTS[vi0][vi1];
}

int Tetrahedron::
get_face_edge_index(const size_t faceInd,
					const size_t faceEdgeInd) const
{
	return tet_rules::FACE_EDGE_INDS[faceInd][2 - faceEdgeInd];
}

bool Tetrahedron::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, Vertex* newVertex,
					std::vector<Vertex*>* pvSubstituteVertices)
{
//	if an edge of a tetrahedron is collapsed, nothing remains.
	vNewVolumesOut.clear();
	return true;
}

std::pair<GridBaseObjectId, int> Tetrahedron::
get_opposing_object(Vertex* vrt) const
{
	using namespace tet_rules;
	for(int i = 0; i < tet_rules::NUM_VERTICES; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(static_cast<GridBaseObjectId>(OPPOSED_OBJECT[i][0]),
							 OPPOSED_OBJECT[i][1]);
		}
	}

	UG_THROW("Specified vertex is not part of this element.");
}

bool Tetrahedron::refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices,
							vector3* corners,
							bool* isSnapPoint)
{
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

	return Refine<TetrahedronClass>(vNewVolumesOut, ppNewVertexOut,
									newEdgeVertices, newFaceVertices,
									newVolumeVertex, prototypeVertex,
									vrts, tet_rules::Refine, corners,
									isSnapPoint);
}


bool Tetrahedron::is_regular_ref_rule(int edgeMarks) const
{
	return tet_rules::IsRegularRefRule(edgeMarks);
}


void Tetrahedron::get_flipped_orientation(VolumeDescriptor& vdOut)  const
{
//	in order to flip a tetrahedrons orientation, we have to invert the order
//	of the base-vertices
	vdOut.set_num_vertices(4);
	vdOut.set_vertex(0, vertex(2));
	vdOut.set_vertex(1, vertex(1));
	vdOut.set_vertex(2, vertex(0));
	vdOut.set_vertex(3, vertex(3));
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

HexahedronDescriptor::HexahedronDescriptor(const VolumeVertices& vv)
{
	assert((vv.num_vertices() == 8) &&	"Bad number of vertices in volume-descriptor. Should be 8.");
	m_vertex[0] = vv.vertex(0);
	m_vertex[1] = vv.vertex(1);
	m_vertex[2] = vv.vertex(2);
	m_vertex[3] = vv.vertex(3);
	m_vertex[4] = vv.vertex(4);
	m_vertex[5] = vv.vertex(5);
	m_vertex[6] = vv.vertex(6);
	m_vertex[7] = vv.vertex(7);
}

HexahedronDescriptor::HexahedronDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
											Vertex* v5, Vertex* v6, Vertex* v7, Vertex* v8)
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
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
	m_vertices[4] = td.vertex(4);
	m_vertices[5] = td.vertex(5);
	m_vertices[6] = td.vertex(6);
	m_vertices[7] = td.vertex(7);
}

Hexahedron::Hexahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
						Vertex* v5, Vertex* v6, Vertex* v7, Vertex* v8)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
	m_vertices[4] = v5;
	m_vertices[5] = v6;
	m_vertices[6] = v7;
	m_vertices[7] = v8;
}

EdgeDescriptor Hexahedron::edge_desc(int index) const
{
	EdgeDescriptor ed;
	edge_desc(index, ed);
	return ed;
}

void Hexahedron::edge_desc(int index, EdgeDescriptor& edOut) const
{
	using namespace hex_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	edOut.set_vertices(m_vertices[e[0]],
					   m_vertices[e[1]]);
}

uint Hexahedron::num_edges() const
{
	return 12;
}

FaceDescriptor Hexahedron::face_desc(int index) const
{
	FaceDescriptor fd;
	face_desc(index, fd);
	return fd;
}

void Hexahedron::face_desc(int index, FaceDescriptor& fdOut) const
{
	using namespace hex_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];

	fdOut.set_num_vertices(4);
	fdOut.set_vertex(0, m_vertices[f[0]]);
	fdOut.set_vertex(1, m_vertices[f[3]]);
	fdOut.set_vertex(2, m_vertices[f[2]]);
	fdOut.set_vertex(3, m_vertices[f[1]]);
}

uint Hexahedron::num_faces() const
{
	return 6;
}

Edge* Hexahedron::create_edge(int index)
{
	using namespace hex_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	return new RegularEdge(m_vertices[e[0]], m_vertices[e[1]]);
}

Face* Hexahedron::create_face(int index)
{
	using namespace hex_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
	return new Quadrilateral(m_vertices[f[0]], m_vertices[f[3]],
							 m_vertices[f[2]], m_vertices[f[1]]);
}

void Hexahedron::
get_vertex_indices_of_edge (size_t& ind1Out,
							size_t& ind2Out,
							size_t edgeInd) const
{
	assert(edgeInd >= 0 && edgeInd < hex_rules::NUM_EDGES);
	ind1Out = hex_rules::EDGE_VRT_INDS[edgeInd][0];
	ind2Out = hex_rules::EDGE_VRT_INDS[edgeInd][1];
}
											  
void Hexahedron::
get_vertex_indices_of_face (std::vector<size_t>& indsOut,
							size_t side) const
{
	assert(side >= 0 && side < hex_rules::NUM_FACES);
	indsOut.resize(4);
	indsOut[0] = hex_rules::FACE_VRT_INDS[side][0];
	indsOut[1] = hex_rules::FACE_VRT_INDS[side][3];
	indsOut[2] = hex_rules::FACE_VRT_INDS[side][2];
	indsOut[3] = hex_rules::FACE_VRT_INDS[side][1];
}

int Hexahedron::
get_edge_index_from_vertices(	const size_t vi0,
								const size_t vi1) const
{
	return hex_rules::EDGE_FROM_VRTS[vi0][vi1];
}

int Hexahedron::
get_face_edge_index(const size_t faceInd,
					const size_t faceEdgeInd) const
{
	return hex_rules::FACE_EDGE_INDS[faceInd][3 - faceEdgeInd];
}

bool Hexahedron::get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const
{
	using namespace hex_rules;
	int localInd = get_local_side_index(f);
	if(localInd == -1)
		return false;

	face_desc(OPPOSED_FACE[localInd], fdOut);
	return true;
}

std::pair<GridBaseObjectId, int> Hexahedron::
get_opposing_object(Vertex* vrt) const
{
	using namespace hex_rules;
	for(int i = 0; i < hex_rules::NUM_VERTICES; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(static_cast<GridBaseObjectId>(OPPOSED_OBJECT[i][0]),
							 OPPOSED_OBJECT[i][1]);
		}
	}

	UG_THROW("Specified vertex is not part of this element.");
}

bool Hexahedron::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, Vertex* newVertex,
					std::vector<Vertex*>* pvSubstituteVertices)
{
//	NOT YET SUPPORTED!
//TODO: implement Hexahedron::collapse_edge
	vNewVolumesOut.clear();
	UG_LOG("edge-collapse for hexahedrons not yet implemented... sorry\n");
	return false;
}

bool Hexahedron::refine(std::vector<Volume*>& vNewVolumesOut,
						Vertex** ppNewVertexOut,
						Vertex** newEdgeVertices,
						Vertex** newFaceVertices,
						Vertex* newVolumeVertex,
						const Vertex& prototypeVertex,
						Vertex** pSubstituteVertices,
						vector3* corners,
						bool* isSnapPoint)
{
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

	return Refine<HexahedronClass>(vNewVolumesOut, ppNewVertexOut,
								   newEdgeVertices, newFaceVertices,
								   newVolumeVertex, prototypeVertex,
								   vrts, hex_rules::Refine, corners,
								   isSnapPoint);
}

bool Hexahedron::is_regular_ref_rule(int edgeMarks) const
{
	return hex_rules::IsRegularRefRule(edgeMarks);
}

void Hexahedron::get_flipped_orientation(VolumeDescriptor& vdOut)  const
{
//	in order to flip a hexahedrons orientation, we have to swap
//	the bottom and top vertices
	vdOut.set_num_vertices(8);
	vdOut.set_vertex(0, vertex(4));
	vdOut.set_vertex(1, vertex(5));
	vdOut.set_vertex(2, vertex(6));
	vdOut.set_vertex(3, vertex(7));
	vdOut.set_vertex(4, vertex(0));
	vdOut.set_vertex(5, vertex(1));
	vdOut.set_vertex(6, vertex(2));
	vdOut.set_vertex(7, vertex(3));
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

PrismDescriptor::PrismDescriptor(const VolumeVertices& vv)
{
	assert((vv.num_vertices() == 6) &&	"Bad number of vertices in volume-descriptor. Should be 6.");
	m_vertex[0] = vv.vertex(0);
	m_vertex[1] = vv.vertex(1);
	m_vertex[2] = vv.vertex(2);
	m_vertex[3] = vv.vertex(3);
	m_vertex[4] = vv.vertex(4);
	m_vertex[5] = vv.vertex(5);
}

PrismDescriptor::PrismDescriptor(Vertex* v1, Vertex* v2, Vertex* v3,
									Vertex* v4, Vertex* v5, Vertex* v6)
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
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
	m_vertices[4] = td.vertex(4);
	m_vertices[5] = td.vertex(5);
}

Prism::Prism(Vertex* v1, Vertex* v2, Vertex* v3,
						Vertex* v4, Vertex* v5, Vertex* v6)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
	m_vertices[4] = v5;
	m_vertices[5] = v6;
}

EdgeDescriptor Prism::edge_desc(int index) const
{
	EdgeDescriptor ed;
	edge_desc(index, ed);
	return ed;
}

void Prism::edge_desc(int index, EdgeDescriptor& edOut) const
{
	using namespace prism_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	edOut.set_vertices(m_vertices[e[0]],
					   m_vertices[e[1]]);
}

uint Prism::num_edges() const
{
	return 9;
}

FaceDescriptor Prism::face_desc(int index) const
{
	FaceDescriptor fd;
	face_desc(index, fd);
	return fd;
}

void Prism::face_desc(int index, FaceDescriptor& fdOut) const
{
	using namespace prism_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
	if(f[3] == -1){
		fdOut.set_num_vertices(3);
		fdOut.set_vertex(0, m_vertices[f[0]]);
		fdOut.set_vertex(1, m_vertices[f[2]]);
		fdOut.set_vertex(2, m_vertices[f[1]]);
	}
	else{
		fdOut.set_num_vertices(4);
		fdOut.set_vertex(0, m_vertices[f[0]]);
		fdOut.set_vertex(1, m_vertices[f[3]]);
		fdOut.set_vertex(2, m_vertices[f[2]]);
		fdOut.set_vertex(3, m_vertices[f[1]]);
	}
}

uint Prism::num_faces() const
{
	return 5;
}

Edge* Prism::create_edge(int index)
{
	using namespace prism_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	return new RegularEdge(m_vertices[e[0]], m_vertices[e[1]]);
}

Face* Prism::create_face(int index)
{
	using namespace prism_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
	if(f[3] == -1){
		return new Triangle(m_vertices[f[0]], m_vertices[f[2]],
							m_vertices[f[1]]);
	}
	else{
		return new Quadrilateral(m_vertices[f[0]], m_vertices[f[3]],
								 m_vertices[f[2]], m_vertices[f[1]]);
	}
}

void Prism::
get_vertex_indices_of_edge (size_t& ind1Out,
							size_t& ind2Out,
							size_t edgeInd) const
{
	assert(edgeInd >= 0 && edgeInd < prism_rules::NUM_EDGES);
	ind1Out = prism_rules::EDGE_VRT_INDS[edgeInd][0];
	ind2Out = prism_rules::EDGE_VRT_INDS[edgeInd][1];
}
											  
void Prism::
get_vertex_indices_of_face (std::vector<size_t>& indsOut,
							size_t side) const
{
	assert(side >= 0 && side < prism_rules::NUM_FACES);

	if(prism_rules::FACE_VRT_INDS[side][3] == -1){
		indsOut.resize(3);
		indsOut[0] = prism_rules::FACE_VRT_INDS[side][0];
		indsOut[1] = prism_rules::FACE_VRT_INDS[side][2];
		indsOut[2] = prism_rules::FACE_VRT_INDS[side][1];
	}
	else{
		indsOut.resize(4);
		indsOut[0] = prism_rules::FACE_VRT_INDS[side][0];
		indsOut[1] = prism_rules::FACE_VRT_INDS[side][3];
		indsOut[2] = prism_rules::FACE_VRT_INDS[side][2];
		indsOut[3] = prism_rules::FACE_VRT_INDS[side][1];
	}
}

int Prism::
get_edge_index_from_vertices(	const size_t vi0,
								const size_t vi1) const
{
	return prism_rules::EDGE_FROM_VRTS[vi0][vi1];
}

int Prism::
get_face_edge_index(const size_t faceInd,
					const size_t faceEdgeInd) const
{
	if(prism_rules::FACE_EDGE_INDS[faceInd][3] == -1)
		return prism_rules::FACE_EDGE_INDS[faceInd][2 - faceEdgeInd];
	else
		return prism_rules::FACE_EDGE_INDS[faceInd][3 - faceEdgeInd];
}

bool Prism::get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const
{
	using namespace prism_rules;
	int localInd = get_local_side_index(f);
	if(localInd == -1)
		return false;

	int opposedInd = OPPOSED_FACE[localInd];
	if(opposedInd == -1)
		return false;

	face_desc(opposedInd, fdOut);
	return true;
}

std::pair<GridBaseObjectId, int> Prism::
get_opposing_object(Vertex* vrt) const
{
	using namespace prism_rules;
	for(int i = 0; i < prism_rules::NUM_VERTICES; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(static_cast<GridBaseObjectId>(OPPOSED_OBJECT[i][0]),
							 OPPOSED_OBJECT[i][1]);
		}
	}
	UG_THROW("Specified vertex is not part of this element.");
}

bool Prism::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, Vertex* newVertex,
					std::vector<Vertex*>* pvSubstituteVertices)
{
	using namespace prism_rules;

	int elemInds[MAX_NUM_COLLAPSE_INDS_OUT];
	int elemIndsSize = CollapseEdge(
							elemInds,
							EDGE_VRT_INDS[edgeIndex][0],
							EDGE_VRT_INDS[edgeIndex][1]);

	if(elemIndsSize > 0){
		Vertex** vrts;
		if(pvSubstituteVertices)
			vrts = &pvSubstituteVertices->front();
		else
			vrts = m_vertices;

		CreateVolumesFromElementIndexList(
			vNewVolumesOut,
			elemInds,
			elemIndsSize,
			vrts);
		return !vNewVolumesOut.empty();
	}
	else{
		vNewVolumesOut.clear();
		return false;
	}
}

bool Prism::refine(std::vector<Volume*>& vNewVolumesOut,
					Vertex** ppNewVertexOut,
					Vertex** newEdgeVertices,
					Vertex** newFaceVertices,
					Vertex* newVolumeVertex,
					const Vertex& prototypeVertex,
					Vertex** pSubstituteVertices,
					vector3* corners,
					bool* isSnapPoint)
{
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

	return Refine<PrismClass>(vNewVolumesOut, ppNewVertexOut,
							  newEdgeVertices, newFaceVertices,
							  newVolumeVertex, prototypeVertex,
							  vrts, prism_rules::Refine, corners,
							  isSnapPoint);
}

bool Prism::is_regular_ref_rule(int edgeMarks) const
{
	return prism_rules::IsRegularRefRule(edgeMarks);
}

void Prism::get_flipped_orientation(VolumeDescriptor& vdOut) const
{
//	in order to flip a prisms orientation, we have to swap
//	the bottom and top vertices
	vdOut.set_num_vertices(6);
	vdOut.set_vertex(0, vertex(3));
	vdOut.set_vertex(1, vertex(4));
	vdOut.set_vertex(2, vertex(5));
	vdOut.set_vertex(3, vertex(0));
	vdOut.set_vertex(4, vertex(1));
	vdOut.set_vertex(5, vertex(2));
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

PyramidDescriptor::PyramidDescriptor(const VolumeVertices& vv)
{
	assert((vv.num_vertices() == 5) &&	"Bad number of vertices in volume-descriptor. Should be 5.");
	m_vertex[0] = vv.vertex(0);
	m_vertex[1] = vv.vertex(1);
	m_vertex[2] = vv.vertex(2);
	m_vertex[3] = vv.vertex(3);
	m_vertex[4] = vv.vertex(4);
}

PyramidDescriptor::PyramidDescriptor(Vertex* v1, Vertex* v2, Vertex* v3,
									Vertex* v4, Vertex* v5)
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
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
	m_vertices[4] = td.vertex(4);
}

Pyramid::Pyramid(Vertex* v1, Vertex* v2, Vertex* v3,
				Vertex* v4, Vertex* v5)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
	m_vertices[4] = v5;
}

EdgeDescriptor Pyramid::edge_desc(int index) const
{
	EdgeDescriptor ed;
	edge_desc(index, ed);
	return ed;
}

void Pyramid::edge_desc(int index, EdgeDescriptor& edOut) const
{
	using namespace pyra_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	edOut.set_vertices(m_vertices[e[0]],
					   m_vertices[e[1]]);
}

uint Pyramid::num_edges() const
{
	return 8;
}

FaceDescriptor Pyramid::face_desc(int index) const
{
	FaceDescriptor fd;
	face_desc(index, fd);
	return fd;
}

void Pyramid::face_desc(int index, FaceDescriptor& fdOut) const
{
	using namespace pyra_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
	if(f[3] == -1){
		fdOut.set_num_vertices(3);
		fdOut.set_vertex(0, m_vertices[f[0]]);
		fdOut.set_vertex(1, m_vertices[f[2]]);
		fdOut.set_vertex(2, m_vertices[f[1]]);
	}
	else{
		fdOut.set_num_vertices(4);
		fdOut.set_vertex(0, m_vertices[f[0]]);
		fdOut.set_vertex(1, m_vertices[f[3]]);
		fdOut.set_vertex(2, m_vertices[f[2]]);
		fdOut.set_vertex(3, m_vertices[f[1]]);
	}
}

uint Pyramid::num_faces() const
{
	return 5;
}

Edge* Pyramid::create_edge(int index)
{
	using namespace pyra_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	return new RegularEdge(m_vertices[e[0]], m_vertices[e[1]]);
}

Face* Pyramid::create_face(int index)
{
	using namespace pyra_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
	if(f[3] == -1){
		return new Triangle(m_vertices[f[0]], m_vertices[f[2]],
							m_vertices[f[1]]);
	}
	else{
		return new Quadrilateral(m_vertices[f[0]], m_vertices[f[3]],
								 m_vertices[f[2]], m_vertices[f[1]]);
	}
}

void Pyramid::
get_vertex_indices_of_edge (size_t& ind1Out,
							size_t& ind2Out,
							size_t edgeInd) const
{
	assert(edgeInd >= 0 && edgeInd < pyra_rules::NUM_EDGES);
	ind1Out = pyra_rules::EDGE_VRT_INDS[edgeInd][0];
	ind2Out = pyra_rules::EDGE_VRT_INDS[edgeInd][1];
}
											  
void Pyramid::
get_vertex_indices_of_face (std::vector<size_t>& indsOut,
							size_t side) const
{
	assert(side >= 0 && side < pyra_rules::NUM_FACES);

	if(pyra_rules::FACE_VRT_INDS[side][3] == -1){
		indsOut.resize(3);
		indsOut[0] = pyra_rules::FACE_VRT_INDS[side][0];
		indsOut[1] = pyra_rules::FACE_VRT_INDS[side][2];
		indsOut[2] = pyra_rules::FACE_VRT_INDS[side][1];
	}
	else{
		indsOut.resize(4);
		indsOut[0] = pyra_rules::FACE_VRT_INDS[side][0];
		indsOut[1] = pyra_rules::FACE_VRT_INDS[side][3];
		indsOut[2] = pyra_rules::FACE_VRT_INDS[side][2];
		indsOut[3] = pyra_rules::FACE_VRT_INDS[side][1];
	}
}

int Pyramid::
get_edge_index_from_vertices(	const size_t vi0,
								const size_t vi1) const
{
	return pyra_rules::EDGE_FROM_VRTS[vi0][vi1];
}

int Pyramid::
get_face_edge_index(const size_t faceInd,
					const size_t faceEdgeInd) const
{
	if(pyra_rules::FACE_EDGE_INDS[faceInd][3] == -1)
		return pyra_rules::FACE_EDGE_INDS[faceInd][2 - faceEdgeInd];
	else
		return pyra_rules::FACE_EDGE_INDS[faceInd][3 - faceEdgeInd];
}

std::pair<GridBaseObjectId, int> Pyramid::
get_opposing_object(Vertex* vrt) const
{
	using namespace pyra_rules;
	for(int i = 0; i < pyra_rules::NUM_VERTICES; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(static_cast<GridBaseObjectId>(OPPOSED_OBJECT[i][0]),
							 OPPOSED_OBJECT[i][1]);
		}
	}
	UG_THROW("Specified vertex is not part of this element.");
}

bool Pyramid::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
					int edgeIndex, Vertex* newVertex,
					std::vector<Vertex*>* pvSubstituteVertices)
{
	using namespace pyra_rules;

	int elemInds[MAX_NUM_COLLAPSE_INDS_OUT];
	int elemIndsSize = CollapseEdge(
							elemInds,
							EDGE_VRT_INDS[edgeIndex][0],
							EDGE_VRT_INDS[edgeIndex][1]);

	if(elemIndsSize > 0){
		Vertex** vrts;
		if(pvSubstituteVertices)
			vrts = &pvSubstituteVertices->front();
		else
			vrts = m_vertices;

		CreateVolumesFromElementIndexList(
			vNewVolumesOut,
			elemInds,
			elemIndsSize,
			vrts);
		return !vNewVolumesOut.empty();
	}
	else{
		vNewVolumesOut.clear();
		return false;
	}
}

bool Pyramid::refine(std::vector<Volume*>& vNewVolumesOut,
						Vertex** ppNewVertexOut,
						Vertex** newEdgeVertices,
						Vertex** newFaceVertices,
						Vertex* newVolumeVertex,
						const Vertex& prototypeVertex,
						Vertex** pSubstituteVertices,
						vector3* corners,
						bool* isSnapPoint)
{
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

	return Refine<PyramidClass>(vNewVolumesOut, ppNewVertexOut,
									newEdgeVertices, newFaceVertices,
									newVolumeVertex, prototypeVertex,
									vrts, pyra_rules::Refine, corners,
									isSnapPoint);
}

bool Pyramid::is_regular_ref_rule(int edgeMarks) const
{
	return pyra_rules::IsRegularRefRule(edgeMarks);
}

void Pyramid::get_flipped_orientation(VolumeDescriptor& vdOut) const
{
//	in order to flip a pyramids orientation, we have to invert the order
//	of the base-vertices
	vdOut.set_num_vertices(5);
	vdOut.set_vertex(0, vertex(3));
	vdOut.set_vertex(1, vertex(2));
	vdOut.set_vertex(2, vertex(1));
	vdOut.set_vertex(3, vertex(0));
	vdOut.set_vertex(4, vertex(4));
}


////////////////////////////////////////////////////////////////////////
//	OctahedronDescriptor
OctahedronDescriptor::OctahedronDescriptor(const OctahedronDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
	m_vertex[3] = td.vertex(3);
   	m_vertex[4] = td.vertex(4);
  	m_vertex[5] = td.vertex(5);
}

OctahedronDescriptor::OctahedronDescriptor(const VolumeVertices& vv)
{
	assert((vv.num_vertices() == 6) &&	"Bad number of vertices in volume-descriptor. Should be 6.");
	m_vertex[0] = vv.vertex(0);
	m_vertex[1] = vv.vertex(1);
	m_vertex[2] = vv.vertex(2);
	m_vertex[3] = vv.vertex(3);
   	m_vertex[4] = vv.vertex(4);
   	m_vertex[5] = vv.vertex(5);
}

OctahedronDescriptor::OctahedronDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, Vertex* v5, Vertex* v6)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
   	m_vertex[4] = v5;
   	m_vertex[5] = v6;
}

////////////////////////////////////////////////////////////////////////
//	Octahedron
Octahedron::Octahedron(const OctahedronDescriptor& td)
{
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
	m_vertices[3] = td.vertex(3);
   	m_vertices[4] = td.vertex(4);
   	m_vertices[5] = td.vertex(5);
}

Octahedron::Octahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, Vertex* v5, Vertex* v6)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
    m_vertices[4] = v5;
	m_vertices[5] = v6;
}

EdgeDescriptor Octahedron::edge_desc(int index) const
{
	EdgeDescriptor ed;
	edge_desc(index, ed);
	return ed;
}

void Octahedron::edge_desc(int index, EdgeDescriptor& edOut) const
{
	using namespace oct_rules;
	assert(index >= 0 && index <= NUM_EDGES);
	edOut.set_vertices(m_vertices[EDGE_VRT_INDS[index][0]],
					   m_vertices[EDGE_VRT_INDS[index][1]]);
}

uint Octahedron::num_edges() const
{
	return 12;
}

FaceDescriptor Octahedron::face_desc(int index) const
{
	FaceDescriptor fd;
	face_desc(index, fd);
	return fd;
}

void Octahedron::face_desc(int index, FaceDescriptor& fdOut) const
{
	using namespace oct_rules;
	assert(index >= 0 && index < NUM_FACES);

	fdOut.set_num_vertices(3);
	fdOut.set_vertex(0, m_vertices[FACE_VRT_INDS[index][0]]);
	fdOut.set_vertex(1, m_vertices[FACE_VRT_INDS[index][2]]);
	fdOut.set_vertex(2, m_vertices[FACE_VRT_INDS[index][1]]);
}

uint Octahedron::num_faces() const
{
	return 8;
}

Edge* Octahedron::create_edge(int index)
{
	using namespace oct_rules;
	assert(index >= 0 && index < NUM_EDGES);
	const int* e = EDGE_VRT_INDS[index];
	return new RegularEdge(m_vertices[e[0]], m_vertices[e[1]]);
}

Face* Octahedron::create_face(int index)
{
	using namespace oct_rules;
	assert(index >= 0 && index < NUM_FACES);

	const int* f = FACE_VRT_INDS[index];
    return new Triangle(m_vertices[f[0]], m_vertices[f[2]], m_vertices[f[1]]);
}

void Octahedron::
get_vertex_indices_of_edge (size_t& ind1Out,
							size_t& ind2Out,
							size_t edgeInd) const
{
	assert(edgeInd >= 0 && edgeInd < oct_rules::NUM_EDGES);
	ind1Out = oct_rules::EDGE_VRT_INDS[edgeInd][0];
	ind2Out = oct_rules::EDGE_VRT_INDS[edgeInd][1];
}
											  
void Octahedron::
get_vertex_indices_of_face (std::vector<size_t>& indsOut,
							size_t side) const
{
	assert(side >= 0 && side < oct_rules::NUM_FACES);

	indsOut.resize(3);
	indsOut[0] = oct_rules::FACE_VRT_INDS[side][0];
	indsOut[1] = oct_rules::FACE_VRT_INDS[side][2];
	indsOut[2] = oct_rules::FACE_VRT_INDS[side][1];
}

int Octahedron::
get_edge_index_from_vertices(	const size_t vi0,
								const size_t vi1) const
{
	return oct_rules::EDGE_FROM_VRTS[vi0][vi1];
}

int Octahedron::
get_face_edge_index(const size_t faceInd,
					const size_t faceEdgeInd) const
{
	return oct_rules::FACE_EDGE_INDS[faceInd][2 - faceEdgeInd];
}

std::pair<GridBaseObjectId, int> Octahedron::
get_opposing_object(Vertex* vrt) const
{
	using namespace oct_rules;
	for(int i = 0; i < oct_rules::NUM_VERTICES; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(static_cast<GridBaseObjectId>(OPPOSED_OBJECT[i][0]),
							 OPPOSED_OBJECT[i][1]);
		}
	}

	UG_THROW("Specified vertex is not part of this element.");
}

bool Octahedron::collapse_edge(std::vector<Volume*>& vNewVolumesOut,
                                int edgeIndex, Vertex* newVertex,
                                std::vector<Vertex*>* pvSubstituteVertices)
{
//	NOT YET SUPPORTED!
//TODO: implement octahedron::collapse_edge
	vNewVolumesOut.clear();
	UG_LOG("edge-collapse for octahedrons not yet implemented... sorry\n");
	return false;
}


bool Octahedron::refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices,
							vector3* corners,
							bool* isSnapPoint)
{
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

	return Refine<OctahedronClass>(vNewVolumesOut, ppNewVertexOut,
									newEdgeVertices, newFaceVertices,
									newVolumeVertex, prototypeVertex,
									vrts, oct_rules::Refine, corners,
									isSnapPoint);
}

bool Octahedron::is_regular_ref_rule(int edgeMarks) const
{
	return oct_rules::IsRegularRefRule(edgeMarks);
}

void Octahedron::get_flipped_orientation(VolumeDescriptor& vdOut)  const
{
//	in order to flip a pyramids orientation, we have to invert the order
//	of the base-vertices
	vdOut.set_num_vertices(6);
	vdOut.set_vertex(0, vertex(0));
	vdOut.set_vertex(1, vertex(4));
	vdOut.set_vertex(2, vertex(3));
	vdOut.set_vertex(3, vertex(2));
	vdOut.set_vertex(4, vertex(1));
	vdOut.set_vertex(5, vertex(5));
}

}//	end of namespace
