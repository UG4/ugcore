//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d04

#include <vector>
#include <algorithm>
#include "grid_objects.h"
#include "common/common.h"
#include "tetrahedron_rules.h"
#include "pyramid_rules.h"
#include "prism_rules.h"
#include "hexahedron_rules.h"

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

/**	This refinement helper is called by the different refine implementations.
 * The last parameter is the actual refinement procedure as defined in
 * ug::tet_rules, ug::pyra_rules, ug::hex_rules or ug::prism_rules.
 */
template <class TElemClass>
static bool Refine(std::vector<Volume*>& vNewVolumesOut,
					Vertex** ppNewVertexOut,
					Vertex** newEdgeVertices,
					Vertex** newFaceVertices,
					Vertex* newVolumeVertex,
					const Vertex& prototypeVertex,
					Vertex** vrts,
					int (*funcRefine)(int*, int*, bool&, vector3*),
					vector3* corners = NULL)
{
	vNewVolumesOut.clear();
	*ppNewVertexOut = NULL;

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
	int numElemInds = funcRefine(newElemInds, newEdgeVrts, centerVrtRequired, corners);

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

	UG_LOG("allVrts:");
	for(int i = 0; i < allVrtsSize; ++i){
		UG_LOG(" " << allVrts[i]);
	}
	UG_LOG(endl);
*/

//	the VolumeDescriptor will be used to create new volumes
	VolumeDescriptor vd;

	for(int i = 0; i < numElemInds;){
		int num = newElemInds[i++];
		vd.set_num_vertices(num);
		for(int j = 0; j < num; ++j){
			assert(allVrts[newElemInds[i]]);
			vd.set_vertex(j, allVrts[newElemInds[i++]]);
		}

		switch(num){
			case 4:	vNewVolumesOut.push_back(new Tetrahedron(vd));	break;
			case 5:	vNewVolumesOut.push_back(new Pyramid(vd));		break;
			case 6:	vNewVolumesOut.push_back(new Prism(vd));		break;
			case 8:	vNewVolumesOut.push_back(new Hexahedron(vd));	break;
		}
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

EdgeBase* Tetrahedron::create_edge(int index)
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
get_local_vertex_indices_of_edge(size_t& ind1Out,
								  size_t& ind2Out,
								  size_t edgeInd) const
{
	assert(edgeInd >= 0 && edgeInd < 6);
	ind1Out = tet_rules::EDGE_VRT_INDS[edgeInd][0];
	ind2Out = tet_rules::EDGE_VRT_INDS[edgeInd][1];
}
											  
void Tetrahedron::
get_local_vertex_indices_of_face(std::vector<size_t>& indsOut,
								 size_t side) const
{
	assert(side >= 0 && side < 4);
	indsOut.resize(3);
	indsOut[0] = tet_rules::FACE_VRT_INDS[side][0];
	indsOut[1] = tet_rules::FACE_VRT_INDS[side][2];
	indsOut[2] = tet_rules::FACE_VRT_INDS[side][1];
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
							vector3* corners)
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
									vrts, tet_rules::Refine, corners);
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

EdgeBase* Hexahedron::create_edge(int index)
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
						vector3*)
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
								   vrts, hex_rules::Refine);
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

EdgeBase* Prism::create_edge(int index)
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
//	NOT YET SUPPORTED!
//TODO: implement prism::collapse_edge
	UG_LOG("edge-collapse for prism not yet implemented... sorry\n");
	return false;
}

bool Prism::refine(std::vector<Volume*>& vNewVolumesOut,
					Vertex** ppNewVertexOut,
					Vertex** newEdgeVertices,
					Vertex** newFaceVertices,
					Vertex* newVolumeVertex,
					const Vertex& prototypeVertex,
					Vertex** pSubstituteVertices,
					vector3*)
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
							  vrts, prism_rules::Refine);
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

EdgeBase* Pyramid::create_edge(int index)
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
//	NOT YET SUPPORTED!
//TODO: implement pyramids::collapse_edge
	UG_LOG("edge-collapse for pyramids not yet implemented... sorry\n");
	return false;
}

bool Pyramid::refine(std::vector<Volume*>& vNewVolumesOut,
						Vertex** ppNewVertexOut,
						Vertex** newEdgeVertices,
						Vertex** newFaceVertices,
						Vertex* newVolumeVertex,
						const Vertex& prototypeVertex,
						Vertex** pSubstituteVertices,
						vector3*)
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
									vrts, pyra_rules::Refine);
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

}//	end of namespace
