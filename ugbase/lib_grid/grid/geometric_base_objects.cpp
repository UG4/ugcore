//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d26

#include "geometric_base_objects.h"
#include "grid_util.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	implementation of edge
bool EdgeBase::get_opposing_side(VertexBase* v, VertexBase** vrtOut)
{
	if(v == m_vertices[0])
		*vrtOut = m_vertices[1];
	else if(v == m_vertices[1])
		*vrtOut = m_vertices[0];
	else
		return false;
	return true;
}

////////////////////////////////////////////////////////////////////////
//	implementation of edge-descriptor
EdgeDescriptor::EdgeDescriptor()
{
}

EdgeDescriptor::EdgeDescriptor(const EdgeDescriptor& ed)
{
	EdgeVertices::assign_edge_vertices(ed);
}

EdgeDescriptor::EdgeDescriptor(VertexBase* vrt1, VertexBase* vrt2)
{
	m_vertices[0] = vrt1;
	m_vertices[1] = vrt2;
}

EdgeDescriptor& EdgeDescriptor::operator = (const EdgeDescriptor& ed)
{
	EdgeVertices::assign_edge_vertices(ed);
	return *this;
}


////////////////////////////////////////////////////////////////////////
//	implementation of face
int Face::get_local_side_index(EdgeVertices* e) const
{
	EdgeDescriptor ed;
	for(size_t i = 0; i < num_sides(); ++i){
		edge_desc(i, ed);
		if(CompareVertices(e, &ed))
			return (int)i;
	}
	return -1;
}


////////////////////////////////////////////////////////////////////////
//	implementation of face-descriptor
FaceDescriptor::FaceDescriptor() :
	m_numVertices(0)
{
}

FaceDescriptor::FaceDescriptor(uint numVertices) :
	m_numVertices(numVertices)
{
}

FaceDescriptor::FaceDescriptor(const FaceDescriptor& fd)
{
	m_numVertices = fd.m_numVertices;
	for(uint i = 0; i < m_numVertices; ++i)
		m_vertices[i] = fd.m_vertices[i];
}

FaceDescriptor& FaceDescriptor::operator = (const FaceDescriptor& fd)
{
	m_numVertices = fd.m_numVertices;
	for(uint i = 0; i < m_numVertices; ++i)
		m_vertices[i] = fd.m_vertices[i];

	return *this;
}


////////////////////////////////////////////////////////////////////////
//	implementation of Volume
int Volume::get_local_side_index(FaceVertices* f) const
{
	FaceDescriptor fd;
	for(size_t i = 0; i < num_sides(); ++i){
		face_desc(i, fd);
		if(CompareVertices(f, &fd))
			return (int)i;
	}
	return -1;
}

void Volume::get_flipped_orientation(VolumeDescriptor& vdOut) const
{
	throw(int(0));
	vdOut = *this;
}

////////////////////////////////////////////////////////////////////////
//	implementation of volume-descriptor
VolumeDescriptor::VolumeDescriptor()
{
}

VolumeDescriptor::VolumeDescriptor(uint numVertices)
{
	set_num_vertices(numVertices);
}

VolumeDescriptor::VolumeDescriptor(const VolumeDescriptor& vd)
{
	m_numVertices = vd.m_numVertices;
	for(uint i = 0; i < m_numVertices; ++i)
		m_vertices[i] = vd.m_vertices[i];
}

VolumeDescriptor& VolumeDescriptor::operator = (const VolumeDescriptor& vd)
{
	m_numVertices = vd.m_numVertices;
	for(uint i = 0; i < m_numVertices; ++i)
		m_vertices[i] = vd.m_vertices[i];

	return *this;
}

VolumeDescriptor& VolumeDescriptor::operator = (const VolumeVertices& vv)
{
	m_numVertices = vv.num_vertices();
	for(uint i = 0; i < m_numVertices; ++i)
		m_vertices[i] = vv.vertex(i);

	return *this;
}

////////////////////////////////////////////////////////////////////////
//	inline implementation of hash-keys
//	this methods are used by the template-specializations of hash_key<...>
///	sums the squared hash-values of associated vertices.
static inline unsigned long HashKey(const EdgeVertices* key)
{
	unsigned long a = key->vertex(0)->get_hash_value();
	unsigned long b = key->vertex(1)->get_hash_value();

	//return (unsigned long)(a + b);
	return (unsigned long)(a * a + b * b);
	//if(b > a) return (unsigned long) ((a+b) * (b-a));
	//else return (unsigned long) ((a+b) * (a-b));
}

///	sums the squared hash-values of associated vertices.
static inline unsigned long HashKey(const FaceVertices* key)
{
	unsigned long retVal = 0;
	size_t numVrts = key->num_vertices();
	FaceVertices::ConstVertexArray vrts = key->vertices();

	for(size_t i = 0; i < numVrts; ++i)
	{
		unsigned long a = vrts[i]->get_hash_value();
		retVal += (a*a);
	}
	return retVal;
}

///	sums the squared hash-values of associated vertices.
static inline unsigned long HashKey(const VolumeVertices* key)
{
	unsigned long retVal = 0;
	size_t numVrts = key->num_vertices();
	VolumeVertices::ConstVertexArray vrts = key->vertices();
	for(size_t i = 0; i < numVrts; ++i)
	{
		unsigned long a = vrts[i]->get_hash_value();
		retVal += (a*a);
	}

	return retVal;
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for vertices
//	returns the hash-value of the vertex.
template <>
size_t hash_key<PVertexBase>(const PVertexBase& key)
{
	return (unsigned long)key->get_hash_value();
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for edges
template <>
size_t hash_key<PEdgeVertices>(const PEdgeVertices& key)
{
	return HashKey(key);
}

template <>
size_t hash_key<PEdgeBase>(const PEdgeBase& key)
{
	return HashKey(key);
}

template <>
size_t hash_key<PEdgeDescriptor>(const PEdgeDescriptor& key)
{
	return HashKey(key);
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for faces
template <>
size_t hash_key<PFaceVertices>(const PFaceVertices& key)
{
	return HashKey(key);
}

template <>
size_t hash_key<PFace>(const PFace& key)
{
	return HashKey(key);
}

template <>
size_t hash_key<PFaceDescriptor>(const PFaceDescriptor& key)
{
	return HashKey(key);
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for volumes
template <>
size_t hash_key<PVolumeVertices>(const PVolumeVertices& key)
{
	return HashKey(key);
}

template <>
size_t hash_key<PVolume>(const PVolume& key)
{
	return HashKey(key);
}

template <>
size_t hash_key<PVolumeDescriptor>(const PVolumeDescriptor& key)
{
	return HashKey(key);
}

}//	end of namespace
