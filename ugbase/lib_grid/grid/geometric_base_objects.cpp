//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d26

#include "geometric_base_objects.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	inline implementation of hash-keys
//	this methods are used by the template-specializations of hash_key<...>
///	sums the squared hash-values of associated vertices.
static inline unsigned long HashKey(const EdgeVertices* key)
{
	uint32 a = key->vertex(0)->get_hash_value();
	uint32 b = key->vertex(1)->get_hash_value();

	//return (unsigned long)(a + b);
	return (unsigned long)(a * a + b * b);
	//if(b > a) return (unsigned long) ((a+b) * (b-a));
	//else return (unsigned long) ((a+b) * (a-b));
}

///	sums the squared hash-values of associated vertices.
static inline unsigned long HashKey(const FaceVertices* key)
{
	unsigned long retVal = key->vertex(0)->get_hash_value();
	retVal *= retVal;
	uint numVrts = key->num_vertices();
	for(uint i = 1; i < numVrts; ++i)
	{
		uint32 a = key->vertex(i)->get_hash_value();
		retVal += (a*a);
	}
	return retVal;
}

///	sums the squared hash-values of associated vertices.
static inline unsigned long HashKey(const VolumeVertices* key)
{
	unsigned long retVal = key->vertex(0)->get_hash_value();
	retVal *= retVal;
	uint numVrts = key->num_vertices();
	for(uint i = 1; i < numVrts; ++i)
	{
		uint32 a = key->vertex(i)->get_hash_value();
		retVal += (a*a);
	}

	return retVal;
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for vertices
//	returns the hash-value of the vertex.
template <>
unsigned long hash_key<PVertexBase>(const PVertexBase& key)
{
	return (unsigned long)key->get_hash_value();
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for edges
template <>
unsigned long hash_key<PEdgeVertices>(const PEdgeVertices& key)
{
	return HashKey(key);
}

template <>
unsigned long hash_key<PEdgeBase>(const PEdgeBase& key)
{
	return HashKey(key);
}

template <>
unsigned long hash_key<PEdgeDescriptor>(const PEdgeDescriptor& key)
{
	return HashKey(key);
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for faces
template <>
unsigned long hash_key<PFaceVertices>(const PFaceVertices& key)
{
	return HashKey(key);
}

template <>
unsigned long hash_key<PFace>(const PFace& key)
{
	return HashKey(key);
}

template <>
unsigned long hash_key<PFaceDescriptor>(const PFaceDescriptor& key)
{
	return HashKey(key);
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for volumes
template <>
unsigned long hash_key<PVolumeVertices>(const PVolumeVertices& key)
{
	return HashKey(key);
}

template <>
unsigned long hash_key<PVolume>(const PVolume& key)
{
	return HashKey(key);
}

template <>
unsigned long hash_key<PVolumeDescriptor>(const PVolumeDescriptor& key)
{
	return HashKey(key);
}

}//	end of namespace
