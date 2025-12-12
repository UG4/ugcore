/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include "grid_base_objects.h"

#include "grid_util.h"

namespace ug {

const char* GRID_BASE_OBJECT_SINGULAR_NAMES[] = {"vertex", "edge", "face", "volume"};
const char* GRID_BASE_OBJECT_PLURAL_NAMES[] = {"vertices", "edges", "faces", "volume"};

////////////////////////////////////////////////////////////////////////
//	implementation of edge
bool Edge::get_opposing_side(Vertex* v, Vertex** vrtOut)
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

EdgeDescriptor::EdgeDescriptor(const EdgeDescriptor& ed)
{
	EdgeVertices::assign_edge_vertices(ed);
}

EdgeDescriptor::EdgeDescriptor(Vertex* vrt1, Vertex* vrt2)
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
size_t hash_key(Vertex* key)
{
	return (unsigned long)key->get_hash_value();
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for edges
size_t hash_key(EdgeVertices* key)
{
	return HashKey(key);
}

size_t hash_key(const EdgeVertices* key)
{
	return HashKey(key);
}

size_t hash_key(Edge* key)
{
	return HashKey(key);
}

size_t hash_key(EdgeDescriptor* key)
{
	return HashKey(key);
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for faces
size_t hash_key(FaceVertices* key)
{
	return HashKey(key);
}

size_t hash_key(const FaceVertices* key)
{
	return HashKey(key);
}

size_t hash_key(Face* key)
{
	return HashKey(key);
}

size_t hash_key(FaceDescriptor* key)
{
	return HashKey(key);
}

////////////////////////////////////////////////////////////////////////
//	hash-funtions for volumes
size_t hash_key(VolumeVertices* key)
{
	return HashKey(key);
}

size_t hash_key(const VolumeVertices* key)
{
	return HashKey(key);
}

size_t hash_key(Volume* key)
{
	return HashKey(key);
}

size_t hash_key(VolumeDescriptor* key)
{
	return HashKey(key);
}

}//	end of namespace
