// created by Sebastian Reiter
// y09 m07 d31
// s.b.reiter@googlemail.com

#include "geometric_object_collection.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GeometricObjectCollection
GeometricObjectCollection::
GeometricObjectCollection(size_t levelEstimate)
{
	m_levels.reserve(levelEstimate);
}

GeometricObjectCollection::
GeometricObjectCollection(ElementStorage<VertexBase>::SectionContainer* vrtCon,
							ElementStorage<EdgeBase>::SectionContainer* edgeCon,
							ElementStorage<Face>::SectionContainer* faceCon,
							ElementStorage<Volume>::SectionContainer* volCon)
{
	m_levels.reserve(1);
	add_level(vrtCon, edgeCon, faceCon, volCon);
}

GeometricObjectCollection::
GeometricObjectCollection(const GeometricObjectCollection& mgoc)
{
	assign(mgoc);
}

GeometricObjectCollection&
GeometricObjectCollection::
operator =(const GeometricObjectCollection& mgoc)
{
	assign(mgoc);
	return *this;
}

void
GeometricObjectCollection::
assign(const GeometricObjectCollection& mgoc)
{
	m_levels.resize(mgoc.num_levels());
	for(size_t i = 0; i < m_levels.size(); ++i)
		m_levels[i] = mgoc.m_levels[i];
}

void
GeometricObjectCollection::
add_level(ElementStorage<VertexBase>::SectionContainer* vrtCon,
			ElementStorage<EdgeBase>::SectionContainer* edgeCon,
			ElementStorage<Face>::SectionContainer* faceCon,
			ElementStorage<Volume>::SectionContainer* volCon)
{
	m_levels.push_back(ContainerCollection(	vrtCon,
											edgeCon,
											faceCon,
											volCon));
}


GeometricObjectCollection::ContainerCollection::
ContainerCollection(ElementStorage<VertexBase>::SectionContainer* vrtCon,
					ElementStorage<EdgeBase>::SectionContainer* edgeCon,
					ElementStorage<Face>::SectionContainer* faceCon,
					ElementStorage<Volume>::SectionContainer* volCon)
{
	vrtContainer = vrtCon;
	edgeContainer = edgeCon;
	faceContainer = faceCon;
	volContainer = volCon;
}

}//	end of namespace
