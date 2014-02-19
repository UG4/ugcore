// created by Sebastian Reiter
// y09 m07 d31
// s.b.reiter@googlemail.com

#include "grid_object_collection.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GridObjectCollection
GridObjectCollection::
GridObjectCollection(size_t levelEstimate)
{
	m_levels.reserve(levelEstimate);
}

GridObjectCollection::
GridObjectCollection(ElementStorage<Vertex>::SectionContainer* vrtCon,
							ElementStorage<EdgeBase>::SectionContainer* edgeCon,
							ElementStorage<Face>::SectionContainer* faceCon,
							ElementStorage<Volume>::SectionContainer* volCon)
{
	m_levels.reserve(1);
	add_level(vrtCon, edgeCon, faceCon, volCon);
}

GridObjectCollection::
GridObjectCollection(const GridObjectCollection& mgoc)
{
	assign(mgoc);
}

GridObjectCollection&
GridObjectCollection::
operator =(const GridObjectCollection& mgoc)
{
	assign(mgoc);
	return *this;
}

void
GridObjectCollection::
assign(const GridObjectCollection& mgoc)
{
	m_levels.resize(mgoc.num_levels());
	for(size_t i = 0; i < m_levels.size(); ++i)
		m_levels[i] = mgoc.m_levels[i];
}

void
GridObjectCollection::
add_level(ElementStorage<Vertex>::SectionContainer* vrtCon,
			ElementStorage<EdgeBase>::SectionContainer* edgeCon,
			ElementStorage<Face>::SectionContainer* faceCon,
			ElementStorage<Volume>::SectionContainer* volCon)
{
	m_levels.push_back(ContainerCollection(	vrtCon,
											edgeCon,
											faceCon,
											volCon));
}


GridObjectCollection::ContainerCollection::
ContainerCollection(ElementStorage<Vertex>::SectionContainer* vrtCon,
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
