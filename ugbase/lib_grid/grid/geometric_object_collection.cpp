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
GeometricObjectCollection(GeometricObjectSectionContainer* pVrtSection,
									GeometricObjectSectionContainer* pEdgeSection,
									GeometricObjectSectionContainer* pFaceSection,
									GeometricObjectSectionContainer* pVolSection)
{
	m_levels.reserve(1);
	add_level(pVrtSection, pEdgeSection, pFaceSection, pVolSection);
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
add_level(	GeometricObjectSectionContainer* pVrtSection,
			GeometricObjectSectionContainer* pEdgeSection,
			GeometricObjectSectionContainer* pFaceSection,
			GeometricObjectSectionContainer* pVolSection)
{
	m_levels.push_back(ContainerCollection(	pVrtSection,
											pEdgeSection,
											pFaceSection,
											pVolSection));
}


GeometricObjectCollection::ContainerCollection::
ContainerCollection(GeometricObjectSectionContainer* vrtCon,
					GeometricObjectSectionContainer* edgeCon,
					GeometricObjectSectionContainer* faceCon,
					GeometricObjectSectionContainer* volCon)
{
	pSectionContainers[0] = vrtCon;
	pSectionContainers[1] = edgeCon;
	pSectionContainers[2] = faceCon;
	pSectionContainers[3] = volCon;
}

}//	end of namespace
