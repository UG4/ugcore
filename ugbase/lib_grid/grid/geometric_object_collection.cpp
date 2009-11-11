// created by Sebastian Reiter
// y09 m07 d31
// s.b.reiter@googlemail.com

#include "geometric_object_collection.h"

namespace ug
{

GeometricObjectCollection::GeometricObjectCollection()
{
	for(uint i = 0; i < 4; ++i)
		m_pSectionContainers[i] = NULL;
}

GeometricObjectCollection::GeometricObjectCollection(
									GeometricObjectSectionContainer* pVrtSection,
									GeometricObjectSectionContainer* pEdgeSection,
									GeometricObjectSectionContainer* pFaceSection,
									GeometricObjectSectionContainer* pVolSection)
{
	m_pSectionContainers[0] = pVrtSection;
	m_pSectionContainers[1] = pEdgeSection;
	m_pSectionContainers[2] = pFaceSection;
	m_pSectionContainers[3] = pVolSection;
}
		
GeometricObjectCollection::GeometricObjectCollection(
									const GeometricObjectCollection& goc)
{
//	simply copy the pointers.
	for(uint i = 0; i < 4; ++i)
		m_pSectionContainers[i] = goc.m_pSectionContainers[i];
}

GeometricObjectCollection&
GeometricObjectCollection::operator =(const GeometricObjectCollection& goc)
{
//	simply copy the pointers.
	for(uint i = 0; i < 4; ++i)
		m_pSectionContainers[i] = goc.m_pSectionContainers[i];
	return *this;
}
								

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	MultiLevelGeometricObjectCollection
MultiLevelGeometricObjectCollection::MultiLevelGeometricObjectCollection()
{
}

MultiLevelGeometricObjectCollection::
MultiLevelGeometricObjectCollection(int levelEstimate)
{
	m_levels.reserve(levelEstimate);
}
		
MultiLevelGeometricObjectCollection::
MultiLevelGeometricObjectCollection(const MultiLevelGeometricObjectCollection& mgoc)
{
	assign(mgoc);
}

MultiLevelGeometricObjectCollection&
MultiLevelGeometricObjectCollection::
operator =(const MultiLevelGeometricObjectCollection& mgoc)
{
	assign(mgoc);
	return *this;
}

void
MultiLevelGeometricObjectCollection::
assign(const MultiLevelGeometricObjectCollection& mgoc)
{
	m_levels.resize(mgoc.num_levels());
	for(size_t i = 0; i < m_levels.size(); ++i)
		m_levels[i] = mgoc.m_levels[i];
}

void
MultiLevelGeometricObjectCollection::
add_level(	GeometricObjectSectionContainer* pVrtSection,
			GeometricObjectSectionContainer* pEdgeSection,
			GeometricObjectSectionContainer* pFaceSection,
			GeometricObjectSectionContainer* pVolSection)
{
	m_levels.push_back(GeometricObjectCollection(	pVrtSection,
													pEdgeSection,
													pFaceSection,
													pVolSection));
}
				
}
