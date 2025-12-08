/*
 * DiamondsEstablish3D.cpp
 *
 *  Created on: 08.12.2025
 *      Author: Markus Knodel
 */

#include <lib_grid/algorithms/extrusion/DiamondsEstablish3D.h>

namespace ug
{

namespace diamonds
{

DiamondsEstablish3D::DiamondsEstablish3D( Grid & grid,
										  SubsetHandler & sh,
										  DiamondsEstablish3D::VecVolManifVrtxCombi const & vecVolManifVrtxC
										 )
		:
			m_grid(grid),
			m_sh(sh),
			m_vecVolManifVrtxCombiToShrink4Diams(vecVolManifVrtxC)
{
}

bool DiamondsEstablish3D::createTheDiamonds()
{
	IndexType sudosVols = m_sh.num_subsets();

	IndexType sudosEdges = sudosVols + 1;

	IndexType sudosFaces = sudosVols + 2;

	for( auto & vmvcd : m_vecVolManifVrtxCombiToShrink4Diams )
	{
		Volume* vol;
		vmvcd.spuckVol(vol);

		m_sh.assign_subset(vol, sudosVols);

		IndexType numLowdimElmsFnd = vmvcd.computeTheLowdimElm(m_grid);

		if( numLowdimElmsFnd != 1 )
		{
			UG_LOG("number of lowdim elems found strange " << numLowdimElmsFnd << std::endl);
			UG_THROW("number of lowdim elems found strange " << numLowdimElmsFnd << std::endl);
		}

		Edge* edge;
		vmvcd.spuckLowdimElem( edge );

		if( edge == nullptr )
		{
			UG_LOG("Edge nicht gefunden " << std::endl);
			UG_THROW("Edge nicht gefunden " << std::endl);
		}
		m_sh.assign_subset( edge, sudosEdges);

		Face * fac;
		vmvcd.spuckManif(fac);
		m_sh.assign_subset( fac, sudosFaces );
	}

	UG_LOG("Established diamonds" << std::endl);


	return true;
}

DiamondsEstablish3D::~DiamondsEstablish3D()
{
	// Auto-generated destructor stub
}

} /* namespace diamonds */


} /* namespace ug */
