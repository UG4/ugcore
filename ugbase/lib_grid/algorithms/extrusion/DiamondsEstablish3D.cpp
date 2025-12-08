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
	IndexType sudos = m_sh.num_subsets();

	for( auto & vmvcd : m_vecVolManifVrtxCombiToShrink4Diams )
	{
		Volume* vol;
		vmvcd.spuckVol(vol);

		m_sh.assign_subset(vol, sudos);
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
