/*
 * DiamondInfo.h
 *
 *  Created on: 05.12.2025
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <type_traits>

namespace ug
{

////////////////////////////////////////////7

namespace diamonds
{

///////////////////////////////////////////

////////////////////////////////////////////////////////////////

template <
typename VOLELEM,
typename MANIFELEM,
typename LOWDIMELM,
typename VERTEXTYP,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<VOLELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class VolManifVrtxCombi
{
public:

	using VrtxPair = std::pair<VERTEXTYP,VERTEXTYP>;

	VolManifVrtxCombi( VOLELEM const & vol,
					   MANIFELEM const & manif,
					   VrtxPair const & oldAndShiftVrtx,
					   INDEXTYP	sudo
					 )
	: m_volElm(vol),
	  m_manifElm(manif),
	  m_oldAndshiftVrtx(oldAndShiftVrtx),
	  m_sudo(sudo),
	  m_lowDimElm(nullptr)
	{
//		INDEXTYP numLowdimElmsFnd = computeTheLowdimElm();
//
//		if( numLowdimElmsFnd != 1 )
//		{
//			UG_LOG("number of lowdim elems found strange " << numLowdimElmsFnd << std::endl);
//			UG_THROW("number of lowdim elems found strange " << numLowdimElmsFnd << std::endl);
//		}
//
//		if( m_lowDimElm == nullptr )
//		{
//			UG_LOG("Edge nicht gefunden " << std::endl);
//			UG_THROW("Edge nicht gefunden " << std::endl);
//		}
	};

	void spuckVol( VOLELEM & vol ) { vol = m_volElm; }

	void spuckManif( MANIFELEM & manif ) { manif = m_manifElm; }

	void spuckOldAndShiftVrtx( VrtxPair & vrtp ) { vrtp = m_oldAndshiftVrtx; }

	void spuckLowdimElem ( LOWDIMELM & lowdimElm ) { lowdimElm = m_lowDimElm; }

	INDEXTYP spuckSudo() { return m_sudo; }

	void changeVol( VOLELEM const & vol ) { m_volElm = vol; }

	void changeManif( MANIFELEM const & manif ) { m_manifElm = manif; }

	// compute automatisch the edge which connects the both vertices!!!!
	template
	<
		typename = std::enable_if<std::is_same<Volume*,VOLELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELM>::value>,
		typename = std::enable_if<std::is_same<Vertex*,VERTEXTYP>::value>
	>
	bool computeTheLowdimElm( Grid & grid )
	{
		INDEXTYP edgeFound = 0;

		for(size_t i_edge = 0; i_edge < m_volElm->num_edges(); ++i_edge)
		{
			LOWDIMELM lowDimElm = grid.get_edge( m_volElm, i_edge );

			if(    EdgeContains(lowDimElm, m_oldAndshiftVrtx.first)
				&& EdgeContains(lowDimElm, m_oldAndshiftVrtx.second )
			   )
			{
				m_lowDimElm = lowDimElm;
				edgeFound++;
			}
		}

		return (edgeFound == 1);
	}


private:

	VOLELEM m_volElm;
	MANIFELEM m_manifElm;
	VrtxPair m_oldAndshiftVrtx;
	INDEXTYP m_sudo;
	LOWDIMELM m_lowDimElm;



};


/////////////////////////////////////////////////////////////////



} // end of namespace diamonds

} // end of namespace ug





#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_ */
