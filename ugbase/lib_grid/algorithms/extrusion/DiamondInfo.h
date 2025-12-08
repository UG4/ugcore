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
typename VERTEXTYP,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<VOLELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
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
	  m_sudo(sudo)
	{};

	void spuckVol( VOLELEM & vol ) { vol = m_volElm; }

	void spuckManif( MANIFELEM & manif ) { manif = m_manifElm; }

	void spuckOldAndShiftVrtx( VrtxPair & vrtp ) { vrtp = m_oldAndshiftVrtx; }

	INDEXTYP spuckSudo() { return m_sudo; }

	void changeVol( VOLELEM const & vol ) { m_volElm = vol; }

	void changeManif( MANIFELEM const & manif ) { m_manifElm = manif; }


private:

	VOLELEM m_volElm;
	MANIFELEM m_manifElm;
	VrtxPair m_oldAndshiftVrtx;
	INDEXTYP m_sudo;

	// TODO FIXME compute automatisch the edge which connects the both vertices!!!!

};


/////////////////////////////////////////////////////////////////



} // end of namespace diamonds

} // end of namespace ug





#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_ */
