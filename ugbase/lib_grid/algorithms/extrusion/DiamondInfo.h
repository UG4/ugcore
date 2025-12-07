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

//template <
//typename VOLELEM,
//typename MANIFELEM,
//typename = std::enable_if< std::is_pointer<VOLELEM>::value>,
//typename = std::enable_if< std::is_pointer<MANIFELEM>::value>
//>
//class VolManifCombi
//{
//
//public:
//
//	using VolPair = std::pair<VOLELEM>;
//
//	VolManifCombi( MANIFELEM const & manifEl, VolPair const & volPr )
//	: m_manifElm(manifEl),
//	  m_volPr(volPr),
//	  m_sudo(-1)
//	{}
//
//	void spuckVolManifCombi( MANIFELEM & manifEl, VolPair & volPr )
//	{
//		manifEl = m_manifElm;
//		volPr = m_volPr;
//	}
//
//private:
//
//	MANIFELEM m_manifElm;
//	VolPair m_volPr;
//	int m_sudo;
//
//};

//template<
//typename VERTEXTYPE
//>
//class CrossingVerticesInfo
//{
//public:
//
//	CrossingVerticesInfo( VERTEXTYPE const & oldVrtx, VERTEXTYPE const & shiftVrtx, bool  )
//
//
//};
//
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

	using VrtxPairVec = std::vector<std::pair<VERTEXTYP,VERTEXTYP>>;

	VolManifVrtxCombi( VOLELEM const & vol,
					   MANIFELEM const & manif,
					   VrtxPairVec const & oldAndShiftVrtcs,
					   INDEXTYP	sudo
					 )
	: m_volElm(vol),
	  m_manifElm(manif),
	  m_oldAndshiftVrtcs(oldAndShiftVrtcs),
	  m_sudo(sudo)
	{};

	void spuckVol( VOLELEM & vol ) { vol = m_volElm; }

	void spuckManif( MANIFELEM & manif ) { manif = m_manifElm; }

	INDEXTYP spuckOldAndShiftVrtcs( VrtxPairVec & vrtpv )
	{
		vrtpv = m_oldAndshiftVrtcs;
		return vrtpv.size();
	}

	INDEXTYP spuckSudo() { return m_sudo; }

private:

	VOLELEM m_volElm;
	MANIFELEM m_manifElm;
	VrtxPairVec m_oldAndshiftVrtcs;
	INDEXTYP m_sudo;

};


/////////////////////////////////////////////////////////////////
//
//template <
//typename VOLELEM,
//typename MANIFELEM,
//typename LOWDIMELM,
//typename VERTEXTYP,
//typename VECTORTYP,
//typename INDEXTYP,
//typename = std::enable_if< std::is_pointer<VOLELEM>::value>,
//typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
//typename = std::enable_if< std::is_pointer<LOWDIMELM>::value>,
//typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
//>
//class DiamondSegmentInfo
//{
//
//public:
//
//	using VolManifCombiVec = std::vector<VolManifCombi<VOLELEM,MANIFELEM>;
//	using VolPair = std::pair<VOLELEM>;
//	using LowDimPair = std::pair<LOWDIMELM>;
//	using VrtxPair = std::pair<VERTEXTYP>;
//
//	DiamondSegmentInfo( // VERTEXTYP const & vertexOld,
//				 VrtxPair const & vertcsShift,
//				 VolManifCombiVec const & volManifCV,
//				 LowDimPair const & lowdimElems,
//				 INDEXTYP sudo )
//	: // m_vertexOld(vertexOld),
//	  m_vertcsShift(vertcsShift),
//	  m_volManifCV(volManifCV),
//	  m_lowdimElems(lowdimElems),
//	  m_sudo(sudo)
//	{};
//
////	void spuckVertexOld( VERTEXTYP & vertexOld )
////	{
////		vertexOld = m_vertexOld;
////	}
//
//	void spuckVrtcsShift( VrtxPair & vertcsShift )
//	{
//		vertcsShift = m_vertcsShift;
//	}
//
//	INDEXTYP spuckvolManifCombiVec( VolManifCombiVec & volManifCV )
//	{
//		volManifCV = m_volManifCV;
//		return volManifCV.size();
//	}
//
//	void spuckLowdimElems( LowDimPair & lowdimElems )
//	{
//		lowdimElems = m_lowdimElems;
//	}
//
//	INDEXTYP spuckSudo()
//	{
//		return m_sudo;
//	}
//
//
//
//private:
//
////	VERTEXTYP m_vertexOld;
//	VrtxPair m_vertcsShift;
//	VolManifCombiVec m_volManifCV;
//	LowDimPair m_lowdimElems;
//	INDEXTYP m_sudo;
//
//};
//
//template <
//typename VOLELEM,
//typename MANIFELEM,
//typename LOWDIMELM,
//typename VERTEXTYP,
//typename VECTORTYP,
//typename INDEXTYP,
//typename = std::enable_if< std::is_pointer<VOLELEM>::value>,
//typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
//typename = std::enable_if< std::is_pointer<LOWDIMELM>::value>,
//typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
//>
//class DiamondInfo
//{
//
//public:
//
//	DiamondInfo( VERTEXTYP const & vertexOld )
//	: m_vertexOld(vertexOld)
//	{};
//
//	void spuckVertexOld( VERTEXTYP & vertexOld )
//	{
//		vertexOld = m_vertexOld;
//	}
//
//private:
//
//	VERTEXTYP m_vertexOld;
//};



} // end of namespace diamonds

} // end of namespace ug





#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_ */
