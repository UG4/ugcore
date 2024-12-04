/*
 * support3D.h
 *
 *  Created on: 31.10.2024
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_

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

#include "support.h"

// TODO FIXME
// verschlanken,
// Verdoppelungen weg,
// Namen generalisieren
// nicht vertex spezifisch zB
// sondern lowest dim oder so

namespace ug
{

namespace support
{
template<typename ELEMTYP, typename INDEX_TYP>
class ElemInfo
{
public:

	ElemInfo( ELEMTYP const & elem, INDEX_TYP sudo )
	: m_elem(elem), m_sudo(sudo)
	{
	}



private:

	ELEMTYP m_elem;
	INDEX_TYP m_sudo;

	ElemInfo() {};


};

template<
typename FACETYP,
typename NORMALTYP,
typename VOLUMETYP,
typename EDGETYP,
typename INDEX_TYP
>
class VertexFractureQuadrupel
{
public:

	using ElemInfoFac = ElemInfo<FACETYP,INDEX_TYP>;

	//face, normal, volume, edge

//	VertexFractureQuadrupel()
//	{};


	VertexFractureQuadrupel( ElemInfoFac const & fracFaceInfo,
							 VOLUMETYP const & attVolume,
							 NORMALTYP const & normal,
							 std::pair<EDGETYP,EDGETYP> const & volCutEdges,
							 std::pair<ElemInfoFac,ElemInfoFac> const & volCutEdgeFaces )
	: m_fracFaceInfo(fracFaceInfo),
	  m_attVolume(attVolume),
	  m_normal(normal),
	  m_volCutEdges(volCutEdges),
	  m_volCutEdgeFaces(volCutEdgeFaces)
	{
	}

	// todo fixme getter und ggf auch setter, aber vermutlich nur getter implementieren!!!

private:

//	FACETYP const getFace() const { return m_full; }
//	NORMALTYP const getNormal() const { return m_normal; }
//	VOLUMETYP const getVolume() const { return m_volume; }
//	EDGETYP const getEdge() const { return m_edge; }

	ElemInfo<FACETYP,INDEX_TYP> m_fracFaceInfo;
	VOLUMETYP m_attVolume;
	NORMALTYP m_normal;
	std::pair<EDGETYP,EDGETYP> m_volCutEdges;
	std::pair<ElemInfoFac,ElemInfoFac> m_volCutEdgeFaces;

//private:
//
//	FACETYP m_face;
//	NORMALTYP m_normal;
//	VOLUMETYP m_volume;
//	EDGETYP m_edge;

	VertexFractureQuadrupel()
	{};
};

// TODO FIXME vertex fracture triplett
// vereinigen mit AttElem!!!!
// doppelte Strukturen!!!


template <
typename MANIFOLDTYP, // 3D: Face
typename INDEXTYP, // int oder unsinged int oder short oder unsigned short etc
typename FULLDIMTYP, // 3D: Volume
typename SENKRECHTENTYP, // 3D und 2D: ug::vector3
typename LOWDIMTYP // 3D: Edge (2D nicht benötigt)
>
class VertexFractureTripleMF
{

public:

	VertexFractureTripleMF( MANIFOLDTYP const & manifElm, INDEXTYP sudo,
							FULLDIMTYP const & fullElm,
							SENKRECHTENTYP const & normal,
							std::pair<LOWDIMTYP,LOWDIMTYP> const & pairLowElm )
	: m_manifElm(manifElm), m_sudo(sudo), m_fullElm(fullElm),
	  m_normal(normal), m_newNormal(normal),
	  m_pairLowElm(pairLowElm)
	{
	};

	MANIFOLDTYP const getManifElm() const { return m_manifElm; }

	INDEXTYP const getSudoElm() const { return m_sudo; }

	FULLDIMTYP const getFullElm() const { return m_fullElm; }

	SENKRECHTENTYP const getNormal() const { return m_normal; }

	void setNewNormal( SENKRECHTENTYP const & chNorml ) { m_newNormal = chNorml; }

	SENKRECHTENTYP const getNewNormal() const { return m_newNormal; }

	std::pair<LOWDIMTYP,LOWDIMTYP> const getPairLowElm() const { return m_pairLowElm; }

private:

	MANIFOLDTYP m_manifElm;
	INDEXTYP m_sudo;
	FULLDIMTYP m_fullElm;
	SENKRECHTENTYP  m_normal;
	SENKRECHTENTYP  m_newNormal;
	std::pair<LOWDIMTYP,LOWDIMTYP> m_pairLowElm;

	VertexFractureTripleMF()
	{};

};


//////////////////////////////////////////////////////////////////

// TODO FIXME das muss angepasst werden, ist noch wie für 2D Fall bisher
enum FracTypVol { SingleFrac = 2, TEnd = 3, XCross = 4 };

template < typename VRT, typename IndTyp >
class CrossingVertexInfoVol
{

public:

	CrossingVertexInfoVol( VRT const & crossVrt, FracTypVol fracTyp )
	: m_crossVrt(crossVrt), m_fracTyp( fracTyp ),
	  m_vecShiftedVrts(std::vector<VRT>())
	,  m_vecShiftedVrtsWithTypInf(std::vector<std::pair<VRT,bool>>())
	, m_numberAtFreeSide(0)
	{
	}

	VRT getCrossVertex() const { return m_crossVrt; }

	FracTypVol getFracTyp() const { return m_fracTyp; }

	void addShiftVrtx( VRT const & vrt, bool isAtFreeSide = false )
	{
		m_vecShiftedVrts.push_back(vrt);

//		if( m_fracTyp == TEnd )
//		{
//			std::pair<VRT, bool > addSVI( vrt, isAtFreeSide );
//			m_vecShiftedVrtsWithTypInf.push_back(addSVI);
//
//			if( isAtFreeSide )
//				m_numberAtFreeSide++;
//
//			if( m_numberAtFreeSide > 1 )
//				UG_THROW("was ist das fuer ein T Ende" << std::endl);
//		}

	}

	void setShiftVrtx( std::vector<VRT> const & vecVrt ) { m_vecShiftedVrts = vecVrt; }

	std::vector<VRT> getVecShiftedVrts() const
	{
		return m_vecShiftedVrts;
	}

	std::vector<std::pair<VRT,bool>> getVecShiftedVrtsWithTypInfo() const
	{
//		if( m_fracTyp != TEnd )
//			UG_THROW("fuer Kreuz nicht erlaubt " << std::endl);

		return m_vecShiftedVrtsWithTypInf;
	}


private:

	VRT m_crossVrt;
	std::vector<VRT> m_vecShiftedVrts;
	std::vector<std::pair<VRT,bool>> m_vecShiftedVrtsWithTypInf;
	FracTypVol m_fracTyp;
	IndTyp m_numberAtFreeSide;
};

///////////////////////////////////////////////////////////////////////////////

// info for a vertex: face and attached edges
template <
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP
>
class AttachedElem
{
public:
	using PairLowEl = std::pair<LOWDIMELM,LOWDIMELM>;

	AttachedElem( MANIFELM const & manifElm,
				  PairLowEl & lowElm,
				  INDEX_TXP sudo )
	:
		m_manifElm(manifElm), m_lowElm(lowElm), m_sudo(sudo)
	{
	};

	MANIFELM const getManifElm() const { return m_manifElm;}
	PairLowEl const getLowElm() const { return m_lowElm; }
	INDEX_TXP const getSudo() const { return m_sudo; };

	bool const testIfEquals( AttachedElem<MANIFELM,LOWDIMELM,INDEX_TXP> const & attElm )
	const
	{
		MANIFELM manifElmOther = attElm.getManifElm();
		PairLowEl lowElmOther = attElm.getLowElm();
		INDEX_TXP sudoOther = attElm.getSudo();

		if(    manifElmOther == this->m_manifElm
			&& lowElmOther == this->m_lowElm
			&& sudoOther == this->m_sudo
		)
		{
			return true;
		}

		return false;
	}

	bool const isNeighboured(  AttachedElem<MANIFELM,LOWDIMELM,INDEX_TXP> const & attElm )
	const
	{
//		MANIFELM manifElmOther = attElm.getManifElm();
		PairLowEl lowElmOther = attElm.getLowElm();
//		INDEX_TXP sudoOther = attElm.getSudo();

		PairLowEl lowElmThis = this->m_lowElm;

		std::vector<bool> test;

		test.push_back( lowElmOther.first  == lowElmThis.first );
		test.push_back( lowElmOther.second  == lowElmThis.first );
		test.push_back( lowElmOther.first  == lowElmThis.second );
		test.push_back( lowElmOther.second  == lowElmThis.second );

		INDEX_TXP countCorr = 0;

		for( auto const t : test )
		{
			if( t )
				countCorr++;
		}

		if( countCorr == 1 )
			return true;

		if( countCorr > 1 )
			UG_THROW("zu viele gleiche Ecken " << std::endl);

		return false;
	}

	bool const isNeighbouredAtSpecificSide( AttachedElem<MANIFELM,LOWDIMELM,INDEX_TXP> const & attElm,
			  	  	  	  	  	  	  	  	LOWDIMELM const & specificLDE )
	const
	{
		PairLowEl lowElmOther = attElm.getLowElm();

		PairLowEl lowElmThis = this->m_lowElm;

		// test if the specific element is part of at least
		// one of the faces

		bool otherFirst = ( lowElmOther.first == specificLDE );
		bool otherSecond = ( lowElmOther.second == specificLDE );

		bool thisFirst = ( lowElmThis.first == specificLDE );
		bool thisSecond = ( lowElmThis.second == specificLDE );

		bool isPartOfThisFace = ( thisFirst || thisSecond );
		bool isPartOfOtherFace = ( otherFirst || otherSecond );

		if( ! isPartOfOtherFace || ! isPartOfThisFace )
		{
			UG_LOG("not part of one of the faces " << std::endl);
			return false;
		}

		if( otherFirst && thisFirst )
		{
			if( lowElmOther.first  == lowElmThis.first )
				return true;
		}
		else if( otherFirst && thisSecond )
		{
			if( lowElmOther.first  == lowElmThis.second )
				return true;
		}
		else if( otherSecond && thisFirst )
		{
			if( lowElmOther.second  == lowElmThis.first )
				return true;
		}
		else if( otherSecond && thisSecond )
		{
			if( lowElmOther.second  == lowElmThis.second )
				return true;
		}

		return false;
	}

private:
	MANIFELM m_manifElm;
	PairLowEl m_lowElm;
	INDEX_TXP m_sudo;


};

//////////////////////////////////////////////////////////////////////////////

// class to help count and store a bool and a number of templete type
// comparable to std::pair<bool,int> but more dedicated to the specific aim
// TODO FIXME adapt for 3D case, figure out if inner end, and number of fracs sourrounding
// CAUTION is also used for edges, but still uses
// vertex as indicator - name should be made more flexible

template<
typename T,
typename ATT_ELEM
>
class VertexFracturePropertiesVol
{
public:

	using pairTB = std::pair<T,bool>;
	using VecPairTB = std::vector<pairTB>;


//	VertexFracturePropertiesVol( bool isBndFracVertex,  T numberCrossingFracsInVertex )
//	: m_isBndFracVertex(isBndFracVertex), m_numberCountedFracsInVertex(numberCrossingFracsInVertex)
//	{
//	};

	enum VrtxFracStatus { noFracSuDoAtt = 0,
						  oneFracSuDoAtt = 1,
						  twoFracSuDoAtt = 2,
						  threeFracSuDoAtt = 3 };

	VertexFracturePropertiesVol()
	: m_isBndFracVertex(false), m_numberCountedFracsInVertex(0),
	  m_status( noFracSuDoAtt ),
	  m_sudoList( std::vector<T>() ),
	  m_sudosClosed(VecPairTB()),
	  m_vecAttElem(std::vector<ATT_ELEM>())
//		VertexFracturePropertiesVol( false, 0 )
	{
	};

	void setIsBndFracVertex( bool iBDV = true )
	{
		m_isBndFracVertex = iBDV;
	}

	void setNumberCrossingFracsInVertex( T const & nCFIV )
	{
		m_numberCountedFracsInVertex = nCFIV;
	}

	bool getIsBndFracVertex()
	{
		return m_isBndFracVertex;
	}

// 	T getCountedNumberFracsInVertex()
// 	{
// 		return m_numberCountedFracsInVertex;
// 	}


	T getNumberFracEdgesInVertex()
	{
		return m_numberCountedFracsInVertex;
	}

// 	T getNumberCrossingFracsInVertex()
//  {
// 		if( m_isBndFracVertex )
// 			return m_numberCountedFracsInVertex;
//
// 		// for inner vertices, each edge passed when
// 		// fractures are counted along their edges
// 		// that the vertizes get hit twice for each fracture run
// 		// only for boundary vertices, this happens only once per fracture
// 		T multipeInnerHits = 2;
//
// 		T rest = m_numberCountedFracsInVertex % multipeInnerHits;
//
// 		if( rest != 0 )
// 		{
//// 			UG_THROW("Expand layers: rest division frac counting not zero " << m_numberCountedFracsInVertex << std::endl);
//
// 			throw std::runtime_error("error");
//
// 			return 0;
// 		}
//
// 		return m_numberCountedFracsInVertex / multipeInnerHits;
//  }

	VertexFracturePropertiesVol & operator++( int a )
	{
		m_numberCountedFracsInVertex++;
		return *this;
	}


	VrtxFracStatus getVrtxFracStatus()
	{
		return m_status;
	}

	// caution: returns false if sudo already known, but no problem, feature, not bug
	// so never stop the program if false returned here, this is good news
	bool addFractSudo( T const & sudo )
	{
		bool alreadyInList = false;

		for( auto const & availSudo : m_sudoList )
		{
			if( sudo == availSudo )
				alreadyInList = true;
		}

		if( ! alreadyInList )
		{
			m_sudoList.push_back( sudo );
		}

		return adaptVrtxFracStatus();
	}

	// setter private to avoid abusive use
	std::vector<T> const getSudoList() const
	{
		return m_sudoList;
	}

	bool setIsAClosedFracture( T sudoNow, bool isClosedNow )
	{

		T alreadyKnownMult = 0;

		for( auto & suSu : m_sudosClosed )
		{
			T & sudoVal = suSu.first;
			bool & isClosedVal = suSu.second;

			if( sudoVal == sudoNow )
			{
				alreadyKnownMult++;

				UG_LOG("Reassign sudo surround " << std::endl);

				if( isClosedVal != isClosedNow )
					UG_THROW("change property sudo surrounded, why?" << std::endl);

				isClosedVal = isClosedNow;
			}
		}

		if( alreadyKnownMult == 0 )
		{
			pairTB infoSudoSurr( sudoNow, isClosedNow );

			m_sudosClosed.push_back( infoSudoSurr );

		}
		else if( alreadyKnownMult > 1 )
		{
			UG_THROW("zu oft bekannt " << std::endl);
			return false;
		}

		// check if now correct

		T testKnownFine = 0;

		for( auto const & suSu : m_sudosClosed )
		{
			T & sudoVal = suSu.first;
			bool & isClosedVal = suSu.second;

			if( sudoVal == sudoNow )
			{
				testKnownFine++;

				if( isClosedVal != isClosedNow )
				{
					UG_THROW("NOT set property sudo surrounded, why?" << std::endl);
					return false;
				}

			}
		}

		if( testKnownFine == 0 || testKnownFine > 1 )
		{
			UG_THROW("immer noch nicht bekannt?" << std::endl);
			return false;

		}

		return true;
	}

	bool getIsAClosedFracture( T sudoNow )
	{
		T foundMultpl = 0;

		bool isClosedReturn = false;

		for( auto const & suSu : m_sudosClosed )
		{
			T const & sudoVal = suSu.first;
			bool const & isClosedVal = suSu.second;

			if( sudoVal == sudoNow )
			{
				foundMultpl++;
				isClosedReturn = isClosedVal;
			}
		}

		if( foundMultpl != 1 )
		{
			UG_THROW("not known status closed or not sudo" << std::endl);
			return false;
		}

		return isClosedReturn;
	}

	bool setInfoAllFractureSudosIfClosed( VecPairTB const & sudosClosed )
	{
		m_sudosClosed = sudosClosed;

		return true;
	}

	VecPairTB const getInfoAllFracSudosIfClosed() const
	{
		return m_sudosClosed;
	}

	// if all open or closed
	template<bool B>
	bool const getInfoAllFracturesSameClosedState() const
	{
		bool allFracsSame = true;

		for( auto const & suSu : m_sudosClosed )
		{
			T const & sudoVal = suSu.first;
			bool const & isClosedVal = suSu.second;

			if( isClosedVal != B )
				allFracsSame = false;
		}

		return allFracsSame;
	}

	bool addAttachedElem( ATT_ELEM const & attElem )
	{
		bool alreadyKnown = false;

		for( auto const & aE : m_vecAttElem )
		{
			if( aE.testIfEquals(attElem) )
				alreadyKnown = true;
		}

		if( ! alreadyKnown )
			m_vecAttElem.push_back(attElem);

		// returns true if ads it, false if no need as known
		return ! alreadyKnown;
	}

	std::vector<ATT_ELEM> const & getAllAttachedElems()
	const
	{
		return m_vecAttElem;
	}

private:
	bool m_isBndFracVertex;
	T m_numberCountedFracsInVertex;

	VrtxFracStatus m_status;

	static VrtxFracStatus constexpr m_maxStatus = VrtxFracStatus::threeFracSuDoAtt;

	std::vector<T> m_sudoList;

	// better private, to avoid confusion
	bool setVrtxFracStatus( VrtxFracStatus status )
	{
		m_status = status;

		if( status < noFracSuDoAtt || status > threeFracSuDoAtt )
			return false;

		return true;
	}

	bool adaptVrtxFracStatus()
	{
		auto sudosNum = m_sudoList.size();
		if( sudosNum > static_cast<T>( m_maxStatus ) )
		{
			UG_THROW("zu viele subdomains crossing in one Punkt" << std::endl);
			return false;
		}

		m_status = static_cast<VrtxFracStatus>(sudosNum);

		return true;
	}

	VecPairTB m_sudosClosed;

	bool setSudoList( std::vector<T> const & sudoList )
	{
		m_sudoList = sudoList;

		return true;
	}

	std::vector<ATT_ELEM> m_vecAttElem;

};


// intention, explained for volume:
template<typename FULLDIM_ELEM,
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP
>
class AttachedFullDimElemInfo
{

public:

	using AttachedManifElemInfo = AttachedElem<MANIFELM,LOWDIMELM,INDEX_TXP>;

	using VecAttachedManifElemInfo = std::vector<AttachedManifElemInfo>;

	AttachedFullDimElemInfo( FULLDIM_ELEM const & fullDimElm )
	: m_fullDimElm(fullDimElm),
	  m_vecManifElm(VecAttachedManifElemInfo()),
	  m_vecManifElmTouchInfo(VecAttManifElmTouchInf()),
	  m_allSidesTouched(false)
	{
	}

	FULLDIM_ELEM const getFulldimElem() const
	{
		return m_fullDimElm;
	}


	bool addManifElem( AttachedManifElemInfo const & manifElm, Grid & grid )
	{
		// Caution: first test if manifold elem is at all part of the fulldim elem manifols
		// if not, return false directly

		if( ! fullDimElmContainsManif( manifElm.getManifElm(), grid ) )
			return false;

		// if manif elem is in principle part of the fulldim elem manifolds,
		// then we need to check if it is already integrated

		bool hasElemAlready = false;


		for( auto const & me : m_vecManifElm )
		{
			if( manifElm.testIfEquals(me) )
			{
				hasElemAlready = true;
				break;
			}
		}

		if( ! hasElemAlready )
		{
			m_vecManifElm.push_back( manifElm );

			AttManifElmTouchInf pamei( manifElm, false );

			m_vecManifElmTouchInfo.push_back(pamei);
		}

		return ! hasElemAlready;
	}

	VecAttachedManifElemInfo const getVecManifElem() const
	{
		return m_vecManifElm;
	}

	bool const tryToTouchManifElem( AttachedManifElemInfo const & manifElemOther ) const
	{
		bool managed2Touch = false;

		for( auto & ameti : m_vecManifElmTouchInfo )
		{
			AttachedManifElemInfo & manifElmTest = ameti.first;
			bool & alreadyTouched = ameti.second;

			if( ! alreadyTouched )
			{
				if( manifElemOther.testIfEquals( manifElmTest ) )
				{
					alreadyTouched = true;
				}
			}

			if( alreadyTouched )
			{
				managed2Touch = true;
				break;
			}
		}

		return managed2Touch;
	}

	VecAttachedManifElemInfo const getAlreadyTouchedManifElems() const
	{
		VecAttachedManifElemInfo alreadyTouchedManifElms;

		for( const auto & ameti : m_vecManifElmTouchInfo )
		{
			if( ameti.second )
				alreadyTouchedManifElms.push_back( ameti );
		}

		return alreadyTouchedManifElms;
	}

	VecAttachedManifElemInfo const getSoFarUnTouchedManifElems() const
	{
		VecAttachedManifElemInfo unTouchedManifElms;

		for( const auto & ameti : m_vecManifElmTouchInfo )
		{
			if( ! ameti.second )
				unTouchedManifElms.push_back( ameti );
		}

		return unTouchedManifElms;
	}

	bool const testIfAllSidesTouched() const
	{
		if( m_allSidesTouched )
			return true;

		bool allSidesTouched = true;

		for( const auto & ameti : m_vecManifElmTouchInfo )
		{
			if( ! ameti.second )
			{
				allSidesTouched = false;
			}
		}

		m_allSidesTouched = allSidesTouched;

		return m_allSidesTouched;
	}

private:

	using AttManifElmTouchInf = std::pair<AttachedManifElemInfo,bool>;
	using VecAttManifElmTouchInf = std::vector<AttManifElmTouchInf>;

	FULLDIM_ELEM m_fullDimElm;
	VecAttachedManifElemInfo m_vecManifElm;
	VecAttManifElmTouchInf m_vecManifElmTouchInfo;

	bool m_allSidesTouched;

//    template<
////	typename std::enable_if<std::is_same<FULLDIM_ELEM,Volume* const &>::value, true>,
////	typename std::enable_if<std::is_same<MANIFELM,Face* const &>::value, true>
////	typename std::enable_if<std::is_same<FULLDIM_ELEM,Volume* const &>::value, true>,
//	typename std::enable_if<std::is_same<MANIFELM,Face* const &>::value,bool>
//	>
	//std::enable_if<std::is_same<MANIFELM,Face* const &>::value>
	template
	<
	typename = std::enable_if<std::is_same<Volume*,FULLDIM_ELEM>::value>,
	typename = std::enable_if<std::is_same<Face*,MANIFELM>::value>
//	typename = std::enable_if<std::is_same<Volume* const &,FULLDIM_ELEM>::value>,
//	typename = std::enable_if<std::is_same<Face* const &,MANIFELM>::value>
	>
	bool fullDimElmContainsManif( MANIFELM const & manifEl, Grid & grid )
	{
		bool contained = false;

		for( INDEX_TXP iFac = 0; iFac < m_fullDimElm->num_faces(); iFac++ )
		{

			static_assert(std::is_same<decltype(m_fullDimElm), Volume *>::value);

			Face * fac = grid.get_face(m_fullDimElm,iFac);

			if( fac == manifEl )
			{
				contained = true;
				break;
			}
		}

		return false;
	}



};



}

}

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_ */
