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

#if 0

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
#endif


template <
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP
>
class AttachedGeneralElem
{
public:
	using PairLowEl = std::pair<LOWDIMELM,LOWDIMELM>;

	using AttGenElm = AttachedGeneralElem<MANIFELM,LOWDIMELM,INDEX_TXP>;

	// for fracture elements
	AttachedGeneralElem( MANIFELM const & manifElm,
				  PairLowEl const & lowElm
				  )
	:
		m_manifElm(manifElm), m_pairLowElm(lowElm)
	{
	};

	MANIFELM const getManifElm() const { return m_manifElm;}
//	PairLowEl const getLowElm() const { return m_lowElm; }
	PairLowEl const getPairLowElm() const { return m_pairLowElm; }

	bool const isNeighboured( AttGenElm const & attElm )
	const
	{
//		MANIFELM manifElmOther = attElm.getManifElm();
		PairLowEl lowElmOther = attElm.getPairLowElm();
//		INDEX_TXP sudoOther = attElm.getSudo();

		PairLowEl lowElmThis = this->m_pairLowElm;

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

	bool const isNeighbouredAtSpecificSide( AttGenElm const & attElm,
			  	  	  	  	  	  	  	  	LOWDIMELM const & specificLDE )
	const
	{
		PairLowEl lowElmOther = attElm.getPairLowElm();

		PairLowEl lowElmThis = this->m_pairLowElm;

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

	bool const testIfEquals( AttGenElm const & attElm )
	const
	{
		MANIFELM manifElmOther = attElm.getManifElm();
		PairLowEl lowElmOther = attElm.getPairLowElm();

		if(    manifElmOther == this->m_manifElm
			&& hasSameEdgePair( lowElmOther )
//				&& lowElmOther == this->m_pairLowElm
		)
		{
			return true;
		}

		if( manifElmOther == this->m_manifElm && ! hasSameEdgePair( lowElmOther ) )
		{
			UG_LOG("gleiches face aber andere Ecken???" << std::endl);
			UG_THROW("gleiches face aber andere Ecken???" << std::endl);
		}

		return false;
	}


protected:

	MANIFELM m_manifElm;
	PairLowEl m_pairLowElm;

	bool const hasSameEdgePair( PairLowEl const & epTwo ) const
	{
		PairLowEl const & epOne = this->m_pairLowElm;

		if(    ( epOne.first == epTwo.first && epOne.second == epTwo.second )
			||	( epOne.first == epTwo.second && epOne.second == epTwo.first )
		)
		{
			return true;
		}

		return false;
	}

};



///////////////////////////////////////////////////////////////////////////////


// TODO FIXME vertex fracture triplett
// vereinigen mit  AttachedGeneralElem !!! davon ableiten!!!
// doppelte Strukturen!!!

#if 0
// [[DEPRECATED]]
// wird abgelöst durch fortschrittlichere Klassen, bald nicht mehr nötig
template <
typename MANIFOLDTYP, // 3D: Face
typename INDEXTYP, // int oder unsinged int oder short oder unsigned short etc
typename FULLDIMTYP, // 3D: Volume
typename SENKRECHTENTYP, // 3D und 2D: ug::vector3
typename LOWDIMTYP // 3D: Edge (2D nicht benötigt)
>
class VertexFractureTripleMF
: public AttachedGeneralElem<MANIFOLDTYP,LOWDIMTYP,INDEXTYP>
{

private:
	using AttGenEl = AttachedGeneralElem<MANIFOLDTYP,LOWDIMTYP,INDEXTYP>;

public:

	using PairLowEl = std::pair<LOWDIMTYP,LOWDIMTYP>;

	VertexFractureTripleMF( MANIFOLDTYP const & manifElm, INDEXTYP sudo,
							FULLDIMTYP const & fullElm,
							SENKRECHTENTYP const & normal,
							PairLowEl const & pairLowElm )
	: //m_manifElm(manifElm),
	  AttGenEl(manifElm,pairLowElm),
	  m_sudo(sudo), m_fullElm(fullElm),
	  m_normal(normal), m_newNormal(normal)
//		,
//	  m_pairLowElm(pairLowElm)
	{
	};

//	MANIFOLDTYP const getManifElm() const { return m_manifElm; }

	INDEXTYP const getSudoElm() const { return m_sudo; }

	FULLDIMTYP const getFullElm() const { return m_fullElm; }

	SENKRECHTENTYP const getNormal() const { return m_normal; }

	// TODO FIXME unklar, ob neue Normale irgendwo gebraucht wird
	// falls notwendig in die Fracture Klasse einführen,
	// die diese Klasse hier mittelfristig ablösen soll vollständig
	void setNewNormal( SENKRECHTENTYP const & chNorml ) { m_newNormal = chNorml; }

	SENKRECHTENTYP const getNewNormal() const { return m_newNormal; }

//	PairLowEl const getPairLowElm() const { return m_pairLowElm; }

private:

//	MANIFOLDTYP m_manifElm;
	INDEXTYP m_sudo;
	FULLDIMTYP m_fullElm;
	SENKRECHTENTYP  m_normal;
	SENKRECHTENTYP  m_newNormal;
//	PairLowEl m_pairLowElm;

	VertexFractureTripleMF()
	{};

};
#endif

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



// info for a vertex: face and attached edges, for fractures only
//  derived from the similar class without fracture property!
template <
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP,
typename NORMAL_VEC
>
class AttachedFractElem
: public AttachedGeneralElem<MANIFELM,LOWDIMELM,INDEX_TXP>
// TODO FIXME derive from AttachedGeneralElem
{
public:
	using PairLowEl = std::pair<LOWDIMELM,LOWDIMELM>;

	using AttFractElm = AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>;

	using AttGenElm = AttachedGeneralElem<MANIFELM,LOWDIMELM,INDEX_TXP>;

	// for fracture elements
	AttachedFractElem( MANIFELM const & manifElm,
				  PairLowEl & lowElm,
				  INDEX_TXP sudo,
				  NORMAL_VEC const & normalVec )
	:
		AttGenElm(manifElm,lowElm),
		//m_manifElm(manifElm), m_lowElm(lowElm),
		m_sudo(sudo),
		m_normalVec(normalVec)
	{
	};


//	MANIFELM const getManifElm() const { return m_manifElm;}
//	PairLowEl const getLowElm() const { return m_lowElm; }
	INDEX_TXP const getSudo() const { return m_sudo; };

	NORMAL_VEC const getNormalVec() const { return m_normalVec; }

	bool const testIfEquals( AttFractElm const & attElm )
	const
	{
		bool geomEqu = AttGenElm::testIfEquals(attElm);

//		MANIFELM manifElmOther = attElm.getManifElm();
//		PairLowEl lowElmOther = attElm.getLowElm();
		INDEX_TXP sudoOther = attElm.getSudo();

//		if(    manifElmOther == this->m_manifElm
//			&& lowElmOther == this->m_lowElm
//			&& sudoOther == this->m_sudo
//		)
//		{
//			return true;
//		}

		if( geomEqu && sudoOther == this->m_sudo )
		{
			return true;
		}

		return false;
	}


//	bool const testIfEquals( AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP> const & attElm )
//	const
//	{
//		MANIFELM manifElmOther = attElm.getManifElm();
//		PairLowEl lowElmOther = attElm.getLowElm();
//		INDEX_TXP sudoOther = attElm.getSudo();
//
//		if(    manifElmOther == this->m_manifElm
//			&& lowElmOther == this->m_lowElm
//			&& sudoOther == this->m_sudo
//		)
//		{
//			return true;
//		}
//
//		return false;
//	}

//	bool const isNeighboured(  AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP> const & attElm )
//	const
//	{
////		MANIFELM manifElmOther = attElm.getManifElm();
//		PairLowEl lowElmOther = attElm.getLowElm();
////		INDEX_TXP sudoOther = attElm.getSudo();
//
//		PairLowEl lowElmThis = this->m_lowElm;
//
//		std::vector<bool> test;
//
//		test.push_back( lowElmOther.first  == lowElmThis.first );
//		test.push_back( lowElmOther.second  == lowElmThis.first );
//		test.push_back( lowElmOther.first  == lowElmThis.second );
//		test.push_back( lowElmOther.second  == lowElmThis.second );
//
//		INDEX_TXP countCorr = 0;
//
//		for( auto const t : test )
//		{
//			if( t )
//				countCorr++;
//		}
//
//		if( countCorr == 1 )
//			return true;
//
//		if( countCorr > 1 )
//			UG_THROW("zu viele gleiche Ecken " << std::endl);
//
//		return false;
//	}
//
//	bool const isNeighbouredAtSpecificSide( AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP> const & attElm,
//			  	  	  	  	  	  	  	  	LOWDIMELM const & specificLDE )
//	const
//	{
//		PairLowEl lowElmOther = attElm.getLowElm();
//
//		PairLowEl lowElmThis = this->m_lowElm;
//
//		// test if the specific element is part of at least
//		// one of the faces
//
//		bool otherFirst = ( lowElmOther.first == specificLDE );
//		bool otherSecond = ( lowElmOther.second == specificLDE );
//
//		bool thisFirst = ( lowElmThis.first == specificLDE );
//		bool thisSecond = ( lowElmThis.second == specificLDE );
//
//		bool isPartOfThisFace = ( thisFirst || thisSecond );
//		bool isPartOfOtherFace = ( otherFirst || otherSecond );
//
//		if( ! isPartOfOtherFace || ! isPartOfThisFace )
//		{
//			UG_LOG("not part of one of the faces " << std::endl);
//			return false;
//		}
//
//		if( otherFirst && thisFirst )
//		{
//			if( lowElmOther.first  == lowElmThis.first )
//				return true;
//		}
//		else if( otherFirst && thisSecond )
//		{
//			if( lowElmOther.first  == lowElmThis.second )
//				return true;
//		}
//		else if( otherSecond && thisFirst )
//		{
//			if( lowElmOther.second  == lowElmThis.first )
//				return true;
//		}
//		else if( otherSecond && thisSecond )
//		{
//			if( lowElmOther.second  == lowElmThis.second )
//				return true;
//		}
//
//		return false;
//	}

private:
//	MANIFELM m_manifElm;
//	PairLowEl m_lowElm;
	INDEX_TXP m_sudo;
	NORMAL_VEC m_normalVec;


};

//////////////////////////////////////////////////////////////////////////////

// a quasi exact double, but only used for boundary faces, to avoid mismatch with frac faces
template <
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP,
typename NORMAL_VEC
>
class AttachedBoundryElem
: public AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>
{
public:
	using PairLowEl = std::pair<LOWDIMELM,LOWDIMELM>;

	using AttBndryElm = AttachedBoundryElem<MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>;

	using AttFractElm = AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>;


	// for boundary elements
	AttachedBoundryElem( MANIFELM const & manifElm,
				  PairLowEl & lowElm,
				  INDEX_TXP sudo,
				  NORMAL_VEC const & normalVec )
	:
		AttFractElm( manifElm, lowElm, sudo, normalVec)
	{
	};
};

////////////////////////////////////////////////////////////////////////////

// class to help count and store a bool and a number of templete type
// comparable to std::pair<bool,int> but more dedicated to the specific aim
// TODO FIXME adapt for 3D case, figure out if inner end, and number of fracs sourrounding
// CAUTION is also used for edges, but still uses
// vertex as indicator - name should be made more flexible

// die meisten Funktionen in dieser Klasse:
// DEPRECATED, to be replaced in near future everywhere, not really useful any more
// due to the Stasi algorithm
// [[deprecated]] ab C++14, leider nicht passend zur Konvention C++11
// die Sudo-Liste wollen wir aber lassen
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

	// DEPRECATED, kann entfernt werden vermutlich, CHECK, TODO FIXME
	enum VrtxFracStatus { noFracSuDoAtt = 0,
						  oneFracSuDoAtt = 1,
						  twoFracSuDoAtt = 2,
						  threeFracSuDoAtt = 3 };

	VertexFracturePropertiesVol()
	: m_isBndFracVertex(false), m_numberCountedFracsInVertex(0),
	  m_status( noFracSuDoAtt ),
	  m_sudoList( std::vector<T>() ) //,
//	  m_sudosClosed(VecPairTB()),
//	  m_vecAttElem(std::vector<ATT_ELEM>())
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

//	bool setIsAClosedFracture( T sudoNow, bool isClosedNow )
//	{
//
//		T alreadyKnownMult = 0;
//
//		for( auto & suSu : m_sudosClosed )
//		{
//			T & sudoVal = suSu.first;
//			bool & isClosedVal = suSu.second;
//
//			if( sudoVal == sudoNow )
//			{
//				alreadyKnownMult++;
//
//				UG_LOG("Reassign sudo surround " << std::endl);
//
//				if( isClosedVal != isClosedNow )
//					UG_THROW("change property sudo surrounded, why?" << std::endl);
//
//				isClosedVal = isClosedNow;
//			}
//		}
//
//		if( alreadyKnownMult == 0 )
//		{
//			pairTB infoSudoSurr( sudoNow, isClosedNow );
//
//			m_sudosClosed.push_back( infoSudoSurr );
//
//		}
//		else if( alreadyKnownMult > 1 )
//		{
//			UG_THROW("zu oft bekannt " << std::endl);
//			return false;
//		}
//
//		// check if now correct
//
//		T testKnownFine = 0;
//
//		for( auto const & suSu : m_sudosClosed )
//		{
//			T & sudoVal = suSu.first;
//			bool & isClosedVal = suSu.second;
//
//			if( sudoVal == sudoNow )
//			{
//				testKnownFine++;
//
//				if( isClosedVal != isClosedNow )
//				{
//					UG_THROW("NOT set property sudo surrounded, why?" << std::endl);
//					return false;
//				}
//
//			}
//		}
//
//		if( testKnownFine == 0 || testKnownFine > 1 )
//		{
//			UG_THROW("immer noch nicht bekannt?" << std::endl);
//			return false;
//
//		}
//
//		return true;
//	}
//
//	bool getIsAClosedFracture( T sudoNow )
//	{
//		T foundMultpl = 0;
//
//		bool isClosedReturn = false;
//
//		for( auto const & suSu : m_sudosClosed )
//		{
//			T const & sudoVal = suSu.first;
//			bool const & isClosedVal = suSu.second;
//
//			if( sudoVal == sudoNow )
//			{
//				foundMultpl++;
//				isClosedReturn = isClosedVal;
//			}
//		}
//
//		if( foundMultpl != 1 )
//		{
//			UG_THROW("not known status closed or not sudo" << std::endl);
//			return false;
//		}
//
//		return isClosedReturn;
//	}
//
//	bool setInfoAllFractureSudosIfClosed( VecPairTB const & sudosClosed )
//	{
//		m_sudosClosed = sudosClosed;
//
//		return true;
//	}
//
//	VecPairTB const getInfoAllFracSudosIfClosed() const
//	{
//		return m_sudosClosed;
//	}
//
//	// if all open or closed
//	template<bool B>
//	bool const getInfoAllFracturesSameClosedState() const
//	{
//		bool allFracsSame = true;
//
//		for( auto const & suSu : m_sudosClosed )
//		{
//			//T const & sudoVal = suSu.first;
//			bool const & isClosedVal = suSu.second;
//
//			if( isClosedVal != B )
//				allFracsSame = false;
//		}
//
//		return allFracsSame;
//	}

	// DEPRECATED; REMOVE
//	bool addAttachedFractElem( ATT_ELEM const & attElem )
//	{
//		bool alreadyKnown = false;
//
//		for( auto const & aE : m_vecAttElem )
//		{
//			if( aE.testIfEquals(attElem) )
//				alreadyKnown = true;
//		}
//
//		if( ! alreadyKnown )
//			m_vecAttElem.push_back(attElem);
//
//		// returns true if ads it, false if no need as known
//		return ! alreadyKnown;
//	}

//	std::vector<ATT_ELEM> const & getAllAttachedFractElems()
//	const
//	{
//		return m_vecAttElem;
//	}

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

//	VecPairTB m_sudosClosed;

	bool setSudoList( std::vector<T> const & sudoList )
	{
		m_sudoList = sudoList;

		return true;
	}

//	std::vector<ATT_ELEM> m_vecAttElem;

};


// intention, explained for volume:
template<
typename FULLDIM_ELEM,
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP,
typename NORMAL_VEC
>
class AttachedFullDimElemInfo
{

public:

	using AttachedFractManifElemInfo = AttachedFractElem<MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>;
	using AttachedGenerManifElemInfo = AttachedGeneralElem<MANIFELM,LOWDIMELM,INDEX_TXP>;
	using AttachedBndryManifElemInfo = AttachedBoundryElem<MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>;

	using VecAttachedFractManifElemInfo = std::vector<AttachedFractManifElemInfo>;
	using VecAttachedGenerManifElemInfo = std::vector<AttachedGenerManifElemInfo>;
	using VecAttachedBndryManifElemInfo = std::vector<AttachedBndryManifElemInfo>;

	using AttFullDimElmInfo = AttachedFullDimElemInfo<FULLDIM_ELEM,MANIFELM,LOWDIMELM,INDEX_TXP,NORMAL_VEC>;

	AttachedFullDimElemInfo( FULLDIM_ELEM const & fullDimElm )
	: m_fullDimElm(fullDimElm),
	  m_elementMarked(false),
	  m_vecFractManifElm(VecAttachedFractManifElemInfo()),
	  m_vecUnclosedFractManifElm(VecAttachedFractManifElemInfo()),
//	  m_vecFractManifElmTouchInfo(VecAttFractManifElmTouchInf()),
//	  m_allSidesTouched(false),
	  m_vecGenerManifElm(VecAttachedGenerManifElemInfo()),
	  m_vecBndryManifElm(VecAttachedBndryManifElemInfo())
//	  m_vecGenerManifElmTouchInfo(VecAttFractManifElmTouchInf())
//	  m_vecInnerSegmentManifElm(VecAttachedGenerManifElemInfo())
	{
	}

	FULLDIM_ELEM const getFulldimElem() const
	{
		return m_fullDimElm;
	}

	bool const hasSameFulldimElem( AttFullDimElmInfo const & otherFullDimElmInf )
	const
	{
		FULLDIM_ELEM otherFullDimElm = otherFullDimElmInf.getFulldimElem();

		return ( otherFullDimElm == m_fullDimElm );
	}

	bool const isMarked() const { return m_elementMarked; }

	void markIt() { m_elementMarked = true; } // to allow in loop to be start element

	bool const hasFracture() const { return ( m_vecFractManifElm.size() > 0 ); }

	bool const hasUnclosedFracture() const { return ( m_vecUnclosedFractManifElm.size() > 0 ); }

	bool addFractManifElem( AttachedFractManifElemInfo const & manifFractElm, Grid & grid )
	{
		return addManifElem( manifFractElm, m_vecFractManifElm, grid );
	}

	bool addGenerManifElem( AttachedGenerManifElemInfo const & manifGenerElm, Grid & grid )
	{
//		static_assert(std::is_same<AttachedGenerManifElemInfo, decltype( manifGenerElm ) >::value);

		return addManifElem( manifGenerElm, m_vecGenerManifElm, grid );
	}

	bool addBndryManifElem( AttachedBndryManifElemInfo const & manifBndryElm, Grid & grid )
	{
		return addManifElem( manifBndryElm, m_vecBndryManifElm, grid );
	}

	// necessary to avoid stupid casting from derived class AttachedFractManifElemInfo
	// else, addGenerManifElem would also eat objects of derived class
	// however not it only accepts explicit base class objects
	template <typename NOGEN>
	bool addGenerManifElem( NOGEN const & noGener, Grid & grid )
	= delete;

	template <typename NOGEN>
	bool addBndryManifElem( NOGEN const & noGener, Grid & grid )
	= delete;


//	// will form new "surface" of next inner round in a segment
//	bool addInnerSegmentManifElem( AttachedGenerManifElemInfo const & manifInnerSegmentElm, Grid & grid )
//	{
//		return addManifElem( manifInnerSegmentElm, m_vecInnerSegmentManifElm, grid );
//	}
//
//	template <typename NOGEN>
//	bool addInnerSegmentElem( NOGEN const & noGener, Grid & grid )
//	= delete;



//	bool addFractureManifElem( AttachedFractManifElemInfo const & manifElm, Grid & grid )
//	{
//		// Caution: first test if manifold elem is at all part of the fulldim elem manifols
//		// if not, return false directly
//
//		if( ! fullDimElmContainsManif( manifElm.getManifElm(), grid ) )
//			return false;
//
//		// if manif elem is in principle part of the fulldim elem manifolds,
//		// then we need to check if it is already integrated
//
//		bool hasElemAlready = false;
//
//
//		for( auto const & me : m_vecFractManifElm )
//		{
//			if( manifElm.testIfEquals(me) )
//			{
//				hasElemAlready = true;
//				break;
//			}
//		}
//
//		if( ! hasElemAlready )
//		{
//			m_vecFractManifElm.push_back( manifElm );
//
////			AttFractManifElmTouchInf pamei( manifElm, false );
////
////			m_vecFractManifElmTouchInfo.push_back(pamei);
//		}
//
//		return ! hasElemAlready;
//	}

	VecAttachedFractManifElemInfo const getVecFractManifElem() const
	{
		return m_vecFractManifElm;
	}

	VecAttachedFractManifElemInfo const getVecUnclosedFractManifElem() const
	{
		return m_vecUnclosedFractManifElm;
	}


	VecAttachedGenerManifElemInfo const getVecGenerManifElem() const
	{
		return m_vecGenerManifElm;
	}

	VecAttachedBndryManifElemInfo const getVecBndryManifElem() const
	{
		return m_vecBndryManifElm;
	}


	bool const searchGenerManifElem( AttachedGenerManifElemInfo const & manifGenerElemOther, bool eraseFound = true )
	{
		bool found = searchManifElem( manifGenerElemOther, m_vecGenerManifElm, eraseFound );

		if( found && eraseFound )
		{
			m_elementMarked = true;
		}

		return found;
	}

	bool const testFullDimElmNeighbour( AttFullDimElmInfo const & attFullDimElmInfOther, bool eraseFoundManif = true )
	{
		VecAttachedGenerManifElemInfo const & vecGenerManifElmOther = attFullDimElmInfOther.getVecGenerManifElem();

		bool manifNeighbored = false;

		for( AttachedGenerManifElemInfo const & generManifElemOther : vecGenerManifElmOther )
		{
			if(	searchManifElem( generManifElemOther, m_vecGenerManifElm, eraseFoundManif ) )
				manifNeighbored = true;

		}

		if( manifNeighbored && eraseFoundManif )
		{
			m_elementMarked = true;
		}

		return manifNeighbored;
	}


	template <typename NOGEN>
	bool searchGenerManifElem( NOGEN const & manifGenerElemOther, bool eraseFound ) = delete;

//	bool const searchFractManifElem( AttachedFractManifElemInfo const & manifFractElemOther, bool eraseFound = true )
	bool const searchFractManifElem( AttachedFractManifElemInfo const & manifFractElemOther, bool shiftToUnclosedFracts = true )
	{
		bool found = searchManifElem( manifFractElemOther, m_vecFractManifElm, shiftToUnclosedFracts );

		if( found && shiftToUnclosedFracts )
		{
			// for the case that a fracture is not closed at a vertex, shift the in principle
			// fracture vertex to the gerneral vertices, as useless for segmente construction
			// and useless for expansion

			m_vecUnclosedFractManifElm.push_back(manifFractElemOther);


//			MANIFELM const & manifel = manifFractElemOther.getManifElm();
//			typename AttachedGenerManifElemInfo::PairLowEl const & pairLowEl = manifFractElemOther.getPairLowElm();
//
//			AttachedGenerManifElemInfo agmei( manifel, pairLowEl );

//			m_vecUnclosedFractManifElm.push_back(agmei);
		}

		return found;

	}

	template <typename NOGEN>
	bool const searchFractManifElem( NOGEN const & manifFractElemOther, bool shiftToGeneral ) = delete;

	bool const searchBndryManifElem( AttachedBndryManifElemInfo const & manifBndryElemOther )
	{
		return searchManifElem( manifBndryElemOther, m_vecBndryManifElm, false );
	}

//	bool const searchInnerSegmentManifElem( AttachedGenerManifElemInfo const & manifInnerSegmElemOther, bool eraseFound = true )
//	{
//		return searchManifElem( manifInnerSegmElemOther, m_vecInnerSegmentManifElm, eraseFound );
//	}

//	bool const tryToTouchManifElem( AttachedFractManifElemInfo const & manifElemOther ) const
//	{
////		bool managed2Touch = false;
//
//		for( typename VecAttachedFractManifElemInfo::iterator afeiIter  = m_vecFractManifElm.begin();
//													 afeiIter != m_vecFractManifElm.end();
//													 afeiIter++
//		)
//		{
//			AttachedFractManifElemInfo & manifElmTest = *afeiIter;
//
//			if( manifElemOther.testIfEquals(manifElmTest) )
//			{
//				managed2Touch = true;
//
//				m_vecFractManifElm.erase(afeiIter);
//
//				return true;
//			}
//		}
//
//		return false;
//
////		return managed2Touch;
//
////		for( auto & ameti : m_vecFractManifElmTouchInfo )
////		{
////			AttachedFractManifElemInfo & manifElmTest = ameti.first;
////			bool & alreadyTouched = ameti.second;
////
////			if( ! alreadyTouched )
////			{
////				if( manifElemOther.testIfEquals( manifElmTest ) )
////				{
////					alreadyTouched = true;
////				}
////			}
////
////			if( alreadyTouched )
////			{
////				managed2Touch = true;
////				break;
////			}
////		}
//
//		return managed2Touch;
//	}

//	VecAttachedFractManifElemInfo const getAlreadyTouchedManifElems() const
//	{
//		VecAttachedFractManifElemInfo alreadyTouchedManifElms;
//
//		for( const auto & ameti : m_vecFractManifElmTouchInfo )
//		{
//			if( ameti.second )
//				alreadyTouchedManifElms.push_back( ameti );
//		}
//
//		return alreadyTouchedManifElms;
//	}

//	VecAttachedFractManifElemInfo const getSoFarUnTouchedManifElems() const
//	{
//		VecAttachedFractManifElemInfo unTouchedManifElms;
//
//		for( const auto & ameti : m_vecFractManifElmTouchInfo )
//		{
//			if( ! ameti.second )
//				unTouchedManifElms.push_back( ameti );
//		}
//
//		return unTouchedManifElms;
//	}

//	bool const testIfAllSidesTouched() const
//	{
//		if( m_allSidesTouched )
//			return true;
//
//		bool allSidesTouched = true;
//
//		for( const auto & ameti : m_vecFractManifElmTouchInfo )
//		{
//			if( ! ameti.second )
//			{
//				allSidesTouched = false;
//			}
//		}
//
//		m_allSidesTouched = allSidesTouched;
//
//		return m_allSidesTouched;
//	}

private:

//	bool m_allSidesTouched;

	FULLDIM_ELEM m_fullDimElm;

	bool m_elementMarked;

	VecAttachedFractManifElemInfo m_vecFractManifElm;

	VecAttachedFractManifElemInfo m_vecUnclosedFractManifElm;


//	using AttFractManifElmTouchInf = std::pair<AttachedFractManifElemInfo,bool>;
//	using VecAttFractManifElmTouchInf = std::vector<AttFractManifElmTouchInf>;
//
//	VecAttFractManifElmTouchInf m_vecFractManifElmTouchInfo;

	VecAttachedGenerManifElemInfo m_vecGenerManifElm;

	VecAttachedBndryManifElemInfo m_vecBndryManifElm;

//	VecAttachedGenerManifElemInfo m_vecInnerSegmentManifElm;

//	using AttGenerManifElmTouchInf = std::pair<AttachedGenerManifElemInfo,bool>;
//	using VecAttGenerManifElmTouchInf = std::vector<AttGenerManifElmTouchInf>;
//
//	VecAttFractManifElmTouchInf m_vecGenerManifElmTouchInfo;


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
//		bool contained = false;

		for( INDEX_TXP iFac = 0; iFac < m_fullDimElm->num_faces(); iFac++ )
		{

//			static_assert(std::is_same<decltype(m_fullDimElm), Volume *>::value);

			Face * fac = grid.get_face(m_fullDimElm,iFac);

			if( fac == manifEl )
			{
				return true;
//				contained = true;
//				break;
			}
		}

//		return contained;
		return false;
	}


	template <typename ATT_MANIF_ELM_INF >
	bool const searchManifElem( ATT_MANIF_ELM_INF const & manifElemOther,
								std::vector<ATT_MANIF_ELM_INF> & memVecManifElem,
								bool eraseFound = true )
	const
	{

		for( typename std::vector<ATT_MANIF_ELM_INF>::iterator afeiIter  = memVecManifElem.begin();
															   afeiIter != memVecManifElem.end();
															   afeiIter++
		)
		{
			ATT_MANIF_ELM_INF & manifElmTest = *afeiIter;

			if( manifElemOther.testIfEquals(manifElmTest) )
			{

				if( eraseFound )
					memVecManifElem.erase(afeiIter);

				return true;
			}
		}

		return false;
	}



	template <typename ATT_MANIF_ELM_INFO >
	bool addManifElem( ATT_MANIF_ELM_INFO const & manifElm,
					   std::vector<ATT_MANIF_ELM_INFO> & memVecManifElm,
					   Grid & grid )
	{
		// Caution: first test if manifold elem is at all part of the fulldim elem manifols
		// if not, return false directly

		if( ! fullDimElmContainsManif( manifElm.getManifElm(), grid ) )
			return false;

		// if manif elem is in principle part of the fulldim elem manifolds,
		// then we need to check if it is already integrated

		for( auto const & me : memVecManifElm )
		{
			if( manifElm.testIfEquals(me) )
			{
				return false;
			}
		}

		// not contained so far, but part of the manifolds of the fulldim elem
		memVecManifElm.push_back(manifElm);

		return true;

	}


};

//////////////////////////////////////////////////////////////////


// Ebenentyp: a x1 + b x2 + c x3 = rhs, normal * ( vecX - baseVect ) = 0
template<typename VECTOR_TYP>
class ManifoldDescriptor
{
public:

	enum ManifoldType { isFracture, isBoundary, isArtificial };

	template<typename = std::enable_if<std::is_same<VECTOR_TYP,vector3>::value>>
//			ManifoldDescriptor( int sudo = -1, number scaleShiftNormal = 0 )
	ManifoldDescriptor()
	: m_normalVect(vector3()),
	  m_baseVect(vector3()),
	  m_rhs(0),
	  m_scaleShiftNormal(0),
	  m_dim(3),
	  m_sudo(-1),
	  m_manifTyp( isArtificial )
	{
	}


	template<typename = std::enable_if<std::is_same<VECTOR_TYP,vector3>::value>>
	ManifoldDescriptor( VECTOR_TYP const & normalVect,
						VECTOR_TYP const & baseVect,
						int sudo = -1,
						ManifoldType manifTyp = isArtificial,
						number scaleShiftNormal = 0
						)
	: m_normalVect(normalVect),
	  m_baseVect(baseVect),
	  m_scaleShiftNormal(scaleShiftNormal),
	  m_dim(3),
	  m_sudo(sudo),
	  m_manifTyp( manifTyp )
	{
		m_rhs = 0;

		UG_LOG("Ebenenkoordinatenform ");

		for( int i = 0; i < m_dim; i++ )
		{
			m_rhs += normalVect[i]*baseVect[i];

			UG_LOG( " + " << normalVect[i] << " x_" << i << " " );
		}

		UG_LOG(" = " << m_rhs << std::endl);

	}

	VECTOR_TYP const & spuckNormalVector() const { return m_normalVect; }
	VECTOR_TYP const & spuckBaseVector() const { return m_baseVect; }

	int const spuckSudo() const { return m_sudo; }

	ManifoldType const spuckManifTyp() const { return m_manifTyp; }

	number const spuckScaleShiftNormal() const { return m_scaleShiftNormal; }

	void schluckSudo( int sudo ) { m_sudo = sudo; }

	void schluckManifTyp( ManifoldType manifTyp )  { m_manifTyp = manifTyp; }

	void schluckScaleShiftNormal( number scaleShiftNormal )
	{
		m_scaleShiftNormal = scaleShiftNormal;
	}

	number const & spuckRHS() const { return m_rhs; }

	template<typename = std::enable_if<    std::is_same<VECTOR_TYP,vector3>::value
//										|| std::is_same<VECTOR_TYP,vector2>::value>
									>
	>
	bool spuckPlaneShiftedAlong( VECTOR_TYP const & shiftVec, ManifoldDescriptor & manifoldDescr )
	{
		VECTOR_TYP shiftedBaseVect;
//		number rhsShiftedPlane = 0;

		VecAdd( shiftedBaseVect, m_baseVect, shiftVec );

		UG_LOG("Ebenenkoordinatenform Shifted Plane " << std::endl);

//		ManifoldDescriptor( VECTOR_TYP const & normalVect,
//							VECTOR_TYP const & baseVect,
//							int sudo = -1,
//							ManifoldType = isArtificial,
//							number scaleShiftNormal = 0
//							)

		manifoldDescr = ManifoldDescriptor( m_normalVect, shiftedBaseVect, m_sudo, m_manifTyp, 0 );

		return true;

	}

	template<typename = std::enable_if<std::is_same<VECTOR_TYP,vector3>::value >>
	bool spuckPlaneShifted( ManifoldDescriptor & manifoldDescr )
	{

		UG_LOG("Ebenenkoordinatenform Shifted Plane " << std::endl);

		manifoldDescr = ManifoldDescriptor( m_normalVect, spuckShiftedBaseVect(), m_sudo, m_manifTyp, 0 );

		return true;

	}


	template<typename = std::enable_if< std::is_same<VECTOR_TYP,vector3>::value>>
	VECTOR_TYP spuckShiftedBaseVect()
	{
		VECTOR_TYP shiftVec;

		VecScale(shiftVec, m_normalVect, m_scaleShiftNormal );

		VECTOR_TYP shiftedBaseVec;

		VecAdd( shiftedBaseVec, m_baseVect, shiftVec );

		return shiftedBaseVec;
	}

private:

	VECTOR_TYP m_normalVect;
	VECTOR_TYP m_baseVect;
	number m_rhs;
	number m_scaleShiftNormal;

	// could be defined as static const variable, but then depending on the template parameter
	// might be an idea for future, but no real use, but would be "nice" from the meta programmig point of view
	int m_dim;
	int m_sudo;
	ManifoldType m_manifTyp;
};


//////////////////////////////////////////////////////////////////


template <
typename FULLDIM_ELEM,
typename MANIFELM,
typename LOWDIMELM,
typename INDEX_TXP,
typename VECTOR_TYP,
typename VRTXTYP
>
class SegmentSides
{
public:

	enum VrtxFracStatus { noFracSuDoAtt = 0,
						  oneFracSuDoAtt = 1,
						  twoFracSuDoAtt = 2,
						  threeFracSuDoAtt = 3 };

	using AttFractElm = AttachedFractElem<MANIFELM, LOWDIMELM, INDEX_TXP, VECTOR_TYP>;
	using AttBndryElm = AttachedBoundryElem<MANIFELM, LOWDIMELM, INDEX_TXP, VECTOR_TYP>;

	using VecAttFractElm = std::vector<AttFractElm>;
	using VecAttBndryElm = std::vector<AttBndryElm>;

	using PairSudoNormlV = std::pair<INDEX_TXP,VECTOR_TYP>;
	using VecPairSudoNormlV = std::vector<PairSudoNormlV>;
	using ManifDescr = ManifoldDescriptor<VECTOR_TYP>;
	using VecManifDescr = std::vector<ManifDescr>;

	// TODO FIXME das soll gleich durch den Manifold Descriptor ersetzt werden
	// oder eine Basisklasse von ihm, die nicht so viele Infos enthält
	// aber mindestens NormalenVektor, Sudo, und ob Boundary oder Fracture
	// kann auch vielleicht einfach eine Klasse sein, die noch einen Parameter enthält
	// der sich abfragen lässt, auch als Template vielleicht, true false, fracture or not
	// also was wie template < index, normal, bool > pairsudonormlbla, oder sowas.....
	// oder man kann einen Parameter setzen für diese Klasse, die extern definiert wird......
	// bool isFracture true false.....

	SegmentSides( VRTXTYP const & vrt, bool isBndry = false )
	: m_vrt(vrt),
	  m_vecAttFractElms(VecAttFractElm()),
	  m_vecAttUnclosedFractElms(VecAttFractElm()),
	  m_vecAttBndryElms(VecAttBndryElm()),
	  m_vecFractSudosNormlV(VecPairSudoNormlV()),
	  m_vecBndrySudosNormlV(VecPairSudoNormlV()),
	  m_isBoundary(isBndry),
	  m_averaged(false),
	  m_contribFulldimElm(std::vector<FULLDIM_ELEM>())
	{};

//	template<typename = std::enable_if< std::is_pointer<VRTXTYP>::value>>
//	SegmentSides()
//	: m_vrt(nullptr),
//	  m_vecAttFractElms(VecAttFractElm()),
//	  m_vecAttUnclosedFractElms(VecAttFractElm()),
//	  m_vecAttBndryElms(VecAttBndryElm()),
//	  m_vecFractSudosNormlV(VecPairSudoNormlV()),
//	  m_vecBndrySudosNormlV(VecPairSudoNormlV()),
//	  m_isBoundary(false),
//	  m_averaged(false),
//	  m_contribFulldimElm(std::vector<FULLDIM_ELEM>())
//	{};


	bool const isBoundary() const { return m_isBoundary; }

	VRTXTYP const spuckVertex() const
	{
		return m_vrt;
	}

	void schluckFulldimElem( FULLDIM_ELEM const & fudielm )
	{
		m_contribFulldimElm.push_back(fudielm);
	}

	void spuckVecFulldimElem( std::vector<FULLDIM_ELEM> & fudielm ) const
	{
		fudielm = m_contribFulldimElm;
	}

	VrtxFracStatus const spuckCrossingTyp() const
	{
		VrtxFracStatus vfsFract =  static_cast<VrtxFracStatus>(m_vecFractSudosNormlV.size());

		if( m_isBoundary )
		{
			VrtxFracStatus vfsBndry = static_cast<VrtxFracStatus>(m_vecBndrySudosNormlV.size());

			return static_cast<VrtxFracStatus>( static_cast<INDEX_TXP>(vfsFract) + static_cast<INDEX_TXP>(vfsBndry) );
		}

		return vfsFract;
	}

	bool schluckVecAttFractElm( std::vector<AttFractElm> const & vecAtFracEl )
	{
		return schluckVecAttElm( vecAtFracEl, m_vecAttFractElms );
	}

	template< typename NOFRACT >
	bool schluckVecAttFractElm( std::vector<NOFRACT> const & vecAtFracEl ) = delete;

	bool schluckAttFractElm( AttFractElm const & afeNew )
	{
		return schluckAttElm( afeNew, m_vecAttFractElms );

//		if( ! isStillUnknown( afeNew, m_vecAttFractElms ) )
//		{
//			return false;
//		}
//
//		m_vecAttFractElms.push_back(afeNew);
//
//		return true;
	}

	template< typename NOFRACT >
	bool schluckAttFractElm( NOFRACT const & afeNew ) = delete;

	// soll auch in der Lage sein, die einzenlen Fracture faces wieder aus zu spucken als Liste
	// analog auch dan	ach die boundary Geschichten
	bool const spuckVecAttFractElm( std::vector<AttFractElm> & vecAttFracEl ) const
	{
		vecAttFracEl = m_vecAttFractElms;
		return true;
	}

	bool schluckVecAttBndryElm( std::vector<AttBndryElm> const & vecAtBndryEl )
	{
//		if( ! checkIfIsAtBndry() )
//			return false;

		return schluckVecAttElm( vecAtBndryEl, m_vecAttBndryElms );
	}

	bool schluckAttBndryElm( AttBndryElm const & afeNew )
	{
//		if( ! checkIfIsAtBndry() )
//			return false;

		return schluckAttElm( afeNew, m_vecAttBndryElms );

//		if( ! isStillUnknown( afeNew, m_vecAttBndryElms ) )
//		{
//			return false;
//		}
//
//		m_vecAttBndryElms.push_back(afeNew);
//
//		return true;
	}

	bool spuckVecAttBndryElm( std::vector<AttBndryElm> & vecAtBndryEl )
	{
		vecAtBndryEl = m_vecAttBndryElms;
		return true;
	}


	bool schluckVecAttUnclosedFractElm( std::vector<AttFractElm> const & vecAtFracEl )
	{
		return schluckVecAttElm( vecAtFracEl, m_vecAttUnclosedFractElms, true );
	}

	template< typename NOFRACT >
	bool schluckVecAttUnclosedFractElm( std::vector<NOFRACT> const & vecAtFracEl ) = delete;

	bool schluckAttUnclosedFractElm( AttFractElm const & afeNew )
	{
		return schluckAttElm( afeNew, m_vecAttUnclosedFractElms, true );
	}

	template< typename NOFRACT >
	bool schluckAttUnclosedFractElm( NOFRACT const & afeNew ) = delete;

	bool const spuckVecAttUnclosedFractElm( std::vector<AttFractElm> & vecAttFracEl ) const
	{
		vecAttFracEl = m_vecAttUnclosedFractElms;
		return true;
	}

	bool const hasUnclosedFaces() const
	{
		return ( m_vecAttUnclosedFractElms.size() > 0 );
	}


	bool averageAll()
	{
		if( m_isBoundary )
		{
			if( ! averageBndryNormals() )
				return false;
		}

		if( ! averageFractNormals() )
			return false;

		m_averaged = true;

		return m_averaged;
	}

	bool const spuckFractSudoNormls( VecPairSudoNormlV & vecFractSudosNormlV ) const
	{
		if( ! m_averaged )
		{
			UG_LOG("please average " << std::endl);
			UG_THROW("please average " << std::endl);
			return false;
		}

		vecFractSudosNormlV = m_vecFractSudosNormlV;

		return true;
	}

	template<typename = std::enable_if< std::is_same<VECTOR_TYP,vector3>::value>>
	bool const spuckFractManifDescr( VecManifDescr & vecManifDesc,
									 Grid::VertexAttachmentAccessor<APosition> const & aaPos,
									 bool clearDescVec = true
	) const
	{
		return spuckManifDescr<ManifDescr::ManifoldType::isFracture>( vecManifDesc, aaPos, m_vecFractSudosNormlV, clearDescVec );
//		return spuckManifDescr<0>( vecManifDesc, aaPos, m_vecFractSudosNormlV );
	}


//	bool spuckFractManifDescr( VecManifDescr & vecManifDesc, Grid::VertexAttachmentAccessor<APosition> const & aaPos )
//	{
//		if( ! m_averaged )
//		{
//			UG_LOG("please average " << std::endl);
//			UG_THROW("please average " << std::endl);
//			return false;
//		}
//
//		VecPairSudoNormlV const & vecFractSudosNormlV = m_vecFractSudosNormlV;
//
//		vecManifDesc.clear();
//
//		for( PairSudoNormlV const & psn : vecFractSudosNormlV )
//		{
//			VECTOR_TYP posVrt = aaPos[m_vrt];
//
//			int sudo = psn.first;
//			VECTOR_TYP normlVec = psn.second;
//
////			ManifoldDescriptor( VECTOR_TYP const & normalVect,
////								VECTOR_TYP const & baseVect,
////								int sudo = -1,
////								ManifoldType = isArtificial,
////								number scaleShiftNormal = 0
////								)
////
//
//			ManifDescr manifDesc( normlVec, posVrt, sudo, ManifDescr::ManifoldType::isFracture );
//
//			vecManifDesc.push_back( manifDesc );
//		}
//
//		return true;
//	}

	bool spuckBndrySudoNormls( VecPairSudoNormlV & vecBndrySudosNormlV )
	{
//		if( ! checkIfIsAtBndry() )
//			return false;

		if( ! m_averaged )
		{
			UG_LOG("please average " << std::endl);
			UG_THROW("please average " << std::endl);
			return false;
		}

		vecBndrySudosNormlV = m_vecBndrySudosNormlV;

		return true;
	}

	template<typename = std::enable_if< std::is_same<VECTOR_TYP,vector3>::value>>
	bool const spuckBndryManifDescr( VecManifDescr & vecManifDesc,
									 Grid::VertexAttachmentAccessor<APosition> const & aaPos,
									 bool clearDescVec = true
	) const
	{
		return spuckManifDescr<ManifDescr::ManifoldType::isBoundary>( vecManifDesc, aaPos, m_vecBndrySudosNormlV, clearDescVec );
		//		return spuckManifDescr<2>( vecManifDesc, aaPos, m_vecFractSudosNormlV );
	}


private:

	VRTXTYP m_vrt;


	VecAttFractElm m_vecAttFractElms;
	VecAttFractElm m_vecAttUnclosedFractElms;

	VecAttBndryElm m_vecAttBndryElms;



	VecPairSudoNormlV m_vecFractSudosNormlV;
	VecPairSudoNormlV m_vecBndrySudosNormlV;

	bool m_isBoundary;
	bool m_averaged;

	std::vector<FULLDIM_ELEM> m_contribFulldimElm;

	template <typename ATT_ELM, typename VEC_ATT_ELM,
			  typename = std::enable_if<std::is_same<std::vector<ATT_ELM>,VEC_ATT_ELM>::value>,
			  typename 	= std::enable_if<std::is_base_of<AttFractElm,ATT_ELM>::value>
			>
	bool isStillUnknown( ATT_ELM const & afeNew, VEC_ATT_ELM const & vecAttELm, bool acceptUnknowns = false )
	{
		for( ATT_ELM const & afeAlt : vecAttELm )
		{
			if( afeAlt.testIfEquals(afeNew) )
			{
				UG_LOG("Strange, already known?" << std::endl);
				if( ! acceptUnknowns )
					UG_THROW("Strange, already known?" << std::endl);
				return false;
			}
		}

		return true;

	}

	template <typename ATT_ELM,
			  typename = std::enable_if<std::is_base_of<AttFractElm,ATT_ELM>::value>
			>
	bool extractSudoList( std::vector<INDEX_TXP> & sudoListSegment, std::vector<ATT_ELM> const & vecAttELm )
	{
		for( AttFractElm const & me : vecAttELm )
		{
			INDEX_TXP sudoNeeded = me.getSudo();

			bool sudoIsKnown = false;

			for( INDEX_TXP sudoInList : sudoListSegment )
			{
				if( sudoInList == sudoNeeded )
				{
					sudoIsKnown = true;
					break;
				}
			}

			if( ! sudoIsKnown )
				sudoListSegment.push_back(sudoNeeded);
		}

		return true;
	}

	template <typename ATT_ELM,
			  typename = std::enable_if<std::is_base_of<AttFractElm,ATT_ELM>::value>
			>
	bool averageNormlForEachSudo( std::vector<ATT_ELM> const & vecAttElm, VecPairSudoNormlV & vecPSudoNrml )
	{
		// first determine appearing sudos

		std::vector<INDEX_TXP> sudoListSegment;

		extractSudoList(sudoListSegment,vecAttElm);

		for( INDEX_TXP sudo : sudoListSegment )
		{
			VECTOR_TYP normlAvrg;

			if( ! averageNormalForSpecificSudo( sudo, vecAttElm, normlAvrg ) )
				return false;

			std::pair<INDEX_TXP, VECTOR_TYP> sudoNorml( sudo, normlAvrg );

			vecPSudoNrml.push_back(sudoNorml);
		}

		return true;
	}

	template <typename ATT_ELM,
			  typename = std::enable_if<std::is_base_of<AttFractElm,ATT_ELM>::value>
			>
	bool averageNormalForSpecificSudo( INDEX_TXP specfcSudo, std::vector<ATT_ELM> const & vecAttElm, VECTOR_TYP & normlAvrg )
	{
		VECTOR_TYP normsSum(0,0,0);
		INDEX_TXP numContrNrmls = 0;

		for( ATT_ELM const & ae : vecAttElm )
		{
			INDEX_TXP sudoElm = ae.getSudo();

			if( specfcSudo == sudoElm )
			{
				VECTOR_TYP normElm = ae.getNormalVec();

				VECTOR_TYP tmpSum = normsSum;

				VecAdd( normsSum, normElm, tmpSum );

				numContrNrmls++;
			}
		}

		if( numContrNrmls == 0 )
		{
			UG_LOG("Kein Beitrag in SUdo? " << std::endl);
			UG_THROW("Kein Beitrag in SUdo? " << std::endl);
			return false;
		}

		VecScale( normlAvrg, normsSum, 1. / static_cast<number>(numContrNrmls) );

		return true;

	}

	bool averageFractNormals()
	{
		return averageNormlForEachSudo( m_vecAttFractElms, m_vecFractSudosNormlV );
	}

	bool averageBndryNormals()
	{
		if( m_isBoundary )
		{
			return averageNormlForEachSudo( m_vecAttBndryElms, m_vecBndrySudosNormlV );
		}
		else
		{
			UG_LOG("no boundary, no averaging");
			return false;
		}
	}

	template
	< typename ATT_ELM,
	  typename = std::enable_if<std::is_base_of<AttFractElm,ATT_ELM>::value>
	>
	bool schluckVecAttElm( std::vector<ATT_ELM> const & vecAttElNew, std::vector<ATT_ELM> & vecAttElmKnown, bool acceptUnknowns = false )
	{
		bool allUnknown = true;

		for( ATT_ELM const & aeN : vecAttElNew )
		{
			if( ! schluckAttElm( aeN, vecAttElmKnown, acceptUnknowns ) )
			{
				allUnknown = false;
				UG_LOG("ist schon bekannt" << std::endl);
				if( ! acceptUnknowns)
					UG_THROW("ist schon bekannt" << std::endl);
				//return false;
			}
		}

		return allUnknown;
	}

	template
	< typename ATT_ELM,
	  typename = std::enable_if<std::is_base_of<AttFractElm,ATT_ELM>::value>
	>
	bool schluckAttElm( ATT_ELM const & attElNew, std::vector<ATT_ELM> & vecAttElmKnown, bool acceptUnknowns = false )
	{
		m_averaged = false;

		if( ! isStillUnknown( attElNew, vecAttElmKnown, acceptUnknowns ) )
		{
			UG_LOG("ist schon bekannt" << std::endl);
			if( ! acceptUnknowns )
				UG_THROW("ist schon bekannt" << std::endl);
			return false;
		}

		vecAttElmKnown.push_back(attElNew);

		return true;
	}

	bool checkIfIsAtBndry()
	{
		if( ! m_isBoundary )
		{
			UG_LOG("gibts keine Bndry " << std::endl);
			UG_THROW("gibts keine Bndry " << std::endl);
			return false;
		}

		return m_isBoundary;
	}

//	template<ManifDescr::ManifoldType manifTyp,
//	template<typename ManifDescr::ManifoldType manifTyp,
	template<typename ManifDescr::ManifoldType manifTyp,
			 typename = std::enable_if< std::is_same<VECTOR_TYP,vector3>::value>
			>
	bool const spuckManifDescr( VecManifDescr & vecManifDesc,
								Grid::VertexAttachmentAccessor<APosition> const & aaPos,
								VecPairSudoNormlV const & vecFractSudosNormlV,
								bool clearDescVec = true
	) const
	{
		if( ! m_averaged )
		{
			UG_LOG("please average " << std::endl);
			UG_THROW("please average " << std::endl);
			return false;
		}

		if( clearDescVec )
			vecManifDesc.clear();

		for( PairSudoNormlV const & psn : vecFractSudosNormlV )
		{
			VECTOR_TYP posVrt = aaPos[m_vrt];

			int sudo = psn.first;
			VECTOR_TYP normlVec = psn.second;

//			ManifoldDescriptor( VECTOR_TYP const & normalVect,
//								VECTOR_TYP const & baseVect,
//								int sudo = -1,
//								ManifoldType = isArtificial,
//								number scaleShiftNormal = 0
//								)
//

			UG_LOG("ASSIGN MANIF TYP " << manifTyp << std::endl);
			ManifDescr manifDesc( normlVec, posVrt, sudo, manifTyp );

			vecManifDesc.push_back( manifDesc );
		}

		return true;
	}


};

//////////////////////////////////////////////////////////////////



#if 0
template <typename VEC_AVEI, typename OPERATION, typename INDX_TYP  >
bool switchFulldimInfo( VEC_AVEI & vecAttVolElemInfoCop,
					    VEC_AVEI const & vecAttVolElemInfo,
					    VEC_AVEI & segmentAVEI,
					    OPERATION opera,
						INDX_TYP switchInd = 0
					 )
{
	auto & startVolInfoThisSegment = vecAttVolElemInfoCop[switchInd];

	auto const & startVol = startVolInfoThisSegment.opera();

	for( auto & possibleOrigVolInfo : vecAttVolElemInfo )
	{
		auto const & possVol = possibleOrigVolInfo.opera();

		if( possVol == startVol )
		{
			segmentAVEI().push_back(possibleOrigVolInfo);
			break;
		}
	}

	if( segmentAVEI().size() != 1 )
	{
		UG_LOG("No start volume reconstructible " << std::endl);
		UG_THROW("No start volume reconstructible " << std::endl);
		return false;
	}

	if( ! vecAttVolElemInfoCop.erase( vecAttVolElemInfoCop.begin() + switchInd ) )
		return false;

	return true;

}
#endif


//////////////////////////////////////////////////////////////////////////

template
<
typename FULLDIMEL,
typename MANIFEL,
typename LOWDIMEL,
typename VRTXTYP
>
class EndingCrossingFractSegmentInfo
{
public:

	EndingCrossingFractSegmentInfo(
									VRTXTYP const & vrt,
									MANIFEL const & endingFractManifCutting,
									MANIFEL const & endingFractManifNotCutting,
									LOWDIMEL const & oldLowDimElCut,
									std::pair<MANIFEL,MANIFEL> const & pairNeighbouredFractClosedManifEl
								  )
	:
		m_unclosedVrtx(vrt),
		m_endingFractManifCutting(endingFractManifCutting),
		m_endingFractManifNotCutting(endingFractManifNotCutting),
		m_pairNeighbouredFractClosedManifEl(pairNeighbouredFractClosedManifEl),
		m_vecClosedFracManifEl(std::vector<MANIFEL>()),
		m_oldLowDimElCut( oldLowDimElCut ),
		m_sudoFractEnding(-1),
		m_sudoFractNotEnding(-1),
		m_vecFulldimEl(std::vector<FULLDIMEL>())
	{
	};

private:

	// TODO FIXME vielleicht statt Manifel die Klasse,  Attached Fracture Objekte? mit richtig geordneten Edges?

	VRTXTYP m_unclosedVrtx;
	MANIFEL m_endingFractManifCutting;
	MANIFEL m_endingFractManifNotCutting;

	std::pair<MANIFEL,MANIFEL> m_pairNeighbouredFractClosedManifEl;

	std::vector<MANIFEL> m_vecClosedFracManifEl;

	LOWDIMEL m_oldLowDimElCut; // common edge between ending frac face with one sudo and durchgehende frac faces with another sudo
//	LOWDIMEL m_newLowDimElCut;

	int m_sudoFractEnding;
	int m_sudoFractNotEnding;

	std::vector<FULLDIMEL> m_vecFulldimEl;
};

//////////////////////////////////////////////////////////////////////////

}

}

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_ */
