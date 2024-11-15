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


template <
typename MANIFOLDTYP,
typename FULLDIMTYP,
typename SENKRECHTENTYP
>
class VertexFractureTripleMF
{

public:

	VertexFractureTripleMF( MANIFOLDTYP const & manif, FULLDIMTYP const & full, SENKRECHTENTYP const & normal   )
	: m_manif(manif), m_full(full), m_normal(normal), m_newNormal(normal)
	{
	};

	MANIFOLDTYP const getManif() const { return m_manif; }

	FULLDIMTYP const getFull() const { return m_full; }

	SENKRECHTENTYP const getNormal() const { return m_normal; }

	void setNewNormal( SENKRECHTENTYP const & chNorml ) { m_newNormal = chNorml; }
	SENKRECHTENTYP const getNewNormal() const { return m_newNormal; }

private:

	MANIFOLDTYP m_manif;
	FULLDIMTYP m_full;
	SENKRECHTENTYP  m_normal;
	SENKRECHTENTYP  m_newNormal;

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


//////////////////////////////////////////////////////////////////////////////

// class to help count and store a bool and a number of templete type
// comparable to std::pair<bool,int> but more dedicated to the specific aim
// TODO FIXME adapt for 3D case, figure out if inner end, and number of fracs sourrounding

template< typename T >
class VertexFracturePropertiesVol
{
public:

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
	  m_sudoList( std::vector<T>() )
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

		if( m_sudoList.size() > static_cast<T>( m_maxStatus ) )
		{
			UG_THROW("zu viele subdomains crossing in one Punkt" << std::endl);
			return false;
		}

		return alreadyInList;
	}

	// TODO FIXME
	bool getIsAClosedFractur( T sudo )
	{
//		static_assert(false);

		return {};
	}

	void setIsAClosedFractur( T sudo, bool isClosed )
	{
		// TODO FIXME
		//static_assert(false);
	}

	bool getInfoAllFracturesClosed()
	{
		//static_assert(false);

		return {};
	}

	bool getInfoNoFracturesClosed()
	{
		//static_assert(false);

		return {};
	}

	void setFracSurroundProts( bool surrounded, T sudo )
	{
		// TODO FIXME
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

};



}

}

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_ */
