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
#include <algorithm>

namespace ug
{

namespace arte
{
////////////////////////////////////////////7

namespace diamonds
{

///////////////////////////////////////////

template
<
typename FULLDIMELM,
typename VERTEXTYP,
typename INDEXTYP
>
class SideElemsOfDiamsCreatedDirectly
{
public:

	SideElemsOfDiamsCreatedDirectly( FULLDIMELM const & fudielm, VERTEXTYP const & connectVrtx, INDEXTYP subsEnd, INDEXTYP subsNotEnd )
	: m_fulldimElm(fudielm), m_connectVrtx(connectVrtx), m_subsEnding(subsEnd), m_subsNotEnding(subsNotEnd)
	{}

	FULLDIMELM spuckFulldimElm() { return m_fulldimElm; }

	VERTEXTYP spuckConnectVrtx() { return m_connectVrtx; }

	INDEXTYP spuckSubsEnding() { return m_subsEnding; }

	INDEXTYP spuckSubsNotEnding() { return m_subsNotEnding; }

private:

	FULLDIMELM m_fulldimElm;
	VERTEXTYP m_connectVrtx; // the vertex at the other end of the cut edge, where the fractures are all closed
	INDEXTYP m_subsEnding, m_subsNotEnding;
};

////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////

template <
typename FULLDIMELEM,
typename MANIFELEM,
typename LOWDIMELM,
typename VERTEXTYP,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class VolManifVrtxCombi
{
public:

	using VrtxPair = std::pair<VERTEXTYP,VERTEXTYP>;

	VolManifVrtxCombi( FULLDIMELEM const & vol,
					   MANIFELEM const & manif,
					   VrtxPair const & oldAndShiftVrtx,
					   INDEXTYP	sudo,
					   bool hasOneEndingCrossingCleft = false,
					   FULLDIMELEM const & partnerVol = nullptr
					 )
	: m_volElm(vol),
	  m_manifElm(manif),
	  m_oldAndshiftVrtx(oldAndShiftVrtx),
	  m_sudo(sudo),
	  m_lowDimElm(nullptr),
	  m_hasOneEndingCrossingCleft(hasOneEndingCrossingCleft),
	  m_volPartnerElm( partnerVol )
	{};

	void spuckFulldimElem( FULLDIMELEM & vol ) { vol = m_volElm; }

	void spuckFulldimPartnerElem( FULLDIMELEM & partVol ) { partVol = m_volPartnerElm; }

	void spuckManif( MANIFELEM & manif ) { manif = m_manifElm; }

	void spuckOldAndShiftVrtx( VrtxPair & vrtp ) { vrtp = m_oldAndshiftVrtx; }

	void spuckLowDimElem ( LOWDIMELM & lowdimElm ) { lowdimElm = m_lowDimElm; }

	INDEXTYP spuckSudo() { return m_sudo; }

	bool spuckHasOneEndingCrossingCleft() { return m_hasOneEndingCrossingCleft; }

	void changeTheElems( FULLDIMELEM const & vol,
						 MANIFELEM const & manif,
						 VERTEXTYP const & newBaseVrtx,
						 FULLDIMELEM const & partnerVol = nullptr
					   )
	{
		m_volElm = vol;
		m_manifElm = manif;
		m_oldAndshiftVrtx.first = newBaseVrtx;
		// only assigned again after integrity check:
		m_lowDimElm = nullptr;
		m_volPartnerElm = partnerVol;
	}

	// compute the edge which connects the both vertices!!!!
	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELM>::value>,
		typename = std::enable_if<std::is_same<Vertex*,VERTEXTYP>::value>
	>
	bool checkIntegrity( Grid & grid )
	{
		INDEXTYP edgeFound = 0;

		FULLDIMELEM testElm = m_volElm;

		if( m_volPartnerElm && m_hasOneEndingCrossingCleft )
			testElm = m_volPartnerElm;

//		for(size_t i_edge = 0; i_edge < m_volElm->num_edges(); ++i_edge)
		for(size_t i_edge = 0; i_edge < testElm->num_edges(); ++i_edge)
		{
//			LOWDIMELM lowDimElm = grid.get_edge( m_volElm, i_edge );
			LOWDIMELM lowDimElm = grid.get_edge( testElm, i_edge );

			if(    EdgeContains(lowDimElm, m_oldAndshiftVrtx.first)
				&& EdgeContains(lowDimElm, m_oldAndshiftVrtx.second )
			   )
			{
				m_lowDimElm = lowDimElm;
				edgeFound++;
			}
		}

		if( edgeFound != 1)
		{
			UG_LOG("edge number found strange " << edgeFound << std::endl);
			return false;
		}

		return true;
	}


private:

	FULLDIMELEM m_volElm;
	MANIFELEM m_manifElm;
	VrtxPair m_oldAndshiftVrtx;
	INDEXTYP m_sudo;
	LOWDIMELM m_lowDimElm;
	bool m_hasOneEndingCrossingCleft;
	FULLDIMELEM m_volPartnerElm;


};


/////////////////////////////////////////////////////////////////

// a pair of a volume and a face, the face should be part of the volume!
template <
typename FULLDIMELEM,
typename LOWDIMELEM,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>
>
class FulldimLowdimTwin
{
public:

	// TODO FIXME Paarvolumen dazu machen, das bei ECC dabei ist ggf

	FulldimLowdimTwin( FULLDIMELEM const & fulldimElem,
					   LOWDIMELEM const & lowdimElem,
					   INDEXTYP sudo,
					   FULLDIMELEM const & fulldimElemPartner = nullptr
				     )
	: m_fullDimElem(fulldimElem), m_lowDimElem(lowdimElem), m_sudo(sudo),
	  m_fullDimElemPartner(fulldimElemPartner),
	  m_relevantFullDimElem(nullptr), m_connectingFullDimElem(nullptr)
	{}

	FulldimLowdimTwin()
	: m_fullDimElem(nullptr), m_lowDimElem(nullptr), m_sudo(0),
	  m_fullDimElemPartner(nullptr),
	  m_relevantFullDimElem(nullptr), m_connectingFullDimElem(nullptr)
	{}

	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool setRelevantFulldimElem( Grid & grid )
	{
		if( ! m_fullDimElem )
		{
			UG_LOG("cannot test an empty volume " << std::endl);
			return false;
		}

		m_connectingFullDimElem = m_fullDimElem;

		if( edgeIsInVolume( m_fullDimElem, m_lowDimElem, grid ) )
		{
			m_relevantFullDimElem = m_fullDimElem;

			return true;
		}

		if( ! m_fullDimElemPartner )
		{
			UG_LOG("edge not in main elem but also no partner " << std::endl);
			return false;
		}

		if( edgeIsInVolume( m_fullDimElemPartner, m_lowDimElem, grid ) )
		{
			m_relevantFullDimElem = m_fullDimElemPartner;

			return true;
		}

		UG_LOG("somehow no Zuordnung volume edge " << std::endl);

		return false;
	}


	void spuckFullDimElem( FULLDIMELEM & fulldimElem )
	{
		fulldimElem = m_fullDimElem;
	}

	void spuckFullDimElemPartner( FULLDIMELEM & fulldimElemPartner )
	{
		fulldimElemPartner = m_fullDimElemPartner;
	}

	bool spuckRelevantFullDimElem( FULLDIMELEM & relevantFulldimElem )
	{
		if( ! m_relevantFullDimElem )
		{
			UG_LOG("relevant element not determined so far " << std::endl);
			return false;
		}

		relevantFulldimElem = m_relevantFullDimElem;

		return true;
	}

	bool spuckConnectingFullDimElem( FULLDIMELEM & connectingFulldimElem )
	{
		if( ! m_connectingFullDimElem )
		{
			UG_LOG("connecting element not determined so far " << std::endl);
			return false;
		}

		connectingFulldimElem = m_connectingFullDimElem;

		return true;
	}


	void spuckLowDimElem( LOWDIMELEM & lowdimElem )
	{
		lowdimElem = m_lowDimElem;
	}

	INDEXTYP spuckSudo() { return m_sudo; }

	bool changeTheElems( FULLDIMELEM const & fulldimElem,
						 LOWDIMELEM const & lowdimElem,
						 INDEXTYP sudo,
						 FULLDIMELEM const & fullDimElemPartner = nullptr )
	{
		m_fullDimElem = fulldimElem;
		m_lowDimElem = lowdimElem;
		m_sudo = sudo;

		if( fullDimElemPartner )
			m_fullDimElemPartner = fullDimElemPartner;

		if( ! setRelevantFulldimElem() )
		{
			UG_LOG("Change but no relevant edge " << std::endl);
			return false;
		}

		return true;
	}

	// template check if volume and edge valid......
	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool checkIntegrity()
	{
		// TODO FIXME den Partner checken
		if( ! VolumeContains(m_fullDimElem, m_lowDimElem))
		{
			UG_LOG("Volume does not contain edge for diams " << std::endl);

			if( m_fullDimElemPartner )
			{
				UG_LOG("einen Partner hat es schonmal " << std::endl);

				if( ! VolumeContains(m_fullDimElemPartner, m_lowDimElem))
				{
					UG_LOG("aber auch der hat die Ecke nicht" << std::endl);
					return false;
				}

				UG_LOG("Partner doch gefunden" << std::endl);
			}
			else
			{
				UG_LOG("aber ohne Partner" << std::endl);
				return false;
			}
		}

		return true;
	}

private:

	FULLDIMELEM m_fullDimElem;
	LOWDIMELEM m_lowDimElem;
	INDEXTYP m_sudo;
	FULLDIMELEM m_fullDimElemPartner;
	FULLDIMELEM m_relevantFullDimElem;
	FULLDIMELEM m_connectingFullDimElem;


	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool edgeIsInVolume( FULLDIMELEM & vol, LOWDIMELEM & edg, Grid & grid )
	{
		for( INDEXTYP iE = 0; iE < vol->num_edges(); iE++ )
		{
			LOWDIMELEM testEdg = grid.get_edge( vol, iE );

			if( testEdg == edg )
				return true;
		}

		UG_LOG("Die relevante Ecke unauffindbar " << std::endl);

		return false;
	}

};

////////////////////////////////////////////////////////////


// exists now several times in several namespaces, should be put in a general file.....
template<typename ELEMTYP>
bool addElem(std::vector<ELEMTYP> & knownElems, ELEMTYP elemToAdd )
{
	bool unknown = true;

	for( ELEMTYP elmKnown : knownElems )
	{
		if( elemToAdd == elmKnown )
		{
			unknown = false;
			break;
		}
	}

	if( unknown )
		knownElems.push_back(elemToAdd);

	return unknown;

}


/////////////////////////////////////////////////////////////////

// the pairs of per manifold connected volumes including the edges between old and shift vertex
template <
typename FULLDIMELEM,
typename MANIFELEM,
typename LOWDIMELEM,
typename VERTEXTYP,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class FullLowDimManifQuintuplet
{
public:

	using FullLowDimTwin = FulldimLowdimTwin<FULLDIMELEM,LOWDIMELEM,INDEXTYP>;
	using PairFullLowDimTwin = std::pair<FullLowDimTwin,FullLowDimTwin>;
	using PairVrtcs = std::pair<VERTEXTYP,VERTEXTYP>;
	using PairLowDimElem = std::pair<LOWDIMELEM,LOWDIMELEM>;

	// TODO FIXME auch die Möglichkeit von Triplett
	// dann durch Vektor
	// oder vielleicht immer Vektor aber muss 2 oder 3 Länge sein

	FullLowDimManifQuintuplet( PairFullLowDimTwin const & fullLowPr, MANIFELEM const & manif )
	: m_pairFullLowDimTwin(fullLowPr),
	  m_manifElem(manif),
	  m_centerVrtx(nullptr),
	  m_shiftVrtcs(PairVrtcs()),
	  m_pairLowDimElem(PairLowDimElem()),
	  m_sudo(0)
	{};

	FullLowDimManifQuintuplet()
	: m_pairFullLowDimTwin(PairFullLowDimTwin()),
	  m_manifElem(nullptr),
	  m_centerVrtx(nullptr),
	  m_shiftVrtcs(PairVrtcs()),
	  m_pairLowDimElem(PairLowDimElem()),
	  m_sudo(0)
	{};

	// TODO FIXME das checkAlsoFace wieder abschaffen, muss immer richtig sein
	// oder lassen aber dann besonderen Umgang für Fall, dass Triplett ist
	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
//	bool checkIntegrity( bool checkAlsoFace = true )
	bool checkIntegrity()
//	bool checkIntegrity( SubsetHandler & sh )
	{
		// TODO FIXME das muss den Partner kennen, damit es gründlich testen kann
		// WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
		if( ! checkIntegrityVols() )
		{
			UG_LOG("Vols not integer " << std::endl);
			return false;
		}

		// TODO FIXME das muss auch wieder gecheckt werden, aber intelligenter
//		if( checkAlsoFace )
//		if() // TODO FIXME test if partner available in one or both volumes
//		{
//		if( ! checkIntegrityFaceInBothVols( sh ) )
		if( ! checkIntegrityFaceInBothVols() )
		{
			UG_LOG("A face in the dark " << std::endl);
			return false;
		}
//		}

		if( ! figureOutMajorVertices() )
		{
			UG_LOG("major vertices not found " << std::endl);
			return false;
		}

		if( ! figureOutSudo() )
		{
			UG_LOG("sudo not unique " << std::endl);
			return false;
		}

		return true;
	}

	void spuckCenterVertex( VERTEXTYP & vrt )
	{
		vrt = m_centerVrtx;
	}

	void spuckShiftVrtcs( PairVrtcs & pv )
	{
		pv = m_shiftVrtcs;
	}

	void spuckPairFullLowDimTwin( PairFullLowDimTwin & pfldt )
	{
		pfldt = m_pairFullLowDimTwin;
	}

	void spuckManifElem( MANIFELEM & m )
	{
		m = m_manifElem;
	}

//	bool swapEntries( SubsetHandler & sh )
	bool swapEntries()
	{
		std::swap( m_pairFullLowDimTwin.first, m_pairFullLowDimTwin.second );

		VERTEXTYP vrtOne = m_shiftVrtcs.first;
		VERTEXTYP vrtTwo = m_shiftVrtcs.second;

		std::swap( m_shiftVrtcs.first, m_shiftVrtcs.second );

//		FullLowDimTwin fldOne = m_pairFullLowDimTwin.first;
//		FullLowDimTwin fldTwo = m_pairFullLowDimTwin.second;


//		if( fldOne != m_pairFullLowDimTwin.second || fldTwo != m_pairFullLowDimTwin.first )
//		{
//			UG_LOG("swappign not worked " << std::endl);
//			return false;
//		}

//		if( fldOne != m_pairFullLowDimTwin.second || fldTwo != m_pairFullLowDimTwin.first )
//		{
//			UG_LOG("swappign not worked " << std::endl);
//			return false;
//		}

		if( vrtOne != m_shiftVrtcs.second || vrtTwo != m_shiftVrtcs.first )
		{
			UG_LOG("swaping vrtx not work " << std::endl);
			return false;
		}

//		if( ! checkIntegrity( sh ) )
		if( ! checkIntegrity() )
		{
			UG_LOG("not integer any more after swap entries" << std::endl);
			return false;
		}

		return true;
	}

	void spuckPairLowDimElem( PairLowDimElem & prLdE )
	{
		prLdE = m_pairLowDimElem;
	}

	INDEXTYP spuckSudo() { return m_sudo; }

private:

	PairFullLowDimTwin m_pairFullLowDimTwin;
	MANIFELEM m_manifElem;
	VERTEXTYP m_centerVrtx;
	PairVrtcs m_shiftVrtcs;
	PairLowDimElem m_pairLowDimElem;
	INDEXTYP m_sudo;

	bool checkIntegrityVols()
	{
		if( ! m_pairFullLowDimTwin.first.checkIntegrity() || ! m_pairFullLowDimTwin.second.checkIntegrity() )
		{
			UG_LOG("Vol integrity not passed " << std::endl);
			return false;
		}

		return true;
	}


	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool checkIntegrityFaceInVol( FullLowDimTwin & fldt )
//	bool checkIntegrityFaceInVol( FullLowDimTwin & fldt, SubsetHandler & sh )
	{
		// TODO FIXME hier muss getestet werden, ob es einen Partner gibt, und wenn ja,
		// dann muss getestet werden, ob er das face enthält!
		// in einem Spezialfall bezieht sich das Ganze dann auf eine Ecke nur
		// EEEEEEEEEEEEEEEEEEEEEEEEEEEEe

		FULLDIMELEM connectingFulldimElem;

		fldt.spuckConnectingFullDimElem( connectingFulldimElem );

		if( ! VolumeContains( connectingFulldimElem, m_manifElem ) )
		{
			UG_LOG("Face not in Vol " << std::endl);
			return false;
		}

		return true;
	}


//		FULLDIMELEM fudielPartner = nullptr;
//
//		fldt.spuckFullDimElemPartner( fudielPartner );
//
//		FULLDIMELEM fudiel;
//
//		fldt.spuckFullDimElem(fudiel);
//
//		// check who contains the edge, the partner volume or the volume itself
//		// that one which contains the edge has to be tested, and NOT the other one
//		// in case of no partner, it is simple
//
//		LOWDIMELEM lodiEl;
//		fldt.spuckLowDimElem(lodiEl);
//
//
//		if( ! fudielPartner || VolumeContains( fudiel, lodiEl ) )
//		{
//			if( ! VolumeContains(fudiel, m_manifElem ) )
//			{
//				UG_LOG("Face not in Vol " << std::endl);
//				return false;
//			}
//		}
//		else
//		{
//			if( ! VolumeContains(fudielPartner, m_manifElem || ! VolumeContains( fudielPartner, lodiEl ) ) )
//			{
//				UG_LOG("Face Partner is not in Vol " << std::endl);
//
//				sh.assign_subset( fudielPartner, sh.num_subsets() );
//				sh.assign_subset( fudiel, sh.num_subsets() );
//				sh.assign_subset( m_manifElem, sh.num_subsets() );
//
//				return false;
//			}
//		}
//
//		return true;



	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool checkIntegrityFaceInBothVols()
//	bool checkIntegrityFaceInBothVols( SubsetHandler & sh )
	{
//		FULLDIMELEM firstV, secondV;
//		m_pairFullLowDimTwin.first.spuckFullDimElem(firstV);
//		m_pairFullLowDimTwin.second.spuckFullDimElem(secondV);

//		if( ! checkIntegrityFaceInVol( firstV ) || ! checkIntegrityFaceInVol( secondV )
//		if(    ! checkIntegrityFaceInVol( m_pairFullLowDimTwin.first, sh)
//			|| ! checkIntegrityFaceInVol( m_pairFullLowDimTwin.second, sh )
//		  )
		if(    ! checkIntegrityFaceInVol( m_pairFullLowDimTwin.first )
			|| ! checkIntegrityFaceInVol( m_pairFullLowDimTwin.second )
		  )
		{
			UG_LOG("face not in one vol at least " << std::endl);
			return false;
		}

		return true;
	}

	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Vertex*,VERTEXTYP>::value>
	>
	bool figureOutMajorVertices()
	{
		Edge * edgeOne;
		Edge * edgeTwo;

		m_pairFullLowDimTwin.first.spuckLowDimElem(edgeOne);
		m_pairFullLowDimTwin.second.spuckLowDimElem(edgeTwo);
		//		edgeOne = m_pairFullLowDimTwin.first.spuckLowDimElem();
		//		edgeTwo = m_pairFullLowDimTwin.first.spuckLowDimElem();

		m_pairLowDimElem = PairLowDimElem( edgeOne, edgeTwo );

		Vertex * centerVrtx;
		Vertex * shiftVrtxOne;
		Vertex * shiftVrtxTwo;
		centerVrtx = nullptr;
		shiftVrtxOne = nullptr;
		shiftVrtxTwo = nullptr;

		if( ! findConnectingAndExtrnlVertex(edgeOne, edgeTwo, centerVrtx, shiftVrtxOne, shiftVrtxTwo) )
		{
			UG_LOG("Vertices not found " << std::endl);
			return false;
		}

		return true;
	}

	template
	<
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Vertex*,VERTEXTYP>::value>
	>
	bool findConnectingAndExtrnlVertex( LOWDIMELEM const & lowDimElemOne,
										LOWDIMELEM const & lowDimElemTwo,
										VERTEXTYP & connctVrtx,
										VERTEXTYP & outerVrtxOne,
										VERTEXTYP & outerVrtxTwo
									   )
	{
		std::vector<Vertex*> verticesFromEdges;

		// all edges have two vertices, so max number two
		for( int iVrt = 0; iVrt < 2 ; iVrt++ )
		{
			Vertex * vrtOne = lowDimElemOne->vertex(iVrt);
			Vertex * vrtTwo = lowDimElemTwo->vertex(iVrt);

			addElem( verticesFromEdges, vrtOne );
			addElem( verticesFromEdges, vrtTwo );
		}

		if( verticesFromEdges.size() != 3 )
		{
			UG_LOG("vertex number for diamonds strange " << verticesFromEdges.size() << std::endl);
		}

		int connVrtFound = 0;
		int extVrtOneFound = 0;
		int extVrtTwoFound = 0;

		for( auto const & v : verticesFromEdges )
		{
			if( EdgeContains(lowDimElemOne, v) && EdgeContains( lowDimElemTwo, v) )
			{
				connctVrtx = v;
				connVrtFound++;
			}
			else if( EdgeContains(lowDimElemOne, v) && ! EdgeContains( lowDimElemTwo, v) )
			{
				outerVrtxOne = v;
				extVrtOneFound++;
			}
			else if( ! EdgeContains(lowDimElemOne, v) && EdgeContains( lowDimElemTwo, v) )
			{
				outerVrtxTwo = v;
				extVrtTwoFound++;
			}
			else
			{
				UG_LOG("strange case of vertices and edges " << std::endl);
				return false;
			}
		}

		if(    connVrtFound != 1
			|| connctVrtx == nullptr
			|| extVrtOneFound != 1
			|| extVrtTwoFound != 1
			|| outerVrtxOne == nullptr
			|| outerVrtxTwo == nullptr
		)
		{
			UG_LOG("connection vertex number " << connVrtFound << std::endl);
			UG_LOG("ext one vertex number " << extVrtOneFound << std::endl);
			UG_LOG("ext two vertex number " << extVrtTwoFound << std::endl);
			return false;
		}

		m_centerVrtx = connctVrtx;
		m_shiftVrtcs = PairVrtcs( outerVrtxOne, outerVrtxTwo );

		return true;
	}

	bool figureOutSudo()
	{
		int sudoFirst = m_pairFullLowDimTwin.first.spuckSudo();
		int sudoSecond = m_pairFullLowDimTwin.second.spuckSudo();

		if( sudoFirst != sudoSecond )
		{
			UG_LOG("Sudos do not coincide " << std::endl );
		}

		m_sudo = sudoFirst;

		return true;
	}

};



/////////////////////////////////////////////////////////////////


template <
typename FULLDIMELEM,
typename MANIFELEM,
typename LOWDIMELEM,
typename VERTEXTYP,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class ElemsToBeQuenched4DiamSpace
{
public:

	using FullLowDimManifQntpl = FullLowDimManifQuintuplet<FULLDIMELEM,MANIFELEM,LOWDIMELEM,VERTEXTYP,INDEXTYP>;
	using VecFullLowDimManifQuintuplet = std::vector<FullLowDimManifQntpl>;
	using PairLowDimElem = std::pair<LOWDIMELEM,LOWDIMELEM>;
	using PairVrtcs = std::pair<VERTEXTYP,VERTEXTYP>;


	ElemsToBeQuenched4DiamSpace( VecFullLowDimManifQuintuplet const & vfldm5 )
	: m_centerVrtx(nullptr),
	  m_originalCenterVrtx(nullptr),
	  m_vecFullLowDimManifQuintpl(vfldm5),
	  m_pairLowDimElem(PairLowDimElem()),
	  m_sudo(0),
	  m_shiftVrtcs(PairVrtcs()),
	  m_midPointOfShiftVrtcs(nullptr)
	{}

	ElemsToBeQuenched4DiamSpace()
	: m_centerVrtx(nullptr),
	  m_originalCenterVrtx(nullptr),
	  m_vecFullLowDimManifQuintpl(VecFullLowDimManifQuintuplet()),
	  m_pairLowDimElem(PairLowDimElem()),
	  m_sudo(0),
	  m_shiftVrtcs(PairVrtcs()),
	  m_midPointOfShiftVrtcs(nullptr)
	{}


//	bool checkIntegrity( SubsetHandler & sh )
	bool checkIntegrity()
	{
		bool centerAssigned = false;

		bool sudoAssigned = false;

		bool pairLowdimElmAssigned = false;

		bool pairShiftVrtcsAssigned = false;

		for( auto & fldmq : m_vecFullLowDimManifQuintpl )
		{
//			if( ! fldmq.checkIntegrity( sh ) )
			if( ! fldmq.checkIntegrity() )
			{
				UG_LOG("fulllowdim manif 5 not integer " << std::endl);
				return false;
			}

			if( ! centerAssigned )
			{
				fldmq.spuckCenterVertex( m_centerVrtx );
				centerAssigned = true;

				if( m_originalCenterVrtx == nullptr )
				{
					// must be first call
					m_originalCenterVrtx = m_centerVrtx;
				}

			}
			else
			{
				Vertex * testVrtx = nullptr;
				fldmq.spuckCenterVertex( testVrtx );

				if( testVrtx != m_centerVrtx )
				{
					UG_LOG("different centers " << std::endl);
					return false;
				}
			}

			if( ! pairLowdimElmAssigned )
			{
				fldmq.spuckPairLowDimElem(m_pairLowDimElem);
				pairLowdimElmAssigned = true;

			}
			else
			{
				PairLowDimElem testPrLDE;
				fldmq.spuckPairLowDimElem(testPrLDE);

				if( m_pairLowDimElem != testPrLDE )
				{
					// check for swaped pair
					std::swap(testPrLDE.first, testPrLDE.second);

					if( testPrLDE !=  m_pairLowDimElem )
					{
						UG_LOG("different lowdim pairs " << std::endl);
						return false;
					}

//					if( ! fldmq.swapEntries( sh ) )
					if( ! fldmq.swapEntries() )
					{
						UG_LOG("entries not swappable" << std::endl);
						return false;
					}

//					if( ! fldmq.checkIntegrity( sh ) )
					if( ! fldmq.checkIntegrity() )
					{
						UG_LOG("quintuplet not integer any more " << std::endl);
						return false;
					}
				}
			}

			if( ! sudoAssigned )
			{
				m_sudo = fldmq.spuckSudo();
				sudoAssigned = true;
			}
			else
			{
				INDEXTYP sudoTest = fldmq.spuckSudo();

				if( sudoTest != m_sudo )
				{
					UG_LOG("sudos not the same " << std::endl);
					return false;
				}
			}

			if( ! pairShiftVrtcsAssigned )
			{
				fldmq.spuckShiftVrtcs(m_shiftVrtcs);
				pairShiftVrtcsAssigned = true;
			}
			else
			{
				PairVrtcs testPrV;

				fldmq.spuckShiftVrtcs(testPrV);

				if( m_shiftVrtcs != testPrV )
				{
					// check for swaped pair
					std::swap(testPrV.first, testPrV.second);

					if( testPrV !=  m_shiftVrtcs )
					{
						UG_LOG("different shift vertex pairs V" << std::endl);
						return false;
					}

//					if( ! fldmq.swapEntries( sh ) )
					if( ! fldmq.swapEntries() )
					{
						UG_LOG("entries not swappable V" << std::endl);
						return false;
					}

//					if( ! fldmq.checkIntegrity( sh ) )
					if( ! fldmq.checkIntegrity() )
					{
						UG_LOG("quintuplet not integer any more V" << std::endl);
						return false;
					}
				}
			}

//			if( fldmq.spuckCenterVertex() != m_centerVrtx )
//			{
//				UG_LOG("Center vertex not identical " << std::endl);
//				return false;
//			}
		}

		if(    ! centerAssigned
			|| m_centerVrtx == nullptr
			|| ! pairLowdimElmAssigned
			|| ! sudoAssigned
		  )
		{
			UG_LOG("CEnter problem " << std::endl);
			return false;
		}

		return true;
	}

	bool changeElems( VecFullLowDimManifQuintuplet const & vfldm5 )
	{
		if( vfldm5.size() != m_vecFullLowDimManifQuintpl.size())
		{
			UG_LOG("change elems but apply different size " << std::endl);
			return false;
		}

		m_vecFullLowDimManifQuintpl = vfldm5;

		if( ! checkIntegrity() )
		{
			UG_LOG("Quintuplet not integer " << std::endl);
			return false;
		}

		return true;
	}

	void spuckCenterVertex( VERTEXTYP & center )
	{
		center = m_centerVrtx;
	}

	void spuckOrigCenterVertex( VERTEXTYP & origCenterVrtx )
	{
		origCenterVrtx = m_originalCenterVrtx;
	}

	void spuckVecFullLowDimManifQuintuplet( VecFullLowDimManifQuintuplet & vfldm5 )
	{
		vfldm5 = m_vecFullLowDimManifQuintpl;
	}

	void spuckPairLowDimElem( PairLowDimElem & plde )
	{
		plde = m_pairLowDimElem;
	}

	INDEXTYP spuckSudo() { return m_sudo; }

	void spuckShiftVrtcs( PairVrtcs & pv )
	{
		pv = m_shiftVrtcs;
	}

	void assignMidPointOfShiftVrtcs( VERTEXTYP const & mp )
	{
		m_midPointOfShiftVrtcs = mp;
	}

	void spuckMidPointOfShiftVrtcs( VERTEXTYP & mp )
	{
		mp = m_midPointOfShiftVrtcs;
	}

private:

	VERTEXTYP m_centerVrtx, m_originalCenterVrtx;
	VecFullLowDimManifQuintuplet m_vecFullLowDimManifQuintpl;
	PairLowDimElem m_pairLowDimElem;
	INDEXTYP m_sudo;
	PairVrtcs m_shiftVrtcs;
	VERTEXTYP m_midPointOfShiftVrtcs;

};


/////////////////////////////////////////////////////////////////

template <
typename FULLDIMELEM,
typename MANIFELEM,
typename LOWDIMELEM,
typename VERTEXTYP,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class ElemGroupVrtxToBeQuenched4DiamSpace
{
public:

	using Elems2BQuenched4Diams = ElemsToBeQuenched4DiamSpace<FULLDIMELEM,MANIFELEM,LOWDIMELEM,VERTEXTYP,INDEXTYP>;

	using VecElems2BQuenched4Diams = std::vector<Elems2BQuenched4Diams>;

	ElemGroupVrtxToBeQuenched4DiamSpace( VecElems2BQuenched4Diams const & ve2bq4d )
	: m_vecElems2bQuenched(ve2bq4d), m_origCenterVrtx(nullptr), m_vecSudos(std::vector<INDEXTYP>())
	{}

	ElemGroupVrtxToBeQuenched4DiamSpace()
	: m_vecElems2bQuenched(VecElems2BQuenched4Diams()), m_origCenterVrtx(nullptr), m_vecSudos(std::vector<INDEXTYP>())
	{}


//	bool checkIntegrity( SubsetHandler & sh )
	bool checkIntegrity()
	{
		bool centerVrtxFound = false;

		for( Elems2BQuenched4Diams & e2bq : m_vecElems2bQuenched )
		{
//			if( ! e2bq.checkIntegrity( sh ) )
			if( ! e2bq.checkIntegrity() )
			{
				UG_LOG("quench elements not integer in vector at one point " << std::endl);
				return false;
			}

			if( ! centerVrtxFound )
			{
				e2bq.spuckOrigCenterVertex(m_origCenterVrtx);
				centerVrtxFound = true;
			}
			else
			{
				VERTEXTYP vrtTest;

				e2bq.spuckOrigCenterVertex(vrtTest);

				if( vrtTest != m_origCenterVrtx )
				{
					UG_LOG("center vertices do not agree in group " << std::endl);
					return false;
				}
			}

			INDEXTYP sudo = e2bq.spuckSudo();

			addElem(m_vecSudos,sudo);
		}

		std::sort(m_vecSudos.begin(), m_vecSudos.end());

		return true;
	}

	void spuckOrigCenterVertex( VERTEXTYP & ocv )
	{
		ocv = m_origCenterVrtx;
	}

	void spuckVecElems2BQuenched4Diams( VecElems2BQuenched4Diams & ve2bq )
	{
		ve2bq = m_vecElems2bQuenched;
	}

	void changeVecElems2BQuenched4Diams( VecElems2BQuenched4Diams const & ve2bq )
	{
		m_vecElems2bQuenched = ve2bq;
	}

	std::vector<INDEXTYP> spuckSudoList() { return m_vecSudos; }

private:

	VecElems2BQuenched4Diams m_vecElems2bQuenched;
	VERTEXTYP m_origCenterVrtx;
	std::vector<INDEXTYP> m_vecSudos;

};


////////////////////////////////////////////////////////


template<
typename VERTEXTYP
>
class CombiPairSingle
{

public:

	using VRTXPAIR = std::pair<VERTEXTYP,VERTEXTYP>;

	CombiPairSingle( VRTXPAIR const & vertPr, VERTEXTYP const & midVrtx )
	: m_vrtxPair(vertPr), m_midVrtx(midVrtx)
	{
	}

	void spuckShiftVrtxPair( VRTXPAIR & vp )
	{
		vp = m_vrtxPair;
	}

	void spuckSinglVrtx( VERTEXTYP & vrt )
	{
		vrt = m_midVrtx;
	}

private:

	VRTXPAIR m_vrtxPair;
	VERTEXTYP m_midVrtx;

};

////////////////////////////////////////////////////////

template<
typename FULLDIMELEM,
typename MANIFELEM,
typename LOWDIMELEM,
typename VRTXTYPE,
typename INDEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<VRTXTYPE>::value>
>
class CombiEntitiesProperties
{
public:

	using VrtxVec = std::vector<VRTXTYPE>;
	using IndxVec = std::vector<INDEXTYP>;
	using VrtxIndxPair = std::pair<VRTXTYPE,IndxVec>;
	using VrtxIndxCombi = std::vector<VrtxIndxPair>;

	using FulldimPr = std::pair<FULLDIMELEM,FULLDIMELEM>;

	CombiEntitiesProperties( FULLDIMELEM const & fulldmelm,
							 MANIFELEM const & manifelm,
							 VrtxIndxCombi const & centerVrtcsSudos,
							 VrtxVec const & shiftVrtcs,
							 VrtxVec const & midPtVrtcs,
							 int isThreeCross = -1
							)
							: m_fulldimElem(fulldmelm),
							  m_manifElem(manifelm),
							  m_lowdimElem(nullptr),
							  m_centerVrtcsSudos(centerVrtcsSudos),
							  m_shiftVrtcs(shiftVrtcs),
							  m_midPtVrtcs(midPtVrtcs),
							  m_prNewFulldimElems(FulldimPr()),
							  m_newSplitVrtx(nullptr),
							  m_isThreeCross(isThreeCross)
							{
							}

	CombiEntitiesProperties()
							: m_fulldimElem(nullptr),
							  m_manifElem(nullptr),
							  m_lowdimElem(nullptr),
							  m_centerVrtcsSudos(VrtxIndxCombi()),
							  m_shiftVrtcs(VrtxVec()),
							  m_midPtVrtcs(VrtxVec()),
							  m_prNewFulldimElems(FulldimPr()),
							  m_newSplitVrtx(nullptr),
							  m_isThreeCross(-1)
							{
							}

	bool schluckLowdimElem( LOWDIMELEM const & ld )
	{
		if( ld != nullptr)
		{
			m_lowdimElem = ld;
			return true;
		}

		return false;
	}

	void spuckFulldimElem( FULLDIMELEM & fde )
	{
		fde = m_fulldimElem;
	}

	void spuckManifElem( MANIFELEM & manif )
	{
		manif = m_manifElem;
	}

	void spuckLowdimElem( LOWDIMELEM & lde )
	{
		lde = m_lowdimElem;
	}

	void spuckCenterVrtcsSudos( VrtxIndxCombi & vic )
	{
		vic = m_centerVrtcsSudos;
	}

	void spuckShiftVrtcs( VrtxVec & sv )
	{
		sv = m_shiftVrtcs;
	}

	void spuckMidPtVrtcs( VrtxVec & mpv )
	{
		mpv = m_midPtVrtcs;
	}

	int spuckThreeCrossIndex() { return m_isThreeCross; }

	void schluckNewFulldimPair( FulldimPr const & fdp )
	{
		m_prNewFulldimElems = fdp;
	}

	void schluckNewSplitVrtx( VRTXTYPE const & vrt )
	{
		m_newSplitVrtx = vrt;
	}

	void spuckNewFulldimPair( FulldimPr & fdp )
	{
		fdp = m_prNewFulldimElems;
	}

	void spuckNewSplitVrtx( VRTXTYPE & vrt )
	{
		vrt = m_newSplitVrtx;
	}

private:
	FULLDIMELEM m_fulldimElem;
	MANIFELEM m_manifElem;
	LOWDIMELEM m_lowdimElem;
	VrtxIndxCombi m_centerVrtcsSudos;
	VrtxVec m_shiftVrtcs;
	VrtxVec m_midPtVrtcs;
	FulldimPr m_prNewFulldimElems;
	VRTXTYPE m_newSplitVrtx;
	int m_isThreeCross;

};

////////////////////////////////////////////////////////

template
<
typename VERTEXTYP,
typename INDEXTYP
>
class CombiCenterVrtxSudo
{
	using VrtxVec = std::vector<VERTEXTYP>;

	using IndxVec = std::vector<INDEXTYP>;

public:

	CombiCenterVrtxSudo( IndxVec const & sdV, INDEXTYP newSudo )
	: m_sudoVec(sdV), m_vrtxVec(VrtxVec()), m_newSudo(newSudo)
	{
	}

	void schluckVertex( VERTEXTYP const & vrt )
	{
		m_vrtxVec.push_back(vrt);
	}

	IndxVec spuckSudoVec() { return m_sudoVec; }

	void spuckVertexVec( VrtxVec & vv )
	{
		vv = m_vrtxVec;
	}

	INDEXTYP spuckNewSudo() { return m_newSudo; }

private:
	IndxVec m_sudoVec;
	VrtxVec m_vrtxVec;
	INDEXTYP m_newSudo;
};


////////////////////////////////////////////////////////




} // end of namespace diamonds

} // end of namespace arte

} // end of namespace ug





#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_ */
