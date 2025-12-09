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
					   INDEXTYP	sudo
					 )
	: m_volElm(vol),
	  m_manifElm(manif),
	  m_oldAndshiftVrtx(oldAndShiftVrtx),
	  m_sudo(sudo),
	  m_lowDimElm(nullptr)
	{};

	void spuckVol( FULLDIMELEM & vol ) { vol = m_volElm; }

	void spuckManif( MANIFELEM & manif ) { manif = m_manifElm; }

	void spuckOldAndShiftVrtx( VrtxPair & vrtp ) { vrtp = m_oldAndshiftVrtx; }

	void spuckLowdimElem ( LOWDIMELM & lowdimElm ) { lowdimElm = m_lowDimElm; }

	INDEXTYP spuckSudo() { return m_sudo; }

	void changeTheElems( FULLDIMELEM const & vol,
						 MANIFELEM const & manif,
						 VERTEXTYP const & newBaseVrtx
					   )
	{
		m_volElm = vol;
		m_manifElm = manif;
		m_oldAndshiftVrtx.first = newBaseVrtx;
		// only assigned again after integrity check:
		m_lowDimElm = nullptr;
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



};


/////////////////////////////////////////////////////////////////

// a pair of a volume and a face, the face should be part of the volume!
template <
typename FULLDIMELEM,
typename LOWDIMELEM
>
class FulldimLowdimTwin
{
public:

	FulldimLowdimTwin( FULLDIMELEM const & fulldimElem, LOWDIMELEM const & lowdimElem  )
	: m_fullDimElem(fulldimElem), m_lowDimElem(lowdimElem)
	{}

	void spuckFullDimElem( FULLDIMELEM & fulldimElem )
	{
		fulldimElem = m_fullDimElem;
	}

	void spuckLowDimElem( LOWDIMELEM & lowdimElem )
	{
		lowdimElem = m_lowDimElem;
	}

	void changeTheElems( FULLDIMELEM const & fulldimElem, LOWDIMELEM const & lowdimElem )
	{
		m_fullDimElem = fulldimElem;
		m_lowDimElem = lowdimElem;
	}

	// template check if volume and edge valid......
	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool checkIntegrety()
	{
		if( ! VolumeContains(m_fullDimElem,m_lowDimElem))
		{
			UG_LOG("Volume does not contain edge for diams " << std::endl);
			return false;
		}

		return true;
	}

private:

	FULLDIMELEM m_fullDimElem;
	LOWDIMELEM m_lowDimElem;
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
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class FullLowDimManifQuintuplet
{
public:

	using FullLowDimTwin = FulldimLowdimTwin<FULLDIMELEM,LOWDIMELEM>;
	using PairFullLowDimTwin = std::pair<FullLowDimTwin,FullLowDimTwin>;
	using PairVrtcs = std::pair<VERTEXTYP,VERTEXTYP>;

	FullLowDimManifQuintuplet( PairFullLowDimTwin const & fullLowPr, MANIFELEM const & manif )
	: m_pairFullLowDimTwin(fullLowPr),
	  m_manifElem(manif),
	  m_centerVrtx(nullptr),
	  m_shiftVrtcs(PairVrtcs())
	{};

	template
	<
		typename = std::enable_if<std::is_same<Volume*,FULLDIMELEM>::value>,
		typename = std::enable_if<std::is_same<Face*,MANIFELEM>::value>,
		typename = std::enable_if<std::is_same<Edge*,LOWDIMELEM>::value>
	>
	bool checkIntegrety()
	{
		if( ! checkIntegrityVols() )
		{
			UG_LOG("Vols not integer " << std::endl);
			return false;
		}

		if( ! checkIntegrityFaceInBothVols())
		{
			UG_LOG("A face in the dark " << std::endl);
			return false;
		}

		if( ! figureOutMajorVertices() )
		{
			UG_LOG("major vertices not found " << std::endl);
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

private:

	PairFullLowDimTwin m_pairFullLowDimTwin;
	MANIFELEM m_manifElem;
	VERTEXTYP m_centerVrtx;
	PairVrtcs m_shiftVrtcs;

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
	bool checkIntegrityFaceInVol( FullLowDimTwin const & fldt )
	{
		if( ! VolumeContains(fldt.spuckFullDimElem(), m_manifElem ) )
		{
			UG_LOG("Face not in Vol " << std::endl);
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
	bool checkIntegrityFaceInBothVols()
	{
		if(    ! checkIntegrityFaceInVol( m_pairFullLowDimTwin.first.spuckVol())
			|| ! checkIntegrityFaceInVol( m_pairFullLowDimTwin.second.spuckVol())
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
		Edge * edgeOne, edgeTwo;

		edgeOne = m_pairFullLowDimTwin.first.spuckLowDimElem();
		edgeTwo = m_pairFullLowDimTwin.first.spuckLowDimElem();

		Vertex * centerVrtx, shiftVrtxOne, shiftVrtxTwo;
		centerVrtx = nullptr;
		shiftVrtxOne = nullptr;
		shiftVrtxTwo = nullptr;

		if( ! findConnectingAndExtrnlVertex() )
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

		return true;
	}

};



/////////////////////////////////////////////////////////////////


template <
typename FULLDIMELEM,
typename MANIFELEM,
typename LOWDIMELEM,
typename VERTEXTYP,
typename = std::enable_if< std::is_pointer<FULLDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<MANIFELEM>::value>,
typename = std::enable_if< std::is_pointer<LOWDIMELEM>::value>,
typename = std::enable_if< std::is_pointer<VERTEXTYP>::value>
>
class ElemsToBeQuenched4DiamSpace
{
public:

	using FullLowDimManifQntpl = FullLowDimManifQuintuplet<FULLDIMELEM,MANIFELEM,LOWDIMELEM,VERTEXTYP>;
	using VecFullLowDimManifQuintuplet = std::vector<FullLowDimManifQntpl>;

	ElemsToBeQuenched4DiamSpace( VecFullLowDimManifQuintuplet const & vfldm5 )
	: m_centerVrtx(nullptr), m_vecFullLowDimManifQuintpl(vfldm5)
	{}

	bool checkIntegrity()
	{
		bool centerAssigned = false;

		for( auto & fldmq : m_vecFullLowDimManifQuintpl )
		{
			if( ! fldmq.checkIntegrity() )
			{
				UG_LOG("fulllowdim manif 5 not integer " << std::endl);
				return false;
			}

			if( ! centerAssigned )
			{
				fldmq.spuckCenterVertex( m_centerVrtx );
				centerAssigned = true;
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

//			if( fldmq.spuckCenterVertex() != m_centerVrtx )
//			{
//				UG_LOG("Center vertex not identical " << std::endl);
//				return false;
//			}
		}

		if( ! centerAssigned || m_centerVrtx == nullptr )
		{
			UG_LOG("CEnter problem " << std::endl);
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

		return true;
	}

	void spuckCenterVertex( VERTEXTYP & center )
	{
		center = m_centerVrtx;
	}

private:

	VERTEXTYP m_centerVrtx;
	VecFullLowDimManifQuintuplet m_vecFullLowDimManifQuintpl;

};


/////////////////////////////////////////////////////////////////


} // end of namespace diamonds

} // end of namespace ug





#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDINFO_H_ */
