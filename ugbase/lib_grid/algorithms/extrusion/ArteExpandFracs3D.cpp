/*
 * ArteExpandFracs3D.cpp
 *
 *  Created on: 06.10.2024
 *      Author: Markus M. Knodel
 *
 *  * expand fractures using the Arte algorithm, 3D case
 *
 * Author: Markus Knodel, inspired by Arte from Fuchs and Sebastian Reiters code for fracture expansion without Arte
 *
 * implementing a class that gives the basic tools for Arte in 3D
 * might be templated at a later stage to fulfill 2D and 3D purposes, if suitable
 * ( so far 2D case one entire function, not so perfect, but running....)
 *
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <boost/function.hpp>

#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"

#include <stack>
#include <utility>
#include <vector>
#include <type_traits>
#include <limits>
#include <atomic>
#include <cstddef>
#include <bitset>
#include <string>
#include <cmath>

#include "support.h"
#include "support3D.h"


#include <lib_grid/algorithms/extrusion/ArteExpandFracs3D.h>

#include "simpleMatrixOps.h"



namespace ug
{

ArteExpandFracs3D::ArteExpandFracs3D(
		Grid & grid, SubsetHandler & sh,
	    std::vector<FractureInfo> const & fracInfos,
		bool useTrianglesInDiamonds, bool establishDiamonds )
	: m_grid(grid),
	  m_sh(sh),
	  m_fracInfos(fracInfos),
	  m_useTrianglesInDiamonds(useTrianglesInDiamonds),
	  m_establishDiamonds(establishDiamonds),
	  m_aaPos(Grid::VertexAttachmentAccessor<APosition>()),
//	  m_facDescr(FaceDescriptor()),
//	  m_volDescr(VolumeDescriptor()),
	  m_fracInfosBySubset(std::vector<FractureInfo>()),
	  //m_sel(Selector()),
	  m_aAdjMarkerVFP(AttVertFracProp()),
	  m_aaMarkVrtVFP( Grid::VertexAttachmentAccessor<AttVertFracProp>()),
//	  m_aAdjMarkerVFP(AttVertFracProp()),
	  m_aaMarkEdgeVFP(Grid::EdgeAttachmentAccessor<AttVertFracProp>()),
	  m_aAdjMarkerFaceIsFracB(ABool()),
	  m_aaMarkFaceIsFracB(Grid::FaceAttachmentAccessor<ABool>()),
//	  m_aaMarkFaceHasUnclosedFracSideB(Grid::FaceAttachmentAccessor<ABool>()),
//	  m_aAdjMarkerVrtxHasUnclosedFracB(ABool()),
//	  m_aaMarkVrtxHasUnclosedFracB(Grid::VertexAttachmentAccessor<ABool>()),
	  m_aAdjMarkerFaceWithEndingCrossingCleft(ABool()),
	  m_aaMarkFaceWithEndingCrossingCleft(Grid::FaceAttachmentAccessor<ABool>()),
	  m_aAdjMarkerVrtxAtEndingCrossingCleft(ABool()),
	  m_aaMarkVrtxAtEndingCrossingCleft(Grid::VertexAttachmentAccessor<ABool>()),
//	  m_aAdjMarkerVrtx2AtInnerEndOfEndingCrossingFract(ABool()),
//	  m_aaMarkVrtx2AtInnerEndOfEndingCrossingFract(Grid::VertexAttachmentAccessor<ABool>()),
	  m_needToSplitEdgesConnectingNeighbrdEndingCrossCleftVrtx(false),
	  m_originalFractureFaces(std::vector<Face*>()),
//	  m_attVrtVec(AttVrtVec()),
//	  m_aaVrtVecVol( Grid::VolumeAttachmentAccessor<AttVrtVec>() ),
	  m_aAdjInfoEdges(AttVecEdge()),
	  m_aAdjInfoFaces(AttVecFace()),
//	  m_aAdjInfoVols(AttVecVol()),
	  m_aaVrtInfoAssoEdges( Grid::VertexAttachmentAccessor<AttVecEdge>()),
	  m_aaVrtInfoAssoFaces( Grid::VertexAttachmentAccessor<AttVecFace>()),
//	  m_aaVrtInfoAssoVols( Grid::VertexAttachmentAccessor<AttVecVol>()),
//	  m_aAdjInfoAVVFT( AttVecVertFracTrip() ),
//	  m_aaVrtInfoFraTri(Grid::VertexAttachmentAccessor<AttVecVertFracTrip>()),
//	  m_vrtxFractrQuadrplVec(VrtxFractrQuadrplArte3DVec())
	  m_attVrtVec(AttVrtVec()),
	  m_aaVrtVecVol( Grid::VolumeAttachmentAccessor<AttVrtVec>() ),
	  m_vecCrossVrtInf(std::vector<CrossVertInf>()),
//	  m_aAdjVolElmInfo(AttVecAttachedVolumeElemInfo()),
//	  m_aaVolElmInfo(Grid::VertexAttachmentAccessor<AttVecAttachedVolumeElemInfo>()),
	  m_attAdjVecSegVolElmInfo( AttVecSegmentVolElmInfo() ),
	  m_accsAttVecSegVolElmInfo( Grid::VertexAttachmentAccessor<AttVecSegmentVolElmInfo>() ),
	  m_attVecSegmLimSid(AttVecSegmLimSid()),
	  m_vrtxAttAccsVecSegmLimSid(Grid::VertexAttachmentAccessor<AttVecSegmLimSid>()),
	  m_vecEdgeDirectConnectingEndingCrossCleftVrtcs( std::vector<Edge* >() ),
	  m_vecEndCrossFractSegmInfo( VecEndingCrossingFractureSegmentInfo() ),
//	  m_vecEndCrossFractSegmInfo(VecEndingCrossingFractureSegmentInfo()),
	  m_attAtVrtxVecEndingCrossFractSegmInfo(AttVecEndingCrossingFractureSegmentInfo()),
	  m_vrtxAttAccsVecEndingCrossFractSegmInfo(Grid::VertexAttachmentAccessor<AttVecEndingCrossingFractureSegmentInfo>()),
//	  m_attAtVrtxIfVrtxIsEndingCrossingCleftVrtx(ABool()),
//	  m_vrtxAttAccsVrtxIsEndingCrossingCleftVrtx(Grid::VertexAttachmentAccessor<ABool>()),
	  m_attAtFaceIfFaceIsSegmLimFaceEndingCrossingCleft(ABool()),
	  m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft(Grid::FaceAttachmentAccessor<ABool>()),
	  m_attAtVolIfVolTouchesEndingCrossingCleft(ABool()),
	  m_volAttAccsVolTouchesEndingCrossingCleft(Grid::VolumeAttachmentAccessor<ABool>()),
	  m_attAtVrtxIfVrtxArisesFromExpandedEndingCrossingCleft(ABool()),
	  m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft(Grid::VertexAttachmentAccessor<ABool>()),
	  m_vrtxArisesFromExpandedEndingCrossingCleft(std::vector<Vertex*>())
{
//	// Notloesung, nicht in die erste Initialisierung vor geschweifter Klammer, da copy constructor privat
	m_sel = Selector();
}


ArteExpandFracs3D::~ArteExpandFracs3D()
{
	//  Auto-generated destructor stub
}

bool ArteExpandFracs3D::run( bool & needToRestart )
{
	if( ! initialize() )
		return false;

	UG_LOG("initialisiert" << std::endl);

	if( ! setSelector() )
		return false;

	UG_LOG("selektiert" << std::endl);

	if( ! attachMarkers() )
		return false;

	UG_LOG("attached" << std::endl);

	int splittedEdges = splitInnerFreeFracEdgs();

	UG_LOG("splitted edges " << splittedEdges << std::endl);
//
//	return false;
//
//	UG_LOG("Splitted inner free frac edges" << splittedEdges << std::endl);

	if( ! countAndSelectFracBaseNums() )
		return false;

	UG_LOG("gezaehlt" << std::endl);

//	constexpr bool assignFracInfosFirst = true;
//
//	if( assignFracInfosFirst )
//	{
	if( ! assignOrigFracInfos() )
		return false;

	UG_LOG("assigniert zuerst " << std::endl);

	if( ! enableVolOptAutoGenFac() )
	{
		UG_LOG("autogen war schon eingestellt" << std::endl);
	}
	else
	{
		UG_LOG("Autogen einstellen" << std::endl);
	}
	// TODO FIXME für was gebraucht????

//	}

//	constexpr bool generVrtInfoFirst = true;
//
//	if( generVrtInfoFirst )
//	{
//	if( ! generateVertexInfos() )
//		return false;
//
//	UG_LOG("generiert zuerst " << std::endl);

//	}

//	if( ! prepareStasi() )
//		return false;
//
//	UG_LOG("Stasi vorbereitet " << std::endl);

	if( ! distinguishSegments() )
		return false;

	UG_LOG("Segmente erzeugt " << std::endl);


	UG_LOG("check unclosed frac faces" << std::endl);

	if( ! establishSegmentLimitingSidesInfo() )
		return false;

	UG_LOG("established Segment limiting sides info" << std::endl);

	if( ! detectEndingCrossingCleftsSegmBased() )
	{
		if( m_needToSplitEdgesConnectingNeighbrdEndingCrossCleftVrtx )
		{
			needToRestart = true;
			UG_LOG("Restart due to splitted edges" << std::endl);
		}
		else
		{
			needToRestart = false;
			UG_LOG("No need to restart due to splitted edges" << std::endl);
		}

		return false;
	}

	UG_LOG("Ending crossing clefts detected" << std::endl);

	if( ! seletForSegmented() )
		return false;

	UG_LOG("Closed Open Vertex Fract untersucht " << std::endl);

//	if( ! assignFracInfosFirst )
//	{
//		if( ! assignOrigFracInfos() )
//			return false;
//
//		UG_LOG("assigniert danach " << std::endl);
//
//	}


	if( ! establishNewVrtBase() )
		return false;

	UG_LOG("etabliert" << std::endl);

//	if( ! generVrtInfoFirst )
//	{
//
//		if( ! generateVertexInfos() )
//			return false;
//
//		UG_LOG("generiert danach " << std::endl);
//	}

	if( ! createConditionForNewVrtcs() )
		return false;

	UG_LOG("kreiert" << std::endl);

	if( ! loop2EstablishNewVertices() )
		return false;

	UG_LOG("loopiert" << std::endl);

	UG_LOG("under construction " << std::endl);

	if( ! createNewElements() )
		return false;

	UG_LOG("new elements created " << std::endl);

	if( ! detachMarkers() )
		return false;

	UG_LOG("detachiert" << std::endl);

//	assignDebugSubsets( false );

	UG_LOG("Debug subsets assigned" << std::endl );

	IndexType numDelEndCrossCleftFacs = deleteEndingCrossingCleftOrigFacs();

	UG_LOG("deleted ending crossing cleft changed faces and direction edges " << numDelEndCrossCleftFacs << std::endl);

	return true;
}

////////////////////////////////////////////////

bool ArteExpandFracs3D::splitEdgesOfNeighboredEndingCrossingFracVrtcs()
{
//	IndexType numberSplittedEdges = 0;

	// TODO FIXME vergleichen, was bei den debug Ecken alles an Randvertizes ist, und was bei
	// den hochstehenden Edges für Randvertizes sind, ob da was verloren gegangen ist
	// oder ob bei den modernen SegmentEnding Infos was fehlt an Vertizes.....
	// angefangen mit dem simplen Beispiel!
	// wenn das klappt, können die debug m_d_ Variablen entfernt werden.....

	UG_LOG("NUMBER EDGES DIRECT " << m_vecEdgeDirectConnectingEndingCrossCleftVrtcs.size() << std::endl );

	UG_LOG("NUMBER VERTICES FRAC SEGMENTS " << m_vecEndCrossFractSegmInfo.size() << std::endl );

	for( Edge * edg : m_vecEdgeDirectConnectingEndingCrossCleftVrtcs )
	{

		UG_LOG("trying to split edge " << CalculateCenter( edg, m_aaPos ) << std::endl);

		UG_LOG("E VERTEX ONE " << m_aaPos[ edg->vertex(0) ] << std::endl);
		UG_LOG("E VERTEX TWO " << m_aaPos[ edg->vertex(1) ] << std::endl);

		for( EndingCrossingFractureSegmentInfo const & ecfsiOne : m_vecEndCrossFractSegmInfo )
		{
			Vertex * vrtOne = ecfsiOne.spuckUnclosedVrtx();

			UG_LOG("V ONE " << m_aaPos[vrtOne] << std::endl);

			for( EndingCrossingFractureSegmentInfo const & ecfsiTwo : m_vecEndCrossFractSegmInfo )
			{
				Vertex * vrtTwo = ecfsiTwo.spuckUnclosedVrtx();

				UG_LOG("V TWO " << m_aaPos[vrtTwo] << std::endl);

				if( vrtOne != vrtTwo )
				{
					if( EdgeContains(edg, vrtOne) && EdgeContains( edg, vrtTwo) )
					{
						UG_LOG("Edge needs to be splitted" << std::endl);

						int suse = m_sh.num_subsets();

						m_sh.assign_subset( edg, suse );

						UG_LOG("need to detach markers" << std::endl);
						detachMarkers();

						vector3 center = CalculateCenter(edg, m_aaPos);
						UG_LOG("splitting ECCV edge at " << center << std::endl);
						RegularVertex* vrtSE = SplitEdge<RegularVertex>(m_grid, edg, false);
						m_aaPos[vrtSE] = center;

						UG_LOG("Edge splitted, please restart process with the thus changed geometry" << std::endl);

						return true;

					}

				}

			}


		}
	}

	return false;

//	if( numberSplittedEdges > 0 )
//		UG_LOG("Edge splitted, please restart process with the thus changed geometry, new edges " << numberSplittedEdges << std::endl);
//
//	return numberSplittedEdges;

	// TODO FIXME das ganze nochmal mit dem Debug Zeug

//	UG_LOG("DEBUG" << std::endl);
//
//	int suse = m_sh.num_subsets();
//
//	for( Edge * edg : m_d_allContributingEdges )
//	{
//		UG_LOG("E VERTEX ONE " << m_aaPos[ edg->vertex(0) ] << std::endl);
//		UG_LOG("E VERTEX TWO " << m_aaPos[ edg->vertex(1) ] << std::endl);
//
//		for( Vertex * vrtOne : m_d_endingCrossingCleftVrtcs )
//		{
//			UG_LOG("V ONE " << m_aaPos[vrtOne] << std::endl);
//
//			for( Vertex * vrtTwo : m_d_endingCrossingCleftVrtcs )
//			{
//
//				UG_LOG("V Two " << m_aaPos[vrtTwo] << std::endl);
//
//				if( vrtOne != vrtTwo )
//				{
//					if( EdgeContains(edg, vrtOne) && EdgeContains( edg, vrtTwo) )
//					{
//						UG_LOG("Edge needs to be splitted" << std::endl);
//
//						m_sh.assign_subset( edg, suse );
//
//						vector3 center = CalculateCenter(edg, m_aaPos);
//						UG_LOG("splitting ECCV edge at " << center << std::endl);
//						RegularVertex* vrtSE = SplitEdge<RegularVertex>(m_grid, edg, false);
//						m_aaPos[vrtSE] = center;
//
//						UG_LOG("Edge splitted, please restart process with the thus changed geometry" << std::endl);
//
//						return 1;
//
//					}
//				}
//			}
//
//		}
//	}
//
//	return 0;
//
//	return false;

}

////////////////////////////////////////////////

void ArteExpandFracs3D::assignDebugSubsets( bool intermediate )
{
//	return;

	std::vector<Face*> d_endingCrossingCleftFaces;
	std::vector<Face*> d_endingCrossingCleftFacesNoCut;
	std::vector<Vertex*> d_endingCrossingCleftVrtcs;
	std::vector<Edge*> d_cuttingEdges;
	std::vector<Edge*> d_shiftEdges;
	std::vector<Face*> d_crossingNeighboredNotEndingFaces;
	std::vector<Face*> d_crossingNeighboredNotEndingFacesCommEdg;
	std::vector<Face*> d_notEndingCrossingFacesNotNeighbour;
	std::vector<Volume*> d_vols;
	std::vector<Vertex*> d_shiftVrtcs;

	UG_LOG("Number ending crossing clefts " << m_vecEndCrossFractSegmInfo.size()  << std::endl);

	for( EndingCrossingFractureSegmentInfo const & ecfsi : m_vecEndCrossFractSegmInfo )
	{
		Vertex * vrt = ecfsi.spuckUnclosedVrtx();

		d_endingCrossingCleftVrtcs.push_back(vrt);

		std::vector<Face*> vecClosFracFacNoNeig = ecfsi.spuckVecClosedFracManifElNoNeighbr();

		for( Face * fac : vecClosFracFacNoNeig )
		{
			d_crossingNeighboredNotEndingFaces.push_back(fac);
		}

		Face * endingFacCut = ecfsi.spuckEndingFractManifCutting();

		d_endingCrossingCleftFaces.push_back(endingFacCut);

		Face * endingFacNoCut = ecfsi.spuckEndingFractManifNotCutting();

		if( endingFacNoCut != nullptr )
		{
			d_endingCrossingCleftFacesNoCut.push_back(endingFacNoCut);
		}

		Edge * cutEdge = ecfsi.spuckOldLowDimElCut();

		d_cuttingEdges.push_back(cutEdge);

		std::pair<Face*,Face*> const & neighrdFacsClos = ecfsi.spuckPairNeighbouredFractClosedManifEl();

		d_crossingNeighboredNotEndingFacesCommEdg.push_back(neighrdFacsClos.first);
		d_crossingNeighboredNotEndingFacesCommEdg.push_back(neighrdFacsClos.second);

		std::vector<Volume*> vols = ecfsi.spuckVecFulldimEl();

		for( Volume * v : vols )
		{
			d_vols.push_back(v);
		}

		Edge * shiftEdge = ecfsi.spuckLowdimElmShiftDirection();
		d_shiftEdges.push_back(shiftEdge);


		if( ! intermediate )
		{
			Vertex * shiVe = ecfsi.spuckShiftVrtx();

			if( shiVe == nullptr )
			{
				m_sh.assign_subset(vrt, m_sh.num_subsets());
				m_sh.assign_subset(shiftEdge, m_sh.num_subsets());

				UG_LOG("no shift vertex for ending crossing cleft segment " << std::endl);
				UG_THROW("no shift vertex for ending crossing cleft segment " << std::endl);

			}
			d_shiftVrtcs.push_back(shiVe);
		}



	}

	int suse = m_sh.num_subsets();

	if( ! intermediate )
	{

		for( Face * fac : d_endingCrossingCleftFaces )
		{
			if( fac != nullptr )
				m_sh.assign_subset( fac, suse );

		}

		suse = m_sh.num_subsets();

		for( Face * fac : d_endingCrossingCleftFacesNoCut )
		{
			if( fac != nullptr )
				m_sh.assign_subset( fac, suse );

		}

		suse = m_sh.num_subsets();

		for( Vertex * vrt : d_endingCrossingCleftVrtcs )
		{
			if( vrt != nullptr )
				m_sh.assign_subset( vrt, suse );
		}

		suse = m_sh.num_subsets();

		for( Edge * edg : d_cuttingEdges )
		{
			if( edg == nullptr )
			{
				UG_LOG("NULL UNERLAUBT" << std::endl);
				UG_THROW("NULL UNERLAUBT" << std::endl);
			}

			m_sh.assign_subset( edg, suse );

		}

		suse = m_sh.num_subsets();

		for( Edge * edg : d_shiftEdges )
		{
			if( edg == nullptr )
			{
				UG_LOG("NULL UNERLAUBT" << std::endl);
				UG_THROW("NULL UNERLAUBT" << std::endl);
			}

			m_sh.assign_subset( edg, suse );

		}


		suse = m_sh.num_subsets();

		for( Face * fac : d_crossingNeighboredNotEndingFaces )
		{
			if( fac != nullptr )
				m_sh.assign_subset( fac, suse );
		}

		suse = m_sh.num_subsets();

		for( Face * fac : d_notEndingCrossingFacesNotNeighbour )
		{
			m_sh.assign_subset( fac, suse );
		}

		suse = m_sh.num_subsets();

		for( Face * fac : d_crossingNeighboredNotEndingFacesCommEdg )
		{
			m_sh.assign_subset( fac, suse );
		}

		suse = m_sh.num_subsets();

		for( Vertex * vrtx : d_shiftVrtcs )
		{
			if( vrtx != nullptr )
				m_sh.assign_subset( vrtx, suse );
		}

//	suse = m_sh.num_subsets();
	}
	else
	{

		for( Volume * v : d_vols )
		{
			if( v != nullptr )
				m_sh.assign_subset( v, suse );
		}
	}

	return;

#if 0
//	if( numEndingCrossingClefts == 0 )
//		return true;

	// debug for ending crossing clefts

	suse = m_sh.num_subsets();

	for( Face * fac : m_d_endingCrossingCleftFaces )
	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
		m_sh.assign_subset( fac, suse );

	}

	suse = m_sh.num_subsets();

	for( Vertex * vrt : m_d_endingCrossingCleftVrtcs )
	{
//		m_sh.assign_subset( vrt, m_sh.num_subsets());
		m_sh.assign_subset( vrt, suse );
	}

	suse = m_sh.num_subsets();

	for( Edge * edg : m_d_cuttingEdges )
	{
		if( edg == nullptr )
		{
			UG_LOG("NULL UNERLAUBT" << std::endl);
			UG_THROW("NULL UNERLAUBT" << std::endl);
		}

//		m_sh.assign_subset( edg, m_sh.num_subsets());
		m_sh.assign_subset( edg, suse );

	}

	suse = m_sh.num_subsets();

	for( Face * fac : m_d_crossingNeighboredNotEndingFaces )
	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
		m_sh.assign_subset( fac, suse );
	}

//	for( Edge * edg : otherEdgeOfCrossingNotEndingFace )
//	{
//		if( edg == nullptr )
//		{
//			UG_LOG("NULL C UNERLAUBT" << std::endl);
//			UG_THROW("NULL C UNERLAUBT" << std::endl);
//		}
//
//		m_sh.assign_subset( edg, m_sh.num_subsets());
//
//	}
	suse = m_sh.num_subsets();

	for( Face * fac : m_d_notEndingCrossingFacesNotNeighbour )
	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
		m_sh.assign_subset( fac, suse );
	}

	suse = m_sh.num_subsets();

	for( Face * fac : m_d_crossingNeighboredNotEndingFacesCommEdg )
	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
		m_sh.assign_subset( fac, suse );
	}

//	suse = m_sh.num_subsets();
//
//	for( Edge * edg : m_d_allContributingEdges )
//	{
//		for( Vertex * vrtOne : m_d_endingCrossingCleftVrtcs )
//		{
//			for( Vertex * vrtTwo : m_d_endingCrossingCleftVrtcs )
//			{
//				if( vrtOne != vrtTwo )
//				{
//					if( EdgeContains(edg, vrtOne) && EdgeContains( edg, vrtTwo) )
//					{
//						UG_LOG("Edge needs to be splitted" << std::endl);
//
//						m_sh.assign_subset( edg, suse );
//
//
////						vector3 center = CalculateCenter(edg, m_aaPos);
////						UG_LOG("splitting ECCV edge at " << center << std::endl);
////						RegularVertex* vrtSE = SplitEdge<RegularVertex>(m_grid, edg, false);
////						m_aaPos[vrtSE] = center;
//
//					}
//				}
//			}
//
//		}
//	}


	return;
#endif

//	return false;
}

////////////////////////////////////////////////

bool ArteExpandFracs3D::initialize()
{
	UG_LOG("initialize " << std::endl);

	//	access position attachment
	if(!m_grid.has_vertex_attachment(aPosition) )
	{
		UG_LOG("Error in ExpandFractures Arte 3D: Missing position attachment");
		return false;
	}

	m_aaPos = Grid::VertexAttachmentAccessor<APosition>(m_grid, aPosition);

	//	make sure that the required options are enabled.

	if( ! m_grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES) )
	{
		UG_LOG("WARNING in Arte 3D init : grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
		m_grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}

	if( ! m_grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES) )
	{
		UG_LOG("WARNING in Arte 3D init: grid option FACEOPT_AUTOGENERATE_EDGES autoenabled.\n");
		m_grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	//	vectors that allow to access fracture properties by subset index
	m_fracInfosBySubset = std::vector<FractureInfo>( m_sh.num_subsets(), FractureInfo(-1, -1, 0) );

	for( size_t i = 0; i < m_fracInfos.size(); ++i)
	{
		if( m_fracInfos[i].subsetIndex >= m_sh.num_subsets())
		{
			throw(UGError("Bad subsetIndex in given fracInfos."));
		}

		m_fracInfosBySubset[ m_fracInfos[i].subsetIndex] = m_fracInfos[i];

	}


	return true;
}


bool ArteExpandFracs3D::setSelector()
{
	//	Collect surrounding volumes, faces and edges of all fractures in a selector
	//	and select fracture faces, edges and vertices too.

//	m_sel = Selector(m_grid);

	m_sel.assign_grid(m_grid);

	m_sel.enable_autoselection(false);
	m_sel.enable_selection_inheritance(true);	//required for select and mark, disabled later
	m_sel.enable_strict_inheritance(false);

//	bool strictInherit = m_sel.strict_inheritance_enabled();
//
//	UG_LOG("strikte Inheritenz ist " << strictInherit << " und false ist " << false << std::endl);

	return true;
}


bool ArteExpandFracs3D::attachMarkers()
{
	// first part

	// attachment pair boundary is fracture, number fractures crossing

	m_aAdjMarkerVFP = AttVertFracProp();

//	support::VertexFracturePropertiesVol<IndexType> vfp0; // false, 0 );
	VertxFracPropts vfp0; // false, 0 );
	// default value: no boundary fracture, no fractures crossing

	m_grid.attach_to_vertices_dv( m_aAdjMarkerVFP, vfp0 );
	m_aaMarkVrtVFP = Grid::VertexAttachmentAccessor<AttVertFracProp> ( m_grid, m_aAdjMarkerVFP );

	//m_aAdjMarkerVFP = AttVertFracProp();
	m_grid.attach_to_edges_dv( m_aAdjMarkerVFP, vfp0 );

	m_aaMarkEdgeVFP = Grid::EdgeAttachmentAccessor<AttVertFracProp>( m_grid, m_aAdjMarkerVFP );

	m_aAdjMarkerFaceIsFracB = ABool(); // used to know if an face is frac face

	m_grid.attach_to_faces_dv( m_aAdjMarkerFaceIsFracB, false );
	m_aaMarkFaceIsFracB = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerFaceIsFracB );

//	m_aAdjMarkerFaceHasUnclosedFracSideB = ABool();
//
//	m_grid.attach_to_faces_dv( m_aAdjMarkerFaceHasUnclosedFracSideB, false );
//	m_aaMarkFaceHasUnclosedFracSideB = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerFaceHasUnclosedFracSideB );
//
//	m_aAdjMarkerVrtxHasUnclosedFracB = ABool();
//
//	m_grid.attach_to_vertices_dv( m_aAdjMarkerVrtxHasUnclosedFracB, false );
//	m_aaMarkVrtxHasUnclosedFracB = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerVrtxHasUnclosedFracB );
//
	m_aAdjMarkerFaceWithEndingCrossingCleft = ABool();

	m_grid.attach_to_faces_dv( m_aAdjMarkerFaceWithEndingCrossingCleft, false );

	m_aaMarkFaceWithEndingCrossingCleft = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerFaceWithEndingCrossingCleft );

	m_aAdjMarkerVrtxAtEndingCrossingCleft = ABool();

	m_grid.attach_to_vertices_dv( m_aAdjMarkerVrtxAtEndingCrossingCleft, false );

	m_aaMarkVrtxAtEndingCrossingCleft = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerVrtxAtEndingCrossingCleft );
//
//	m_aAdjMarkerVrtx2AtInnerEndOfEndingCrossingFract = ABool();
//
//	m_grid.attach_to_vertices_dv( m_aAdjMarkerVrtx2AtInnerEndOfEndingCrossingFract, false );
//
//	m_aaMarkVrtx2AtInnerEndOfEndingCrossingFract = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_aAdjMarkerVrtx2AtInnerEndOfEndingCrossingFract );

	// second part

//	m_grid.attach_to_volumes(m_attVrtVec);
//	m_aaVrtVecVol = Grid::VolumeAttachmentAccessor<AttVrtVec>( m_grid, m_attVrtVec);

	std::vector<Edge*> noEdge;
	std::vector<Face*> noFace;
//	std::vector<Volume*> noVol;

	m_aAdjInfoEdges = AttVecEdge();
	m_aAdjInfoFaces = AttVecFace();
//	m_aAdjInfoVols = AttVecVol();

	m_grid.attach_to_vertices_dv( m_aAdjInfoEdges, noEdge );
	m_aaVrtInfoAssoEdges = Grid::VertexAttachmentAccessor<AttVecEdge>( m_grid, m_aAdjInfoEdges );

	m_grid.attach_to_vertices_dv( m_aAdjInfoFaces, noFace );
	m_aaVrtInfoAssoFaces = Grid::VertexAttachmentAccessor<AttVecFace>( m_grid, m_aAdjInfoFaces );

//	m_grid.attach_to_vertices_dv( m_aAdjInfoVols, noVol );
//	m_aaVrtInfoAssoVols = Grid::VertexAttachmentAccessor<AttVecVol>( m_grid, m_aAdjInfoVols );


	//  TODO FIXME
	//  das fehlt hier , Analogon 2D Fall!!!!!!!!!!!!! der geht hier eigentlich weiter
	// die Vertizes, Faces und Edges, die mit einer Kluft zu tun haben
	//	using VertFracTrip = VertexFractureTriple<Edge*, Face*, vector3>;
	//	using VecVertFracTrip = std::vector<VertFracTrip>;
	//	VecVertFracTrip vertexNoInfo;

	// AttVecVertFracTrip m_aAdjInfoAVVFT;

//	VecVertFracTrip vertexNoInfo;
//
//	m_aAdjInfoAVVFT = AttVecVertFracTrip();
//
//	m_grid.attach_to_vertices_dv( m_aAdjInfoAVVFT, vertexNoInfo );
//
//	m_aaVrtInfoFraTri = Grid::VertexAttachmentAccessor<AttVecVertFracTrip>(m_grid,  m_aAdjInfoAVVFT );


	//	associate a vector of vertices for each volume adjacent to the frac.
	//	An entry later will contain the new vertex, if the
	//	corresponding vertex is an inner fracture vertex, and nullptr if not.

	m_attVrtVec = AttVrtVec();

	m_grid.attach_to_volumes(m_attVrtVec);

	m_aaVrtVecVol = Grid::VolumeAttachmentAccessor<AttVrtVec>(m_grid, m_attVrtVec);


//	VecAttachedVolumeElemInfo noVolInfo;

//	m_aAdjVolElmInfo = AttVecAttachedVolumeElemInfo();
//
//	m_grid.attach_to_vertices_dv(m_aAdjVolElmInfo,noVolInfo);
//
//	m_aaVolElmInfo = Grid::VertexAttachmentAccessor<AttVecAttachedVolumeElemInfo>(m_grid, m_aAdjVolElmInfo);


	VecSegmentVolElmInfo noSegmts;

	m_attAdjVecSegVolElmInfo = AttVecSegmentVolElmInfo();

	m_grid.attach_to_vertices_dv( m_attAdjVecSegVolElmInfo, noSegmts );

	m_accsAttVecSegVolElmInfo = Grid::VertexAttachmentAccessor<AttVecSegmentVolElmInfo>( m_grid, m_attAdjVecSegVolElmInfo );

	VecSegmentLimitingSides emptyVecSegmLimSid;

	AttVecSegmLimSid m_attVecSegmLimSid;

	m_grid.attach_to_vertices_dv( m_attVecSegmLimSid, emptyVecSegmLimSid );

	m_vrtxAttAccsVecSegmLimSid = Grid::VertexAttachmentAccessor<AttVecSegmLimSid>( m_grid, m_attVecSegmLimSid );

	VecEndingCrossingFractureSegmentInfo emptyVECFSI;

	m_attAtVrtxVecEndingCrossFractSegmInfo = AttVecEndingCrossingFractureSegmentInfo();

	m_grid.attach_to_vertices_dv( m_attAtVrtxVecEndingCrossFractSegmInfo, emptyVECFSI );

	m_vrtxAttAccsVecEndingCrossFractSegmInfo = Grid::VertexAttachmentAccessor<AttVecEndingCrossingFractureSegmentInfo>( m_grid, m_attAtVrtxVecEndingCrossFractSegmInfo );

//	m_attAtVrtxIfVrtxIsEndingCrossingCleftVrtx = ABool();
//
//	m_grid.attach_to_vertices_dv( m_attAtVrtxIfVrtxIsEndingCrossingCleftVrtx, false );
//
//	m_vrtxAttAccsVrtxIsEndingCrossingCleftVrtx = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_attAtVrtxIfVrtxIsEndingCrossingCleftVrtx );

	m_attAtFaceIfFaceIsSegmLimFaceEndingCrossingCleft = ABool();

	m_grid.attach_to_faces_dv( m_attAtFaceIfFaceIsSegmLimFaceEndingCrossingCleft, false );

	m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_attAtFaceIfFaceIsSegmLimFaceEndingCrossingCleft );

	m_attAtVolIfVolTouchesEndingCrossingCleft = ABool();

	m_grid.attach_to_volumes_dv( m_attAtVolIfVolTouchesEndingCrossingCleft, false );

	m_volAttAccsVolTouchesEndingCrossingCleft = Grid::VolumeAttachmentAccessor<ABool>( m_grid, m_attAtVolIfVolTouchesEndingCrossingCleft );

	m_attAtVrtxIfVrtxArisesFromExpandedEndingCrossingCleft = ABool();

	m_grid.attach_to_vertices_dv( m_attAtVrtxIfVrtxArisesFromExpandedEndingCrossingCleft, false );

	m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_attAtVrtxIfVrtxArisesFromExpandedEndingCrossingCleft );

	return true;
}


bool ArteExpandFracs3D::detachMarkers()
{
	m_grid.detach_from_vertices( m_aAdjMarkerVFP );
	m_grid.detach_from_edges( m_aAdjMarkerVFP );
	m_grid.detach_from_faces( m_aAdjMarkerFaceIsFracB );
//	m_grid.detach_from_faces( m_aAdjMarkerFaceHasUnclosedFracSideB );
//
//	m_grid.detach_from_vertices( m_aAdjMarkerVrtxHasUnclosedFracB );
//
	m_grid.detach_from_faces( m_aAdjMarkerFaceWithEndingCrossingCleft );
	m_grid.detach_from_vertices( m_aAdjMarkerVrtxAtEndingCrossingCleft );
//	m_grid.detach_from_vertices( m_aAdjMarkerVrtx2AtInnerEndOfEndingCrossingFract );

	m_grid.detach_from_vertices( m_aAdjInfoEdges );
	m_grid.detach_from_vertices( m_aAdjInfoFaces );
//	m_grid.detach_from_vertices( m_aAdjInfoVols );

//	m_grid.detach_from_vertices( m_aAdjInfoAVVFT  );

	m_grid.detach_from_volumes( m_attVrtVec );

//	m_grid.detach_from_vertices(m_aAdjVolElmInfo);

	m_grid.detach_from_vertices( m_attAdjVecSegVolElmInfo );

	m_grid.detach_from_vertices( m_attVecSegmLimSid );

	m_grid.detach_from_vertices( m_attAtVrtxVecEndingCrossFractSegmInfo );

//	m_grid.detach_from_vertices( m_attAtVrtxIfVrtxIsEndingCrossingCleftVrtx );

	m_grid.detach_from_faces( m_attAtFaceIfFaceIsSegmLimFaceEndingCrossingCleft );

	m_grid.detach_from_volumes( m_attAtVolIfVolTouchesEndingCrossingCleft );

	m_grid.detach_from_vertices( m_attAtVrtxIfVrtxArisesFromExpandedEndingCrossingCleft );

	return true;
}

bool ArteExpandFracs3D::enableVolOptAutoGenFac()
{
	// brauchen wir das? für was? von SR irgendwie übernommen, wo dort was entfernt ähnliches gemacht wird....
	if(! m_grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES) )
	{
		UG_LOG("WARNING grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
		m_grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
		return true;
	}

	return false;

}

//////////////////////////////////////////////////////////////////

int ArteExpandFracs3D::splitInnerFreeFracEdgs()
{
	int splittedEdges = 0;

	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		UG_LOG("sudo ind " << fracIndSudo << std::endl);

		AttVertFracProp func_aAdjMarkerVFP;
		VertxFracPropts vfp0;
		m_grid.attach_to_edges_dv( func_aAdjMarkerVFP, vfp0 );
		Grid::EdgeAttachmentAccessor<AttVertFracProp> func_aaMarkEdgeVFP( m_grid, func_aAdjMarkerVFP );

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

			std::vector<Edge*> facEdges;

			CollectEdges( facEdges, m_grid, fac );

			for( auto const & edg : facEdges )
			{
				if( func_aaMarkEdgeVFP[edg].getNumberFracEdgesInVertex() != 0 )
					UG_THROW("Attachment nicht auf default " << std::endl);
//				func_aaMarkEdgeVFP[edg].setNumberCrossingFracsInVertex(0);
			}
		}


		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

			std::vector<Edge*> facEdges;

			CollectEdges( facEdges, m_grid, fac );

			for( auto const & edg : facEdges )
			{
				func_aaMarkEdgeVFP[edg]++;
			}
		}

		std::vector<Edge *> edgesToBeSplitted;

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

			std::vector<Edge*> facEdges;

			CollectEdges( facEdges, m_grid, fac );

			IndexType openEdges = 0;

			for( auto const & edg : facEdges )
			{
				VertxFracPropts & edgeFracPrps = func_aaMarkEdgeVFP[edg];

				IndexType fracEdgesOverlap = edgeFracPrps.getNumberFracEdgesInVertex();

				if( fracEdgesOverlap == 1 )
				{
					openEdges++;
				}
				else if( fracEdgesOverlap == 2 )
				{
					; // fine, inner edge of fracture
				}
				else
				{
					m_sh.assign_subset(edg, m_sh.num_subsets());

					UG_LOG("how many fractures at this edge " << fracEdgesOverlap << std::endl);
					UG_LOG("sudo " << fracIndSudo << std::endl);
					UG_THROW("how many fractures at this edge " << fracEdgesOverlap << std::endl);
					UG_THROW("sudo " << fracIndSudo << std::endl);
//					return splittedEdges;
//					UG_THROW("how many fractures at this edge " << fracEdgesOverlap << std::endl);
				}
			}

			if( openEdges == 2 )
			{
				// figure out that edge that is not open, this must be splitted

				IndexType innerEdges = 0;

				for( auto const & edg : facEdges )
				{
					VertxFracPropts & edgeFracPrps = func_aaMarkEdgeVFP[edg];

					IndexType fracEdgesOverlap = edgeFracPrps.getNumberFracEdgesInVertex();

					if( fracEdgesOverlap == 2 )
					{
						edgesToBeSplitted.push_back(edg);
						innerEdges++;
					}
				}

				if( innerEdges != 1 )
				{
					UG_LOG("inner edge number strange " << innerEdges << std::endl);
					UG_THROW("inner edge number strange " << innerEdges << std::endl);
				}


			}
		}

		for( Edge * edg : edgesToBeSplitted )
		{
			vector3 center = CalculateCenter(edg, m_aaPos);
			UG_LOG("splitting edge at " << center << std::endl);
			RegularVertex* vrt = SplitEdge<RegularVertex>(m_grid, edg, false);
			m_aaPos[vrt] = center;
			splittedEdges++;
		}

		m_grid.detach_from_edges( func_aAdjMarkerVFP );
	}

	return splittedEdges;
}

#if 0
int ArteExpandFracs3D::splitInnerFreeFracEdgs()
{
	int splittedEdges = 0;

	UG_LOG("search inner frac edges with two inner boundary edges" << std::endl);

	// TODO FIXME


	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		UG_LOG("sudo ind " << fracIndSudo << std::endl);

//		std::vector<Edge *> edgesToBeSplitted;

		std::vector<Face*> facesAtInnerBoundary;

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

//			std::vector<Vertex*> assoVrt;
//
//			CollectVertices( assoVrt, m_grid, fac );
//
//			// check number of edges which have only one fracture face around
//
//			for( Vertex * vrt : assoVrt )
//			{
			std::vector<Edge*> assoEdg;

			CollectEdges( assoEdg, m_grid, fac );

			for( Edge * edg : assoEdg )
			{
				// check if edge is inner boundary edge, i.e. if has at one side non-fracture face
				// only check for those edges which have one fracture face on one side

				// split only inner edges at the moment in case of problems
				if( IsBoundaryEdge3D(m_grid, edg) )
				{
					continue;
				}

				std::vector<Face *> assoFac;

				CollectFaces( assoFac, m_grid, edg );

				IndexType numFacsFromFracSudo = 0;

				for( Face * testFac : assoFac )
				{
					IndexType testSudo = m_sh.get_subset_index(testFac);

					if( testSudo == fracIndSudo )
						numFacsFromFracSudo++;
				}

				if( numFacsFromFracSudo == 0 )
				{
//					UG_LOG("no facs at edg " << std::endl );
						; // nothing to do, belongs from an edge not relevant for the fracture
				}
				else if( numFacsFromFracSudo == 1 )
				{
					// relevant! edge needs to be split

//					UG_LOG("one fac at edg" << std::endl);

					addElem( facesAtInnerBoundary, fac );
				}
				else if( numFacsFromFracSudo == 2 )
				{
//						UG_LOG("two fac sides at edge " << std::endl);
						; // nothing to do, edge at both sides surrounded by fracture face
				}
				else
				{
					UG_LOG("komische Ecke" << std::endl);
					UG_THROW("komische Ecke" << std::endl);
				}
			}

		}

		std::vector<Edge *> edgesToBeSplitted;

		for( Face * fac : facesAtInnerBoundary )
		{
			std::vector<Vertex *> assoVrt;

			CollectVertices(assoVrt, m_grid, fac);

			// check if at a vertex of a frac face is associated with two inner boundary edges

			for( Vertex * vrt : assoVrt )
			{
				std::vector<Edge*> assoEdg;

				CollectEdges( assoEdg, m_grid, vrt );

				IndexType innerBoundaryEdges = 0;

				// count number of edges of the vertex which have one free side

				for( Edge * edg : assoEdg )
				{
					// check if edge is inner boundary edge, i.e. if has at one side non-fracture face
					// only check for those edges which have one fracture face on one side

					// split only inner edges at the moment in case of problems
					if( IsBoundaryEdge3D(m_grid, edg) )
					{
						continue;
					}

					std::vector<Face *> assoFac;

					CollectFaces( assoFac, m_grid, edg );

					IndexType numFacsFromFracSudo = 0;

					for( Face * testFac : assoFac )
					{
						IndexType testSudo = m_sh.get_subset_index(testFac);

						if( testSudo == fracIndSudo )
							numFacsFromFracSudo++;
					}

					if( numFacsFromFracSudo == 0 )
					{
//						UG_LOG("no facs at edg " << std::endl );
						; // nothing to do, belongs from an edge not relevant for the fracture
					}
					else if( numFacsFromFracSudo == 1 )
					{
						// relevant! edge needs to be split

//						UG_LOG("one fac at edg" << std::endl);

//						addElemToSplit( edgesToBeSplitted, edg );
						innerBoundaryEdges++;
					}
					else if( numFacsFromFracSudo == 2 )
					{
//						UG_LOG("two fac sides at edge " << std::endl);
						; // nothing to do, edge at both sides surrounded by fracture face
					}
					else
					{
						UG_LOG("komische Ecke" << std::endl);
						UG_THROW("komische Ecke" << std::endl);
					}

				}

				if( innerBoundaryEdges == 0  || 1 || 3 )
				{
//					UG_LOG("komische innere Grenze ohne Grenze " << std::endl);
//					UG_THROW("komische innere Grenze ohne Grenze " << std::endl);
//					; // nothing to do
//				}
//				else if( innerBoundaryEdges == 1 )
//				{
					; // nothing to do, entire boundary edge with no problem, or external edge
				}
				else if( innerBoundaryEdges == 2 )
				{
					// figure out that edge that is NOT the boundary edge

					UG_LOG("we have two inner boundary edges" << std::endl);

					std::vector<Edge*> assoEdg;

					CollectEdges( assoEdg, m_grid, fac );

					for( Edge * edg : assoEdg )
					{
						// check if edge is inner boundary edge, i.e. if has at one side non-fracture face
						// only check for those edges which have one fracture face on one side

						// split only inner edges at the moment in case of problems
						if( IsBoundaryEdge3D(m_grid, edg) )
						{
							continue;
						}

						std::vector<Face *> assoFac;

						CollectFaces( assoFac, m_grid, edg );

						IndexType numFacsFromFracSudo = 0;

						for( Face * testFac : assoFac )
						{
							IndexType testSudo = m_sh.get_subset_index(testFac);

							if( testSudo == fracIndSudo )
								numFacsFromFracSudo++;
						}

						if( numFacsFromFracSudo == 0 ||  numFacsFromFracSudo == 1 )
						{
							// relevant! edge needs to be split

		//					UG_LOG("one fac at edg" << std::endl);
							;
//							addElem( facesAtInnerBoundary, fac );
						}
						else if( numFacsFromFracSudo == 2 )
						{
		//						UG_LOG("two fac sides at edge " << std::endl);
								; // nothing to do, edge at both sides surrounded by fracture face
								addElem(edgesToBeSplitted, edg);
						}
						else
						{
							UG_LOG("komische Ecke" << std::endl);
							UG_THROW("komische Ecke" << std::endl);
						}
					}

				}
				else
				{
					UG_LOG("how many inner boundary edges at a fracture face???" << std::endl) ;
					UG_THROW("how many inner boundary edges at a fracture face???" << std::endl) ;
				}
			}

			// figure out that edge that needs to be splitted - it is that one which is in touch with the fracture on two sides

//			for( Face * fac : facesAtInnerBoundary )



		}

		for( Edge * edg : edgesToBeSplitted )
		{
			vector3 center = CalculateCenter(edg, m_aaPos);
			UG_LOG("splitting edge at " << center << std::endl);
			RegularVertex* vrt = SplitEdge<RegularVertex>(m_grid, edg, false);
			m_aaPos[vrt] = center;
		}

	}


	return splittedEdges;

}

#endif

//////////////////////////////////////////////////////////////////

template<typename ELEMTYP>
bool ArteExpandFracs3D::addElem(std::vector<ELEMTYP> & knownElems, ELEMTYP elemToAdd )
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

//////////////////////////////////////////////////////////////////


bool ArteExpandFracs3D::countAndSelectFracBaseNums()
{
	UG_LOG("countandselect" << std::endl);

	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		UG_LOG("sudo ind " << fracIndSudo << std::endl);

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

//			UG_LOG("Gesicht " << m_aaPos[fac] << std::endl);

			vector3 facCenter = CalculateCenter( fac, m_aaPos );
			UG_LOG("fac center " << facCenter << std::endl);

			for( IndexType i = 0; i < fac->num_vertices(); i++ )
				UG_LOG("Vertex " << i << " -> " << m_aaPos[fac->vertex(i)] << std::endl);

			UG_LOG("alle Vertizes" << std::endl);

			m_sel.select(fac);

			UG_LOG("selektiert msel " << fac << std::endl);

			m_aaMarkFaceIsFracB[fac] = true;

			UG_LOG("mark bool " << fac << std::endl);

//			return true;

			std::vector<Edge*> facEdges;

			UG_LOG("kollektiere ecken" << std::endl);

			CollectEdges( facEdges, m_grid, fac );

			IndexType d_anzahlEcken = facEdges.size();

			if( d_anzahlEcken == 0 )
			{
				UG_LOG("keine Ecken " << std::endl);
				return true;
			}

			UG_LOG("Anzahl Ecken " << d_anzahlEcken << std::endl);

			for( auto const & edg : facEdges )
			{
				m_aaMarkEdgeVFP[edg]++;

				if( IsBoundaryEdge3D( m_grid, edg ) )
					m_aaMarkEdgeVFP[edg].setIsBndFracVertex();

				m_sel.select(edg);
			}

			UG_LOG("Ecken gesammelt " << std::endl);

			for(size_t i = 0; i < fac->num_vertices(); ++i)
			{
				Vertex* vrt = fac->vertex(i);
				m_sel.select(vrt);
				// TODO FIXME hier anpassen, herausfinden ob fracture geschlossen
				// oder inneres Ende, und Anzahl der umgebenden fractures bestimmen!!!

				auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

//				m_aaMarkVrtVFP[vrt]++;
				vrtxFracPrps++;
				// vielleicht auch in nachfolgendem Loop über alle selektierten Vertizes,
				// wo dann die attached faces abgelaufen und dabei die Subdomain Nummer ermittelt wird
				// kann eventuell sogar im Hauptloop gemacht werden, dann könnte man praktisch
				// das alte vertex fracture properties von 2D weiter verwenden, noch zu klären

//				m_aaMarkVrtVFP[vrt].addFractSudo(fracIndSudo);
				vrtxFracPrps.addFractSudo(fracIndSudo);

				if( IsBoundaryVertex3D(m_grid, vrt) )
				{
//					m_aaMarkVrtVFP[vrt].setIsBndFracVertex();
					vrtxFracPrps.setIsBndFracVertex();
				}

				// nicht mehr notwendig
//				// die Ecken heraus filtern, die mit diesem Vertex assoziert sind
//				std::vector<Edge*> attEdg;
//
//				for( auto const & ae: facEdges )
//				{
//					if( EdgeContains(ae,vrt) )
//						attEdg.push_back(ae);
//				}
//
//				if( attEdg.size() == 2 )
//				{
//					EdgePair edgPr( attEdg[0], attEdg[1] );
//
//					AttachedFractFaceEdgeSudo afes( fac, edgPr, fracIndSudo );
//
//					vrtxFracPrps.addAttachedFractElem(afes);
//				}
//				else
//				{
//					UG_THROW("number of attached edges wrong " << std::endl);
//				}

			}
		}
	}


	//	now make sure that no inner edge is associated with two
	//	boundary vertices (referring to the selection)

	constexpr bool splitEdgesTwoBdryVrt = false;

	if( splitEdgesTwoBdryVrt )
	{
		std::vector<Edge*> tmpEdges;

		for(EdgeIterator iterEdg = m_sel.begin<Edge>(); iterEdg != m_sel.end<Edge>(); iterEdg++ )
		{
			Edge* edg = *iterEdg;

			Vertex * vrtZer = edg->vertex(0);
			Vertex * vrtOne = edg->vertex(1);

			VertxFracPropts & vrtxFracPrpsVrtZer = m_aaMarkVrtVFP[ vrtZer ];
			VertxFracPropts & vrtxFracPrpsVrtOne = m_aaMarkVrtVFP[ vrtOne ];
			VertxFracPropts & edgeFracPrps = m_aaMarkEdgeVFP[ edg ];

			if(      vrtxFracPrpsVrtZer.getIsBndFracVertex()
				&&   vrtxFracPrpsVrtOne.getIsBndFracVertex()
				&& ! edgeFracPrps.getIsBndFracVertex()
			)
			{
				tmpEdges.push_back(edg);
			}

		}

		for( Edge * edg : tmpEdges )
		{
			vector3 center = CalculateCenter(edg, m_aaPos);
			RegularVertex* vrt = SplitEdge<RegularVertex>(m_grid, edg, false);
			m_aaPos[vrt] = center;
			m_sel.select(vrt);
			auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];
			vrtxFracPrps++;

			vrtxFracPrps.setIsBndFracVertex(false);

			//	assign adjacency values for associated selected edges (2 to each)
			for(Grid::AssociatedEdgeIterator iterEdg  = m_grid.associated_edges_begin(vrt);
											 iterEdg != m_grid.associated_edges_end(vrt);
											 iterEdg++
			)
			{
				Edge * assoEdg = *iterEdg;

				if( m_sel.is_selected(assoEdg) )
				{
					auto & edgFracPrps = m_aaMarkEdgeVFP[assoEdg];
					edgFracPrps.setIsBndFracVertex(false);
				}
			}
		}

		// TODO FIXME unsicher, ob das hier richtig übertragen von Prof. Reiter......
	}

#if 0

	// TODO FIXME das ist was komisches, was von Prof. Reiter da ist, es werden edges gesplittet, für was?

	edges.clear();
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(aaMarkVRT[e->vertex(0)] != 2 &&
			aaMarkVRT[e->vertex(1)] != 2 &&
			aaMarkEDGE[e] > 1)
		{
			edges.push_back(e);
		}
	}

	for(size_t i = 0; i < edges.size(); ++i){
		vector3 center = CalculateCenter(edges[i], aaPos);
		RegularVertex* v =	SplitEdge<RegularVertex>(grid, edges[i], false);
		aaPos[v] = center;
		aaMarkVRT[v] = 2;
		sel.select(v);
	//	assign adjacency values for associated selected edges (2 to each)
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
			iter != grid.associated_edges_end(v); ++iter)
		{
			if(sel.is_selected(*iter))
				aaMarkEDGE[*iter] = 2;
		}
	}


#endif

	return true;

}

int ArteExpandFracs3D::prepareStasi( Vertex * const & vrt, AttachedVolumeElemInfo & attVolElmInfo )
{
	// NOTE returns number of boundary faces

	// Voraussetzung  FÜR StammiBene Aufrufung
	// Stammi-Bene-Vorbereitung
//	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
//	{
//		Vertex* vrt = *iter;

//		std::vector<Volume*> & attVol = m_aaVrtInfoAssoVols[vrt];

	auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

//	bool isBndryVrtx = vrtxFracPrps.getIsBndFracVertex();

	IndexType numBndryFacs = 0;

		// TODO FIXME das hier soll wegfallen, und an dieser Stelle direkt berechnet werden
		// damit es einheitlich für echte fracture faces und auch für boundary faces gilt,
		// die eine Art Fracture sind, an denen mit dem Durchmesser 0 der neue Vertex verschoben wird

		// die erzeugen wir an Ort und Stelle neu und verteilen nicht alles auf hundert Plätze
//		VecAttachedFractFaceEdgeSudo vecAttFacEdgSudo = vrtxFracPrps.getAllAttachedFractElems();

//		auto & vecVolElmInfo = m_aaVolElmInfo[vrt];

		// TODO FIXME eigentlich eine Dummheit, das auf zu teilen in ein VolElemInfo und ein VrtInfoAssoVols
		// denn die InfoAssoVols haben als Info nur die Volumen, die VolElmInfos haben die Volumen
		// und noch viel mehr Infos zu den Faces und den Edges....
		// mittelfristig die m_aaVrtInfoAssoVols abschaffen und alles auf die AttachedFullDimElemInfo
		// übertragen, dann geht der folgende Loop gleich über den Vektor darüber, bzw. gleichbedeutend
		// über m_aaVolElmInfo
//		for( auto & vol : attVol )
//		for( AttachedVolumeElemInfo & attVolElmInfo : vecVolElmInfo )
//		{
//			AttachedVolumeElemInfo attVolElmInfo( vol );
	Volume * vol = attVolElmInfo.getFulldimElem();

			// add those faces which are fracture faces
			// TODO FIXME hier müssen die fracture faces neu erzeugt und addiert werden, oder weiter unten.....
//			for( auto & afes : vecAttFacEdgSudo )
//			{
//				attVolElmInfo.addFractManifElem(afes, m_grid);
//			}

			// add those faces which are NOT fracture faces, assign them arbitraryly subdomain  -1
			// to indicate that they are not from the manifold, independent of their subdomain

			// collect all volume faces which incorporate the vertex

	std::vector<Face*> volFacesContainingVrtx;

	for( IndexType iFac = 0; iFac < vol->num_faces(); iFac++ )
	{
		Face * fac = m_grid.get_face(vol,iFac);

		if( FaceContains( fac, vrt ) )
		{
			volFacesContainingVrtx.push_back( fac );
		}
	}

	for( auto const & fac : volFacesContainingVrtx )
	{
				// get the edges of the face connected to the vertex

		std::vector<Edge*> vecEdgesFaceVrtx;

				// need to be two edges always, check

		for( IndexType iEdge = 0; iEdge < fac->num_edges(); iEdge++ )
		{
			Edge * edg = m_grid.get_edge(fac,iEdge);

			if( EdgeContains(edg,vrt) )
			{
				vecEdgesFaceVrtx.push_back(edg);
			}
		}

		if( vecEdgesFaceVrtx.size() != 2 )
		{
			UG_LOG("edge number Unsinn " << vecEdgesFaceVrtx.size() << std::endl);
			UG_THROW("edge number Unsinn " << vecEdgesFaceVrtx.size() << std::endl);
			return false;
		}

		EdgePair edgesFaceVrtx( vecEdgesFaceVrtx[0], vecEdgesFaceVrtx[1] );

				// test the subdomain first, if from the subdomains of the cleft manifolds

		IndexType sudoThisFace = m_sh.get_subset_index(fac);

		std::vector<IndexType> const & sudoList = vrtxFracPrps.getSudoList();

				// test if sudo of face belongs to the fracture face subdom list

		bool faceBelongsToFracSudo = false;

		for( auto const & sudoFrac : sudoList )
		{
			if( sudoFrac == sudoThisFace )
			{
				faceBelongsToFracSudo = true;
				break;
			}
		}

				// TODO FIXME für Boundary faces und für Fracture faces soll die Normale berechnet werden ins VOlumen hinein
				// hier die KuhVol-Prozedur

		bool isBoundaryFace = false;

		if( ! faceBelongsToFracSudo )
		{
			if( IsBoundaryFace3D( m_grid, fac ) )
			{
				isBoundaryFace = true;
			}
		}

		if( isBoundaryFace && faceBelongsToFracSudo )
			UG_THROW("not allowed to be fracture at boundary" << std::endl);

		// will get a nonzero value in case that fracture face or boundary face
		NormalVectorFacIntoVol normalIntoVol(0,0,0);

		if( isBoundaryFace || faceBelongsToFracSudo )
		{
				// TODO FIXME compute normal into volume, needed to know!!!
					// kuhVol Procedure
			if( ! computeNormalKuhVolProcedure(vol,fac,normalIntoVol) )
			{
				UG_LOG("Kuh VOl schief gegangen " << std::endl);
				UG_THROW("Kuh VOl schief gegangen " << std::endl);
				return numBndryFacs;
			}
		}

		if( faceBelongsToFracSudo )
		{
					// if it belongs, construct it again and test if it already belongs to the fracture faces
					// MUST be already part of the list, else major error appeared!

			UG_LOG("Kuh Normale " << normalIntoVol << std::endl);

			AttachedFractFaceEdgeSudo afesFract( fac, edgesFaceVrtx, sudoThisFace, normalIntoVol );
					//  hier soll auch gleich die Normale relativ zum VOlumen dazu gespeichert werden!!

			attVolElmInfo.addFractManifElem( afesFract, m_grid );

					// vorher soll es nicht mehr gespeichert werden, da die Infos vorher weg fallen sollen
//					if( attVolElmInfo.addFractManifElem( afesFract, m_grid ) )
//					{
//						UG_LOG("manifold element already contained!" << std::endl);
//						UG_THROW("manifold element already contained!" << std::endl);
//						return false;
//					}

					// nothing to do, already added before hoffentlich

		}
		else
		{
					// TODO FIXME hier muss entschieden werden, ob es eine boundary face ist
			// die kommt dann auch nochmal in ein anderes Konstrukt, was aber
					// dem fracture face ähnelt
					// am Ende müssen die Frac Faces, die nicht geschlossen sind, und doppelt in einem Segment
					// in die general faces verschoben werden, da stören sie nicht, und sind egal

					// zeitweilig fehlten alle Faces, die keine fractures sind
					// die müssen noch irgendwie als nicht-fracture-faces auch dazu gefügt werden
					// die sind in den attached volumes schon enthalten,
					// Frage: wie prüfen, dass die gar keine fractures sind, die Infos sollten bekannt sein.....
					// wichtig auch, dass nur die faces dazu kommen, die den Vertex enthalten!!!
					// irgendwas von der Art "nonFractureFaceInfos" oder sowas dazu hängen, mit Info
					// ob schon getouched oder noch nicht.....


					// we construct the attached manifold, given that it is NOT a fracture manifold

					// notwendig auch, dass es eine Markierungsmöglichkeit gibt dafür, ob
					// ein face bei der nächsten weiter inneren Runde quasi äussere Begrenzung sein muss
					// gilt sowohl für fracture faces, die können das in erster Runde auch sein am Ende
					// der Runde, sowie danach nur noch für nicht-fracture-faces
//					if( IsBoundaryFace3D( m_grid, fac ) )
			if( isBoundaryFace )
			{
				AttachedBndryFaceEdgeSudo afesBndry( fac, edgesFaceVrtx, sudoThisFace, normalIntoVol );

				attVolElmInfo.addBndryManifElem( afesBndry, m_grid );

				UG_LOG("Boundary element added" << m_aaPos[vrt] << " -> " << normalIntoVol << std::endl);

				numBndryFacs++;

//							UG_THROW("da ist es " << std::endl);

				// TODO FIXME sehr wichtig, die Boundary faces brauchen auch
						// noch ihre Normale in das Volumen hinein, analog zu den Fracture faces!!!
						// die haben die allerdings dummerweise in die boundary vertizes gesteckt
						// vielleicht boundaries wie fractures behandeln? keine gute Idee......
						// mithin steckt die Normale NICHT in diesen Informationen drin
						// sondern extra in den vertex fracture tripeln
						// vielleicht dort die boundaries irgendwie dazu würgen.....

			}
			else // normal face, no fracture, no boundary
			{
				AttachedGenerFaceEdgeSudo afesAdd( fac, edgesFaceVrtx );

				attVolElmInfo.addGenerManifElem( afesAdd, m_grid );
			}


		}

	}

//		}
//	}

//	if( isBndryVrtx && numBndryFacs == 0 )
//	{
//		UG_LOG("boundary vertex but no boundary faces ")
//	}

	return numBndryFacs;
}

bool ArteExpandFracs3D::computeNormalKuhVolProcedure( Volume * const & kuhVol, Face * const & fac, NormalVectorFacIntoVol & normalIntoVol )
{
	IndexType numFd = 0;

	for( IndexType iSide = 0; iSide < kuhVol->num_sides(); iSide++ )
	{
		Face * kuhFac = m_grid.get_side(kuhVol, iSide);

//					UG_LOG("Center Kuh Face " << CalculateCenter(kuhFac,m_aaPos) << std::endl);

		//  eigentliches Ziel ist es, den face descriptor des Volumens zu finden,
		// das mit dem face übereinstimmt, alle anderen Seiten des Volumens sind egal
		// Funktion suchen, die ausgehend von einem Face, das ein Volumen begrenzt,
		// den zum Volumen passenden FaceDescriptor findet, also auch richtige Orientierung
		// der Vertices beinhaltet
		// FRAGE: ist ein face descriptor von der Orientierung zum Volumen abhängig
		// oder hängt der nur vom Face ab, das eine vorgegebene Oriertierung hat?

		bool checkCoincide = checkIfFacesVerticesCoincide( kuhFac, fac );

		if( checkCoincide )
		{
			numFd++;

			if( kuhFac != fac )
				UG_LOG("Kuh Fac ist nicht fac " << std::endl);

			FaceDescriptor facDescr;

			// testen, ob der Face Descriptor von der Orientierung abhängt
			// also testen, ob sich der face descriptor ändert, wenn das Volumen
			// auf der einen und auf der anderen Seite des faces hängt
			// deswegen auch die ganze Prozedur mit den kuhFacs, die hoffentlich
			// je nach Volumen anders orientiert sind als das eigentliche Face,
			// aber dieselben Vertices haben, also geometrisch gleich sind, aber anders orientiert!!!!

			//  andere Hergehensweise vielleicht:
			// von m_aaVrtInfoAssoVols ausgehen, darüber loopen, oder die in einen Vektor stecken,
			// wo die Vertices dabei sind, dann kann man sich vielelicht ein paar Klimmzüge sparen,
			// vielleicht aber auch nicht.....

			kuhVol->face_desc( iSide, facDescr );

//			vector3 normal;

			CalculateNormal( normalIntoVol, & facDescr, m_aaPos );

			vector3 facCenter = CalculateCenter( kuhFac, m_aaPos );
			vector3 kuhCenter = CalculateCenter( fac, m_aaPos );
			vector3 kuhVolCenter = CalculateCenter( kuhVol, m_aaPos);

//						UG_LOG("Normale zum face descriptor " << normal << " , " << facCenter << std::endl);
//						UG_LOG("Normale zum Kuhh descriptor " << normal << " , " << kuhCenter << std::endl);
//						UG_LOG("Zentrum des Vol")

//						UG_LOG("fac " << fac << std::endl );
//						UG_LOG("kuh " << kuhFac << std::endl );

			UG_LOG("Normale ist " << normalIntoVol << " fac " << facCenter
					<< " vol " << kuhVolCenter << std::endl);


//						VolNormPair normalsAwayVol( kuhVol, normal );
//
//						vecNormalsAwayVol.push_back( normalsAwayVol );

		}
	}

	if( numFd != 1 )
	{
		UG_THROW("Kein Kuh Volumen gefunden" << std::endl);
		return false;
	}


	return true;
}

//////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::distinguishSegments()
{

	UG_LOG("neuer Beginn Segmente" << std::endl);
//	return true;

	//	StammiBene Algorithmus erklärt und in Worten Plan erstellt
	//	Search the adjacent merged manifold iteratively between each necessary element
	// TODO FIXME hier Loop über alle selektierten Vertizes
	// darin für jeden Vertex die adjungierten Volumen bestimmen ohne Vorbedingung
	// dann den Loop zur Vorbereitung des StammiBene Algorithmus aufrufen
	// für jeden Vertex wieder, der hat dann schon die Basisinfo von der
	// VecAttachedVolumeElemInfo Klasse an jedem Vertex
	// auf die aufbauend fügt er Fracture  Manifolds, general manifolds,
	// und NEU auch boundary manifolds in einer eigenen Member dazu
	// danach wird für jeden Vertex der StammiBene Algorithmus aufgerufen
	// später werden Boundary Faces wie eine eigene Subdomain Ebene
	// mit Expansion null behandelt
	// was an Knicken aussen an boundary zu tun ist, noch zu überlegen

	UG_LOG("Stasi Algo alle Vrt Start << std::endl");

	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

		if( ! stasiAlgo( vrt ) )
		{
			UG_LOG("Stasi schief gegangen " << std::endl);
			UG_THROW("Stasi schief gegangen" << std::endl);
			return false;
		}
	}

	UG_LOG("Stasi Algo alle Vrt End << std::endl");

	// TODO FIXME
	// es müssen in den Segmenten noch die fracture faces zu generellen faces verschoben werden,
	// die nicht abgeschlossen sind, und umgangen werden können, und so in die Landschaft ragen
	// also die fracture faces, die für die Segmenteigenschaft unerheblich sind

	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

		IndexType shiFraFac = shiftUnclosedFracFacesToUnclosedFractFaces( vrt );

		UG_LOG("shifted frac faces at " << m_aaPos[vrt] << shiFraFac << std::endl);
	}


	return true;

}

//////////////////////////////////////////////////////////////////


bool ArteExpandFracs3D::detectEndingCrossingCleftsSegmBased()
{

	m_vecEdgeDirectConnectingEndingCrossCleftVrtcs = std::vector<Edge*>();

	// TODO FIXME eigene Klasse erzeugen, die die ganzen Infos an endenden Klüften speichert, wo sich
	// zwei Klüfte kreuzen
	// Vertex, auslaufende faces, wo anders aufgeblasen werden muss,
	// Edges der Kreuzung, die anhängenden faces der anderen Kluft, die
	// die Basis sind für die Aufblasungen, die auch anders funktionieren
	// von dieser Klasse dann ein attachment an die vertizes und vielleicht die faces,
	// die problematisch sind
	// (noch zu klären, ob vielleicht ein paar der auskommentierten attachments doch sinnvoll)

	IndexType numEndingCrossingClefts = 0;

//	std::vector<Face*> m_d_endingCrossingCleftFaces;
//
//	std::vector<Vertex*> m_d_endingCrossingCleftVrtcs;
//
//	std::vector<Edge*> m_d_cuttingEdges;
//	std::vector<Face*> crossingNotEndingFaces;
//
//	std::vector<Edge*> otherEdgeOfCrossingNotEndingFace;
//	std::vector<Face*> nextFaceOfCrossingNotEndingFaces;

	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

		VecSegmentLimitingSides & vecSegmLimSid = m_vrtxAttAccsVecSegmLimSid[vrt];

		if(  vecSegmLimSid.size() == 0 )
		{
			UG_LOG("keine verschiedenen Segmente ECCV" << std::endl);
//			return true;
		}

		if( vecSegmLimSid.size() == 1 )
		{
			UG_LOG("sonderbarer Fall von nur einem Segment ECCV" << std::endl );
			UG_THROW("sonderbarer Fall von nur einem Segment ECCV" << std::endl );
			return false;
		}

		for( SegmentLimitingSides & segLimSids : vecSegmLimSid )
		{
			if( segLimSids.hasUnclosedFaces() )
			{
				UG_LOG("im distinguisch ungeschlossene Gesichter " << std::endl );

				numEndingCrossingClefts++;

				std::vector<Volume*> volsOfUnclosedSegm;

				if( ! segLimSids.spuckVecFulldimElem(volsOfUnclosedSegm) )
				{
					UG_LOG("No volumes at segment with ending cleft" << std::endl);
					UG_THROW("No volumes at segment with ending cleft" << std::endl);
					return false;
				}

				for( Volume * vol : volsOfUnclosedSegm )
				{
					m_volAttAccsVolTouchesEndingCrossingCleft[vol] = true;
				}

				Vertex * endingCrossingCleftVrtx = segLimSids.spuckVertex();

				if( endingCrossingCleftVrtx != vrt )
				{
					UG_THROW("System von attachments und vektoren durcheinander " << std::endl);
				}

				//  commented out the debug members:  m_d_endingCrossingCleftVrtcs.push_back( endingCrossingCleftVrtx );

				UG_LOG("vertex gefunden an ending crossing cleft " << std::endl);

				m_aaMarkVrtxAtEndingCrossingCleft[vrt] = true;

				// the unclosed ending crossing faces
				VecSegLimSidesFractFace vecSegmLimSiFFUnclosed;

				if( ! segLimSids.spuckVecAttUnclosedFractElm( vecSegmLimSiFFUnclosed ) )
					return false;

				// figure out the faces that are in touch with the ending faces, and the corresponding edges

				VecSegLimSidesFractFace vecSegmLimSiFFClosed;

				if( ! segLimSids.spuckVecAttFractElm( vecSegmLimSiFFClosed ) )
					return false;

				UG_LOG("Faces suchen an ending crossing cleft " << std::endl);

//				Face * endingFractManifCutting = nullptr;
//				Face * endingFractManifNotCutting = nullptr;
//				Edge * oldLowDimElCut = nullptr;
//				std::pair<Face*,Face*> pairNeighbouredFractClosedManifEl(nullptr,nullptr);
//				IndexType sudoFractEnding = std::numeric_limits<IndexType>::max();
//				IndexType sudoFractNotEnding = std::numeric_limits<IndexType>::max();

				std::vector<Face*> vecEndingFractFaceCutting; // must have size two after been filled but both same content
				std::vector<Face*> vecEndingFractFaceNotCutting; // must have size maximum two same content both, or zero
				std::vector<Face*> vecNeighbouredFacesClosedFract; // must have size 2

				std::vector<IndexType> vecSudoFractFacsEnding; // must be filled with same number
				std::vector<IndexType> vecSudoFractFacsNotEnding; // must be filled with same number

				std::vector<Edge*> vecCuttingEdges; // must contain only same element, twice

				std::vector<Edge*> vecUncuttingEdge; // must get same element twice

				UG_LOG("size of vec segm unclos " << vecSegmLimSiFFUnclosed.size() << std::endl);

				for( SegLimSidesFractFace const & slsffUncl : vecSegmLimSiFFUnclosed )
				{
					bool isCuttingTheClosed = false;

					Face * facUncl = slsffUncl.getManifElm();
					//  commented out the debug members:  m_d_endingCrossingCleftFaces.push_back(facUncl);
					m_aaMarkFaceWithEndingCrossingCleft[facUncl] = true;

					IndexType sudoUncl = slsffUncl.getSudo();
					vecSudoFractFacsEnding.push_back(sudoUncl);

					EdgePair const & epU = slsffUncl.getPairLowElm();

					//  commented out the debug members:  m_d_allContributingEdges.push_back(epU.first);
					//  commented out the debug members:  m_d_allContributingEdges.push_back(epU.second);

					addElem( m_vecEdgeDirectConnectingEndingCrossCleftVrtcs, epU.first );
					addElem( m_vecEdgeDirectConnectingEndingCrossCleftVrtcs, epU.second );

					for( SegLimSidesFractFace const & slsffClos : vecSegmLimSiFFClosed )
					{
						// we need different subdoms, as want to have the not ending crossing fracture neighbours
						// that have a common edge
						Edge * commonEdge = nullptr;

						if( fractFacesAreNeighboured<false>( slsffClos, slsffUncl, commonEdge ) )
						{
							isCuttingTheClosed = true;

							Face * facClos = slsffClos.getManifElm();
							//  commented out the debug members:  m_d_crossingNeighboredNotEndingFaces.push_back( facClos );

							if( commonEdge == nullptr )
							{
								UG_LOG("NULL COMMON" << std::endl);
								UG_THROW("NULL COMMON" << std::endl);
							}

							//  commented out the debug members:  m_d_cuttingEdges.push_back(commonEdge);
							//  commented out the debug members:  m_d_crossingNeighboredNotEndingFacesCommEdg.push_back(facUncl);
//							// search the durchgehende fracture face which is neighboured to the durchgehende one
//							// same vertex, other edge

							vecNeighbouredFacesClosedFract.push_back(facClos);
							vecCuttingEdges.push_back(commonEdge);
							vecEndingFractFaceCutting.push_back(facUncl);

							Edge * uncuttingEdgeUnclosedFac = nullptr;

							if( ! assignNotCuttingEdgeUnclosedCrossingCleft( slsffUncl, commonEdge, uncuttingEdgeUnclosedFac ) )
							{
								UG_LOG("NO shift edge" << std::endl);
								UG_THROW("NO shift edge" << std::endl);
							}

							if( uncuttingEdgeUnclosedFac == nullptr )
							{
								UG_LOG("NO shift edge Z" << std::endl);
								UG_THROW("NO shift edge Z" << std::endl);
							}

							vecUncuttingEdge.push_back(uncuttingEdgeUnclosedFac);

							//							EdgePair const & edgesClosedFracFace = slsffClos.getPairLowElm();
//
//							Edge * notCommonEdgeOfNotEndingCleft = nullptr;
//
//							if(  edgesClosedFracFace.first == commonEdge )
//							{
//								notCommonEdgeOfNotEndingCleft = edgesClosedFracFace.second;
//							}
//							else if( edgesClosedFracFace.second == commonEdge )
//							{
//								notCommonEdgeOfNotEndingCleft = edgesClosedFracFace.first;
//							}
//
//							if( notCommonEdgeOfNotEndingCleft == nullptr )
//							{
//								UG_LOG("no not common edge?" << std::endl);
//								UG_THROW("no not common edge?" << std::endl);
//							}
//
////							otherEdgeOfCrossingNotEndingFace.push_back( notCommonEdgeOfNotEndingCleft );
//
//							Edge * closedCommonEdge = nullptr;
//
//							for( SegLimSidesFractFace const & slsffClosOther : vecSegmLimSiFFClosed )
//							{
//								if( fractFacesAreNeighboured<true>( slsffClos, slsffClosOther, closedCommonEdge ) )
//								{
//
//									if( closedCommonEdge == nullptr )
//									{
//										UG_LOG("NULL NOT COMMON" << std::endl);
//										UG_THROW("NULL NOT COMMON" << std::endl);
//									}
//
//									if( closedCommonEdge != notCommonEdgeOfNotEndingCleft )
//									{
//										UG_LOG("Ecke passt nicht " << CalculateCenter(closedCommonEdge, m_aaPos) << " - " << CalculateCenter( notCommonEdgeOfNotEndingCleft, m_aaPos) << std::endl);
////										UG_THROW("Ecke verloren " << std::endl);
//									}
//									else
//									{
//										Face * facClosOther = slsffClosOther.getManifElm();
//
//										nextFaceOfCrossingNotEndingFaces.push_back( facClosOther );
//										otherEdgeOfCrossingNotEndingFace.push_back(notCommonEdgeOfNotEndingCleft);
//									}
//
//
//								}
//							}


						}
					}

					if( ! isCuttingTheClosed )
					{
						vecEndingFractFaceNotCutting.push_back(facUncl);
					}

				}


				std::vector<Face*> vecClosedFracFacNoNeighbr;

				for( SegLimSidesFractFace const & slsffClos : vecSegmLimSiFFClosed )
				{
					Face * facClos = slsffClos.getManifElm();

					EdgePair const & epC = slsffClos.getPairLowElm();

					//  commented out the debug members:  m_d_allContributingEdges.push_back(epC.first);
					//  commented out the debug members:  m_d_allContributingEdges.push_back(epC.second);

					addElem( m_vecEdgeDirectConnectingEndingCrossCleftVrtcs, epC.first );
					addElem( m_vecEdgeDirectConnectingEndingCrossCleftVrtcs, epC.second );


					IndexType sudoClos = slsffClos.getSudo();

					vecSudoFractFacsNotEnding.push_back(sudoClos);

					// valid for all closed faces
					m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft[facClos] = true;

					bool facGiven = false;

					// ATTENTION hier ausgetauscht erst später, Gefahr von Fehler
					// m_d_ Variable war Teil von echtem Algorithmus, gefährlich
					// hoffentlich richtig getauscht
					// im Falle von bugs hier nochmal testen, ob die m_d Variable
					// richtig getauscht wurde!!!
//					for( Face * facClosNoNei : m_d_crossingNeighboredNotEndingFaces )
					for( Face * facClosNoNei : vecNeighbouredFacesClosedFract )
					{
						if( facClosNoNei == facClos )
						{
							facGiven = true;
							break;
						}
					}

					if( ! facGiven )
					{
						//  commented out the debug members:  m_d_notEndingCrossingFacesNotNeighbour.push_back( facClos );
						vecClosedFracFacNoNeighbr.push_back( facClos );
					}
				}


				UG_LOG("Faces vielleicht gefunden an ending crossing cleft " << std::endl);

				//				Face * endingFractManifCutting = nullptr;
				//				Face * endingFractManifNotCutting = nullptr;
				//				Edge * oldLowDimElCut = nullptr;
				//				std::pair<Face*,Face*> pairNeighbouredFractClosedManifEl(nullptr,nullptr);
				//				IndexType sudoFractEnding = std::numeric_limits<IndexType>::max();
				//				IndexType sudoFractNotEnding = std::numeric_limits<IndexType>::max();

				std::vector<Face*> vecContEndingFractFaceCutting; // must have size two after been filled but both same content
				std::vector<Face*> vecContEndingFractFaceNotCutting; // must have size maximum two same content both, or zero
				std::vector<Face*> vecContNeighbouredFacesClosedFract; // must have size 2


				std::vector<IndexType> vecContSudoFractFacsEnding; // must be filled with same number
				std::vector<IndexType> vecContSudoFractFacsNotEnding; // must be filled with same number

				std::vector<Edge*> vecContCuttingEdges; // must contain only same element, twice

				std::vector<Edge*> vecContUncuttingEdge; // must get same element twice

				if( vecEndingFractFaceCutting.size() == 2 )
				{
					if( ! checkIfContentUnique( vecEndingFractFaceCutting, vecContEndingFractFaceCutting, 1 ) )
					{
						UG_THROW("Problem with cutting fac " << std::endl);
					}
				}
				else
				{
					UG_THROW("Problem with cutting fac size " << std::endl);

				}

				Face * endingFractManifCutting = vecContEndingFractFaceCutting[0];

				Face * endingFractManifNotCutting = nullptr;

				if( vecEndingFractFaceNotCutting.size() != 0 )
				{
					if( vecEndingFractFaceNotCutting.size() == 1 )
					{
						if( ! checkIfContentUnique( vecEndingFractFaceNotCutting, vecContEndingFractFaceNotCutting, 1 ) )
						{
							UG_THROW("Problem with not cutting fac " << std::endl);
						}

						endingFractManifNotCutting = vecContEndingFractFaceNotCutting[0];

					}
					else
					{
						for( Face * fac : vecEndingFractFaceNotCutting )
						{
							m_sh.assign_subset( fac, m_sh.num_subsets() );
						}

						UG_LOG("schief gegangen " << std::endl);

						return false;

						UG_THROW("Problem with not cutting fac size " << vecEndingFractFaceNotCutting.size() << std::endl);
					}

				}


				if( vecNeighbouredFacesClosedFract.size() == 2  )
				{
					if( ! checkIfContentUnique( vecNeighbouredFacesClosedFract, vecContNeighbouredFacesClosedFract, 2 ) )
					{
						UG_THROW("Problem with neigh clos fac " << std::endl);
					}

				}
				else
				{
					UG_THROW("Problem with neigh clos fac size " << std::endl);
				}

				std::pair<Face*,Face*> pairNeighbouredFractClosedManifEl( vecContNeighbouredFacesClosedFract[0],
																		  vecContNeighbouredFacesClosedFract[1]
																		 );


				if(  ! checkIfContentUnique( vecSudoFractFacsEnding, vecContSudoFractFacsEnding, 1 ) )
				{
					UG_THROW("Problem with sudo ending " << std::endl);
				}

				IndexType sudoFractEnding = vecContSudoFractFacsEnding[0];

				if(  ! checkIfContentUnique( vecSudoFractFacsNotEnding, vecContSudoFractFacsNotEnding, 1 ) )
				{
					UG_THROW("Problem with sudo not ending " << std::endl);
				}


				IndexType sudoFractNotEnding = vecContSudoFractFacsNotEnding[0];

				if( vecCuttingEdges.size() == 2 )
				{
					if(  ! checkIfContentUnique( vecCuttingEdges, vecContCuttingEdges, 1 ) )
					{
						UG_THROW("Problem with edges unique " << std::endl);
					}

				}
				else
				{
					UG_THROW("Problem with edges unique size " << std::endl);

				}

				Edge * oldLowDimElCut = vecContCuttingEdges[0];

				UG_LOG("Pushing back Edge" << std::endl);

				if( vecUncuttingEdge.size() == 2 )
				{
					if(  ! checkIfContentUnique( vecUncuttingEdge, vecContUncuttingEdge, 1 ) )
					{
						UG_THROW("Problem with edges shift unique " << std::endl);
					}

				}
				else
				{
					UG_THROW("Problem with edges shift unique size " << std::endl);

				}

				Edge * uncuttingLowDimEl = vecContUncuttingEdge[0];

				Edge * shiftLowDimEl = nullptr;

				if( ! findShiftEdgeUncuttingCrossingCleft( vecSegmLimSiFFUnclosed,
														   endingFractManifCutting,
													       endingFractManifNotCutting,
														   uncuttingLowDimEl,
														   shiftLowDimEl
														  )
					)
				{
					UG_LOG("Kein shift Element gefunden" << std::endl);
					UG_THROW("Kein shift Element gefunden" << std::endl);
				}

				if( shiftLowDimEl == nullptr )
				{
					UG_LOG("SHift El null " << std::endl);
					UG_THROW("SHift El null " << std::endl);
				}

				if(  vecClosedFracFacNoNeighbr.size() + vecNeighbouredFacesClosedFract.size()
					!=  vecSegmLimSiFFClosed.size()
				)
				{
					UG_THROW("size sum clos not fit " << std::endl);
				}

				EndingCrossingFractureSegmentInfo ecfsi( vrt,
													     endingFractManifCutting,
													     endingFractManifNotCutting,
														 oldLowDimElCut,
														 pairNeighbouredFractClosedManifEl,
														 shiftLowDimEl,
														 sudoFractEnding,
														 sudoFractNotEnding
													    );

				ecfsi.schluckVecClosedFracManifElNoNeighbr( vecClosedFracFacNoNeighbr );

				ecfsi.schluckVecFulldimElm( volsOfUnclosedSegm );

				m_vrtxAttAccsVecEndingCrossFractSegmInfo[vrt].push_back( ecfsi );

				m_vecEndCrossFractSegmInfo.push_back( ecfsi );


				segLimSids.schluckLowdimElmShiftDirectionIfUnclosedFractPresent(shiftLowDimEl);

				// schlechte Idee: die segment limiting sides müssen die ganze EndingCrossingFractureSegmentInfo zu wissen bekommen
				// klar alternativ sollte die shift direction reichen, um die ending crossing fracture segment info rekonstruieren zu können
				// wäre aber ein Umweg!!!
				// besser aber der Umweg, ansonsten entsteht das Henne Ei Problem:
				// wer ist von was der Ursprung? wer schluckt was,
				// die Segment Seite die Ending crossing fracture segment info, oder umgekehrt?
				// und wenn es dann als Kopie geschluckt wird, was passiert wenn was verändert wird,
				// es müssten also Pointer sein, und es entsteht eine unsinnige Kaskade.....

			}
		}

	}

	UG_LOG("detected ending crossing cleft faces " << numEndingCrossingClefts << std::endl);

	UG_LOG("EDGES number " << m_vecEdgeDirectConnectingEndingCrossCleftVrtcs.size() << std::endl);

//	if( m_d_endingCrossingCleftVrtcs.size() == 0 )
//		return true;

	if( splitEdgesOfNeighboredEndingCrossingFracVrtcs() )
	{
		UG_LOG("needed to split, please restart with this geometry" << std::endl);
		m_needToSplitEdgesConnectingNeighbrdEndingCrossCleftVrtx = true;
		return false;
	}
	else
	{
		m_needToSplitEdgesConnectingNeighbrdEndingCrossCleftVrtx = false;
	}

//	return true;


	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

		VecSegmentLimitingSides & vecSegmLimSid = m_vrtxAttAccsVecSegmLimSid[vrt];

		for( SegmentLimitingSides const & sls : vecSegmLimSid )
		{
			if( sls.hasUnclosedFaces())
				UG_LOG("Hier nach detect hat das Gesicht noch ungeschlossene " << std::endl);
		}

	}

	if( numEndingCrossingClefts == 0 )
			return true;

//	assignDebugSubsets( true );
//
//	UG_LOG("problematic elements highlighted " << std::endl);

	return true;



//	if( numEndingCrossingClefts == 0 )
//		return true;
//
//	// debug for ending crossing clefts
//
//	for( Face * fac : m_d_endingCrossingCleftFaces )
//	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
//
//	}
//
//	for( Vertex * vrt : m_d_endingCrossingCleftVrtcs )
//	{
//		m_sh.assign_subset( vrt, m_sh.num_subsets());
//	}
//
//	for( Edge * edg : m_d_cuttingEdges )
//	{
//		if( edg == nullptr )
//		{
//			UG_LOG("NULL UNERLAUBT" << std::endl);
//			UG_THROW("NULL UNERLAUBT" << std::endl);
//		}
//
//		m_sh.assign_subset( edg, m_sh.num_subsets());
//
//	}
//
//	for( Face * fac : crossingNotEndingFaces )
//	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
//	}
//
//	for( Edge * edg : otherEdgeOfCrossingNotEndingFace )
//	{
//		if( edg == nullptr )
//		{
//			UG_LOG("NULL C UNERLAUBT" << std::endl);
//			UG_THROW("NULL C UNERLAUBT" << std::endl);
//		}
//
//		m_sh.assign_subset( edg, m_sh.num_subsets());
//
//	}
//
//	for( Face * fac : nextFaceOfCrossingNotEndingFaces )
//	{
//		m_sh.assign_subset( fac, m_sh.num_subsets());
//	}
//
//
//	return false;
}

//////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::assignNotCuttingEdgeUnclosedCrossingCleft( SegLimSidesFractFace const & slsffUncl, Edge * const & commonEdge, Edge * & shiftEdge )
{
	EdgePair const & edgesFac = slsffUncl.getPairLowElm();

	Edge * edgeFrst = edgesFac.first;
	Edge * edgeScnd = edgesFac.second;

	if( edgeFrst == commonEdge )
	{
		shiftEdge = edgeScnd;
	}
	else if( edgeScnd == commonEdge )
	{
		shiftEdge = edgeFrst;
	}
	else
	{
		UG_LOG("common and shift edge different face " << std::endl);
		UG_THROW("common and shift edge different face " << std::endl);
	}

	if( shiftEdge != nullptr )
		return true;

	UG_LOG("shift edge null" << std::endl);

	return false;

}

//////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::findShiftEdgeUncuttingCrossingCleft( VecSegLimSidesFractFace const & vecSegmLimSiFFUnclosed,
															 Face * const & endingFractManifCutting,
		  	  	  	  	  	  	  	  	  	  	  	  	     Face * const & endingFractManifNotCutting,
															 Edge * const & uncuttingLowDimEl,
															 Edge * & shiftLowDimEl
															)
{
	if( uncuttingLowDimEl == nullptr )
	{
		UG_LOG("null uncut" << std::endl );
		UG_THROW("null uncut" << std::endl );
	}

	if( endingFractManifNotCutting == nullptr )
	{
		shiftLowDimEl = uncuttingLowDimEl;
		return true;
	}

	IndexType edgeWasFound = 0;

	for( SegLimSidesFractFace const & slsffUncl : vecSegmLimSiFFUnclosed )
	{

		Face * facUncl = slsffUncl.getManifElm();

		if( facUncl == endingFractManifNotCutting )
		{
			EdgePair const & epU = slsffUncl.getPairLowElm();

			Edge * edgeOne = epU.first;
			Edge * edgeTwo = epU.second;

			if( edgeOne == uncuttingLowDimEl )
			{
				shiftLowDimEl = edgeTwo;
			}
			else if( edgeTwo == uncuttingLowDimEl )
			{
				shiftLowDimEl = edgeOne;
			}
			else
			{
				UG_LOG("No edge found uncutting shift?" << std::endl);
				UG_THROW("No edge found uncutting shift?" << std::endl);
			}

			edgeWasFound++;

		}
	}

	if( edgeWasFound == 0 )
	{
		UG_LOG("NO shift edge " << std::endl);
		return false;
	}
	else if( edgeWasFound == 1 )
	{
		return true;
	}
	if( edgeWasFound > 1 )
	{
		UG_LOG("Many shift edges " << edgeWasFound << std::endl);
		return false;
	}

	UG_LOG("Cannot come here at this edge" << std::endl);

	return false;

}

//////////////////////////////////////////////////////////////////

template<typename ELMTYP>
bool ArteExpandFracs3D::checkIfContentUnique( std::vector<ELMTYP> const & vecTest, std::vector<ELMTYP> & content, IndexType mandatoryDifferentElems )
{
	if( vecTest.size() == 0 )
	{
		UG_LOG("zero vector testing nonsense" << std::endl);
		UG_THROW("zero vector testing nonsense" << std::endl);
	}

	std::vector<ELMTYP> vecTestCop = vecTest; // copy

	while( vecTestCop.size() != 0 )
	{
		ELMTYP elmCont = vecTestCop[0];

		content.push_back(elmCont);

		for( typename std::vector<ELMTYP>::iterator itEl = vecTestCop.begin(); itEl != vecTestCop.end(); )
		{
			ELMTYP elmTest = *itEl;

			if( elmCont == elmTest )
			{
				itEl = vecTestCop.erase(itEl);
			}
			else
				++itEl;
		}

	}

	return ( mandatoryDifferentElems == content.size() );
}

//////////////////////////////////////////////////////////////////

template< bool FACES_HAVE_SAME_SUDO >
bool ArteExpandFracs3D::fractFacesAreNeighboured( SegLimSidesFractFace const & fractFaceOne,
							   	   	   	   	   	  SegLimSidesFractFace const & fractFaceTwo,
												  Edge * & commonEdge
												)
{
	commonEdge = nullptr; // general case

	bool facesHaveSameSudo = FACES_HAVE_SAME_SUDO;

	// exclude that the elements are the same

	Face * facOne = fractFaceOne.getManifElm();
	Face * facTwo = fractFaceTwo.getManifElm();

	if( facOne == facTwo )
	{
		UG_LOG("Faces coincide" << std::endl);
		return false;
	}

	EdgePair const & edgesFacOne = fractFaceOne.getPairLowElm();
	EdgePair const & edgesFacTwo = fractFaceTwo.getPairLowElm();

	Edge * edgeOneFrst = edgesFacOne.first;
	Edge * edgeOneScnd = edgesFacOne.second;

	Edge * edgeTwoFrst = edgesFacTwo.first;
	Edge * edgeTwoScnd = edgesFacTwo.second;

	bool edgeCondA = ( edgeOneFrst == edgeTwoFrst );
	bool edgeCondB = ( edgeOneFrst == edgeTwoScnd );
	bool edgeCondC = ( edgeOneScnd == edgeTwoFrst );
	bool edgeCondD = ( edgeOneScnd == edgeTwoScnd );


	bool fractFacsNeighbr = (    edgeCondA
							  || edgeCondB
							  || edgeCondC
							  || edgeCondD
							);

	IndexType sudoOne = fractFaceOne.getSudo();
	IndexType sudoTwo = fractFaceTwo.getSudo();

	bool sudosCoincide = ( sudoOne == sudoTwo );

	bool necessarySudoProperty = ( facesHaveSameSudo == sudosCoincide );

	bool neighbrsWithRequestedSudo = ( fractFacsNeighbr && necessarySudoProperty );

	if( neighbrsWithRequestedSudo )
	{
		if( edgeCondA || edgeCondB )
		{
			commonEdge = edgeOneFrst;
		}
		else if( edgeCondC || edgeCondD )
		{
			commonEdge = edgeOneScnd;
		}

		if( commonEdge == nullptr )
		{
			UG_LOG("COMMON NULL " << std::endl);
			UG_THROW("COMMON NULL " << std::endl);
		}
		else
		{
			UG_LOG("COORDINATES NOT NULL " << CalculateCenter( commonEdge, m_aaPos ) << std::endl );
		}
	}

	return neighbrsWithRequestedSudo;
}


//////////////////////////////////////////////////////////////////


#if 0
bool ArteExpandFracs3D::detectEndingCrossingClefts()
{
//	bool unclosedFracFacesPresent = false;
//
//	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
//	{
//		Vertex* vrt = *iter;
//
//		IndexType shiFraFac = shiftUnclosedFracFacesToUnclosedFractFaces( vrt );
//
//		UG_LOG("shifted frac faces at " << m_aaPos[vrt] << " -> " << shiFraFac << std::endl);
//
//		if( shiFraFac > 0 )
//			unclosedFracFacesPresent = true;
//
////		constexpr bool d_highlightVrtcsWithShifts = false;
////
////		if( d_highlightVrtcsWithShifts )
////		{
////			if( shiFraFac > 0 )
////			{
////				IndexType sudoNum = m_sh.num_subsets();
////
////				m_sh.assign_subset(vrt, sudoNum );
////			}
////		}
//	}

	IndexType numEndingCrossingClefts = 0;

	std::vector<Face*> m_d_endingCrossingCleftFaces;

	std::vector<Vertex*> m_d_endingCrossingCleftVrtcs;


//	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
//	{
//		Vertex* vrt = *iter;
//
//
//
//	}



	// TODO FIXME unterscheide Faces entlang der expandierten Kreuzungskluft,
	// mit 2 Knoten an Kreuzungspunkten, ein Knoten voll expandiert für beide subsets
	// von solchen, die in der Luft hängen am Ende der Kluft, wo nur die durchgehende Kluft expandiert

	// TODO FIXME need to detect attached faces from the durchgehende cleft which gets expanded
	// as specific structures needed also there

	// TODO FIXME auch die Kante, die an dem auslaufenden Face entlang geht im Schnitt mit dem ausgedehnten,
	// per edge attachment markieren! bool
	// und die anhängenden faces der kreuzenden durchgehenden Kluft ebenso markieren irgendwie per face attachment bool




#if 1

	IndexType unclosedFracFacesFound = 0;

	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{

		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); ++iter )
		{
			Face* fac = *iter;

			UG_LOG("Check for unclosed frac faces " << CalculateCenter(fac, m_aaPos) << std::endl);

			if( m_aaMarkFaceHasUnclosedFracSideB[fac] )
			{
				unclosedFracFacesFound++;

//				m_sh.assign_subset( fac, m_sh.num_subsets());
			}

		}

	}

	if( unclosedFracFacesFound > 0 ) // || unclosedFracFacesPresent )
	{
		UG_LOG("unclosed Frac faces " << unclosedFracFacesFound << std::endl);
//		return false;
	}

//	return true;


	//  wieso die ending crossing cleft faces nicht aus den Segmenten raus holen,
	// und zwar aus den ungeschlossenen Faces, die im Segment liegen?
	// müssen doch dieselben sein, dann sollte die komische Prozedur hier unnötig werden!!!


//	UG_THROW("KÄSE" << std::endl);

	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); iter++ )
		{
			Face* fac = *iter;

			UG_LOG("Detect for unclosed ending cleft frac faces " << CalculateCenter(fac, m_aaPos) << std::endl);

			if( m_aaMarkFaceHasUnclosedFracSideB[fac] )
			{

				// figure out the vertice(s) where the face is unclosed
				// check if at one of these vertices, two fractures are crossing and one ending
				// by checking if there are at least two segments as precondition, else no matter

				std::vector<Vertex *> assoVrt;

				CollectVertices(assoVrt, m_grid, fac);

				for( Vertex * vrt : assoVrt )
				{
					if( m_aaMarkVrtxHasUnclosedFracB[vrt] )
					{
						VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[vrt];

						if( vecSegVolElmInf.size() > 1 ) // expansion takes place here
						{
//							m_sh.assign_subset(fac,m_sh.num_subsets());
							numEndingCrossingClefts++;
							m_d_endingCrossingCleftFaces.push_back(fac);
							m_d_endingCrossingCleftVrtcs.push_back(vrt);
							// TODO FIXME vielleicht eine Klasse, wo die Faces und Vertices einer
							// jeweiligen ending Crossing Cleft Stelle zusammengefasst werden?
							// statt unabhängige Listen? gemeinsame Liste?
						}

					}
				}


			}

		}

	}


	UG_LOG("detected ending crossing cleft faces " << numEndingCrossingClefts << std::endl);

	if( numEndingCrossingClefts == 0 )
		return true;

	// debug for ending crossing clefts

	for( Face * fac : m_d_endingCrossingCleftFaces )
	{
		m_sh.assign_subset( fac, m_sh.num_subsets());

	}

	for( Vertex * vrt : m_d_endingCrossingCleftVrtcs )
	{
		m_sh.assign_subset( vrt, m_sh.num_subsets());
	}

#endif


	return false;



}

#endif

////////////////////////////////////////////////////////////////////

#if 0
bool ArteExpandFracs3D::detectEndingCrossingClefts()
{

	//  wieso die ending crossing cleft faces nicht aus den Segmenten raus holen,
	// und zwar aus den ungeschlossenen Faces, die im Segment liegen?
	// müssen doch dieselben sein, dann sollte die komische Prozedur hier unnötig werden!!!

	IndexType numEndingCrossingClefts = 0;

	std::vector<Face*> m_d_endingCrossingCleftFaces;

	std::vector<Vertex*> m_d_endingCrossingCleftVrtcs;

//	UG_THROW("KÄSE" << std::endl);

	for(size_t i_fi = 0; i_fi < m_fracInfos.size(); ++i_fi )
	{
		int fracIndSudo = m_fracInfos[i_fi].subsetIndex;

		for( FaceIterator iter = m_sh.begin<Face>(fracIndSudo); iter != m_sh.end<Face>(fracIndSudo); iter++ )
		{
			Face* fac = *iter;

			UG_LOG("Detect for unclosed ending cleft frac faces " << CalculateCenter(fac, m_aaPos) << std::endl);

			if( m_aaMarkFaceHasUnclosedFracSideB[fac] )
			{

				// figure out the vertice(s) where the face is unclosed
				// check if at one of these vertices, two fractures are crossing and one ending
				// by checking if there are at least two segments as precondition, else no matter

				std::vector<Vertex *> assoVrt;

				CollectVertices(assoVrt, m_grid, fac);

				for( Vertex * vrt : assoVrt )
				{
					if( m_aaMarkVrtxHasUnclosedFracB[vrt] )
					{
						VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[vrt];

						if( vecSegVolElmInf.size() > 1 ) // expansion takes place here
						{
//							m_sh.assign_subset(fac,m_sh.num_subsets());
							numEndingCrossingClefts++;
							m_d_endingCrossingCleftFaces.push_back(fac);
							m_d_endingCrossingCleftVrtcs.push_back(vrt);
							// TODO FIXME vielleicht eine Klasse, wo die Faces und Vertices einer
							// jeweiligen ending Crossing Cleft Stelle zusammengefasst werden?
							// statt unabhängige Listen? gemeinsame Liste?
						}

					}
				}


			}

		}

	}

	// TODO FIXME m_d_endingCrossingCleftVrtcs 	m_d_endingCrossingCleftFaces globale member vielleicht?

#if 0
//	for( VertexIterator iterV = m_sel.begin<Vertex>(); iterV != m_sel.end<Vertex>(); iterV++ )
	for( Vertex * vrt : m_d_endingCrossingCleftVrtcs )
	{
//		Vertex * vrt = *iterV;

		// search neigbours of ending crossing cleft faces

		if( m_aaMarkVrtxHasUnclosedFracB[vrt] )
		{
			VecSegmentVolElmInfo const & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[vrt];

			for( SegmentVolElmInfo const & svei : vecSegVolElmInf )
			{
				for( AttachedVolumeElemInfo const & attVolEI : svei )
				{
					VecAttachedFractFaceEdgeSudo vecAttFracFace = attVolEI.getVecFractManifElem();
					// die sind alle von den echten das Segment begrenzenden Faces

					for( AttachedFractFaceEdgeSudo & afes : vecAttFracFace )
					{
						Face * faceSegmLim = afes.getManifElm();

						IndexType sudoSegmLim = afes.getSudo();

						std::pair<Edge*,Edge*> const & edgePair = afes.getPairLowElm();

						// figure out neighbours of ending crossing cleft faces at ending node

						// rule out free edge

						Edge * freeEdge = nullptr;
						Edge * boundedEdge = nullptr;

						for( Face * eccf : m_d_endingCrossingCleftFaces )
						{
							if( FaceContains( eccf, vrt ) )
							{
								if( eccf != fac )
								{
									if( FaceContains( eccf, edgePair.first ) || FaceContains( eccf, edgePair.second ) )
									{

									}
								}
							}
						}
//						if( )
					}

				}


			}
		}
		else
		{
			UG_LOG("only ending crossing cleft vertices allowed here" << std::endl);
			UG_THROW("only ending crossing cleft vertices allowed here" << std::endl);
			return false;
		}

	}

#endif


	UG_LOG("detected ending crossing cleft faces " << numEndingCrossingClefts << std::endl);

	if( numEndingCrossingClefts == 0 )
		return true;

	// debug for ending crossing clefts

	for( Face * fac : m_d_endingCrossingCleftFaces )
	{
		m_sh.assign_subset( fac, m_sh.num_subsets());

	}

	for( Vertex * vrt : m_d_endingCrossingCleftVrtcs )
	{
		m_sh.assign_subset( vrt, m_sh.num_subsets());
	}

	// TODO FIXME unterscheide Faces entlang der expandierten Kreuzungskluft,
	// mit 2 Knoten an Kreuzungspunkten, ein Knoten voll expandiert für beide subsets
	// von solchen, die in der Luft hängen am Ende der Kluft, wo nur die durchgehende Kluft expandiert

	// TODO FIXME need to detect attached faces from the durchgehende cleft which gets expanded
	// as specific structures needed also there

	// TODO FIXME auch die Kante, die an dem auslaufenden Face entlang geht im Schnitt mit dem ausgedehnten,
	// per edge attachment markieren! bool
	// und die anhängenden faces der kreuzenden durchgehenden Kluft ebenso markieren irgendwie per face attachment bool

	return false;
}

#endif

//////////////////////////////////////////////////////////////////

// TODO FIXME diese Funktion macht nur in Teilen das richtige
// die ungeschlossenen Faces zeigen an, dass eine Kluft hier endet
// wird bisher nicht berücksichtigt
// irgendwie muss das markiert werden, damit die Kluft, die zu Ende geht.
// im Durchstich trotzdem berücksichtigt wird
ArteExpandFracs3D::IndexType ArteExpandFracs3D::shiftUnclosedFracFacesToUnclosedFractFaces( Vertex * const & vrt )
{
	IndexType shiftedFracFaces = 0;

	// TODO FIXME still to be implemented - shift those fracture faces which appear
	// twice in the segment volumes, i.e. which are part of two volumes of the
	// segment, i.e. which touch each other, but are not closed,
	// shift them to the general faces, to avoid that they are taken into account
	// for the creation of new vertices, they must not have any influence

	UG_LOG("SHIFT FRAC 2 GENER" << std::endl);

	VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[vrt];

	IndexType segmentNumber = 0;



	for( SegmentVolElmInfo & svei : vecSegVolElmInf )
	{
//		VecAttachedFractFaceEdgeSudo vecAttFractList;

		for( AttachedVolumeElemInfo & attVolEIOne : svei )
		{
			VecAttachedFractFaceEdgeSudo vecAttFracFaceOne = attVolEIOne.getVecFractManifElem();

			for( AttachedVolumeElemInfo & attVolEITwo : svei )
			{
				if( ! attVolEIOne.hasSameFulldimElem( attVolEITwo ) )
				{
					VecAttachedFractFaceEdgeSudo vecAttFracFaceTwo = attVolEITwo.getVecFractManifElem();

					for( AttachedFractFaceEdgeSudo & afesOne : vecAttFracFaceOne )
					{
						for( AttachedFractFaceEdgeSudo & afesTwo : vecAttFracFaceTwo )
						{
							if( afesOne.testIfEquals( afesTwo ) )
							{
								// beide in die generellen Faces verschieben!

								if( ! attVolEIOne.searchFractManifElem( afesOne, true ) )
								{
									UG_THROW("im einen nicht gefunden "<< std::endl);
								}

								if( ! attVolEITwo.searchFractManifElem( afesTwo, true ) )
								{
									UG_THROW("im anderen nicht gefunden "<< std::endl);
								}

								shiftedFracFaces++;

//								Face * unclosedFace = afesOne.getManifElm();
//								// tested if same as that one from afesTwo
//
////								m_aaMarkFaceHasUnclosedFracSideB[ unclosedFace ] = true;
//
//								UG_LOG("unclosed face " << CalculateCenter( unclosedFace, m_aaPos ) << std::endl);

//								m_sh.assign_subset( unclosedFace, m_sh.num_subsets());

//								m_aaMarkVrtxHasUnclosedFracB[vrt] = true;

								//  added vertex attachment that knows if at vertex there is an unclosed fracture


							}
						}
					}
				}
			}


//			for( AttachedFractFaceEdgeSudo & affe : nextVolFacs )
//			{
//				vecAttFractList.push_back( affe );
//			}
			// create a full chain of all fracture faces, and if one appears twice, shift it for both
			// volumes where it is in to the general faces, i.e. count after establishing, and when twice,
			// then shift
		}
	}

	UG_LOG("SHIFT FRAC 2 GENER" << std::endl);



	return shiftedFracFaces;
}

//////////////////////////////////////////////////////////////////



bool ArteExpandFracs3D::seletForSegmented()
{
	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

#if 0

		bool wahl = true;

		// TODO FIXME
		// hier den Stammi-bene Algorithmus einfügen
		// Search the merged manifold interatively between each basic element
		// dazu noch für boundary faces ein eigenes member einführen,
		// das vom Typ General statt fracture manifold ist, aber sonst wie general
		// später wirken die boundary faces wie eine fracture, die aber
		// um den Wert null nur verschoben wird
		// die Anzahl der Segmente bestimmt, ob der Vertex gewählt wird
		// damit er gewählt wird, muss mehr als ein Segment vorhanden sein

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

//		bool isBnd = m_aaMarkVrtVFP[ vrt ].getIsBndFracVertex();
		bool isBnd = vrtxFracPrps.getIsBndFracVertex();
//		auto numCrosFrac = m_aaMarkVrtVFP[ vrt ].getNumberFracEdgesInVertex();

		VertxFracPropts::VrtxFracStatus vfpStatus =  vrtxFracPrps.getVrtxFracStatus();

		if( vfpStatus == VertxFracPropts::noFracSuDoAtt )
			UG_THROW("vertex selected and no frac " << std::endl);

		// TODO FIXME das ist richtig für den 2D Fall, aber passt das auch im 3D Fall???? nein, da SudoFrac Zahl nötig
//		if( ! isBnd && numCrosFrac == 1 )
//		if( ! isBnd && vfpStatus == VertxFracPropts::oneFracSuDoAtt && )
		// TODO FIXME was, wenn ein Teil geschlossen ist der fractures, und ein anderer nicht???
		//static_assert(false);


		// bestimmen, ob die vertizes der fracture vertizes von ihrer subdomain komplett umzingelt werden
		// muss vor  hier geklärt werden!!!

//		VecPairSudoBool sudoIsSourrounded;



		bool allClosed = false;

		allClosed = isVrtxSurroundedByFracFaces( vrt, vrtxFracPrps ); //, sudoIsSourrounded );
//			// case boundary: figure out if the vertex is surrounded by frac faces of which two end in
//			// boundary edges, similar the case when the boundary face has itself two
//			// boundary edges, where the vertex is connected to both of them, then it is easy

//		if( ! allClosed )
//			return false;

//		return true;

//		if( ! isBnd )
//		{
//			UG_LOG("test if closed" << std::endl);
//
////			m_sh.assign_subset(vrt,m_sh.num_subsets());
//
////			UG_LOG("test if closed assign" << std::endl);
//
//
//		}
//		else
//		{
//
//
//
////			allClosed =
//		}
//		// TODO FIXME auch bei einer boundary muss das gemacht werden!

		UG_LOG("getestet if closed " << m_aaPos[vrt] << std::endl);


//		return true;

		if( allClosed == vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
			UG_THROW("da ist was schief gegangen " << std::endl);

		// TODO FIXME ist das so richtig? kann sein, dass das zu vereinfacht ist!!!!!
		// sudo is suourrounded muss übertragen werden TODO FIXME
//		if( ! isBnd && vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
		// das !isBnd im 2D Fall wichtig, hier bestimmt auch, wie verallgemeinern?
		// bei mehreren subdoms eventuell komplizierter, kommt aber hoffentlich am Rand nicht vor......
		if( vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() )
		{
			wahl = false;
		}

#else

		VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[vrt];

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

//		bool isBnd = m_aaMarkVrtVFP[ vrt ].getIsBndFracVertex();
		bool isBnd = vrtxFracPrps.getIsBndFracVertex();

//		bool wahl = ( vecSegVolElmInf.size() > 1 );
//		bool wahl = ( vecSegVolElmInf.size() > 1 && ! isBnd );
		bool wahl = ( vecSegVolElmInf.size() > 1 );
		// TODO FIXME danach vielleicht auch fuer boundary wieder selektieren
		// sobald die Boundary Behandlung wie Pseudo Fracture mit Expansion null
		// aber solange die bounaries nicht behandelt werden, führt
		// die Selektion der Boundary Vertizes zu Problemen

#endif


//		UG_LOG("SELEKTIERE " << m_aaPos[vrt] << " -> " << vrtxFracPrps.getInfoAllFracturesSameClosedState<false>() << std::endl);

		UG_LOG("SELEKTIERE " << m_aaPos[vrt] << " -> " << wahl << std::endl);


		// was, wenn numCrossFrac == 0 ist?
		// wieso werden die boundary vrt ausgeschlossen, oder sollen die nicht ausgeschlossen werden?
		// schon im 2D Fall unklar, hier noch wirrer!!! TODO FIXME

		if( wahl )
//		if( m_aaMarkVrtVFP[vrt].getNumberFracEdgesInVertex() > 1 ) // TODO FIXME stimmt das so?
		{
			//	select all associated edges, faces and volumes
			m_sel.select( m_grid.associated_edges_begin(vrt), m_grid.associated_edges_end(vrt) );
			m_sel.select( m_grid.associated_faces_begin(vrt), m_grid.associated_faces_end(vrt) );
			m_sel.select( m_grid.associated_volumes_begin(vrt), m_grid.associated_volumes_end(vrt) );

			std::vector<Edge*> assoEdg;
			std::vector<Face*> assoFac;
//			std::vector<Volume*> assoVol;
//			VecAttachedVolumeElemInfo assoVolElemInfo;

			for( std::vector<Edge *>::iterator iterEdg = m_grid.associated_edges_begin(vrt);
											   iterEdg != m_grid.associated_edges_end(vrt);
											   iterEdg++ )
			{
				assoEdg.push_back(*iterEdg);
			}

			for( std::vector<Face *>::iterator iterFac = m_grid.associated_faces_begin(vrt);
											   iterFac != m_grid.associated_faces_end(vrt);
											   iterFac++ )
			{
				assoFac.push_back(*iterFac);
			}

			// TODO FIXME das nach oben verschieben, wo der Stammi Bene Algo sein soll
			// die asso edges und faces brauchen wir vielleicht gar nicht
			// bzw asso edges und asso faces können hier bleiben wo gewählt wird
			// die assoVolElemInfo wird schon oben erzeugt vor der Wahl
			// und dann wird die danach folgende Loop Info
//			for( std::vector<Volume *>::iterator iterVol = m_grid.associated_volumes_begin(vrt);
//											   	 iterVol != m_grid.associated_volumes_end(vrt);
//											   	 iterVol++ )
//			{
//				assoVol.push_back(*iterVol);
//
//				AttachedVolumeElemInfo avei(*iterVol);
//
//				assoVolElemInfo.push_back(avei);
//			}

			m_aaVrtInfoAssoEdges[vrt] = assoEdg;
			m_aaVrtInfoAssoFaces[vrt] = assoFac;
//			m_aaVrtInfoAssoVols[vrt] = assoVol;
//			m_aaVolElmInfo[vrt] = assoVolElemInfo;

		}
	}

	UG_LOG("vertex Infos Runde eins fertig " << std::endl);

#if 0
	// Voraussetzung  FÜR StammiBene Aufrufung
	// Stammi-Bene-Vorbereitung
	for( VertexIterator iter = m_sel.begin<Vertex>(); iter != m_sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;

//		std::vector<Volume*> & attVol = m_aaVrtInfoAssoVols[vrt];

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ vrt ];

		VecAttachedFractFaceEdgeSudo vecAttFacEdgSudo = vrtxFracPrps.getAllAttachedFractElems();

		auto & vecVolElmInfo = m_aaVolElmInfo[vrt];

		// TODO FIXME eigentlich eine Dummheit, das auf zu teilen in ein VolElemInfo und ein VrtInfoAssoVols
		// denn die InfoAssoVols haben als Info nur die Volumen, die VolElmInfos haben die Volumen
		// und noch viel mehr Infos zu den Faces und den Edges....
		// mittelfristig die m_aaVrtInfoAssoVols abschaffen und alles auf die AttachedFullDimElemInfo
		// übertragen, dann geht der folgende Loop gleich über den Vektor darüber, bzw. gleichbedeutend
		// über m_aaVolElmInfo
//		for( auto & vol : attVol )
		for( AttachedVolumeElemInfo & attVolElmInfo : vecVolElmInfo )
		{
//			AttachedVolumeElemInfo attVolElmInfo( vol );
			Volume * vol = attVolElmInfo.getFulldimElem();

			// add those faces which are fracture faces
			for( auto & afes : vecAttFacEdgSudo )
			{
				attVolElmInfo.addFractManifElem(afes, m_grid);
			}

			// add those faces which are NOT fracture faces, assign them arbitraryly subdomain  -1
			// to indicate that they are not from the manifold, independent of their subdomain

			// collect all volume faces which incorporate the vertex

			std::vector<Face*> volFacesContainingVrtx;

			for( IndexType iFac = 0; iFac < vol->num_faces(); iFac++ )
			{
				Face * fac = m_grid.get_face(vol,iFac);

				if( FaceContains( fac, vrt ) )
				{
					volFacesContainingVrtx.push_back( fac );
				}
			}

			for( auto const & fac : volFacesContainingVrtx )
			{
				// get the edges of the face connected to the vertex

				std::vector<Edge*> vecEdgesFaceVrtx;

				// need to be two edges always, check

				for( IndexType iEdge = 0; iEdge < fac->num_edges(); iEdge++ )
				{
					Edge * edg = m_grid.get_edge(fac,iEdge);

					if( EdgeContains(edg,vrt) )
					{
						vecEdgesFaceVrtx.push_back(edg);
					}
				}

				if( vecEdgesFaceVrtx.size() != 2 )
				{
					UG_LOG("edge number Unsinn " << vecEdgesFaceVrtx.size() << std::endl);
					UG_THROW("edge number Unsinn " << vecEdgesFaceVrtx.size() << std::endl);
					return false;
				}

				EdgePair edgesFaceVrtx( vecEdgesFaceVrtx[0], vecEdgesFaceVrtx[1] );

				// test the subdomain first, if from the subdomains of the cleft manifolds

				IndexType sudoThisFace = m_sh.get_subset_index(fac);

				std::vector<IndexType> const & sudoList = vrtxFracPrps.getSudoList();

				// test if sudo of face belongs to the fracture face subdom list

				bool belongsToFracFaceSudo = false;

				for( auto const & sudoFrac : sudoList )
				{
					if( sudoFrac == sudoThisFace )
					{
						belongsToFracFaceSudo = true;
						break;
					}
				}


				if( belongsToFracFaceSudo )
				{
					// if it belongs, construct it again and test if it already belongs to the fracture faces
					// MUST be already part of the list, else major error appeared!

					AttachedFractFaceEdgeSudo afesTest( fac, edgesFaceVrtx, sudoThisFace );

					if( attVolElmInfo.addFractManifElem( afesTest, m_grid ) )
					{
						UG_LOG("manifold element already contained!" << std::endl);
						UG_THROW("manifold element already contained!" << std::endl);
						return false;
					}

					// nothing to do, already added before hoffentlich

				}
				else
				{
					// zeitweilig fehlten alle Faces, die keine fractures sind
					// die müssen noch irgendwie als nicht-fracture-faces auch dazu gefügt werden
					// die sind in den attached volumes schon enthalten,
					// Frage: wie prüfen, dass die gar keine fractures sind, die Infos sollten bekannt sein.....
					// wichtig auch, dass nur die faces dazu kommen, die den Vertex enthalten!!!
					// irgendwas von der Art "nonFractureFaceInfos" oder sowas dazu hängen, mit Info
					// ob schon getouched oder noch nicht.....


					// we construct the attached manifold, given that it is NOT a fracture manifold

					// notwendig auch, dass es eine Markierungsmöglichkeit gibt dafür, ob
					// ein face bei der nächsten weiter inneren Runde quasi äussere Begrenzung sein muss
					// gilt sowohl für fracture faces, die können das in erster Runde auch sein am Ende
					// der Runde, sowie danach nur noch für nicht-fracture-faces

					AttachedGenerFaceEdgeSudo afesAdd( fac, edgesFaceVrtx );

					attVolElmInfo.addGenerManifElem( afesAdd, m_grid );

				}

			}

		}
	}
#endif

	UG_LOG("vertex Infos Runde eins fertig Volumen auch" << std::endl);

	return true;
}

bool ArteExpandFracs3D::stasiAlgo( Vertex * const & oldVrt )
{

	UG_LOG("Stasi start " << m_aaPos[oldVrt] << std::endl);
	// TODO FIXME übernehmen von loop2EstablishNewVertices und establishNewVertices
	// plus Boundary faces ähnlich Fracture
	// am Ende die in Segment verbundenen offenen Fracture Faces verschwinden lassen
	// Segmente erstellen Ziel hier, aber auch raus finden, ob es mehr als eines gibt
	// oder ob es offen ist, wegen wahl

	vector3 posOldVrt = m_aaPos[oldVrt];

	UG_LOG("vertex at " << posOldVrt << std::endl);

	VecSegmentVolElmInfo & vecSegVolElmInfo = m_accsAttVecSegVolElmInfo[oldVrt];

	auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];

	bool isBndryVrtx = vrtxFracPrps.getIsBndFracVertex();

	UG_LOG("under construction Tetrahedra limited Stasi Algo" << std::endl);

//	VecVertFracTrip const & vecVertFracTrip = m_aaVrtInfoFraTri[oldVrt];

	VecAttachedVolumeElemInfo assoVolElemInfo;

	int bndryFacsFnd = 0;

	for( std::vector<Volume *>::iterator iterVol = m_grid.associated_volumes_begin(oldVrt);
									   	 iterVol != m_grid.associated_volumes_end(oldVrt);
									   	 iterVol++ )
	{
		Volume * vol = *iterVol;

		AttachedVolumeElemInfo avei(vol);

		bndryFacsFnd += prepareStasi(oldVrt, avei);

		assoVolElemInfo.push_back(avei);
	}

	if( isBndryVrtx && bndryFacsFnd == 0 )
	{
		UG_LOG("Boundary vertex with no boundary faces adjoint" << std::endl);
		UG_THROW("Boundary vertex with no boundary faces adjoint" << std::endl);
		return false;
	}


//	if( vrtxFracPrps.getIsBndFracVertex() )
//	{
//		for( AttachedVolumeElemInfo ave : assoVolElemInfo )
//		{
//			UG_LOG("BONDVERT " << m_aaPos[oldVrt] << " -> " << ave.getVecBndryManifElem().size()  << std::endl);
//			if(  ave.getVecBndryManifElem().size() == 0 )
//			{
//				UG_THROW("VERLUST" << std::endl);
//			}
//		}
//	}

	// von copy und paste angepasst, die unsinnige Verdopplung von vecAttVolElemInfo und assoVolElemInfo
	// vielleicht noch entfernen
	VecAttachedVolumeElemInfo const & vecAttVolElemInfo = assoVolElemInfo; // m_aaVolElmInfo[oldVrt];

	VecAttachedVolumeElemInfo vecAttVolElemInfoCop = vecAttVolElemInfo; // echte KOPIE

	VecAttachedVolumeElemInfo reconstructedVecAttVolElmInf;

		/*
		 * While Schleifen aufbauen für den
		 * Search the adjacent surface interatively - Algorithmus
		 * (Stasi Algorithmus)
		 *
		 */

	IndexType d_segmenteErledigt = 0;

	while( vecAttVolElemInfoCop.size() != 0 )
	{
		SegmentVolElmInfo segmentAVEI;

		AttachedVolumeElemInfo & startVolInfoThisSegment = vecAttVolElemInfoCop[0];

		startVolInfoThisSegment.markIt();

		Volume * volSta = startVolInfoThisSegment.getFulldimElem();

		vector3 center;

		if( volSta != nullptr )
			center = CalculateCenter(volSta,m_aaPos);

//		UG_LOG("volume center " << center << std::endl );

		int d_loopsDone = 0;

		while( vecAttVolElemInfoCop.size() != 0 )
		{
			// count number of marked elements
			IndexType numMarkedElems = 0;
			IndexType markPoint = 0;

//			IndexType lastMarkPt = 0;
			IndexType startIndexInner = 0;

			for( AttachedVolumeElemInfo const & volElInfCop : vecAttVolElemInfoCop )
			{
				if( volElInfCop.isMarked() )
				{
					Volume * vol = volElInfCop.getFulldimElem();
//					m_sh.assign_subset(vol, m_sh.num_subsets());

					vector3 center = CalculateCenter(vol,m_aaPos);

//					UG_LOG("DAS ZENTRUM " << numMarkedElems << " -> " << center << std::endl);

					startIndexInner = markPoint;
					numMarkedElems++;

				}

				markPoint++;
			}

			UG_LOG("LOOPS DONE " << numMarkedElems << std::endl);

			if( numMarkedElems == 0 )
				break;

//			int startIndexInner = -1;
//
//			for( int i = 0; i < vecAttVolElemInfoCop.size(); i++ )
//			{
//				AttachedVolumeElemInfo vi = vecAttVolElemInfoCop[i];
//
//				Volume * vo = vi.getFulldimElem();
//
//				vector3 center = CalculateCenter(vo,m_aaPos);
//
//				UG_LOG("DAS ZENTRUM ZAHL VOR " << i << " -> " <<  center << std::endl);
//
//				if( vi.isMarked() )
//					startIndexInner = i;
//			}
//
//			if( startIndexInner < 0 )
//			{
//				UG_THROW("kein Anfang gefunden " << std::endl);
//			}
//
//#if 0
//			IndexType startIndexInner = markPoint - 1;
//#endif
			AttachedVolumeElemInfo startVolInfoMarkLoop = vecAttVolElemInfoCop[startIndexInner];

			Volume * stattVoll = startVolInfoMarkLoop.getFulldimElem();

			vector3 centerX = CalculateCenter(stattVoll,m_aaPos);

			UG_LOG("DAS ZENTRUM DANACH STASI" << startIndexInner << " -> " <<  centerX << std::endl);

//			m_sh.assign_subset(stattVoll, m_sh.num_subsets());
#if 0
			for( int i = 0; i < vecAttVolElemInfoCop.size(); i++ )
			{
				AttachedVolumeElemInfo vi = vecAttVolElemInfoCop[i];

				Volume * vo = vi.getFulldimElem();

				vector3 center = CalculateCenter(vo,m_aaPos);

				UG_LOG("DAS ZENTRUM ZAHL " << i << " -> " <<  center << std::endl);

			}
#endif
			for( AttachedVolumeElemInfo const & possibleOrigVolInfo : vecAttVolElemInfo )
			{
				if( possibleOrigVolInfo.hasSameFulldimElem( startVolInfoMarkLoop ) )
				{
					segmentAVEI.push_back(possibleOrigVolInfo);
					reconstructedVecAttVolElmInf.push_back(possibleOrigVolInfo);
					break;
				}
			}

			vecAttVolElemInfoCop.erase( vecAttVolElemInfoCop.begin() + startIndexInner );

//			if( d_loopsDone == 1 )
//				return false;

			for( VecAttachedVolumeElemInfo::iterator aveiIt = vecAttVolElemInfoCop.begin();
													 aveiIt < vecAttVolElemInfoCop.end();
													 aveiIt++
			)
			{
				AttachedVolumeElemInfo & possibleNeighbour = *aveiIt;

				if( possibleNeighbour.hasSameFulldimElem( startVolInfoMarkLoop ) )
				{
					continue;
				}
				else
				{
					bool neighbourFound = possibleNeighbour.testFullDimElmNeighbour( startVolInfoMarkLoop );

					if( neighbourFound )
					{
						Volume * vol = possibleNeighbour.getFulldimElem();

//						m_sh.assign_subset(vol, m_sh.num_subsets());

					}
				}
			}


			d_loopsDone++;


		}

		vecSegVolElmInfo.push_back(segmentAVEI);

//		d_segmenteErledigt++;
//
//		if( d_segmenteErledigt == 1 )
//		return false;
	}

	if( reconstructedVecAttVolElmInf.size() != vecAttVolElemInfo.size() )
	{
		UG_LOG("Rekonstruktion schief gegangen " << std::endl);
		UG_THROW("Rekonstruktion schief gegangen " << std::endl);
		return false;
	}

// for debug purposes

	constexpr bool d_assignSudos2Segments = false;

	if( d_assignSudos2Segments )
	{
		if( vecSegVolElmInfo.size() > 1 )
		{
			for( SegmentVolElmInfo const & svei : vecSegVolElmInfo )
			{
				// TODO FIXME das hier wieder entfernen, die Subdomain Zuweisung, nur für debug Zwecke
				IndexType sudoMax = m_sh.num_subsets();

				for( AttachedVolumeElemInfo const & vei : svei )
				{
					Volume * vol = vei.getFulldimElem();

					m_sh.assign_subset( vol, sudoMax );
				}
			}
		}
	}

	UG_LOG("Stasi END " << m_aaPos[oldVrt] << std::endl);

	return true;
}

#if 0
// Deprecated due to Stasi Algo
// herausfinden für Sudo der frac, ob bezüglich dieser sudo die faces geschlossen sind, oder ob ein Fracture End vorliegt
bool ArteExpandFracs3D::isVrtxSurroundedByFracFaces( Vertex * const & vrt, VertxFracPropts & vrtxFracPrps )
//, VecPairSudoBool & sudoSurrounded )
{
	// TODO FIXME wenn an Boundary, dann auch auf closed open unterscheiden - sowohl wenn nur edge an
	// boundary, aber auch wenn ein ganzes face dort, noch unklar, was das bedeutet

	// ganz ähnlich wie im 2D Fall, Loopen, im Kreis drehen, kucken, ob wir vom Anfang ans Ende kommen,
	// und ob das Ende der edges wieder der Anfang der edges ist, da wir uns auf faces drehen

	// case boundary: figure out if the vertex is surrounded by frac faces of which two end in
	// boundary edges, similar the case when the boundary face has itself two
	// boundary edges, where the vertex is connected to both of them, then it is easy


	VecAttachedFractFaceEdgeSudo vafes = vrtxFracPrps.getAllAttachedFractElems();

	std::vector<IndexType> sudoList = vrtxFracPrps.getSudoList();

//	for( auto const & sudo : sudoList )
//	{
//		std::vector<Face*> tmpFaces;
//
//		CollectFaces( tmpFace, m_grid, vrt );
//
//	}

	// first compare sudo list, if equal

	std::vector<IndexType> sudosTestList;

	std::vector<bool> sudoFound( sudoList.size(), false );

	UG_LOG("sudo list size " << sudoList.size() << std::endl );

	UG_LOG("vafes list size VA " << vafes.size() << std::endl );


	for( auto const & af : vafes )
	{
		bool found = false;

		IndexType sudo = af.getSudo();

		UG_LOG("sudo af " << sudo << std::endl);

		for( IndexType i = 0; i < sudoList.size(); i++ )
		{
			UG_LOG("sudo list dusso is " << sudoList[i] << std::endl);

			if( sudo == sudoList[i] )
			{
				sudoFound[i] = true;
				found = true;
			}
		}


		if( ! found )
			UG_THROW("sudo nicht gefunden muss aber da sein " << std::endl);
	}

	UG_LOG("alles gefunden " << std::endl);


	for( auto const & sf: sudoFound )
	{
		if( sf == false )
		{
			UG_LOG("Falsch" << std::endl);
			UG_THROW("sudo not found but must be there " << std::endl);
		}
	}

	UG_LOG("alles gefunden Test " << std::endl);


	// sort faces with respect to subdomain - macht das wirklich Sinn, wie umständlich das gemacht wird jetzt?

	// if we arrive here, all sudos found, the entire circle closing can start

	VecPairSudoBool sudoSurrounded;

	for( auto const & sudo : sudoList )
	{
		VecAttachedFractFaceEdgeSudo vecAttFacSudo;

		for( auto const & attFac : vafes )
		{
			if( sudo == attFac.getSudo() )
			{
				UG_LOG("die sudo gefunden " << sudo << std::endl);
				vecAttFacSudo.push_back(attFac);
			}
		}

		VecAttachedFractFaceEdgeSudo vecAttFacSudoSort;



		bool isClosed = false;

		bool isBndVrtx = IsBoundaryVertex3D(m_grid,vrt);

		if( ! isBndVrtx )
		{
			UG_LOG("No boundary vertex test closed " << m_aaPos[vrt] << std::endl);

			if( vecAttFacSudo.size() == 1 )
			{
				// no need for further investigations for inner faces
				isClosed = false;
				vecAttFacSudoSort = vecAttFacSudo;
			}
			else
			{
				isClosed = sortElemCircleIsClosed( vecAttFacSudo, vecAttFacSudoSort );
			}

		}
		else // different treatment boundary vertex
		{
			// figure out start face with a boundary edge, and figure out an eventual additional boundary edge
			// if only one face, then check if two attached boundary edges, then true, else false

			if( vecAttFacSudo.size() == 1 )
			{
				AttachedFractFaceEdgeSudo & singleEntry = vecAttFacSudo[0];

				EdgePair faceEdgs = singleEntry.getPairLowElm();

				if( IsBoundaryEdge3D(m_grid, faceEdgs.first) && IsBoundaryEdge3D(m_grid, faceEdgs.second ) )
					isClosed = true;

				// very simple
				vecAttFacSudoSort = vecAttFacSudo;

			}
			else
			{
				// figure out a begin face with a boundary edge and figure out another face with a boundary edge

				Edge * beginEdge = nullptr;
				Edge * endEdge = nullptr;

				IndexType startFaceIndx = 0;
				IndexType endFaceIndx = 0;

				for( auto const & afs : vecAttFacSudo )
				{
					Face * fac = afs.getManifElm();

					EdgePair edgs = afs.getPairLowElm();

					Edge * edgOne = edgs.first;
					Edge * edgTwo = edgs.second;

					if( beginEdge == nullptr )
					{
						if( IsBoundaryEdge3D(m_grid, edgOne) )
						{
							beginEdge = edgTwo;
						}
						else if( IsBoundaryEdge3D(m_grid, edgTwo) )
						{
							beginEdge = edgOne;
						}
						else
						{
							startFaceIndx++;
						}
					}
					else
					{
						if( IsBoundaryEdge3D(m_grid, edgOne) )
						{
							endEdge = edgOne;
						}
						else if( IsBoundaryEdge3D(m_grid, edgTwo) )
						{
							endEdge = edgTwo;
						}
					}

					if( endEdge != nullptr )
						break;

					endFaceIndx++;
				}

				if( beginEdge == nullptr && endFaceIndx == vecAttFacSudo.size() )
//					|| beginEdge == nullptr || endEdge == nullptr )
				{
					UG_LOG("keine boundary edges" << std::endl);

					startFaceIndx = -1;

					if( endEdge != nullptr )
						UG_THROW("Ende nicht null, Anfang null " << std::endl);
				}

				UG_LOG("Boundary vertex test closed " << m_aaPos[vrt] << std::endl);

//				int d_num = 0;
//				for( auto const & afs : vecAttFacSudo )
//				{
//					d_num++;
//					Face * fac = afs.getManifElm();
//					m_sh.assign_subset(fac, m_sh.num_subsets());
//				}
//				UG_LOG("number of surr faces " << d_num << std::endl );

//				m_sh.assign_subset(beginEdge,m_sh.num_subsets());
//				m_sh.assign_subset(endEdge,m_sh.num_subsets());
//				m_sh.assign_subset(vecAttFacSudo[startFaceIndx].getManifElm(),m_sh.num_subsets());

				isClosed = sortElemCircleIsClosed( vecAttFacSudo, vecAttFacSudoSort, startFaceIndx, beginEdge, endEdge );

			}
		}

		UG_LOG("-------------------------------" << std::endl);
		UG_LOG("is closed " << isClosed << " at " << m_aaPos[vrt] << std::endl);
		UG_LOG("-------------------------------" << std::endl);

//		if( isClosed )
//		{
//			m_sh.assign_subset(vrt,3);
//		}
//		else
//		{
//			m_sh.assign_subset(vrt,4);
//		}

		if( vecAttFacSudo.size() != vecAttFacSudoSort.size() )
		{
//			return false;

			UG_THROW("Die Sortierung ist komplett schief gegangen " << std::endl);
		}

		// DEBUG Zeug, später entfernen!!!!!!
//		for( auto const & afss : vecAttFacSudoSort )
//		{
//			Face * fac = afss.getManifElm();
//
//			m_sh.assign_subset(fac, m_sh.num_subsets());
//		}

		PairSudoBool ic( sudo, isClosed );

		sudoSurrounded.push_back( ic );

	}

	vrtxFracPrps.setInfoAllFractureSudosIfClosed(sudoSurrounded);

	bool allClosed = true;

	for( auto const & ic : sudoSurrounded )
	{
		if( ! ic.second )
			allClosed = false;
	}

	return allClosed;
}

bool ArteExpandFracs3D::sortElemCircleIsClosed( VecAttachedFractFaceEdgeSudo const & vecAttFac,
												VecAttachedFractFaceEdgeSudo & vecSortedFac,
												int startFacIndexUser,
//												int endFacIndexUser,
//												IndexType startEdgeIndexUser,
//												IndexType endEdgeIndexUser
//												Face * const & startFacUser,
//												Face * const & endFacUser,
												Edge * const & startEdgUser,
												Edge * const & endEdgUser
												)
{

	UG_LOG("Schliesstest" << std::endl);

	IndexType originalSize = vecAttFac.size();

	if( originalSize == 0 )
	{
		UG_THROW("zu klein zum sortieren " << std::endl);
		return false;
	}
	else if ( originalSize == 1 )
	{
		UG_THROW("should not happen size 1, should have been mentioned before " << std::endl);
	}

	UG_LOG("Kopieren zwecks sortieren " << std::endl);

	for( auto const & af : vecAttFac )
	{
		UG_LOG("die sudos innen sind " << af.getSudo() << std::endl);
	}

	VecAttachedFractFaceEdgeSudo copyVecAttFac = vecAttFac;

	for( auto const & af : copyVecAttFac )
	{
		UG_LOG("die sudos kopiert sind " << af.getSudo() << std::endl);
	}

	IndexType beginIndx = 0;

	UG_LOG("begin Index " << beginIndx << std::endl);

	bool simpleConnectionTest = false;

	Edge * startEdgeForced = nullptr;

	if( startFacIndexUser >= 0 )
	{
		UG_LOG("Veränderung " << startFacIndexUser << std::endl);
		beginIndx = startFacIndexUser;
		simpleConnectionTest = true; // user hopefully did it
	}
	else
	{
		// we need to ensure to start at a fracture which is not in between, in case that circle not closed
		// so ensure that all fracture faces have a fracture face at both edges as neighbour,
		// in principle this is already sufficient to answer the question which we want to know
		// all the rest here is useless playing in principle

		bool broken = false;

		for( IndexType i = 0; i < vecAttFac.size(); i++ )
		{
			IndexType firstSideConnected = 0;
			IndexType secondSideConnected = 0;

			AttachedFractFaceEdgeSudo afBase = vecAttFac[i];

			Face * faceBase = afBase.getManifElm();
			EdgePair edgPairBase = afBase.getPairLowElm();

			Edge * edgeBaseOne = edgPairBase.first;
			Edge * edgeBaseTwo = edgPairBase.second;

			for( IndexType j = 0; j < vecAttFac.size(); j++ )
			{
				if( i != j )
				{
					AttachedFractFaceEdgeSudo afCompr = vecAttFac[j];

					Face * faceCompr = afCompr.getManifElm();
					EdgePair edgPairCompr = afCompr.getPairLowElm();

					Edge * edgeComprOne = edgPairCompr.first;
					Edge * edgeComprTwo = edgPairCompr.second;

					if( edgeComprOne == edgeBaseOne || edgeComprTwo == edgeBaseOne )
						firstSideConnected++;

					if( edgeComprOne == edgeBaseTwo || edgeComprTwo == edgeBaseTwo )
						secondSideConnected++;
				}

			}

			if( vecAttFac.size() > 2 && ( firstSideConnected > 1 || secondSideConnected > 1 ) )
			{
				UG_THROW("zu oft verbunden " << std::endl );
			}
			else if( vecAttFac.size() == 2 && ( firstSideConnected > 2 || secondSideConnected > 2 ) )
			{
				UG_THROW("zu oft verbunden " << std::endl );
			}
			else if( firstSideConnected == 0 )
			{
				// face is open into one direction, already clear that not closed!!

				beginIndx = i;
				simpleConnectionTest = false;
				startEdgeForced = edgeBaseTwo;
				UG_LOG("forcieren 1 " << std::endl);
				broken = true;
				break;
			}
			else if( secondSideConnected == 0 )
			{
				// face is open into one direction, already clear that not closed!!

				beginIndx = i;
				simpleConnectionTest = false;
				startEdgeForced = edgeBaseOne;
				UG_LOG("forcieren 2 " << std::endl);
				broken = true;
				break;
			}
			else if( firstSideConnected == 1 && secondSideConnected == 1 )
			{
				simpleConnectionTest = true; // as long as a look
			}
			else
			{
				UG_THROW("komischer Verbindungsfall " << std::endl);
			}

			if( broken )
				break;

		}
	}

	UG_LOG("begin Index X " << beginIndx << std::endl);


//	AttachedFractFaceEdgeSudo initialAFES = *(copyAttFac.begin());
	AttachedFractFaceEdgeSudo initialAFES = copyVecAttFac[beginIndx];

	IndexType sudo = initialAFES.getSudo();

	UG_LOG("sudo " << beginIndx << " ist " << sudo << std::endl);

	Face * beginFacLoop = initialAFES.getManifElm();

	// TODO FIXME wird das irgendwo verwendet? wieso nicht?
//	Face * endFacLoop = nullptr;

//	if( endFacIndexUser != -1 )
//	{
//		endFacLoop = copyVecAttFac[endFacIndexUser].getManifElm();
//	}

	EdgePair beginEdges = initialAFES.getPairLowElm();

	Edge * beginEdgeLoop = beginEdges.second;
	Edge * targetedEndEdgeLoop = beginEdges.first; // should be closed! should end at same edge as it begins!

//
//	return true;

//	if( startEdgeIndexUser != -1 )
//	{
//		beginEdgeLoop =
//	}

//	if( startFacUser != nullptr )
//	{
//		beginFacLoop = startFacUser;
//
//	}
//

//	if( ! FaceContains(beginFacLoop,beginEdgeLoop->vertex(0)) )
//		UG_THROW("Sortierung Gesicht hat nicht die gewünschte Ecke " << std::endl);
//
//	if( endFacUser != nullptr )
//	{
//		endFacLoop = endFacUser;
//	}
//

	if( startEdgUser != nullptr )
	{
		beginEdgeLoop = startEdgUser;

		// check if part of the begin face!

		if( ! FaceContains( beginFacLoop, beginEdgeLoop ) )
			UG_THROW("Anfangsgesicht hat keine Anfangsecke " << std::endl);

	}

	if( startEdgeForced != nullptr )
	{
		beginEdgeLoop = startEdgeForced;

		// check if part of the begin face!

		if( ! FaceContains( beginFacLoop, beginEdgeLoop ) )
			UG_THROW("Anfangsgesicht hat keine Anfangsecke forced " << std::endl);

//		m_sh.assign_subset(startEdgeForced, m_sh.num_subsets());
		UG_LOG("forciert " << std::endl);

//		return false;

	}


	if( endEdgUser != nullptr )
	{
		// check if part of end face at end of loop somehow
		targetedEndEdgeLoop = endEdgUser;

		if( startEdgeForced != nullptr )
			UG_THROW("das muss schief gehen, Chaos " << std::endl);

//		if( endFacLoop != nullptr )
//		{
//			if( ! FaceContains( endFacLoop, targetedEndEdgeLoop ) )
//				UG_THROW("Endgesicht hat keine Endecke " << std::endl);
//		}
	}
//

	// DEBUG
//	m_sh.assign_subset(beginFacLoop, m_sh.num_subsets());
//	m_sh.assign_subset(beginEdgeLoop, m_sh.num_subsets());
//	m_sh.assign_subset(targetedEndEdgeLoop, m_sh.num_subsets());


	// Du musst sortieren. Du musst einen Zeitplan machen. Das kann man lernen. Du kannst das selber machen.

	IndexType countedCrossedFaces = 1;

	vecSortedFac.push_back( initialAFES );

	copyVecAttFac.erase( copyVecAttFac.begin() + beginIndx );

//	Face * face2Append = beginFacLoop;
	Edge * startEdge2Append = beginEdgeLoop;

//	m_sh.assign_subset(startEdge2Append,m_sh.num_subsets());

	UG_LOG("while loop anfangen " << std::endl);

	IndexType d_whi = 0;

	while( copyVecAttFac.size() != 0 )
	{

		UG_LOG("in while loop " << d_whi << std::endl);
		d_whi++;

		IndexType foundCommEdg = 0;

		Edge * nextEdge = nullptr;

//		for( auto const & caf : copyVecAttFac )
		for( VecAttachedFractFaceEdgeSudo::iterator itAttFES  = copyVecAttFac.begin();
											   itAttFES != copyVecAttFac.end();
											   itAttFES++
		)
		{
			AttachedFractFaceEdgeSudo caf = *itAttFES;

			Face * d_Fac = caf.getManifElm();

//			m_sh.assign_subset(d_Fac,m_sh.num_subsets());

			EdgePair edgPr = caf.getPairLowElm();

			Edge * edgOne = edgPr.first;
			Edge * edgTwo = edgPr.second;

//			m_sh.assign_subset(edgOne,m_sh.num_subsets());
//			m_sh.assign_subset(edgTwo,m_sh.num_subsets());

//			return true;

			IndexType hasEdge = 0;

			Edge * overNextEdge = nullptr;

			if( edgOne == startEdge2Append )
			{
				nextEdge = edgOne;
				overNextEdge = edgTwo;
				hasEdge++;
			}

			if( edgTwo == startEdge2Append )
			{
				nextEdge = edgTwo;
				overNextEdge = edgOne;
				hasEdge++;
			}

			if( hasEdge > 1 )
				UG_THROW("zu viele Ecken und Kanten " << std::endl);

			if( hasEdge == 1 )
			{
				Face * fac2App = caf.getManifElm();
//				m_sh.assign_subset(fac2App,m_sh.num_subsets());
				EdgePair edgesNextFace( edgOne, edgTwo );
				AttachedFractFaceEdgeSudo nextAttFES( fac2App, edgesNextFace, sudo );

				vecSortedFac.push_back(nextAttFES);

				copyVecAttFac.erase(itAttFES);
				foundCommEdg++;
				startEdge2Append = overNextEdge;

				break;
			}


		}

		if( foundCommEdg > 1 )
			UG_THROW("Kein Anschluss unter dieser Nummer " << std::endl);

		if( foundCommEdg == 0 )
		{
			UG_LOG("Kreislauf nicht geschlossen " << std::endl);

			if( nextEdge != nullptr )
				UG_THROW("nicht konsistent, null und null " << std::endl);

			break;
		}

		if( nextEdge == nullptr )
		{
			if( foundCommEdg != 0 )
				UG_THROW("nicht konsistent, null und null v2 " << foundCommEdg << " -> " << nextEdge << std::endl);

	//			return false;
		}

	}

	if( originalSize != vecSortedFac.size() )
	{
		UG_THROW("Sortierung hat nicht funktioniert " << std::endl);
		return false;
	}

	if( startEdge2Append != targetedEndEdgeLoop )
	{
		if( simpleConnectionTest )
			UG_THROW("obwohl offen oder vorgegeben trotzdem Ziel nicht erreicht?" << std::endl);

		UG_LOG("Ende nicht erreicht, Kreis nicht geschlossen, innerer Rand vermutlich" << std::endl);
		return false;
	}


	UG_LOG("Kreislauf Faces 3D Test ob Knoten umrandet geschlossen " << std::endl);

	return true;
}

#endif

bool ArteExpandFracs3D::assignOrigFracInfos()
{
	m_originalFractureFaces.clear();

	for( FaceIterator iter = m_sel.begin<Face>(); iter != m_sel.end<Face>(); ++iter)
	{
		if( m_aaMarkFaceIsFracB[*iter] == true )
			m_originalFractureFaces.push_back(*iter);
	}

	m_fracInfosBySubset = std::vector<FractureInfo>( m_sh.num_subsets(), FractureInfo(-1, -1, 0) );

	for(size_t i = 0; i < m_fracInfos.size(); ++i)
	{
		if( m_fracInfos[i].subsetIndex >= m_sh.num_subsets())
		{
			throw(UGError("Bad subsetIndex in given fracInfos."));
		}

		m_fracInfosBySubset[ m_fracInfos[i].subsetIndex ] = m_fracInfos[i];
	}

//	disable selection inheritance to avoid infinite recursion.
	m_sel.enable_selection_inheritance(false);

	return true;
}

bool ArteExpandFracs3D::establishNewVrtBase()
{
	//	iterate over all surrounding volumes to enable shifted vertices, this loop taken from SR but shortened

	for( VolumeIterator iterSurrVol = m_sel.volumes_begin(); iterSurrVol != m_sel.volumes_end(); ++ iterSurrVol )
	{
		Volume* sv = *iterSurrVol;

		std::vector<Vertex*> & newVrts = m_aaVrtVecVol[sv];
		newVrts.resize(sv->num_vertices());

		for(size_t iVrt = 0; iVrt < sv->num_vertices(); ++ iVrt )
		{
			newVrts[iVrt] = nullptr;
		}

		// erstmal so tun, als ob keine neuen Vertizes erzeugt werden an den alten Vertizes

	}

	return true;
}

#if 0
// Analogon zu VertrexFractureInfo in 2D, wo jeder Vertex eine Liste bekommt, wo alle die ihm angehängten
// Ecken, Faces und Volumen gespeichert werden; dazu die Normalen, und vielleicht noch weitere Infos
bool ArteExpandFracs3D::generateVertexInfos()
{
	UG_LOG("Starte Generierung" << std::endl);

	// TODO FIXME das wird benötigt

	// sowas von der Art als attachement bei den attachments, und dann mit Leben füllen für jeden Vertex
	// in dieser Funktion;
	// vielleicht braucht es auch Edge Infos, oder nur Edge Infos?
//	VecVertFracTrip vertexNoInfo;
//	using AttVecVertFracTrip = Attachment<VecVertFracTrip>;
//	AttVecVertFracTrip aAdjInfoAVVFT;
//	grid.attach_to_vertices_dv( aAdjInfoAVVFT, vertexNoInfo );
//	Grid::VertexAttachmentAccessor<AttVecVertFracTrip> aaVrtInfoFraTri(grid,  aAdjInfoAVVFT );

	// Lebendigmachung in:
	// 	for( auto & fsfpmv : fracSubdom_facePerpendMinVal ) .....
	// von dort lernen!!!!!

	// notwendig: face, normal, volume, edge

	// TODO FIXME das wollen wir nicht, sondern das alte Vertex Fracture Triple
//	m_vrtxFractrQuadrplVec = VrtxFractrQuadrplArte3DVec();
	// TODO FIXME diese Einträge erzeugen

	// TODO FIXME kann vielleicht vereinfacht werden durch einen Loop über alle Vertizes,
	// und das Abfragen der dort schon gespeicherten vertex property Geschichten, sonst
	// wird manches doppelt gemoppelt
	
	for( auto const & fracInf : m_fracInfos )
	{
		IndexType fracSudo = fracInf.subsetIndex;


//		for(EdgeIterator iterEdg = m_sh.begin<Edge>(fracIndSudo); iterEdg != m_sh.end<Edge>(fracIndSudo); iterEdg++ )
//
		for( FaceIterator iterFac = m_sh.begin<Face>(fracSudo); iterFac != m_sh.end<Face>(fracSudo); iterFac++ )
		{

//			VrtxFractrQuadrplArte3D vrtxFractrQuadrpl;
			
			// TODO FIXME die Innereien dieses Loops irgendwie für boundary faces hinbiegen,
			// vermutlich nicht viel verändern, dafür eigene Routine, die dann für jeweils ein face den
			// Müll umsetzt

			Face* fac = *iterFac;
			
			auto sudoFacInnerLoop = m_sh.get_subset_index(fac);
			
			if( sudoFacInnerLoop != fracSudo )
				UG_THROW("Subdomain Index Fehler 3D " << std::endl);

//			std::vector<Vertex*> verticesFac;
//
//			for( IndexType i = 0; i < fac->num_vertices(); i++ )
//			{
//				verticesFac.push_back( fac->vertex(i) );
//			}


			std::vector<Volume*> assoVols;

			// wo kam denn der Käse her?
//			if( ! m_grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES) )
//				UG_THROW("How to collect asso vols?" << std::endl);

			// brauchen wir das? für was? von SR irgendwie übernommen, wo dort was entfernt ähnliches gemacht wird....
			if(! m_grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES) )
			{
				UG_LOG("WARNING grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
				m_grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
			}


//			for( Grid::AssociatedVolumeIterator volIt  = m_grid.associated_volumes_begin(fac);
//												volIt != m_grid.associated_volumes_end(fac);
//												volIt++
//			)
//			{
//				assoVols.push_back(*volIt);
//			}

			CollectVolumes(assoVols, m_grid, fac );
//
//			for( auto const & av : assoVols )
//			{
//				m_sh.assign_subset(av, m_sh.num_subsets());
//			}
//
//			return true;

//			using VolNormPair = std::pair< Volume*, vector3 >;
//
//			using VecVolNormPair = std::vector<VolNormPair>;
//
//			VecVolNormPair vecNormalsAwayVol;

//			UG_LOG("Center Face " << CalculateCenter(fac,m_aaPos) << std::endl);

			for( auto const & kuhVol : assoVols )
			{
				bool facFound = false;
				IndexType numFd = 0;

				for( IndexType iSide = 0; iSide < kuhVol->num_sides(); iSide++ )
				{
					Face * kuhFac = m_grid.get_side(kuhVol, iSide);

//					UG_LOG("Center Kuh Face " << CalculateCenter(kuhFac,m_aaPos) << std::endl);

					// TODO FIXME eigentliches Ziel ist es, den face descriptor des Volumens zu finden,
					// das mit dem face übereinstimmt, alle anderen Seiten des Volumens sind egal
					// Funktion suchen, die ausgehend von einem Face, das ein Volumen begrenzt,
					// den zum Volumen passenden FaceDescriptor findet, also auch richtige Orientierung
					// der Vertices beinhaltet
					// FRAGE TODO FIXME ist ein face descriptor von der Orientierung zum Volumen abhängig
					// oder hängt der nur vom Face ab, das eine vorgegebene Oriertierung hat?

					bool checkCoincide = checkIfFacesVerticesCoincide( kuhFac, fac );

					if( checkCoincide )
					{
						facFound = true;
						numFd++;

						if( kuhFac != fac )
							UG_LOG("Kuh Fac ist nicht fac " << std::endl);

						FaceDescriptor facDescr;

						// TODO FIXME testen, ob der Face Descriptor von der Orientierung abhängt
						// also testen, ob sich der face descriptor ändert, wenn das Volumen
						// auf der einen und auf der anderen Seite des faces hängt
						// deswegen auch die ganze Prozedur mit den kuhFacs, die hoffentlich
						// je nach Volumen anders orientiert sind als das eigentliche Face,
						// aber dieselben Vertices haben, also geometrisch gleich sind, aber anders orientiert!!!!

						// TODO FIXME andere Hergehensweise vielleicht:
						// von m_aaVrtInfoAssoVols ausgehen, darüber loopen, oder die in einen Vektor stecken,
						// wo die Vertices dabei sind, dann kann man sich vielelicht ein paar Klimmzüge sparen,
						// vielleicht aber auch nicht.....

						kuhVol->face_desc( iSide, facDescr );

						vector3 normal;

						CalculateNormal( normal, & facDescr, m_aaPos );

						vector3 facCenter = CalculateCenter( kuhFac, m_aaPos );
						vector3 kuhCenter = CalculateCenter( fac, m_aaPos );
						vector3 kuhVolCenter = CalculateCenter( kuhVol, m_aaPos);

//						UG_LOG("Normale zum face descriptor " << normal << " , " << facCenter << std::endl);
//						UG_LOG("Normale zum Kuhh descriptor " << normal << " , " << kuhCenter << std::endl);
//						UG_LOG("Zentrum des Vol")

//						UG_LOG("fac " << fac << std::endl );
//						UG_LOG("kuh " << kuhFac << std::endl );

						UG_LOG("Normale ist " << normal << " fac " << facCenter
								<< " vol " << kuhVolCenter << std::endl);


//						VolNormPair normalsAwayVol( kuhVol, normal );
//
//						vecNormalsAwayVol.push_back( normalsAwayVol );

						std::vector<Edge*> facEdgs;

						CollectEdges( facEdgs, m_grid, fac);

						for( IndexType iF = 0; iF < fac->num_vertices(); iF++ )
						{
							Vertex * vrt = fac->vertex(iF);

			//				verticesFac.push_back(vrt);

							std::vector<Edge*> edgOfVrtx;

							for( auto const & edg : facEdgs )
							{
								if( EdgeContains(edg, vrt) )
								{
									edgOfVrtx.push_back(edg);
								}
							}

							if( edgOfVrtx.size() == 2 )
							{
								EdgePair commonEdges(edgOfVrtx[0], edgOfVrtx[1]); //  fill values
								// edges commun between face and volume, with the vertex included as well, i.e. two possibilities

								// TODO FIXME diese Info muss woanders hin, in der AttachedFractFaceEdgeSudo Klasse speichern!
								VertFracTrip infoVerticesThisFace( fac, fracSudo, kuhVol, normal, commonEdges );

								// TODO FIXME hier irgendwie graphische Ausgabe von irgendwas

								m_aaVrtInfoFraTri[vrt].push_back( infoVerticesThisFace );

								// DEBUG OUTPUT; REMOVE LATER
//								m_sh.assign_subset(kuhVol,m_sh.num_subsets());
							}
							else
							{
								UG_THROW("Mein Face das hat keine Ecken, keine Ecken hat mein Face" << std::endl);
							}
						}

					}
				}

				if( ! facFound || numFd != 1 )
				{
					UG_THROW("Kein Kuh Volumen gefunden" << std::endl);
					return false;
				}

			}

		}


	}

	UG_LOG("GEnerierung gelungen " << std::endl);

	return true;
}
#endif

bool ArteExpandFracs3D::checkIfFacesVerticesCoincide( Face * const & facOne, Face * const & facTwo )
{

	if( facOne->size() != facTwo->size() )
		return false;

	std::vector<Vertex* > facOneVrtcs, facTwoVrtcs;

	collectFaceVertices( facOneVrtcs, facOne );
	collectFaceVertices( facTwoVrtcs, facTwo );

	for( auto const & vrtOne : facOneVrtcs )
	{
		bool found = false;

		IndexType numFd = 0;

		for( auto const & vrtTwo : facTwoVrtcs )
		{
			if( vrtOne == vrtTwo )
			{
				found = true;
				numFd++;
			}
		}

		if( ! found || numFd != 1 )
			return false;
	}

	return true;
}

bool ArteExpandFracs3D::collectFaceVertices( std::vector<Vertex*> & facVrt, Face * const & fac )
{
	if( fac == nullptr )
		return false;

	facVrt.clear();

	for( IndexType iF = 0; iF < fac->num_vertices(); iF++ )
	{
		Vertex * vrt = fac->vertex(iF);

		facVrt.push_back( vrt );
	}

	return true;
}




// major function of new grid generation, in Keil Style, but functional grid, only the diamonds have to be
// established in additional functionalities independent of this function
bool ArteExpandFracs3D::loop2EstablishNewVertices()
{
	m_vecCrossVrtInf = std::vector<CrossVertInf>();


	// TODO FIXME sowas von der Art wird nötig sein als Vorbereitung für die Diamanten,
	// Infos darin speichern, vielleicht auch noch notwendig, die Kanten zu speichern oder die faces,
	// zu klären im Laufe der Implementation
	// std::vector<CrossVertInf > vecCrossVrtInf;

	// zentraler Loop

	// TODO FIXME vorher noch den attVrtVec und den aaVrtVecFac analog implementieren, das fehlt noch!!!
	// ebenso seine Befüllung, braucht noch eine Funktion dazwischen, die attachments selber in die
	// attach Funktion natürlich
	
	int untilVrt = 0;

	for( VertexIterator iterV = m_sel.begin<Vertex>(); iterV != m_sel.end<Vertex>(); ++ iterV )
	{
		Vertex * oldVrt = *iterV;

		// Position dieses Vertex
		vector3 posOldVrt = m_aaPos[oldVrt];

		// TODO FIXME diese Funktion mit Leben und Analytischer Geometrie 13. Klasse füllen

//		VecVertFracTrip & vecVertFracTrip = m_aaVrtInfoFraTri[oldVrt];
//
//		std::vector<Edge*> & allAssoEdges = m_aaVrtInfoAssoEdges[oldVrt];
//		std::vector<Face*> & allAssoFaces = m_aaVrtInfoAssoFaces[oldVrt];
//		std::vector<Volume*> & allAssoVolumes = m_aaVrtInfoAssoVols[oldVrt];

		UG_LOG("vertex at " << posOldVrt << std::endl);

		VecSegmentLimitingSides & vecSegmLimSid = m_vrtxAttAccsVecSegmLimSid[oldVrt];

		for( SegmentLimitingSides const & sls : vecSegmLimSid )
		{
			if( sls.hasUnclosedFaces())
				UG_LOG("Hier hat das Gesicht noch ungeschlossene " << std::endl);
		}

//		if( ! vrtxIsBndVrt )
//		{
		if( ! establishNewVertizesStasiBased(oldVrt) )
		{
			UG_LOG("Vertex Erzeugung schief gegangen " << std::endl);
			return false;
		}

//		}


		

	}

	// for debugging
//	return false;

	return true;
}

////////////////////////////////////////////////////////////////////

#if 0
bool ArteExpandFracs3D::establishNewVertizesStasiBased( Vertex * const & oldVrt)
{
	// anfangs nur für innere Vertizes mit einer fracture

	VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[oldVrt];

	if( vecSegVolElmInf.size() < 2 )
	{
		UG_LOG("nur ein Segment, aber will frische Vertizes?" << std::endl);
//		UG_THROW("nur ein Segment, aber will frische Vertizes?" << std::endl);
//		return false;
		return true;
	}

	for( SegmentVolElmInfo const & svei : vecSegVolElmInf )
	{
		// count the number of the fracture subdomains surrounding the segment
		// in case of boundary vertex, also count the number of boundary subdomains

		std::vector<IndexType> sudosInSegment;

		if( ! extracFractSudosOfSegment( svei, sudosInSegment ) )
		{
			UG_LOG("kann sudos nicht extrahieren " << std::endl);
			UG_THROW("kann sudos nicht extrahieren " << std::endl);
			return false;
		}

		IndexType sudoNumInSeg = sudosInSegment.size();

		// check if is boundary vertex

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];

		bool vrtxIsBndVrt = vrtxFracPrps.getIsBndFracVertex();

		UG_LOG("is bndry " << vrtxIsBndVrt << std::endl);

		if( ! vrtxIsBndVrt )
		{
			// count number of fracture subdomains in the Segment, should have been done before......

			// TODO FIXME ersetze das altmodische VrtxFracProptsStatus durch ein Zählen
			// der subdomains, die pro Segment wirklich vorkommen, nachdem die
			// auslaufenden fracture faces den generellen Faces zugeschlagen worden sind
			// diese Zählweise ist die erste Folgeaufgabe, damit der komische vrtxFractureProperties Status
			// weg geworfen werden kann, der ist nämlich nutzlos inzwischen, da er auch auslaufende
			// fracture faces zählt, was Unsinn ist.......
			// später folgt auch die Verallgemeinerung auf boundary vertizes, dann muss deren Wahl
			// vielleicht wieder dazu genommen werden.....

			//Vorbereitung Ersetzung VertexFractureProperties Status
			// durch Zählen der Fracture Subdomains von jedem Segment
			// die Zählerei kann dann gemacht werden, wenn die doppelten offenen fracture faces
			// in die allgemeinen faces verschoben werden, dann geht das in einem Abwasch,
			// gehe dazu in die entsprechende Funktion

//			if( vrtxFracPrps.getVrtxFracStatus() == VrtxFracProptsStatus::oneFracSuDoAtt )
			if( sudoNumInSeg == 1 )
			{
				// standard case, one fracture
				if( ! expandWithinTheSegment<1,false>( oldVrt, svei ) )
				{
					UG_LOG("Expandierung einfachster Fall schief gegangen " << std::endl);
					return false;
				}
			}
			else // unterscheiden zwei und drei vermutlich..... aber wichtiger Segmentzahl..... und dann kommt es darauf an, wieviele subdoms im einzelnen Segment sind
			{

			}

		}
		else
		{
			// boundary faces counted as fracture faces, but no expansion
			// zuerst aber die inneren Schnittvertizes, weil bei denen die Nicht-Null Expansion
			// gelernt werden kann, dann anzuwenden auf Null-Expansion senkrecht zu den boundaries......
		}

	}

	return true;
}
#endif

///////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::establishSegmentLimitingSidesInfo()
{

	UG_LOG("ESTABLISH SGEMENT LIM SIDE INFO" << std::endl );

	for( VertexIterator iterV = m_sel.begin<Vertex>(); iterV != m_sel.end<Vertex>(); ++ iterV )
	{
		Vertex * oldVrt = *iterV;

		UG_LOG("establish segm lim sides for " << m_aaPos[oldVrt] << std::endl);

		VecSegmentLimitingSides & vecSegmLimSid = m_vrtxAttAccsVecSegmLimSid[oldVrt];

		UG_LOG("Segment limiting sides got " << std::endl);

		auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];

		bool vrtxIsBndVrt = vrtxFracPrps.getIsBndFracVertex();

		UG_LOG("is bndry " << vrtxIsBndVrt << std::endl);

		VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[oldVrt];

		UG_LOG("CHECK SIZE" << std::endl);

		if( vecSegVolElmInf.size() < 2 )
		{
			UG_LOG("nur ein Segment, aber will frische Vertizes?" << std::endl);
	//		UG_THROW("nur ein Segment, aber will frische Vertizes?" << std::endl);
	//		return false;
//			return true;
			continue;
		}

		UG_LOG("CONTINUED" << std::endl);

		for( SegmentVolElmInfo const & segVolsElInf : vecSegVolElmInf )
		{
			// count the number of the fracture subdomains surrounding the segment
			// in case of boundary vertex, also count the number of boundary subdomains

	//		std::vector<IndexType> sudosInSegment;

			// hier die neue Klasse und ihr averaing einführen, oder im Hauptloop der neuen Elemente.....
			// SegmentSides<....> in .h file für die richtigen template parameter Kurznamen geben und hier
			// Objekt erstellen und je nach ob mit oder ohne Bndry entsprechend aufladen und averagen.....
			// und zwar für jedes Segment einzeln, und abhängig von boundary oder nicht die boundary Geschichten
			// dazu oder auch nicht...... also einen VecSegmentSides erstellen auch, das einzelne Objekt
			// weiss dann, ob es eine Boundary ist..... vielleicht noch den Vektor dazu als übergebener Parameter
			// damit man den da drin weiter reichen kann?
			// wenn jedes Objekt des Vektors von SegmentSIdes sowohl seinen Vertex kennt als auch weiss,
			// ob es boundary oder nicht ist, kann danach darüber geloopt werden, ohne nochmal
			// aussen den Vertex mit geben zu müssen, oder die Info, ob bndry oder nicht.....
			// danach gibts dann Funktionen, die alleine ein Objekt von der Sorte SegmentSides schlucken brauchen
			// und hier nur ein Loop über den ganzen Vektor davon, wo für jedes Element dann die Fkt aufgerufen wird....

			SegmentLimitingSides segLimSids( oldVrt, vrtxIsBndVrt );

	//		std::vector<Volume*> vecVolsOfSegment;

			IndexType boundarySites = 0;

			for( AttachedVolumeElemInfo const & volElmInf : segVolsElInf )
			{
				if( ! segLimSids.schluckVecAttFractElm( volElmInf.getVecFractManifElem() ) )
				{
					UG_LOG("schlucken schief gegangen " << std::endl);
					UG_THROW("schlucken schief gegangen " << std::endl);
					return false;
				}


				if( volElmInf.hasUnclosedFracture() )
				{
					if( ! segLimSids.schluckVecAttUnclosedFractElm( volElmInf.getVecUnclosedFractManifElem() ) )
					{
						UG_LOG("Schlucken von unclosed schief gegangen " << std::endl);
//						UG_THROW("Schlucken von unclosed schief gegangen " << std::endl);
//						return false;
					}
					else
					{
						UG_LOG("ungeschlossene fracture geschluckt " << std::endl);
					}

					if( ! segLimSids.hasUnclosedFaces() )
					{
						UG_LOG("keine ungeschlossenen Gesichter " << std::endl);
						UG_THROW("keine ungeschlossenen Gesichter " << std::endl);
					}
				}

				if( vrtxIsBndVrt )
				{
					auto vecBndryManifelm = volElmInf.getVecBndryManifElem();

					IndexType sizeVecBndryManifElm = volElmInf.getVecBndryManifElem().size();

					boundarySites += sizeVecBndryManifElm;

	//				if( ! segLimSids.schluckVecAttBndryElm( volElmInf.getVecBndryManifElem() ) )
					if( ! segLimSids.schluckVecAttBndryElm( vecBndryManifelm ) )
					{
						UG_LOG("schlucken B schief gegangen " << std::endl);
						UG_THROW("schlucken B schief gegangen " << std::endl);
						return false;
					}

					// TODO FIXME es muss abgefangen werden, wenn bei einem boundary vertex gar keine boundary Seiten da sind
	//				if( (volElmInf.getVecBndryManifElem()).size() == 0  )
					if( sizeVecBndryManifElm == 0  )
						UG_LOG("Grenze verloren gegangen " << m_aaPos[oldVrt] << std::endl);
				}

				Volume * vol2Add = volElmInf.getFulldimElem();
				segLimSids.schluckFulldimElem(vol2Add);
			}

			if( vrtxIsBndVrt && boundarySites == 0 )
			{
				UG_LOG("No boundary sites at " << m_aaPos[oldVrt] << std::endl);
				UG_THROW("No boundary sites at " << m_aaPos[oldVrt] << std::endl);
			}

			if( ! segLimSids.averageAll() )
			{
				UG_LOG("keine Mittelung " << std::endl);
				UG_THROW("keine Mittelung " << std::endl);
				return false;
			}

			vecSegmLimSid.push_back(segLimSids);
		}


	}

	UG_LOG("END ESTABLISH SGEMENT LIM SIDE INFO" << std::endl );


	return true;
}

///////////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::establishNewVertizesStasiBased( Vertex * const & oldVrt)
{

	auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];

	bool vrtxIsBndVrt = vrtxFracPrps.getIsBndFracVertex();

	UG_LOG("is bndry " << vrtxIsBndVrt << std::endl);


//	VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[oldVrt];
//
//	if( vecSegVolElmInf.size() < 2 )
//	{
//		UG_LOG("nur ein Segment, aber will frische Vertizes?" << std::endl);
////		UG_THROW("nur ein Segment, aber will frische Vertizes?" << std::endl);
////		return false;
//		return true;
//	}

	VecSegmentLimitingSides & vecSegmLimSid = m_vrtxAttAccsVecSegmLimSid[oldVrt];

	if(  vecSegmLimSid.size() == 0 )
	{
		UG_LOG("keine verschiedenen Segmente" << std::endl);
		return true;
	}

	if( vecSegmLimSid.size() == 1 )
	{
		UG_LOG("sonderbarer Fall von nur einem Segment " << std::endl );
		UG_THROW("sonderbarer Fall von nur einem Segment " << std::endl );
		return false;
	}

	for( SegmentLimitingSides const & segLimSids : vecSegmLimSid )
	{

		if( segLimSids.hasUnclosedFaces() )
		{
			UG_LOG("beim ausdehnen ungeschlossene Gesichter " << std::endl);
		}

		if( ! expandWithinTheSegment(segLimSids) )
		{
			UG_LOG("schief gegangen Vertex Erzeugung " << std::endl);
			UG_THROW("schief gegangen Vertex Erzeugung " << std::endl);
			return false;
		}
	}

	UG_LOG("Vertex creation hat funktioniert " << std::endl);
	//	UG_THROW("Vertex creation failed " << std::endl);

	return true;

}

/////////////////////////////////////////////////////////////////

#if 0

bool ArteExpandFracs3D::establishNewVertizesStasiBased( Vertex * const & oldVrt)
{
	// anfangs nur für innere Vertizes mit einer fracture

	// testweise, später verallgemeinert mit den Boundary faces

	auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];

	bool vrtxIsBndVrt = vrtxFracPrps.getIsBndFracVertex();

	UG_LOG("is bndry " << vrtxIsBndVrt << std::endl);

//	if( vrtxIsBndVrt )
//	{
//		UG_LOG("boundary noch zu lösen, bisher nix machen" << std::endl);
//		return true;
//	}

	VecSegmentVolElmInfo & vecSegVolElmInf = m_accsAttVecSegVolElmInfo[oldVrt];

	if( vecSegVolElmInf.size() < 2 )
	{
		UG_LOG("nur ein Segment, aber will frische Vertizes?" << std::endl);
//		UG_THROW("nur ein Segment, aber will frische Vertizes?" << std::endl);
//		return false;
		return true;
	}


	for( SegmentVolElmInfo const & segVolsElInf : vecSegVolElmInf )
	{
		// count the number of the fracture subdomains surrounding the segment
		// in case of boundary vertex, also count the number of boundary subdomains

//		std::vector<IndexType> sudosInSegment;

		// hier die neue Klasse und ihr averaing einführen, oder im Hauptloop der neuen Elemente.....
		// SegmentSides<....> in .h file für die richtigen template parameter Kurznamen geben und hier
		// Objekt erstellen und je nach ob mit oder ohne Bndry entsprechend aufladen und averagen.....
		// und zwar für jedes Segment einzeln, und abhängig von boundary oder nicht die boundary Geschichten
		// dazu oder auch nicht...... also einen VecSegmentSides erstellen auch, das einzelne Objekt
		// weiss dann, ob es eine Boundary ist..... vielleicht noch den Vektor dazu als übergebener Parameter
		// damit man den da drin weiter reichen kann?
		// wenn jedes Objekt des Vektors von SegmentSIdes sowohl seinen Vertex kennt als auch weiss,
		// ob es boundary oder nicht ist, kann danach darüber geloopt werden, ohne nochmal
		// aussen den Vertex mit geben zu müssen, oder die Info, ob bndry oder nicht.....
		// danach gibts dann Funktionen, die alleine ein Objekt von der Sorte SegmentSides schlucken brauchen
		// und hier nur ein Loop über den ganzen Vektor davon, wo für jedes Element dann die Fkt aufgerufen wird....

		SegmentLimitingSides segLimSids( oldVrt, vrtxIsBndVrt );

//		std::vector<Volume*> vecVolsOfSegment;

		IndexType boundarySites = 0;

		for( AttachedVolumeElemInfo const & volElmInf : segVolsElInf )
		{
			if( ! segLimSids.schluckVecAttFractElm( volElmInf.getVecFractManifElem() ) )
			{
				UG_LOG("schlucken schief gegangen " << std::endl);
				UG_THROW("schlucken schief gegangen " << std::endl);
				return false;
			}

			if( vrtxIsBndVrt )
			{
				auto vecBndryManifelm = volElmInf.getVecBndryManifElem();

				IndexType sizeVecBndryManifElm = volElmInf.getVecBndryManifElem().size();

				boundarySites += sizeVecBndryManifElm;

//				if( ! segLimSids.schluckVecAttBndryElm( volElmInf.getVecBndryManifElem() ) )
				if( ! segLimSids.schluckVecAttBndryElm( vecBndryManifelm ) )
				{
					UG_LOG("schlucken B schief gegangen " << std::endl);
					UG_THROW("schlucken B schief gegangen " << std::endl);
					return false;
				}

				// TODO FIXME es muss abgefangen werden, wenn bei einem boundary vertex gar keine boundary Seiten da sind
//				if( (volElmInf.getVecBndryManifElem()).size() == 0  )
				if( sizeVecBndryManifElm == 0  )
					UG_LOG("Grenze verloren gegangen " << m_aaPos[oldVrt] << std::endl);
			}

			Volume * vol2Add = volElmInf.getFulldimElem();
			segLimSids.schluckFulldimElem(vol2Add);
		}

		if( vrtxIsBndVrt && boundarySites == 0 )
		{
			UG_LOG("No boundary sites at " << m_aaPos[oldVrt] << std::endl);
			UG_THROW("No boundary sites at " << m_aaPos[oldVrt] << std::endl);
		}

		if( ! segLimSids.averageAll() )
		{
			UG_LOG("keine Mittelung " << std::endl);
			UG_THROW("keine Mittelung " << std::endl);
			return false;
		}

//		if( vrtxIsBndVrt )
//		{
//			if( segLimSids. )
//		}

//		for( SegmentLimitSidesPairSudoNorml segLimSiPSN : vecSegmLimSidPrSudoNorml )
//		{
//
//		}

//		expandWithinTheSegment<SegmentVrtxFracStatus::oneFracSuDoAtt>(segLimSids

		if( ! expandWithinTheSegment(segLimSids) )
		{
			UG_LOG("schief gegangen Vertex Erzeugung " << std::endl);
			UG_THROW("schief gegangen Vertex Erzeugung " << std::endl);
			return false;
		}

//		SegmentVrtxFracStatus segStatFract = segLimSids.spuckCrossingTyp();
//
//		switch( segStatFract )
//		{
//			case SegmentVrtxFracStatus::noFracSuDoAtt :
//			{
//				UG_LOG("is not a fracture, but a segment?" << std::endl);
//				UG_THROW("is not a fracture, but a segment?" << std::endl);
//				return false;
//			}
//			case SegmentVrtxFracStatus::oneFracSuDoAtt :
//			{
//				if( ! expandWithinTheSegment<SegmentVrtxFracStatus::oneFracSuDoAtt,false,false>(segLimSids, vecVolsOfSegment) )
//				{
//					UG_LOG("Expandierung 1 schief gegangen " << std::endl);
//					UG_THROW("Expandierung 1 schief gegangen " << std::endl);
//					return false;
//				}
//				break;
//			}
//			case SegmentVrtxFracStatus::twoFracSuDoAtt :
//			{
//				break; // TODO FIXME implementieren
//			}
//			case SegmentVrtxFracStatus::threeFracSuDoAtt :
//			{
//				if( ! expandWithinTheSegment<SegmentVrtxFracStatus::threeFracSuDoAtt,false,false>(segLimSids, vecVolsOfSegment) )
//				{
//					UG_LOG("Expandierung 3 schief gegangen " << std::endl);
//					UG_THROW("Expandierung 3 schief gegangen " << std::endl);
//					return false;
//				}
//				break;
//			}
//			default :
//			{
//				UG_LOG("strange fracture crossing" << std::endl);
//				UG_THROW("strange fracture crossing?" << std::endl);
//				return false;
//			}
//		}
	}

	UG_LOG("Vertex creation hat funktioniert " << std::endl);
	//	UG_THROW("Vertex creation failed " << std::endl);

	return true;
}

#endif

//		if( ! extracFractSudosOfSegment( svei, sudosInSegment ) )
//		{
//			UG_LOG("kann sudos nicht extrahieren " << std::endl);
//			UG_THROW("kann sudos nicht extrahieren " << std::endl);
//			return false;
//		}
//
//		IndexType sudoNumInSeg = sudosInSegment.size();
//
//		// check if is boundary vertex
//
//		auto & vrtxFracPrps = m_aaMarkVrtVFP[ oldVrt ];
//
//
//		if( ! vrtxIsBndVrt )
//		{
//			// count number of fracture subdomains in the Segment, should have been done before......
//
//			// TODO FIXME ersetze das altmodische VrtxFracProptsStatus durch ein Zählen
//			// der subdomains, die pro Segment wirklich vorkommen, nachdem die
//			// auslaufenden fracture faces den generellen Faces zugeschlagen worden sind
//			// diese Zählweise ist die erste Folgeaufgabe, damit der komische vrtxFractureProperties Status
//			// weg geworfen werden kann, der ist nämlich nutzlos inzwischen, da er auch auslaufende
//			// fracture faces zählt, was Unsinn ist.......
//			// später folgt auch die Verallgemeinerung auf boundary vertizes, dann muss deren Wahl
//			// vielleicht wieder dazu genommen werden.....
//
//			//Vorbereitung Ersetzung VertexFractureProperties Status
//			// durch Zählen der Fracture Subdomains von jedem Segment
//			// die Zählerei kann dann gemacht werden, wenn die doppelten offenen fracture faces
//			// in die allgemeinen faces verschoben werden, dann geht das in einem Abwasch,
//			// gehe dazu in die entsprechende Funktion
//
////			if( vrtxFracPrps.getVrtxFracStatus() == VrtxFracProptsStatus::oneFracSuDoAtt )
//			if( sudoNumInSeg == 1 )
//			{
//				// standard case, one fracture
//				if( ! expandWithinTheSegment<1,false>( oldVrt, svei ) )
//				{
//					UG_LOG("Expandierung einfachster Fall schief gegangen " << std::endl);
//					return false;
//				}
//			}
//			else // unterscheiden zwei und drei vermutlich..... aber wichtiger Segmentzahl..... und dann kommt es darauf an, wieviele subdoms im einzelnen Segment sind
//			{
//
//			}
//
//		}
//		else
//		{
//			// boundary faces counted as fracture faces, but no expansion
//			// zuerst aber die inneren Schnittvertizes, weil bei denen die Nicht-Null Expansion
//			// gelernt werden kann, dann anzuwenden auf Null-Expansion senkrecht zu den boundaries......
//		}
//
//	}

//	UG_LOG("Vertex creation failed " << std::endl);
//	UG_THROW("Vertex creation failed " << std::endl);
//
//
//	return false;
//


////////////////////////////////////////////////////////////////////

#if 0
bool ArteExpandFracs3D::extracFractSudosOfSegment( SegmentVolElmInfo const & segmVolElmInfo, std::vector<ArteExpandFracs3D::IndexType> & sudosInSegment )
{


	VecAttachedFractFaceEdgeSudo vecAttFractFaces;

	for( AttachedVolumeElemInfo const & avei : segmVolElmInfo )
	{
		VecAttachedFractFaceEdgeSudo const & vecAttFractVol = avei.getVecFractManifElem();

		for( AttachedFractFaceEdgeSudo const & affe : vecAttFractVol )
		{
			vecAttFractFaces.push_back(affe);
		}
	}

	IndexType numbContrFracFaces = vecAttFractFaces.size();

	if( numbContrFracFaces < 1 )
	{
		UG_LOG("Kein Affe da " << std::endl);
		UG_THROW("Kein Affe da " << std::endl);
		return false;
	}

	IndexType sudoBase = vecAttFractFaces[0].getSudo();

	sudosInSegment.push_back(sudoBase);

	// add sudos different from the base one
	for( AttachedFractFaceEdgeSudo const & affe : vecAttFractFaces )
	{
		IndexType sudoNeeded = affe.getSudo();

		bool sudoIsKnown = false;

		for( IndexType sudoInList : sudosInSegment )
		{
			if( sudoNeeded == sudoInList )
			{
				sudoIsKnown = true;
			}
		}

		if( ! sudoIsKnown )
		{
			sudosInSegment.push_back(sudoNeeded);
		}

	}

	return true;
}
#endif


////////////////////////////////////////////////////////////////////

// for only one surrounding subdom around the segment, for example only one fracture, or T End like ending side
// TODO FIXME alles in eine einzige Funktion, die verschiedene Unterfunktionen aufruft für verschiedene Zahlen
// von Seiten innen und aussen!!!!
//template<>
//bool ArteExpandFracs3D::expandWithinTheSegment<ArteExpandFracs3D::SegmentVrtxFracStatus::oneFracSuDoAtt>( SegmentLimitingSides const & segmLimSides )
bool ArteExpandFracs3D::expandWithinTheSegment( ArteExpandFracs3D::SegmentLimitingSides const & segmLimSides )
{
	// should not be called for boundary vertices

	bool isBndry = segmLimSides.isBoundary();


//	if( isBndry )
//		return true;

//	if( isBndry )
//	{
//		UG_LOG("boundary noch zu behandeln " << std::endl);
//		UG_THROW("boundary noch zu behandeln " << std::endl);
//		return false;
//	}

//	VecSegmentLimitSidesPairSudoNorml vecSegmLimSidPrSudoNrml;
	VecPlaneDescriptor vecPlaneFracDescr;

	//	if( ! segmLimSides.spuckFractSudoNormls( vecSegmLimSidPrSudoNrml ) )
	if( ! segmLimSides.spuckFractManifDescr( vecPlaneFracDescr, m_aaPos ) )
	{
		UG_LOG("Spucken schief gegangen  " << segmLimSides.spuckCrossingTyp() << std::endl);
		UG_THROW("Spucken schief gegangen  " << segmLimSides.spuckCrossingTyp() << std::endl);
		return false;
	}

	// falls Boundary, hier noch spuckBndryManifDescr aufrufen

	Vertex * oldVrt = segmLimSides.spuckVertex();
	vector3 posOldVrt = m_aaPos[oldVrt];

	VecPlaneDescriptor vecShiftedPlaneDescript;



//	for( SegmentLimitSidesPairSudoNorml const & segLimSidPrSN : vecSegmLimSidPrSudoNrml )
	for( PlaneDescriptor & planeDescr : vecPlaneFracDescr )
	{
		UG_LOG("GOT MANIF TYP " << planeDescr.spuckManifTyp() << std::endl );
//		IndexType const & sudoSide = segLimSidPrSN.first;
//		vector3 const & normlAvrg = segLimSidPrSN.second;
//
//		// normal computed standard mässig directs of of the volume, we need that one into the volume
//		vector3 normalOutsideVol;

		// hier testen, ob die Subdomain von der Liste der fracture subdomains ist,
		// damit für den Fall von sich zwei schneidenden Ebenen für die dritte,
		// die als senkrecht zu den beiden anderen gesetzt werden soll, mit Verschiebung null,
		// analog Boundary sides, die künstliche Weite null erhalten kann hier

		if( planeDescr.spuckManifTyp() != PlaneDescriptorType::isFracture )
			UG_THROW("muss fracture sein " << std::endl);
//
		int sudoSide = planeDescr.spuckSudo();
		number width = m_fracInfosBySubset[sudoSide].width;;

		// ensure that the width is nonzero only for real fractures, not for pseudo vectors or such stuff
//		if( ! isBndry )
//		{
//			if( planeDescr.spuckManifTyp() == PlaneDescriptorType::isFracture )
//			{
//				width = m_fracInfosBySubset[sudoSide].width;
//			}
//			else
//			{
//				UG_LOG("Manif Typ is " << planeDescr.spuckManifTyp() << std::endl );
//				UG_LOG("Fract wäre " << PlaneDescriptorType::isFracture << std::endl );
//				UG_LOG("Bndry wäre " << PlaneDescriptorType:: << std::endl );
//				UG_LOG("Artif wäre " << PlaneDescriptorType::isArtificial << std::endl );
//
//				UG_THROW("keine Boundary, aber will Boundary Manif" << std::endl);
//			}
//		}
//		else
//		{
//			UG_THROW("hier nur keine boundary" << std::endl);
//		}


//		else
//		{
//	TODO FIXME		bei boundary auch noch für die entsprechenden boundaries
//		}

		number shiftScal = width / 2.;
		planeDescr.schluckScaleShiftNormal( - shiftScal );

//		vector3 shiftVec4Plane;
//
//		// as normal outside vol, but we need to go inside the volumes / segment side
//		VecScale( shiftVec4Plane, normalOutsideVol, - shiftScal );

//		vecShiftVec4Plane.push_back( shiftVec4Plane );
//		vecSudosSides.push_back( sudoSide );
//		vecNormalsAveraged.push_back( normalIntoVol );

//		PlaneDescriptor planeDescr( normalOutsideVol, posOldVrt, shiftScal );

		PlaneDescriptor shiftedPlaneDescr;
		planeDescr.spuckPlaneShifted( shiftedPlaneDescr );

		//		planeDescr.spuckPlaneShiftedAlong( shiftVec4Plane, shiftedPlaneDescr );

		vecShiftedPlaneDescript.push_back( shiftedPlaneDescr );
	}

	// will get the sudo of the shifted vertex
	IndexType sudoExample = (vecPlaneFracDescr[0]).spuckSudo();
	vector3 posNewVrt; // to be determined depending on the segment properties

	bool hasUnclosedFaces = segmLimSides.hasUnclosedFaces();

	if( hasUnclosedFaces )
	{
		UG_LOG("hat ungeschlossene Seiten " << std::endl);
	}

	if( ! isBndry && vecShiftedPlaneDescript.size() == 1 )
	{
		// TODO FIXME die unclosed crossing faces behandeln!!!

		// XXXXX muss wieder eingeschaltet werden hier
		// TODO FIXME use the hasUnclosedFaces property, but for visual debugging, better set as false
		bool distinguishUnclosedFaces = segmLimSides.hasUnclosedFaces();
//		constexpr bool distinguishUnclosedFaces = false; // segmLimSides.hasUnclosedFaces();
		UG_LOG("Please use switch unclosed faces in final version, please remove debugging set false" << std::endl);

		if( ! distinguishUnclosedFaces ) // standard case, just shift perpendicular to the plane, no SLE to solve
		{
			//		// one single inner fracture, Testfall Anfang
			//		computeShiftVector( vecShiftedPlaneDescript[0] );
			posNewVrt = (vecPlaneFracDescr[0]).spuckShiftedBaseVect();

		}
		else
		{
			// NOTE: TODO FIXME funktioniert vermutlich auch, wenn der Vertex ein Boundary Vertex ist
			// ob der Spezialfall aber auch praktisch vorkommt, fragwürdig
			// Test der Methode erstmal innen drin, später ggf. Ausweitung auf Grenzvertizes
			// dann würde das hier raus wandern aus dem Frage ob Boundary oder nicht, isBndry wäre egal
			// relevant wäre dass size des vecShiftPlaneDescr 1 ist und withUnclosedFaces true


			Edge * shiftDirectionEdg = nullptr;

			if( ! segmLimSides.spuckLowdimElmShiftDirectionIfUnclosedFractPresent(shiftDirectionEdg) )
			{
				UG_LOG("no shift direction " << std::endl);
				UG_THROW("no shift direction " << std::endl);

			}

			if( shiftDirectionEdg == nullptr )
			{
				UG_LOG("Null Shift cross end?" << std::endl);
				UG_THROW("Null Shift cross end?" << std::endl);
			}

			// jetzt sollte also die shift direction bekannt sein, muss normalisiert werden,
			// vielleicht kann man das im PlaneDescriptor drin auch machen

			PlaneDescriptor const & shiftedPlane = vecShiftedPlaneDescript[0];


			if( ! computeCrossPointOfPlaneWithLine( shiftedPlane, shiftDirectionEdg, oldVrt, posNewVrt ) )
			{
				UG_LOG("Not possible to compute crossing point plane with line" << std::endl);
				UG_THROW("Not possible to compute crossing point plane with line" << std::endl);
			}
		}
	}
	else
	{
		if( ! isBndry && vecShiftedPlaneDescript.size() > 1 )
		{

			if( vecShiftedPlaneDescript.size() == 2 )
			{
				// we need to add an artificial plane which gets extended with scale 0

				vector3 artificialNormal;
				VecCross( artificialNormal, vecShiftedPlaneDescript[0].spuckNormalVector(), vecShiftedPlaneDescript[1].spuckNormalVector() );

				PlaneDescriptor artificialPlane( artificialNormal, posOldVrt );

				vecShiftedPlaneDescript.push_back(artificialPlane);
			}

			if( vecShiftedPlaneDescript.size() > 3  )
				UG_THROW("too much fractures" << std::endl);

			// now we should have three planes to cross

		}
		else if( isBndry )
		{
			VecPlaneDescriptor vecPlaneBndryDescr;

			//	if( ! segmLimSides.spuckFractSudoNormls( vecSegmLimSidPrSudoNrml ) )
			if( ! segmLimSides.spuckBndryManifDescr( vecPlaneBndryDescr, m_aaPos ) )
			{
				UG_LOG("Spucken schief gegangen bndry " << segmLimSides.spuckCrossingTyp() << std::endl);
				UG_THROW("Spucken schief gegangen bndry " << segmLimSides.spuckCrossingTyp() << std::endl);
				return false;
			}

			if( vecPlaneBndryDescr.size() < 1 )
			{
				UG_LOG("at point " << m_aaPos[oldVrt] << std::endl);

				UG_LOG("BOUNRARY PPPPP vertex BOUNDARY size problem " << m_aaPos[oldVrt] << std::endl);

				UG_LOG("Boundaries" << std::endl);
				for( PlaneDescriptor pd : vecPlaneBndryDescr  )
				{
					UG_LOG("Subdom " << pd.spuckSudo() << std::endl);
					UG_LOG("Base " << pd.spuckBaseVector() << std::endl);
					UG_LOG("Normal " << pd.spuckNormalVector() <<std::endl);
				}

				UG_LOG("Shifted" << std::endl);
				for( PlaneDescriptor pd : vecShiftedPlaneDescript  )
				{
					UG_LOG("Subdom S " << pd.spuckSudo() << std::endl);
					UG_LOG("Base S " << pd.spuckBaseVector() << std::endl);
					UG_LOG("Normal S " << pd.spuckNormalVector() <<std::endl);
				}

				UG_LOG("Fracs " << std::endl);
				for( PlaneDescriptor pd : vecPlaneFracDescr  )
				{
					UG_LOG("Subdom F " << pd.spuckSudo() << std::endl);
					UG_LOG("Base F " << pd.spuckBaseVector() << std::endl);
					UG_LOG("Normal F " << pd.spuckNormalVector() <<std::endl);
				}


				// PROBLEM PPPPPPPPPPPPPPPPP
				UG_THROW("NO boundary sudos " << std::endl);
			}

			if( vecPlaneBndryDescr.size() > 2 )
			{
				UG_LOG("dudos" << std::endl);
				for( PlaneDescriptor pd : vecPlaneBndryDescr )
				{
					int sudo = pd.spuckSudo();
					UG_LOG("sudos " << sudo << std::endl);

				}
				UG_LOG("at point " << m_aaPos[oldVrt] << std::endl);
				UG_THROW("too much boundary sudos " << vecPlaneBndryDescr.size() << std::endl);
			}


			if( vecShiftedPlaneDescript.size() + vecPlaneBndryDescr.size() > 3 )
				UG_THROW("too much crossing stuff at boundary"<<std::endl);

			for( PlaneDescriptor const & pbd : vecPlaneBndryDescr )
			{
				vecShiftedPlaneDescript.push_back( pbd );
			}

			if( vecPlaneBndryDescr.size() == 2 )
			{
				// PROBLEM PPPPPPPPPPPPPPPP
				if( vecShiftedPlaneDescript.size() != 3 ||  vecPlaneFracDescr.size() != 1 )
				{
					UG_LOG("BOUNRARY PPPPP vertex TWO problem " << m_aaPos[oldVrt] << std::endl);

					UG_LOG("Boundaries" << std::endl);
					for( PlaneDescriptor pd : vecPlaneBndryDescr  )
					{
						UG_LOG("Subdom " << pd.spuckSudo() << std::endl);
						UG_LOG("Base " << pd.spuckBaseVector() << std::endl);
						UG_LOG("Normal " << pd.spuckNormalVector() <<std::endl);
					}

					UG_LOG("Shifted" << std::endl);
					for( PlaneDescriptor pd : vecShiftedPlaneDescript  )
					{
						UG_LOG("Subdom S " << pd.spuckSudo() << std::endl);
						UG_LOG("Base S " << pd.spuckBaseVector() << std::endl);
						UG_LOG("Normal S " << pd.spuckNormalVector() <<std::endl);
					}

					UG_LOG("Fracs " << std::endl);
					for( PlaneDescriptor pd : vecPlaneFracDescr  )
					{
						UG_LOG("Subdom F " << pd.spuckSudo() << std::endl);
						UG_LOG("Base F " << pd.spuckBaseVector() << std::endl);
						UG_LOG("Normal F " << pd.spuckNormalVector() <<std::endl);
					}

					UG_THROW("noch nicht genug Grenzen " << std::endl);

				}
			}
			else if( vecPlaneBndryDescr.size() == 1 )
			{
				if( vecShiftedPlaneDescript.size() < 3 )
					UG_LOG("noch nicht genug Grenzen " << std::endl);

				if( vecShiftedPlaneDescript.size() > 3 )
					UG_LOG("too much Grenzen " << std::endl);

				if( vecPlaneFracDescr.size() != 1 && vecPlaneFracDescr.size() != 2 )
					UG_THROW("Nicht passend Fract plus Bndry" << std::endl);

//				if( vecPlaneFracDescr.size() == 1 && vecPlaneFracDescr.size() == 2 )
//					if( vecShiftedPlaneDescript.size() != 3 )
//						UG_THROW"SO viele unsinnige Kombinationen " << std::endl );

				if( vecPlaneFracDescr.size() == 2 )
				{
					if( vecShiftedPlaneDescript.size() != 3 )
						UG_THROW("da passt nicht zusammen was 3 sein sollte" << std::endl);

					// sonst nix zu tun
				}
				else if( vecPlaneFracDescr.size() == 1 )
				{
					if( vecShiftedPlaneDescript.size() != 2 )
						UG_LOG("1+2 nicht 3" << std::endl);

					// add an artificial perpendicular plane with move vector zero

					vector3 artificialNormal;
					VecCross( artificialNormal, vecShiftedPlaneDescript[0].spuckNormalVector(), vecShiftedPlaneDescript[1].spuckNormalVector() );

					PlaneDescriptor artificialPlane( artificialNormal, posOldVrt );

					vecShiftedPlaneDescript.push_back(artificialPlane);

				}



			}

		}

		// vector shifted plane descriptor kriegt was dazu, wenn es keine drei sind ohne boundary

		if( vecShiftedPlaneDescript.size() != 3 )
			UG_THROW("wie gross soll es denn sein" << std::endl);

		computeCrossingPointOf3Planes( vecShiftedPlaneDescript, posNewVrt );

	}



//	if( vecSegmLimSidPrSudoNrml.size() != 1 )
//	{
//		UG_LOG("only one fracture, but not one normal and sudo?"  <<  std::endl);
//		UG_THROW("only one fracture, but not one normal and sudo?" << std::endl);
//		return false;
//	}

	// FALL eins extra Funktion, unterscheiden nach Länge

//	if(  )


//	SegmentLimitSidesPairSudoNorml & sudoAndNormal = vecSegmLimSidPrSudoNrml[0];
//
//	IndexType sudoBase = sudoAndNormal.first;
//	vector3 normalsAveraged = sudoAndNormal.second;
//
//	Vertex * oldVrt = segmLimSides.spuckVertex();
//
//	number width = m_fracInfosBySubset[sudoBase].width;
//
//	number scal = width / 2.;
//
//	NormalVectorFacIntoVol scaledNormal;
//
//	VecScale( scaledNormal, normalsAveraged, - scal );
//
//	vector3 posOldVrt = m_aaPos[oldVrt];
//
//	UG_LOG("NORMAL OLD VRTX " << posOldVrt << " -> " << normalsAveraged << std::endl );


//	VecAdd( posNewVrt, posOldVrt, scaledNormal );

	Vertex * newShiftVrtx = *m_grid.create<RegularVertex>();

	if( newShiftVrtx == nullptr )
	{
		UG_LOG("Nullen erzeugt" << std::endl);
		UG_THROW("Nullen erzeugt" << std::endl);
		return false;
	}

	m_aaPos[newShiftVrtx] = posNewVrt;

	UG_LOG("Created new vertex at " << m_aaPos[newShiftVrtx] << std::endl );

	m_sh.assign_subset(newShiftVrtx, sudoExample);
//	m_sh.assign_subset(newShiftVrtx, m_sh.num_subsets());

	// wenn der neue Vertex mit einer ungeschlossenen endenden Kluft zu tun hat
	// muss am Ende vom Kernvertex, old vertex, die edge, die nicht mit der
	// kreuzenden nicht endenden KLuft geht, gelöscht werden
	m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft[newShiftVrtx] = hasUnclosedFaces;

	// TODO FIXME das auch verwenden

	if( hasUnclosedFaces )
	{
		// TEST AUF ANZAHL der kreuzenden Klüfte, zusätzlich, ob es eine ist nur
		// zur Sicherheit
		if( isBndry || vecShiftedPlaneDescript.size() != 1 )
		{
			UG_LOG("Shift vertex but boundary or more than one shift direction " << std::endl);
			UG_THROW("Shift vertex but boundary or more than one shift direction " << std::endl);
		}

		if( ! addElem( m_vrtxArisesFromExpandedEndingCrossingCleft, newShiftVrtx ) )
		{
			UG_LOG("Shifted ecf vertex schon bekannt " << std::endl);
			UG_THROW("Shifted ecf vertex schon bekannt " << std::endl);
		}

		// TODO FIXME figure out segment limiting side of the member vector and add the shift vertex

		Edge * shiftDirection = nullptr;

		if( ! segmLimSides.spuckLowdimElmShiftDirectionIfUnclosedFractPresent( shiftDirection ) )
		{
			UG_LOG("No shift direction in ending crossing cleft segment " << std::endl);
			UG_THROW("No shift direction in ending crossing cleft segment " << std::endl);
		}

		IndexType numFoundECFSI = 0;

		// figure out the Ending crossing cleft segment info corresponding to this shift edge
		for( EndingCrossingFractureSegmentInfo const & ecfs : m_vecEndCrossFractSegmInfo )
		{
			if( ecfs.spuckLowdimElmShiftDirection() == shiftDirection )
				numFoundECFSI++;
		}

		if( numFoundECFSI != 1 )
		{
			UG_LOG("not one cefsi found " << numFoundECFSI << std::endl);
			UG_THROW("not one cefsi found " << numFoundECFSI << std::endl);
		}

		for( EndingCrossingFractureSegmentInfo & ecfs : m_vecEndCrossFractSegmInfo )
		{
			if( ecfs.spuckLowdimElmShiftDirection() == shiftDirection )
			{
				if( ! ecfs.schluckShiftVrtx( newShiftVrtx ) )
				{
					UG_LOG("Shift vertex wird nicht geschluckt" << m_aaPos[newShiftVrtx] << std::endl );
					UG_THROW("Shift vertex wird nicht geschluckt" << m_aaPos[newShiftVrtx] << std::endl );
					return false;
				}

				break;
			}
		}

		IndexType numFoundAssigned = 0;

		for( EndingCrossingFractureSegmentInfo const & ecfs : m_vecEndCrossFractSegmInfo )
		{
			if( ecfs.spuckShiftVrtx() != nullptr && ecfs.spuckLowdimElmShiftDirection() == shiftDirection )
				numFoundAssigned++;
		}

		if( numFoundAssigned != 1 )
		{
			m_sh.assign_subset( oldVrt, m_sh.num_subsets());
			m_sh.assign_subset( shiftDirection, m_sh.num_subsets());


			UG_LOG( "num found assigend wrong " << numFoundAssigned );
			UG_THROW( "num found assigend wrong " << numFoundAssigned );
		}

		UG_LOG("ending crossing cleft shift vertex generated " << std::endl);

	}

	std::vector<Volume*> volsInSegm;

	segmLimSides.spuckVecFulldimElem( volsInSegm );

	for( Volume * const & vol : volsInSegm )
	{
		std::vector<Vertex*> & newVrts4Fac = m_aaVrtVecVol[ vol ];

		for(size_t indVrt = 0; indVrt < (vol)->num_vertices();  indVrt++ )
		{
			Vertex* volVrt = (vol)->vertex(indVrt);

			if(  volVrt == oldVrt )
			{
				newVrts4Fac[ indVrt ] = newShiftVrtx;
			}
		}
	}

	return true;


}

//////////////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::computeCrossPointOfPlaneWithLine( PlaneDescriptor const & shiftedPlane, Edge * const & shiftDirectionEdg, Vertex * const & oldVrt, vector3 & posCrossingPt )
{
	vector3 const & normalShiftedPlane = shiftedPlane.spuckNormalVector();
	vector3 const & baseShiftedPlane = shiftedPlane.spuckBaseVector();
	number const & rhs = shiftedPlane.spuckRHS();

	// vector3 posOldVrt Aufpunkt der Edge, entland der verschoben werden soll
	// establish

	if( ! ( EdgeContains(shiftDirectionEdg, oldVrt ) ) )
	{
		UG_LOG("Shift edge ohne alten Vertex?" << std::endl);
		UG_THROW("Shift edge ohne alten Vertex?" << std::endl);
	}

	vector3 posOldVrt = m_aaPos[oldVrt];

	vector3 const & shiftVecBase = posOldVrt;

	Vertex * const & shiftEdgeVrtOne = shiftDirectionEdg->vertex(0);
	Vertex * const & shiftEdgeVrtTwo = shiftDirectionEdg->vertex(1);

	vector3 shiftVecEnd = shiftVecBase; // needs to be changed still

	if( shiftEdgeVrtOne == oldVrt )
	{
		shiftVecEnd = m_aaPos[shiftEdgeVrtTwo];
	}
	else if( shiftEdgeVrtTwo == oldVrt )
	{
		shiftVecEnd = m_aaPos[shiftEdgeVrtOne];
	}
	else
	{
		UG_LOG("no end vertex" << std::endl);
		UG_THROW("no end vertex" << std::endl);
	}

	vector3 directionVec;

	VecSubtract(directionVec, shiftVecEnd, posOldVrt );

	// x = aufpunkt + s * richtung einsetzen in ( x - x_0 ) * n = 0, auflösen nach r
	// x*n = rhs = x_0*n
	// ( a + s . r -x0 ) * n = 0
	// a * n + s . n * r - x0*n = 0
	// s . ( n * r ) = x0 * n - a * n
	// s = ( x0*n - a*n ) / ( n*r )

	number planeBaseMultNormal = VecDot(baseShiftedPlane,normalShiftedPlane);
	number lineBaseMultNormal = VecDot(posOldVrt,normalShiftedPlane);
	number lineDirMultNormal = VecDot( directionVec, normalShiftedPlane );

	UG_LOG("COMPARE RHS AND LBMN " << rhs << " - " << planeBaseMultNormal << std::endl);

	number fbs = std::fabs( ( rhs - planeBaseMultNormal ) / ( rhs + planeBaseMultNormal) );

	UG_LOG("RELATIVE DIFF " << fbs << std::endl);

	number const d_toler = 1e-4;

	if( fbs > d_toler )
	{
		UG_LOG("RELATIVE DIFF TOO BIG " << std::endl);
		UG_THROW("RELATIVE DIFF TOO BIG " << std::endl);
	}

	if( lineDirMultNormal == 0. )
	{
		UG_THROW("denominator zero division " << std::endl);
		UG_LOG("denominator zero division " << std::endl);
	}

	number moveAlongLine = ( planeBaseMultNormal - lineBaseMultNormal ) / lineDirMultNormal;

	vector3 shiftAlongLine;

	VecScale( shiftAlongLine, directionVec, moveAlongLine );

	VecAdd( posCrossingPt, posOldVrt, shiftAlongLine );

	return true;

}
//////////////////////////////////////////////////////////////////////////


#if 0
template<>
bool ArteExpandFracs3D::expandWithinTheSegment<1,false>( Vertex * const & oldVrt, SegmentVolElmInfo const & segmVolElmInfo )
{

	// get all fracture faces, check if they belong to the same subdomain, must be the case here!

	VecAttachedFractFaceEdgeSudo vecAttFractFaces;

	for( AttachedVolumeElemInfo const & avei : segmVolElmInfo )
	{
		VecAttachedFractFaceEdgeSudo const & vecAttFractVol = avei.getVecFractManifElem();

		for( AttachedFractFaceEdgeSudo const & affe : vecAttFractVol )
		{
			vecAttFractFaces.push_back(affe);
		}
	}

	IndexType numbContrFracFaces = vecAttFractFaces.size();

	if( numbContrFracFaces < 1 )
	{
		UG_LOG("Kein Affe da " << std::endl);
		UG_THROW("Kein Affe da " << std::endl);
		return false;
	}

	IndexType sudoBase = vecAttFractFaces[0].getSudo();

	// check if all sudos equal
	for( AttachedFractFaceEdgeSudo const & affe : vecAttFractFaces )
	{
		IndexType sudoTest = affe.getSudo();

		if( sudoTest != sudoBase )
		{
			UG_LOG("unterschiedliche Sudos an einer einzelnen Fracture?" << std::endl);
			UG_THROW("unterschiedliche Sudos an einer einzelnen Fracture?" << std::endl);
			return false;
		}
	}

	// now we are sure we have the same sudo, now we average the normals

	NormalVectorFacIntoVol normalsFacInVolSummed(0,0,0);

	for( AttachedFractFaceEdgeSudo const & affe : vecAttFractFaces )
	{
		NormalVectorFacIntoVol tmpVec = normalsFacInVolSummed;
		UG_LOG("Normal Vec " << tmpVec << std::endl);
		VecAdd(normalsFacInVolSummed, tmpVec, affe.getNormalVec() );
	}

	NormalVectorFacIntoVol normalsAveraged;

	VecScale( normalsAveraged, normalsFacInVolSummed, 1. / ( static_cast<number>(numbContrFracFaces) ) );

	number width = m_fracInfosBySubset[sudoBase].width;

	number scal = width / 2.;

	NormalVectorFacIntoVol scaledNormal;

	VecScale( scaledNormal, normalsAveraged, - scal );

	vector3 posOldVrt = m_aaPos[oldVrt];

	vector3 posNewVrt;

	VecAdd( posNewVrt, posOldVrt, scaledNormal );

	Vertex * newShiftVrtx = *m_grid.create<RegularVertex>();

	if( newShiftVrtx == nullptr )
	{
		UG_LOG("Nullen erzeugt" << std::endl);
		UG_THROW("Nullen erzeugt" << std::endl);
		return false;
	}

	m_aaPos[newShiftVrtx] = posNewVrt;

	UG_LOG("Created new vertex at " << m_aaPos[newShiftVrtx] << std::endl );

	m_sh.assign_subset(newShiftVrtx, sudoBase);

	for( AttachedVolumeElemInfo const & avei : segmVolElmInfo )
	{
		Volume * vol = avei.getFulldimElem();

		std::vector<Vertex*> & newVrts4Fac = m_aaVrtVecVol[ vol ];

		for(size_t indVrt = 0; indVrt < (vol)->num_vertices();  indVrt++ )
		{
			Vertex* volVrt = (vol)->vertex(indVrt);

			if(  volVrt == oldVrt )
			{
				newVrts4Fac[ indVrt ] = newShiftVrtx;
			}
		}
	}

	return true;
}
#endif


////////////////////////////////////////////////////////////////////

//template<ArteExpandFracs3D::SegmentVrtxFracStatus::oneFracSuDoAtt>
//bool ArteExpandFracs3D::computeShiftVector( ArteExpandFracs3D::VecSegmentLimitSidesPairSudoNorml const & vecSegmLimSidPrSudoNrml )
//{
//	return {};
//}


////////////////////////////////////////////////////////////////////

// TODO FIXME muss verschwinden!
//template<>
//bool ArteExpandFracs3D::computeShiftVector( VecSegmentLimitSidesPairSudoNorml const & vecSegmLimSidPrSudoNrml )
////bool ArteExpandFracs3D::expandWithinTheSegment<ArteExpandFracs3D::SegmentVrtxFracStatus::threeFracSuDoAtt>( SegmentLimitingSides const & segmLimSides )
//{
//	// wenn es nur zwei Segmentnormalen und sudos sind, dann wird die dritte auf Normale senkrecht zu
//	// den beiden anderen gesetzt, und um Null verschoben, ihr
//
////	if( segmLimSides.isBoundary() )
////	{
////		UG_LOG("three fracture at boundary need implementation " << std::endl);
////		UG_THROW("three fracture at boundary need implementation " << std::endl);
////		return false;
////	}
//
////	VecSegmentLimitSidesPairSudoNorml vecSegmLimSidPrSudoNrml;
//
//	Vertex * oldVrt = segmLimSides.spuckVertex();
//
//	if( ! segmLimSides.spuckFractSudoNormls( vecSegmLimSidPrSudoNrml ) )
//	{
//		UG_LOG("Spucken schief gegangen three " << std::endl);
//		UG_THROW("Spucken schief gegangen three " << std::endl);
//		return false;
//	}
//
//	if( vecSegmLimSidPrSudoNrml.size() != 3 )
//	{
//		UG_LOG("three fractures, but not three normals and sudo?"  <<  std::endl);
//		UG_THROW("three fracture, but not three normal and sudo?" << std::endl);
//		return false;
//	}
//
//	// select one of the sudos of the sourrounding sudos of the segment, to which we want to assign the new vertex
//	IndexType sudoExample = (vecSegmLimSidPrSudoNrml[0]).first;
//
////	std::vector<IndexType> vecSudosSides;
////	std::vector<vector3> vecNormalsAveraged;
////	std::vector<vector3> vecShiftVec4Plane;
//
//	vector3 posOldVrt = m_aaPos[oldVrt];
//
//	VecPlaneDescriptor vecShiftedPlaneDescript;
//
////	IndexType side = 0;
//
//	for( SegmentLimitSidesPairSudoNorml const & segLimSidPrSN : vecSegmLimSidPrSudoNrml )
//	{
//		IndexType const & sudoSide = segLimSidPrSN.first;
//		vector3 const & normlAvrg = segLimSidPrSN.second;
//
//		// normal computed standard mässig directs of of the volume, we need that one into the volume
//		vector3 normalOutsideVol;
//
//		// hier testen, ob die Subdomain von der Liste der fracture subdomains ist,
//		// damit für den Fall von sich zwei schneidenden Ebenen für die dritte,
//		// die als senkrecht zu den beiden anderen gesetzt werden soll, mit Verschiebung null,
//		// analog Boundary sides, die künstliche Weite null erhalten kann hier
//		number width = 0;
//
//		// ensure that the width is nonzero only for real fractures, not for pseudo vectors or such stuff
//		if( ( side == 1 && ! artificialNormalTwo ) ||  ( side == 2 && ! artificialNormalThree ) )
//			width = m_fracInfosBySubset[sudoSide].width;
//
//		side++;
//
//		number shiftScal = width / 2.;
//
//		vector3 shiftVec4Plane;
//
//		// as normal outside vol, but we need to go inside the volumes / segment side
//		VecScale( shiftVec4Plane, normalOutsideVol, - shiftScal );
//
////		vecShiftVec4Plane.push_back( shiftVec4Plane );
////		vecSudosSides.push_back( sudoSide );
////		vecNormalsAveraged.push_back( normalIntoVol );
//
//		PlaneDescriptor planeDescr( normalOutsideVol, posOldVrt );
//
//		PlaneDescriptor shiftedPlaneDescr;
//
//		planeDescr.spuckPlaneShiftedAlong( shiftVec4Plane, shiftedPlaneDescr );
//
//		vecShiftedPlaneDescript.push_back( shiftedPlaneDescr );
//	}
//
//	// KÄSE alles!!! wenn es eine boundary ist oder eine Kreuzung von zwei Ebenen innendrin,
//	// dann muss hier DANACH noch was dazu kommen, das heisst, alles, was hier gemacht wird,
//	// soll für den Fall sein, dass es keine einzelne Fracture innen drin ist,
//	// und alle bool template Parameter sind MÜLL!!!!
//	// das für zwei Fractures muss hier auch dazu
//	// vielleicht am einfachsten, wenn man den template parameter auf true oder false stellt
//	// für eine fracture innen oder allen anderen Rest......
//
//	// TODO FIXME compute crossing point of the three new planes!!!! need functions from lib_disc, cf FiniteStrainMechanics, solve LGS exact
//
//	vector3 shiftedCrossingPoint;
//
//	computeCrossingPointOf3Planes( vecShiftedPlaneDescript, shiftedCrossingPoint );
//
//	return true;
//}

////////////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::computeCrossingPointOf3Planes(  VecPlaneDescriptor const & vecPlaneDescr, vector3 & crossingPoint )
{
	// n_i * vecX_i = a_i * n_i , i=1,2,3, ohne Summenkonvention, und n_i, vecX_i, a_i jeweils Vektoren Grösse 3
	// es entsteht also ein LGS für drei Unbekannte und drei Gleichungen

//	 vector<vector<double>> coefficients = {{2, 1}, {1, -3}};
//	 vector<double> constants = {5, -1};
//	 vector<double> solutions = cramer_rule(coefficients, constants);
//	 cout << "Solution for x: " << solutions[0] << endl;
//	 cout << "Solution for y: " << solutions[1] << endl;

	PlaneDescriptor const & planeDescrZero = vecPlaneDescr[0];
	PlaneDescriptor const & planeDescrOne = vecPlaneDescr[1];
	PlaneDescriptor const & planeDescrTwo = vecPlaneDescr[2];

	vector3 const & normalZero = planeDescrZero.spuckNormalVector();
	vector3 const & normalOne = planeDescrOne.spuckNormalVector();
	vector3 const & normalTwo = planeDescrTwo.spuckNormalVector();

	number const & rhsZero = planeDescrZero.spuckRHS();
	number const & rhsOne = planeDescrOne.spuckRHS();
	number const & rhsTwo = planeDescrTwo.spuckRHS();



	std::vector<std::vector<double>> coefficients = { {normalZero[0], normalZero[1], normalZero[2]},
													  {normalOne[0], normalOne[1], normalOne[2]},
													  {normalTwo[0], normalTwo[1], normalTwo[2]},
													};

	std::vector<double> constants = { rhsZero, rhsOne, rhsTwo };

	std::vector<double> solutions = ug::simpleMatrOps::cramerRule(coefficients, constants);

	crossingPoint = vector3( solutions[0], solutions[1], solutions[2] );

	return true;
}

////////////////////////////////////////////////////////////////////



#if 0
//template <>
//bool ArteExpandFracs3D::establishNewVertices< Tetrahedron,
//											  ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt
//											>( Vertex * const & oldVrt )
//template< bool APPLY_GENERAL_SEGMENT_ORDERING,
//		  ArteExpandFracs3D::VrtxFracProptsStatus vfp,
//		  typename std::enable_if< std::integral_constant<bool,APPLY_GENERAL_SEGMENT_ORDERING>>
//		>
template <>
bool ArteExpandFracs3D::establishNewVertices< true,
											  ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt
											>( Vertex * const & oldVrt )
{
	UG_LOG("under construction Tetrahedra limited" << std::endl);

	VecVertFracTrip const & vecVertFracTrip = m_aaVrtInfoFraTri[oldVrt];

	VecAttachedVolumeElemInfo const & vecAttVolElemInfo = m_aaVolElmInfo[oldVrt];

	VecAttachedVolumeElemInfo vecAttVolElemInfoCop = vecAttVolElemInfo; // echte KOPIE

	VecAttachedVolumeElemInfo reconstructedVecAttVolElmInf;

	VecSegmentVolElmInfo vecSegVolElmInfo;

		/*
		 * While Schleifen aufbauen für den
		 * Search the adjacent surface interatively - Algorithmus
		 * (Stasi Algorithmus)
		 *
		 */

	IndexType d_segmenteErledigt = 0;

	while( vecAttVolElemInfoCop.size() != 0 )
	{
		SegmentVolElmInfo segmentAVEI;

		AttachedVolumeElemInfo & startVolInfoThisSegment = vecAttVolElemInfoCop[0];

		startVolInfoThisSegment.markIt();

		Volume * volSta = startVolInfoThisSegment.getFulldimElem();

		vector3 center;

		if( volSta != nullptr )
			center = CalculateCenter(volSta,m_aaPos);

//		UG_LOG("volume center " << center << std::endl );

		int d_loopsDone = 0;

		while( vecAttVolElemInfoCop.size() != 0 )
		{
			// count number of marked elements
			IndexType numMarkedElems = 0;
			IndexType markPoint = 0;

//			IndexType lastMarkPt = 0;
			IndexType startIndexInner = 0;

			for( AttachedVolumeElemInfo const & volElInfCop : vecAttVolElemInfoCop )
			{
				if( volElInfCop.isMarked() )
				{
					Volume * vol = volElInfCop.getFulldimElem();
//					m_sh.assign_subset(vol, m_sh.num_subsets());

					vector3 center = CalculateCenter(vol,m_aaPos);

//					UG_LOG("DAS ZENTRUM " << numMarkedElems << " -> " << center << std::endl);

					startIndexInner = markPoint;
					numMarkedElems++;

				}

				markPoint++;
			}

			UG_LOG("LOOPS DONE " << numMarkedElems << std::endl);

			if( numMarkedElems == 0 )
				break;

//			int startIndexInner = -1;
//
//			for( int i = 0; i < vecAttVolElemInfoCop.size(); i++ )
//			{
//				AttachedVolumeElemInfo vi = vecAttVolElemInfoCop[i];
//
//				Volume * vo = vi.getFulldimElem();
//
//				vector3 center = CalculateCenter(vo,m_aaPos);
//
//				UG_LOG("DAS ZENTRUM ZAHL VOR " << i << " -> " <<  center << std::endl);
//
//				if( vi.isMarked() )
//					startIndexInner = i;
//			}
//
//			if( startIndexInner < 0 )
//			{
//				UG_THROW("kein Anfang gefunden " << std::endl);
//			}
//
//#if 0
//			IndexType startIndexInner = markPoint - 1;
//#endif
			AttachedVolumeElemInfo startVolInfoMarkLoop = vecAttVolElemInfoCop[startIndexInner];

			Volume * stattVoll = startVolInfoMarkLoop.getFulldimElem();

			vector3 centerX = CalculateCenter(stattVoll,m_aaPos);

			UG_LOG("DAS ZENTRUM DANACH " << startIndexInner << " -> " <<  centerX << std::endl);

//			m_sh.assign_subset(stattVoll, m_sh.num_subsets());
#if 0
			for( int i = 0; i < vecAttVolElemInfoCop.size(); i++ )
			{
				AttachedVolumeElemInfo vi = vecAttVolElemInfoCop[i];

				Volume * vo = vi.getFulldimElem();

				vector3 center = CalculateCenter(vo,m_aaPos);

				UG_LOG("DAS ZENTRUM ZAHL " << i << " -> " <<  center << std::endl);

			}
#endif
			for( AttachedVolumeElemInfo const & possibleOrigVolInfo : vecAttVolElemInfo )
			{
				if( possibleOrigVolInfo.hasSameFulldimElem( startVolInfoMarkLoop ) )
				{
					segmentAVEI.push_back(possibleOrigVolInfo);
					reconstructedVecAttVolElmInf.push_back(possibleOrigVolInfo);
					break;
				}
			}

			vecAttVolElemInfoCop.erase( vecAttVolElemInfoCop.begin() + startIndexInner );

//			if( d_loopsDone == 1 )
//				return false;

			for( VecAttachedVolumeElemInfo::iterator aveiIt = vecAttVolElemInfoCop.begin();
													 aveiIt < vecAttVolElemInfoCop.end();
													 aveiIt++
			)
			{
				AttachedVolumeElemInfo & possibleNeighbour = *aveiIt;

				if( possibleNeighbour.hasSameFulldimElem( startVolInfoMarkLoop ) )
				{
					continue;
				}
				else
				{
					bool neighbourFound = possibleNeighbour.testFullDimElmNeighbour( startVolInfoMarkLoop );

					if( neighbourFound )
					{
						Volume * vol = possibleNeighbour.getFulldimElem();

//						m_sh.assign_subset(vol, m_sh.num_subsets());

					}
				}
			}


			d_loopsDone++;


		}

		vecSegVolElmInfo.push_back(segmentAVEI);

//		d_segmenteErledigt++;
//
//		if( d_segmenteErledigt == 1 )
//		return false;
	}

	if( reconstructedVecAttVolElmInf.size() != vecAttVolElemInfo.size() )
	{
		UG_LOG("Rekonstruktion schief gegangen " << std::endl);
		UG_THROW("Rekonstruktion schief gegangen " << std::endl);
		return false;
	}

	for( SegmentVolElmInfo const & svei : vecSegVolElmInfo )
	{
		// TODO FIXME das hier wieder entfernen, die Subdomain Zuweisung, nur für debug Zwecke
		IndexType sudoMax = m_sh.num_subsets();

		for( AttachedVolumeElemInfo const & vei : svei )
		{
			Volume * vol = vei.getFulldimElem();

			m_sh.assign_subset( vol, sudoMax );
		}
	}

	return true;
}


////////////////////////////////////////////////////////////////////



//template <>
//bool ArteExpandFracs3D::establishNewVertices< Hexahedron,
//											  ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt
//											>( Vertex * const & oldVrt )
//template< bool APPLY_GENERAL_SEGMENT_ORDERING,
//		  ArteExpandFracs3D::VrtxFracProptsStatus vfp,
//		  typename std::enable_if< std::integral_constant<bool,!APPLY_GENERAL_SEGMENT_ORDERING>>
//		>
template <>
bool ArteExpandFracs3D::establishNewVertices< false,
											  ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt
											>( Vertex * const & oldVrt )
{
	VecVertFracTrip & vecVertFracTrip = m_aaVrtInfoFraTri[oldVrt];

//	std::vector<Edge*> & allAssoEdges = m_aaVrtInfoAssoEdges[oldVrt];
//	std::vector<Face*> & allAssoFaces = m_aaVrtInfoAssoFaces[oldVrt];
//	std::vector<Volume*> & allAssoVolumes = m_aaVrtInfoAssoVols[oldVrt];

	// TODO FIXME works if at all only for very simple geometries
	// and works only in particular if all asso volumes are part of the triple, i.e. have all a common face with the fracture
	// this is not selbstverständlich at all!!!!

	VecVertFracTrip  vecVertFracTripCopy = m_aaVrtInfoFraTri[oldVrt]; // copy, not reference!


	VecVertFracTrip firstSegment;
	VecVertFracTrip secondSegment;

	IndexType beginIndex = 0;

	VertFracTrip startTrip = vecVertFracTripCopy[beginIndex];

	firstSegment.push_back( startTrip );

	vector3 normalOne = startTrip.getNormal();
	vector3 normalTwo;
	// vermute einfach umgekehrte Normale, ein Trick, der im einfachsten Fall geht......
	VecScale(normalTwo,normalOne,-1);

	vecVertFracTripCopy.erase( vecVertFracTripCopy.begin() + beginIndex );

//	for( VecVertFracTrip::iterator itVFT  = vecVertFracTripCopy.begin();
//									   itVFT != vecVertFracTripCopy.end();
//									   itVFT++
//		)
//		{
//			VertFracTrip actualVFT = *itVFT;
//
//			vector3 volCenter = CalculateCenter(actualVFT.getFullElm(),m_aaPos);
//			UG_LOG("Vol center" << volCenter << std::endl);
//			vecVertFracTripCopy.erase(itVFT);
//		}
//
//	return true;

//	while( vecVertFracTripCopy.size() != 0 )
//	{
//		for( VecVertFracTrip::iterator itVFT  = vecVertFracTripCopy.begin();
//									   itVFT != vecVertFracTripCopy.end();
//									   itVFT++
//		)
	for( auto const & actualVFT : vecVertFracTripCopy )
	{
//			VertFracTrip actualVFT = *itVFT;

		vector3 normalActual = actualVFT.getNormal();

		vector3 volCenter = CalculateCenter(actualVFT.getFullElm(),m_aaPos);
		UG_LOG("Vol center" << volCenter << std::endl);

		// test direction

		number cosinus2One = VecDot(normalOne,normalActual);
		number cosinus2Two = VecDot(normalTwo,normalActual);

		UG_LOG("normal eins " << normalOne << std::endl);
		UG_LOG("normal zwei " << normalTwo << std::endl);
		UG_LOG("normal actu " << normalActual << std::endl);

		UG_LOG("cosi one " << cosinus2One << std::endl );
		UG_LOG("cosi two " << cosinus2Two << std::endl );


		// if cosinus > 0, assume same side

		if( ( cosinus2One >= 0 && cosinus2Two >= 0 ) || ( cosinus2One <= 0 && cosinus2Two <= 0 ) )
			UG_THROW("kann nicht auf zwei Seiten hinken" << std::endl);

		if( cosinus2One >= 0 )
		{
			firstSegment.push_back( actualVFT );
		}
		else if( cosinus2Two >= 0 )
		{
			secondSegment.push_back( actualVFT );
		}
		else
		{
			UG_THROW("muss wo dazu gehoeren wohl" << std::endl);
		}
	}
//			vecVertFracTripCopy.erase(itVFT);
//		}
//	}

	// computer averaged normal

	vector3 normalsOneSummed(0,0,0);
	vector3 normalsTwoSummed(0,0,0);


	for( auto const & seg: firstSegment )
	{
		vector3 tmpVec = normalsOneSummed;
		VecAdd(normalsOneSummed,tmpVec,seg.getNormal());
	}

	for( auto const & seg: secondSegment )
	{
		vector3 tmpVec = normalsTwoSummed;
		VecAdd(normalsTwoSummed,tmpVec,seg.getNormal());
	}

	vector3 normalsOneAveraged;
	vector3 normalsTwoAveraged;

	if( firstSegment.size() != 0 )
	{
		VecScale(normalsOneAveraged, normalsOneSummed, 1./firstSegment.size());
	}

	if( secondSegment.size() != 0 )
	{
		VecScale(normalsTwoAveraged, normalsTwoSummed, 1./secondSegment.size());
	}

	// get to know width of fracture

	IndexType suse = startTrip.getSudoElm();

	number width = m_fracInfosBySubset[suse].width;

	number scal = width / 2.;

	vector3 scaledNormalOne, scaledNormalTwo;
	VecScale( scaledNormalOne, normalOne, - scal );
	VecScale( scaledNormalTwo, normalTwo, - scal );
	// Minuszeichen wichtig, sonst wird in die falsche Richtung gedrückt, und die Volumen gehen über die fracture
	// raus und werden grösser, anstatt kleiner zu werden.....

	vector3 posOldVrt = m_aaPos[oldVrt];

	vector3 posNewVrtOne, posNewVrtTwo;

	VecAdd( posNewVrtOne, posOldVrt, scaledNormalOne);
	VecAdd( posNewVrtTwo, posOldVrt, scaledNormalTwo);

	Vertex * newShiftVrtxOne = *m_grid.create<RegularVertex>();
	Vertex * newShiftVrtxTwo = *m_grid.create<RegularVertex>();

	m_aaPos[newShiftVrtxOne] = posNewVrtOne;
	m_aaPos[newShiftVrtxTwo] = posNewVrtTwo;

	UG_LOG("Created new vertex 1 at " <<m_aaPos[newShiftVrtxOne] << std::endl );
	UG_LOG("Created new vertex 2 at " <<m_aaPos[newShiftVrtxTwo] << std::endl );

	m_sh.assign_subset(newShiftVrtxOne, suse);
	m_sh.assign_subset(newShiftVrtxTwo, suse);

//	m_sh.assign_subset(newShiftVrtxOne, 3);
//	m_sh.assign_subset(newShiftVrtxTwo, 3);

	for( auto const & fs : firstSegment )
	{
		Volume * vol = fs.getFullElm();

		std::vector<Vertex*>& newVrts4Fac = m_aaVrtVecVol[ vol ];

		for(size_t indVrt = 0; indVrt < (vol)->num_vertices();  indVrt++ )
		{
			Vertex* volVrt = (vol)->vertex(indVrt);

			if(  volVrt == oldVrt )
			{
				newVrts4Fac[ indVrt ] = newShiftVrtxOne;
			}
		}
	}

	for( auto const & ses : secondSegment )
	{
		Volume * vol = ses.getFullElm();

		std::vector<Vertex*>& newVrts4Fac = m_aaVrtVecVol[ vol ];

		for(size_t indVrt = 0; indVrt < (vol)->num_vertices();  indVrt++ )
		{
			Vertex* volVrt = (vol)->vertex(indVrt);

			if(  volVrt == oldVrt )
			{
				newVrts4Fac[ indVrt ] = newShiftVrtxTwo;
			}
		}
	}


	return true;
}

#endif

////////////////////////////////////////////////////////////////////


bool ArteExpandFracs3D::createConditionForNewVrtcs()
{

	//	iterate over all surrounding volumes to enable volume changes, this loop taken from SR but shortened
	for(VolumeIterator iterSurrVol = m_sel.volumes_begin(); iterSurrVol != m_sel.volumes_end(); iterSurrVol++ )
	{
		Volume * sv = *iterSurrVol;

		std::vector<Vertex*>& newVrts = m_aaVrtVecVol[sv];
		newVrts.resize(sv->num_vertices());

		for(size_t iVrt = 0; iVrt < sv->num_vertices(); iVrt++ )
		{
			newVrts[iVrt] = nullptr;
		}
			// erstmal so tun, als ob keine neuen Vertizes erzeugt werden an den alten Vertizes
	}


	return true;
}

/////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::createNewElements()
{
	// originally practically copied from Sebastian, as this concept was fine for 2D, but adapted in 3D

	//	create new elements

	UG_LOG("want to create new elems" << std::endl );

	//	holds local side vertex indices
	std::vector<size_t>	locVrtInds;

	//	first we create new edges from selected ones which are connected to
	//	inner vertices. This allows to preserve old subsets.
	//	Since we have to make sure that we use the right vertices,
	//	we have to iterate over the selected volumes and perform all actions on the edges
	//	of those volumes.
	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end(); ++iter_sv)
	{
		Volume* sv = *iter_sv;

		UG_LOG("entering volume to create new elems " << CalculateCenter(sv, m_aaPos) << std::endl);

		//	check for each edge whether it has to be copied.
		for(size_t i_edge = 0; i_edge < sv->num_edges(); ++i_edge)
		{
			Edge* e = m_grid.get_edge(sv, i_edge);

			if(m_sel.is_selected(e))
			{
				//	check the associated vertices through the volumes aaVrtVecVol attachment.
				//	If at least one has an associated new vertex and if no edge between the
				//	new vertices already exists, we'll create the new edge.
				size_t ind0, ind1;
				sv->get_vertex_indices_of_edge(ind0, ind1, i_edge);
				Vertex* nv0 = (m_aaVrtVecVol[sv])[ind0];
				Vertex* nv1 = (m_aaVrtVecVol[sv])[ind1];

				if(nv0 || nv1)
				{
					//	if one vertex has no associated new one, then we use the vertex itself
					if(!nv0)
						nv0 = sv->vertex(ind0);
					if(!nv1)
						nv1 = sv->vertex(ind1);

#if 1
					//	create the new edge if it not already exists.
					if( ! m_grid.get_edge(nv0, nv1))
						m_grid.create_by_cloning(e, EdgeDescriptor(nv0, nv1), e);
#else // only for debugging purposes

					//					if( ! ( m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft[nv0] || m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft[nv0] ) )
					{
						//	create the new edge if it not already exists.
						if( ! m_grid.get_edge(nv0, nv1))
						{
							Edge * cloneEdge = * ( m_grid.create_by_cloning(e, EdgeDescriptor(nv0, nv1), e) );

							if( ( m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft[nv0] || m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft[nv1] ) )
							{
								m_sh.assign_subset(cloneEdge, m_sh.num_subsets());
								UG_LOG("Edge clone from ending crossing cleft new shift vertex" << std::endl);

							}
						}
					}
					//					else
					//					{
					//						UG_LOG("Edge clone from ending crossing cleft new shift vertex" << std::endl);
					//					}

					// Ergebnis dieses debug Teils: alle edges, die geclont werden, sind notwendig
					// auch bei ending crossing clefts
#endif
				}
			}
		}
	}

	UG_LOG("Vol enter clone finished " << std::endl);

	//	now we create new faces from selected ones which are connected to
	//	inner vertices. This allows to preserve old subsets.
	//	Since we have to make sure that we use the right vertices,
	//	we have to iterate over the selected volumes and perform all actions on the side-faces
	//	of those volumes.

	FaceDescriptor fd;


	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end(); ++iter_sv)
	{
		Volume* sv = *iter_sv;
		//	check for each face whether it has to be copied.

		UG_LOG("Face descriptor for vol " << CalculateCenter(sv, m_aaPos) << std::endl);

		for(size_t i_face = 0; i_face < sv->num_faces(); ++i_face)
		{
			Face* sf = m_grid.get_face(sv, i_face);

			if( m_sel.is_selected(sf))
			{
				//	check the associated vertices through the volumes aaVrtVecVol attachment.
				//	If no face between the new vertices already exists, we'll create the new face.
				sv->get_vertex_indices_of_face(locVrtInds, i_face);
				fd.set_num_vertices(sf->num_vertices());

				bool containsEndingCrossingCleftVrtx = false;

				for(size_t i = 0; i < sf->num_vertices(); ++i)
				{
					Vertex* nVrt = (m_aaVrtVecVol[sv])[locVrtInds[i]];

					if(nVrt)
						fd.set_vertex(i, nVrt);
					else
						fd.set_vertex(i, sv->vertex(locVrtInds[i]));
				}

#if 1
				//	if the new face does not already exist, we'll create it
				if(!m_grid.get_face(fd))
					m_grid.create_by_cloning(sf, fd, sf);
#else
				//	if the new face does not already exist, we'll create it
				if( ! m_grid.get_face(fd) )
				{
					Face * fac = *( m_grid.create_by_cloning(sf, fd, sf) );

					bool faceAtEndingCrossingCleft = m_aaMarkFaceWithEndingCrossingCleft[sf];

					bool faceTouchingEndingCrossingCleft = m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft[sf];

					if( faceAtEndingCrossingCleft || faceTouchingEndingCrossingCleft )
					{
						m_sh.assign_subset( fac, m_sh.num_subsets() );
						UG_LOG("Face at ending crossing cleft" << std::endl);
					}
				}

				// Ergebnis dieses debug Teils: alle faces, die geclont werden, sind notwendig
				// auch bei ending crossing clefts
#endif

			}
		}
	}

	UG_LOG("Face descriptor left" << std::endl);

	//	Expand all faces.
	//	Since volumes are replaced on the fly, we have to take care with the iterator.
	//	record all new volumes in a vector. This will help to adjust positions later on.

	std::vector<Volume*> newFractureVolumes;
	std::vector<IndexType> subsOfNewVolumes;

	// create alternative volumes where there are ending crossing clefts

//	etablishVolumesAtEndingCrossingClefts( std::vector<Volume*> & newFractureVolumes, std::vector<IndexType> & subsOfNewVolumes );
	if( ! etablishVolumesAtEndingCrossingClefts( newFractureVolumes, subsOfNewVolumes ) )
	{
		UG_LOG("unable to establish volumes at ending crossing clefts" << std::endl);
		return false;
	}

//	return false;

	VolumeDescriptor vd;

	// beim Volumen fängt das Abfangen an, es muss abgefragt werden, ob das Volumen aus
	// einem Segment ist, wo bisher falsch expandiert wird
	// erst beim Face ab zu fangen ist viel zu spät TODO FIXME
	// nicht zu vergessen, die Edges ordentlich sortiert zu sammeln, die zwei endende Vertizes enthalten,
	// und dann zu splitten, wenn alle Attachements entfernt sind, dann Neustart anfordern mit neuer Geometrie

	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end();)
	{
		Volume* sv = *iter_sv;
		++iter_sv;

		UG_LOG("Volume new creation try at " << CalculateCenter(sv, m_aaPos) << std::endl);

		bool volHasEndingCrossingCleftFace = m_volAttAccsVolTouchesEndingCrossingCleft[sv];

		//	now expand the fracture faces of sv to volumes.
		for(size_t i_side = 0; i_side < sv->num_sides(); ++i_side)
		{
			//	get the local vertex indices of the side of the volume
			sv->get_vertex_indices_of_face(locVrtInds, i_side);

			Face* tFace = m_grid.get_side(sv, i_side);

			if(tFace)
			{
				if( m_aaMarkFaceIsFracB[tFace] ) // && ! m_aaMarkFaceHasUnclosedFracSideB[tFace] )
				{
					bool faceIsSegmLimEndCrossCleft = m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft[tFace];
					bool faceIsEndingCleftCrossFace = m_aaMarkFaceWithEndingCrossingCleft[tFace];

					bool avoidFace = ( volHasEndingCrossingCleftFace && ( faceIsSegmLimEndCrossCleft || faceIsEndingCleftCrossFace ) );

//					for( Face * testFac :  m_d_endingCrossingCleftFaces )
//					{
//						if( testFac == tFace )
//						{
//							avoidFace = true;
//						}
//					}
//
//					for( Face * testFac :  m_d_crossingNeighboredNotEndingFaces )
//					{
//						if( testFac == tFace )
//						{
//							avoidFace = true;
//						}
//					}
//
//					for( Face * testFac :  m_d_crossingNeighboredNotEndingFacesCommEdg )
//					{
//						if( testFac == tFace )
//						{
//							avoidFace = true;
//						}
//					}
//
//					for( Face * testFac :  m_d_notEndingCrossingFacesNotNeighbour )
//					{
//						if( testFac == tFace )
//						{
//							avoidFace = true;
//						}
//					}


					constexpr bool debugTest = true;

					if( avoidFace && debugTest )
						continue;

					Volume* expVol = nullptr;

					if(locVrtInds.size() == 3)
					{
						size_t iv0 = locVrtInds[0];
						size_t iv1 = locVrtInds[1];
						size_t iv2 = locVrtInds[2];

						if(    ( m_aaVrtVecVol[sv] )[iv0]
							&& ( m_aaVrtVecVol[sv] )[iv1]
							&& ( m_aaVrtVecVol[sv] )[iv2]
						)
						{
							//	create a new prism
							expVol = *m_grid.create<Prism>(
											PrismDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
															(m_aaVrtVecVol[sv])[iv2],
															(m_aaVrtVecVol[sv])[iv1],
															(m_aaVrtVecVol[sv])[iv0]));
						}
						else if(    ( m_aaVrtVecVol[sv] )[iv0]
								 && ( m_aaVrtVecVol[sv] )[iv1]
						)
						{
							//	create a new Pyramid
							expVol = *m_grid.create<Pyramid>(
											PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
												(m_aaVrtVecVol[sv])[iv1],
												(m_aaVrtVecVol[sv])[iv0],
												sv->vertex(iv2)));
						}
						else if(    ( m_aaVrtVecVol[sv] )[iv1]
								 && ( m_aaVrtVecVol[sv] )[iv2]
						)
						{
							//	create a new Pyramid
							expVol = *m_grid.create<Pyramid>(
											PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
												(m_aaVrtVecVol[sv])[iv2],
												(m_aaVrtVecVol[sv])[iv1],
												sv->vertex(iv0)));
						}
						else if(    (m_aaVrtVecVol[sv])[iv0]
								 && (m_aaVrtVecVol[sv])[iv2]
						)
						{
							//	create a new Pyramid
							expVol = *m_grid.create<Pyramid>(
											PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
												(m_aaVrtVecVol[sv])[iv0],
												(m_aaVrtVecVol[sv])[iv2],
												sv->vertex(iv1)));
						}
						else if( ( m_aaVrtVecVol[sv])[iv0] )
						{
							//	create a new Tetrahedron
							expVol = *m_grid.create<Tetrahedron>(
											TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																 (m_aaVrtVecVol[sv])[iv0]));
						}
						else if( ( m_aaVrtVecVol[sv])[iv1] )
						{
							//	create a new Tetrahedron
							expVol = *m_grid.create<Tetrahedron>(
											TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																 (m_aaVrtVecVol[sv])[iv1]));
						}
						else if( ( m_aaVrtVecVol[sv])[iv2] )
						{
							//	create a new Tetrahedron
							expVol = *m_grid.create<Tetrahedron>(
											TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																 (m_aaVrtVecVol[sv])[iv2]));
						}
						else
						{
							//	Text from SR, still similar:
							//  this code-block should never be entered. If it is entered then
							//	we either selected the wrong faces (this shouldn't happen), or there
							//	are selected faces, which have fracture-boundary-vertices only.
							//	This is the same is if inner fracture edges exists, which are
							//	connected to two boundary vertices.
							//	Since we tried to remove those edges above, something went wrong.
							//	remove the temporary attachments and throw an error

							UG_LOG("Tetraeder Fehlt eine Loesung " << std::endl);
#if 1
							detachMarkers();
							throw(UGError("Error in ExpandFractures3d Arte Stasi. Implementation Error."));
							return false;
#endif
						}
					}
					else if ( locVrtInds.size() == 4 )
					{
						// newly implemented by Markus to test with Hexahedrons

						size_t iv0 = locVrtInds[0];
						size_t iv1 = locVrtInds[1];
						size_t iv2 = locVrtInds[2];
						size_t iv3 = locVrtInds[3];

						if(    ( m_aaVrtVecVol[sv] )[iv0]
							&& ( m_aaVrtVecVol[sv] )[iv1]
							&& ( m_aaVrtVecVol[sv] )[iv2]
							&& ( m_aaVrtVecVol[sv] )[iv3]
						)
						{
							//	create a new prism
							expVol = *m_grid.create<Hexahedron>(
												HexahedronDescriptor(
															sv->vertex(iv3), sv->vertex(iv2),
															sv->vertex(iv1), sv->vertex(iv0),
															(m_aaVrtVecVol[sv])[iv3],
															(m_aaVrtVecVol[sv])[iv2],
															(m_aaVrtVecVol[sv])[iv1],
															(m_aaVrtVecVol[sv])[iv0]
															)
														);

//							m_sh.assign_subset(expVol, m_fracInfosBySubset.at(m_sh.get_subset_index(tFace)).newSubsetIndex);
//
//							return true;
						}
						else if(    ( m_aaVrtVecVol[sv] )[iv0]
								 && ( m_aaVrtVecVol[sv] )[iv1]

						)
						{
							//	create a new prism
							//	create a new prism
//							expVol = *m_grid.create<Prism>(
//											PrismDescriptor(sv->vertex(iv3),sv->vertex(iv2), sv->vertex(iv1),
//															sv->vertex(iv0),
//															(m_aaVrtVecVol[sv])[iv1],
//															(m_aaVrtVecVol[sv])[iv0])
//															);
							//	create a new Prism
							///	only used to initialize a prism. for all other tasks you should use VolumeDescripor.
							/**
							 * please be sure to pass the vertices in the correct order:
							 * v1, v2, v3: bottom-vertices in counterclockwise order (if viewed from the top).
							 * v4, v5, v6: top-vertices in counterclockwise order (if viewed from the top).
							 * 		PrismDescriptor(Vertex* v1, Vertex* v2, Vertex* v3,
						Vertex* v4, Vertex* v5, Vertex* v6);
							 *
							 */
							expVol = *m_grid.create<Prism>(
											PrismDescriptor( (m_aaVrtVecVol[sv])[iv0],
													sv->vertex(iv0), sv->vertex(iv3),
													(m_aaVrtVecVol[sv])[iv1], sv->vertex(iv1), sv->vertex(iv2)
												)
												);

							UG_LOG("PRISM 0 1 " << std::endl);
						}
						else if(    ( m_aaVrtVecVol[sv] )[iv0]
								 && ( m_aaVrtVecVol[sv] )[iv2]
						)
						{
							UG_LOG("PRISM 0 2 " << std::endl);
							expVol = *m_grid.create<Prism>(
											PrismDescriptor( (m_aaVrtVecVol[sv])[iv0],
													sv->vertex(iv0), sv->vertex(iv3),
													sv->vertex(iv2), (m_aaVrtVecVol[sv])[iv2], sv->vertex(iv1)
												)
												);

//							m_sh.assign_subset(expVol, m_sh.num_subsets());
//							m_sh.assign_subset( sv->vertex(iv0), m_sh.num_subsets());
//							m_sh.assign_subset( sv->vertex(iv1), m_sh.num_subsets());
//							m_sh.assign_subset( sv->vertex(iv2), m_sh.num_subsets());
//							m_sh.assign_subset( sv->vertex(iv3), m_sh.num_subsets());
//							m_sh.assign_subset( (m_aaVrtVecVol[sv])[iv0], m_sh.num_subsets());
//							m_sh.assign_subset( (m_aaVrtVecVol[sv])[iv2], m_sh.num_subsets());


						}
						else if(    ( m_aaVrtVecVol[sv] )[iv0]
								 && ( m_aaVrtVecVol[sv] )[iv3]
						)
						{
							UG_LOG("PRISM 0 3 " << std::endl);

						}
						else if(    ( m_aaVrtVecVol[sv] )[iv1]
								 && ( m_aaVrtVecVol[sv] )[iv2]
						)
						{

							UG_LOG("PRISM 1 2 " << std::endl);

						}
						else if(    ( m_aaVrtVecVol[sv] )[iv2]
								 && ( m_aaVrtVecVol[sv] )[iv3]
						)
						{
							UG_LOG("PRISM 2 3 " << std::endl);

						}



					}
					else
					{
						//	traditionally only tetrahedrons are supported. This section thus raises an error
						// Markus tries to implement also Hexahedra
//							grid.detach_from_vertices(aVrtVec);
//							grid.detach_from_volumes(aVrtVec);
//							grid.detach_from_vertices(aAdjMarker);
//							grid.detach_from_edges(aAdjMarker);
						detachMarkers();
						throw(UGError("Incomplete implementation error in ExpandFractures3d Arte: Only tetrahedrons are supported in the current implementation, and hexahedra are in development."));
						return false;
					}

					if(expVol)
					{

						IndexType newSubs = m_fracInfosBySubset.at(m_sh.get_subset_index(tFace)).newSubsetIndex;

						subsOfNewVolumes.push_back( newSubs );

//						m_sh.assign_subset(expVol, m_fracInfosBySubset.at(m_sh.get_subset_index(tFace)).newSubsetIndex);
						m_sh.assign_subset(expVol, newSubs);

						newFractureVolumes.push_back(expVol);
					}
				}
			}
		}

		//	now set up a new volume descriptor and replace the volume.
		if(vd.num_vertices() != sv->num_vertices())
			vd.set_num_vertices(sv->num_vertices());

		for(size_t i_vrt = 0; i_vrt < sv->num_vertices(); ++i_vrt)
		{
			if( (m_aaVrtVecVol[sv])[i_vrt] )
				vd.set_vertex(i_vrt, (m_aaVrtVecVol[sv])[i_vrt]);
			else
				vd.set_vertex(i_vrt, sv->vertex(i_vrt));
		}

		m_grid.create_by_cloning(sv, vd, sv);
		m_grid.erase(sv);
	}

	UG_LOG("Volumes erzeugt " << std::endl);

//	return false;




	//	we have to clean up unused faces and edges.
	//	note that all selected edges with mark 0 may safley be deleted. - warum?
	for(EdgeIterator iter = m_sel.begin<Edge>(); iter!= m_sel.end<Edge>();)
	{
		//	take care of the iterator
		Edge* e = *iter;
		++iter;

		if( m_aaMarkEdgeVFP[e].getNumberFracEdgesInVertex() == 0 )
			m_grid.erase(e);
	}

	//	make sure that no unused faces linger around (This should never happen!)
	bool foundUnusedFaces = false;

	for(FaceIterator iter = m_sel.begin<Face>(); iter != m_sel.end<Face>();)
	{
		Face* f = *iter;
		++iter;

//		if( m_aaMarkFaceHasUnclosedFracSideB[f] )
//		{
//			UG_LOG("want to delete unclosed frac face " << std::endl);
//			// todo fixme XXXXXXXXXXXXXXXXXX
//			//return false;
//		}

		if( ! m_aaMarkFaceIsFracB[f] )
		{
			foundUnusedFaces = true;
			m_grid.erase(f);
		}
	}

	if(foundUnusedFaces)
	{
		UG_LOG("WARNING in ExpandFractures3D Arte: Unused faces encountered during cleanup. Removing them...\n");
	}

	if( subsOfNewVolumes.size() != newFractureVolumes.size() )
	{
		UG_THROW("andere zahl neue volumes als subdoms " << std::endl);
	}


//	return true;


	IndexType nfn = 0;

	for( auto const & nf : newFractureVolumes )
	{
		for(size_t iFace = 0; iFace < nf->num_faces(); ++iFace)
		{
			Face * fac = m_grid.get_face(nf, iFace);

			m_sh.assign_subset( fac, subsOfNewVolumes[nfn] );

		}

		for(size_t iEdge = 0; iEdge < nf->num_edges(); ++iEdge)
		{
			Edge* edg = m_grid.get_edge(nf, iEdge);

			m_sh.assign_subset( edg, subsOfNewVolumes[nfn] );

		}

		for( size_t iVrt = 0; iVrt < nf->num_vertices(); iVrt++ )
		{
			Vertex * vrt = nf->vertex(iVrt);

			m_sh.assign_subset( vrt, subsOfNewVolumes[nfn] );
		}

		nfn++;
	}

	for( EndingCrossingFractureSegmentInfo const & ecfsi : m_vecEndCrossFractSegmInfo )
	{
		Face * hiddenCutFracFace = ecfsi.spuckHiddenCutFractManifEl();

		IndexType subsECC = ecfsi.spuckSudoFractEnding();

		m_sh.assign_subset( hiddenCutFracFace, subsECC );

		for(size_t iEdge = 0; iEdge < hiddenCutFracFace->num_edges(); ++iEdge)
		{
			Edge* edg = m_grid.get_edge(hiddenCutFracFace, iEdge);

			m_sh.assign_subset( edg, subsECC );

		}

		for( size_t iVrt = 0; iVrt < hiddenCutFracFace->num_vertices(); iVrt++ )
		{
			Vertex * vrt = hiddenCutFracFace->vertex(iVrt);

			m_sh.assign_subset( vrt, subsECC );
		}

	}


	return true;
}

ArteExpandFracs3D::IndexType ArteExpandFracs3D::deleteEndingCrossingCleftOrigFacs()
{
	IndexType numOfDelSegs = 0;

	for( EndingCrossingFractureSegmentInfo const & ecfsi : m_vecEndCrossFractSegmInfo )
	{

		Edge * directionEdge = ecfsi.spuckLowdimElmShiftDirection();

		Face * oldFractFacCut = ecfsi.spuckEndingFractManifCutting();

		Face * oldFractFacNotCut = ecfsi.spuckEndingFractManifNotCutting();

		Edge * commonEdge = nullptr;

		if( oldFractFacNotCut )
		{
			for(size_t iEdge = 0; iEdge < oldFractFacCut->num_edges(); ++iEdge)
			{
				Edge* edgC = m_grid.get_edge(oldFractFacCut, iEdge);

				for(size_t iEdge = 0; iEdge < oldFractFacNotCut->num_edges(); ++iEdge)
				{
					Edge* edgN = m_grid.get_edge(oldFractFacNotCut, iEdge);

					if( edgC == edgN )
					{
						commonEdge = edgC;
						break;
					}

				}

				if( commonEdge )
					break;

			}
		}

		if( oldFractFacCut )
			m_grid.erase(oldFractFacCut);

		if( oldFractFacNotCut )
			m_grid.erase(oldFractFacNotCut);

		if( directionEdge )
			m_grid.erase(directionEdge);

		if( commonEdge )
			m_grid.erase(commonEdge);

		numOfDelSegs++;
	}

	return numOfDelSegs;
}

////////////////////////////////////////////////////////////

bool ArteExpandFracs3D::etablishVolumesAtEndingCrossingClefts( std::vector<Volume*> & newFractureVolumes, std::vector<IndexType> & subsOfNewVolumes )
{

	for( EndingCrossingFractureSegmentInfo & ecfsi : m_vecEndCrossFractSegmInfo )
	{
		Vertex * baseVrtx = ecfsi.spuckUnclosedVrtx();

		Vertex * shiftVrtx = ecfsi.spuckShiftVrtx();

		Edge * cutEdge = ecfsi.spuckOldLowDimElCut();

		Vertex * secondVrtxCutEdge = nullptr;

		if( cutEdge->vertex(0) == baseVrtx )
		{
			secondVrtxCutEdge = cutEdge->vertex(1);
		}
		else if( cutEdge->vertex(1) == baseVrtx )
		{
			secondVrtxCutEdge = cutEdge->vertex(0);
		}
		else
		{
			UG_LOG("no second vertex of cut edge " << std::endl);
			UG_THROW("no second vertex of cut edge " << std::endl);
			return false;
		}

		Edge * divisionEdge = *m_grid.create<RegularEdge>( EdgeDescriptor( secondVrtxCutEdge, shiftVrtx ) );

		IndexType subsetECC = ecfsi.spuckSudoFractEnding();

//		m_sh.assign_subset(divisionEdge, m_sh.num_subsets());
		m_sh.assign_subset(divisionEdge, subsetECC);

		Face * hiddenCutFracFace = *m_grid.create<Triangle>(TriangleDescriptor( baseVrtx, shiftVrtx, secondVrtxCutEdge ));

//		m_sh.assign_subset(hiddenCutFracFace, m_sh.num_subsets());
		m_sh.assign_subset(hiddenCutFracFace, subsetECC);

		IndexType subsNewFacesEdges = m_sh.get_subset_index(hiddenCutFracFace);

		for(size_t iEdge = 0; iEdge < hiddenCutFracFace->num_edges(); ++iEdge)
		{
			Edge* edg = m_grid.get_edge(hiddenCutFracFace, iEdge);

			m_sh.assign_subset( edg, subsNewFacesEdges );
			UG_LOG("HEDGE CENTER " << CalculateCenter( edg, m_aaPos ) << std::endl );

			UG_LOG("HEdge subdom " << m_sh.get_subset_index(edg) << std::endl );

		}

		if( ! ecfsi.schluckHiddenCutFractManifEl( hiddenCutFracFace ))
		{
			UG_LOG("hidden cut face Problem" << std::endl);
			UG_THROW("hidden cut face Problem" << std::endl);
		}

		Face * endingFractFacCutting = ecfsi.spuckEndingFractManifCutting();

		IndexType const triangVrtxNum = 3;

		if( endingFractFacCutting->num_vertices() != triangVrtxNum )
		{
			UG_LOG("only triangles allowed " << std::endl);
			UG_THROW("only triangles allowed " << std::endl);
		}


		IndexType foundNotTouchVrtx = 0;

		Vertex * notTouchingVrtx = nullptr;


		// figure out that vertex of the ending fracture faces that is not touching the crossing not ending cleft
		for( IndexType vrtIndx = 0; vrtIndx < triangVrtxNum; vrtIndx++ )
		{
			Vertex * testVrt = endingFractFacCutting->vertex(vrtIndx);

			if( ! EdgeContains( cutEdge, testVrt ) )
			{
				foundNotTouchVrtx++;
				notTouchingVrtx = testVrt;
			}
		}

		if( foundNotTouchVrtx != 1 || ! notTouchingVrtx  )
		{
			UG_LOG("not touching vertex not found " << std::endl);
			UG_THROW("not touching vertex not found " << std::endl);
		}

		// replace the face that touches with an edge
		Face * replaceEndingFractCutFac = *m_grid.create<Triangle>(TriangleDescriptor( shiftVrtx, secondVrtxCutEdge, notTouchingVrtx ));

//		m_sh.assign_subset( replaceEndingFractCutFac, m_sh.num_subsets() );
		m_sh.assign_subset( replaceEndingFractCutFac, subsetECC );

		IndexType subsNewFacesEdgesC = m_sh.get_subset_index(replaceEndingFractCutFac);

		UG_LOG("EDGE NUMBER CCS CC " << replaceEndingFractCutFac->num_edges() << std::endl);

		for(size_t iEdge = 0; iEdge < replaceEndingFractCutFac->num_edges(); ++iEdge)
		{
			Edge* edg = m_grid.get_edge(replaceEndingFractCutFac, iEdge);

			m_sh.assign_subset( edg, subsNewFacesEdgesC );
			UG_LOG("EDGE CENTER " << CalculateCenter( edg, m_aaPos ) << std::endl );

			UG_LOG("Edge subdom " << m_sh.get_subset_index(edg) << std::endl );

		}

		// replace the face that has only the base vertex common, if existing
		Face * endingFractFacNotCutting = ecfsi.spuckEndingFractManifNotCutting();

		if( endingFractFacNotCutting )
		{
			if( endingFractFacNotCutting->num_vertices() != triangVrtxNum )
			{
				UG_LOG("only triangles allowed NC" << std::endl);
				UG_THROW("only triangles allowed NC" << std::endl);
			}

			std::vector<Vertex*> vrtcsNotBase;

			// figure out those vertices which are not the basis vertex
			for( IndexType vrtIndx = 0; vrtIndx < triangVrtxNum; vrtIndx++ )
			{
				Vertex * testVrt = endingFractFacNotCutting->vertex(vrtIndx);

				if( testVrt != baseVrtx )
					vrtcsNotBase.push_back(testVrt);
			}

			if( vrtcsNotBase.size() != 2 )
			{
				UG_LOG("strange number vertices " << vrtcsNotBase.size() << std::endl);
				UG_THROW("strange number vertices " << vrtcsNotBase.size() << std::endl);
			}

			Face * replaceEndingFractNotCutFac = *m_grid.create<Triangle>(TriangleDescriptor( shiftVrtx, vrtcsNotBase[0], vrtcsNotBase[1] ));

//			m_sh.assign_subset( replaceEndingFractNotCutFac, m_sh.num_subsets() );
			m_sh.assign_subset( replaceEndingFractNotCutFac, subsetECC );

			IndexType subsNewFacesEdges = m_sh.get_subset_index(replaceEndingFractNotCutFac);

			UG_LOG("EDGE NUMBER CCS " << replaceEndingFractNotCutFac->num_edges() << std::endl);

			for(size_t iEdge = 0; iEdge < replaceEndingFractNotCutFac->num_edges(); ++iEdge)
			{
				Edge* edg = m_grid.get_edge(replaceEndingFractNotCutFac, iEdge);

				m_sh.assign_subset( edg, subsNewFacesEdges );

				UG_LOG("EDGE CENTER " << CalculateCenter( edg, m_aaPos ) << std::endl );

				UG_LOG("Edge subdom " << m_sh.get_subset_index(edg) << std::endl );

			}


		}

//		for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end();)
		for( Volume * const & sv : ecfsi.spuckVecFulldimEl() )
		{

			UG_LOG("CEC Volume new creation try at " << CalculateCenter(sv, m_aaPos) << std::endl);

			bool volHasEndingCrossingCleftFace = m_volAttAccsVolTouchesEndingCrossingCleft[sv];

			// should be true!

			if( ! volHasEndingCrossingCleftFace )
			{
				UG_LOG("Was ist da los, ending cleft und doch nicht " << std::endl);
				UG_THROW("Was ist da los, ending cleft und doch nicht " << std::endl);
			}

			// check if volume contains the base vertex

			if( ! VolumeContains(sv, baseVrtx))
			{
				UG_LOG("VOlume ausdehenen ECC ohne base vertex " << std::endl);
				UG_THROW("VOlume ausdehenen ECC ohne base vertex " << std::endl);
			}

			//	now expand the fracture faces of sv to volumes.
			for(size_t i_side = 0; i_side < sv->num_sides(); ++i_side)
			{
				//	holds local side vertex indices
				std::vector<size_t>	locVrtInds;

				//	get the local vertex indices of the side of the volume
				sv->get_vertex_indices_of_face(locVrtInds, i_side);

				Face* tFace = m_grid.get_side(sv, i_side);

				if(tFace)
				{
					if( m_aaMarkFaceIsFracB[tFace] )
					{
						bool faceIsSegmLimEndCrossCleft = m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft[tFace];
						bool faceIsEndingCleftCrossFace = m_aaMarkFaceWithEndingCrossingCleft[tFace];

						UG_LOG("Volumenerzeugung Versuch ECC Anfang" << std::endl);

						Volume* expVol = nullptr;

						Volume* expVolTwo = nullptr;


						if(locVrtInds.size() == 3)
						{
							size_t iv0 = locVrtInds[0];
							size_t iv1 = locVrtInds[1];
							size_t iv2 = locVrtInds[2];

//							Vertex * vrtxOne = ( m_aaVrtVecVol[sv] )[iv0];
//							Vertex * vrtxTwo = ( m_aaVrtVecVol[sv] )[iv1];
//							Vertex * vrtxThree = ( m_aaVrtVecVol[sv] )[iv2];

							// figure out which vertex is the base vertex

							int indBasVrtx = -1;

							// wie kann das nur stimmen, das Volumen hat doch mehr als 3 Vertizes....
							// zufällig für das Testbeispiel vielleicht richtig?
//							for( IndexType vrtxInd = 0; vrtxInd < triangVrtxNum; vrtxInd++ )
							for( IndexType vrtxInd = 0; vrtxInd < sv->num_vertices(); vrtxInd++ )
							{
								if( sv->vertex(vrtxInd) == baseVrtx )
									indBasVrtx = vrtxInd;
							}
							// muss der Loop nicht über die Vertizes des volumens gehen?????

							// der Index des Basisvertex ist tabu für die Ausehnung der endenden fracture

							if( endingFractFacCutting == tFace || endingFractFacNotCutting == tFace )
							{
								if( ! faceIsEndingCleftCrossFace )
								{
									UG_LOG("Widerspruch ending but not ending ECC" << std::endl);
								}

								if( iv0 == indBasVrtx )
								{
									if(    ( m_aaVrtVecVol[sv] )[iv0]
										&& ( m_aaVrtVecVol[sv] )[iv1]
										&& ( m_aaVrtVecVol[sv] )[iv2]
									)
									{
										//	create a new pyramid
//										expVol = *m_grid.create<Pyramid>(
//														PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv1),
//																		(m_aaVrtVecVol[sv])[iv2],
//																		(m_aaVrtVecVol[sv])[iv1],
//																		(m_aaVrtVecVol[sv])[iv0]));

										expVol = *m_grid.create<Pyramid>(
														PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
															(m_aaVrtVecVol[sv])[iv2],
															(m_aaVrtVecVol[sv])[iv1],
															(m_aaVrtVecVol[sv])[iv0]));

									}
									else if(    ( m_aaVrtVecVol[sv] )[iv0]
											 && ( m_aaVrtVecVol[sv] )[iv1]
									)
									{
										//	create a new tetrahedron
										expVol = *m_grid.create<Tetrahedron>(
														TetrahedronDescriptor(sv->vertex(iv1),
															(m_aaVrtVecVol[sv])[iv1],
															(m_aaVrtVecVol[sv])[iv0],
															sv->vertex(iv2)));
									}
//									else if(    ( m_aaVrtVecVol[sv] )[iv1]
//											 && ( m_aaVrtVecVol[sv] )[iv2]
//									)
//									{
//										//	create a new Pyramid
//										expVol = *m_grid.create<Pyramid>(
//														PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
//															(m_aaVrtVecVol[sv])[iv2],
//															(m_aaVrtVecVol[sv])[iv1],
//															sv->vertex(iv0)));
//									}
									else if(    (m_aaVrtVecVol[sv])[iv0]
											 && (m_aaVrtVecVol[sv])[iv2]
									)
									{
										//	create a new tetrahedron
										expVol = *m_grid.create<Tetrahedron>(
														TetrahedronDescriptor(sv->vertex(iv2),
															(m_aaVrtVecVol[sv])[iv0],
															(m_aaVrtVecVol[sv])[iv2],
															sv->vertex(iv1)));
									}
//									else if( ( m_aaVrtVecVol[sv])[iv0] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv0]));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv1] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv1]));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv2] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv2]));
//									}
									else
									{
										//	Text from SR, still similar:
										//  this code-block should never be entered. If it is entered then
										//	we either selected the wrong faces (this shouldn't happen), or there
										//	are selected faces, which have fracture-boundary-vertices only.
										//	This is the same is if inner fracture edges exists, which are
										//	connected to two boundary vertices.
										//	Since we tried to remove those edges above, something went wrong.
										//	remove the temporary attachments and throw an error

										UG_LOG("Tetraeder ECC Fehlt eine Loesung " << std::endl);
										detachMarkers();
										throw(UGError("Error in ExpandFractures3d Arte Stasi. Implementation Error."));
										return false;
									}

								}
								else if( iv1 == indBasVrtx )
								{
									if(    ( m_aaVrtVecVol[sv] )[iv0]
										&& ( m_aaVrtVecVol[sv] )[iv1]
										&& ( m_aaVrtVecVol[sv] )[iv2]
									)
									{
										//	create a new prism
//										expVol = *m_grid.create<Pyramid>(
//														PyramidDescriptor( sv->vertex(iv2), sv->vertex(iv0),
//																				(m_aaVrtVecVol[sv])[iv2],
//																				(m_aaVrtVecVol[sv])[iv1],
//																				(m_aaVrtVecVol[sv])[iv0]));


										expVol = *m_grid.create<Pyramid>(
														PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
															(m_aaVrtVecVol[sv])[iv0],
															(m_aaVrtVecVol[sv])[iv2],
															(m_aaVrtVecVol[sv])[iv1]));

									}
									else if(    ( m_aaVrtVecVol[sv] )[iv0]
											 && ( m_aaVrtVecVol[sv] )[iv1]
									)
									{
										//	create a new Pyramid
										expVol = *m_grid.create<Tetrahedron>(
														TetrahedronDescriptor(sv->vertex(iv0),
															(m_aaVrtVecVol[sv])[iv1],
															(m_aaVrtVecVol[sv])[iv0],
															sv->vertex(iv2)));
									}
									else if(    ( m_aaVrtVecVol[sv] )[iv1]
											 && ( m_aaVrtVecVol[sv] )[iv2]
									)
									{
										//	create a new Pyramid
										expVol = *m_grid.create<Tetrahedron>(
														TetrahedronDescriptor(sv->vertex(iv2),
															(m_aaVrtVecVol[sv])[iv2],
															(m_aaVrtVecVol[sv])[iv1],
															sv->vertex(iv0)));
									}
//									else if(    (m_aaVrtVecVol[sv])[iv0]
//											 && (m_aaVrtVecVol[sv])[iv2]
//									)
//									{
//										//	create a new Pyramid
//										expVol = *m_grid.create<Pyramid>(
//														PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
//															(m_aaVrtVecVol[sv])[iv0],
//															(m_aaVrtVecVol[sv])[iv2],
//															sv->vertex(iv1)));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv0] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv0]));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv1] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv1]));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv2] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv2]));
//									}
									else
									{
										//	Text from SR, still similar:
										//  this code-block should never be entered. If it is entered then
										//	we either selected the wrong faces (this shouldn't happen), or there
										//	are selected faces, which have fracture-boundary-vertices only.
										//	This is the same is if inner fracture edges exists, which are
										//	connected to two boundary vertices.
										//	Since we tried to remove those edges above, something went wrong.
										//	remove the temporary attachments and throw an error

										UG_LOG("Tetraeder ECC Fehlt eine Loesung " << std::endl);
										detachMarkers();
										throw(UGError("Error in ExpandFractures3d Arte Stasi. Implementation Error."));
										return false;
									}
								}
								else if( iv2 == indBasVrtx )
								{
									if(    ( m_aaVrtVecVol[sv] )[iv0]
										&& ( m_aaVrtVecVol[sv] )[iv1]
										&& ( m_aaVrtVecVol[sv] )[iv2]
									)
									{
										//	create a new prism
//										expVol = *m_grid.create<Pyramid>(
//														PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv0),
//																		(m_aaVrtVecVol[sv])[iv2],
//																		(m_aaVrtVecVol[sv])[iv1],
//																		(m_aaVrtVecVol[sv])[iv0]));

										expVol = *m_grid.create<Pyramid>(
														PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
															(m_aaVrtVecVol[sv])[iv1],
															(m_aaVrtVecVol[sv])[iv0],
															(m_aaVrtVecVol[sv])[iv2]));

									}
//									else if(    ( m_aaVrtVecVol[sv] )[iv0]
//											 && ( m_aaVrtVecVol[sv] )[iv1]
//									)
//									{
//										//	create a new Pyramid
//										expVol = *m_grid.create<Pyramid>(
//														PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
//															(m_aaVrtVecVol[sv])[iv1],
//															(m_aaVrtVecVol[sv])[iv0],
//															sv->vertex(iv2)));
//									}
									else if(    ( m_aaVrtVecVol[sv] )[iv1]
											 && ( m_aaVrtVecVol[sv] )[iv2]
									)
									{
										//	create a new Pyramid
										expVol = *m_grid.create<Tetrahedron>(
														TetrahedronDescriptor(sv->vertex(iv1),
															(m_aaVrtVecVol[sv])[iv2],
															(m_aaVrtVecVol[sv])[iv1],
															sv->vertex(iv0)));
									}
									else if(    (m_aaVrtVecVol[sv])[iv0]
											 && (m_aaVrtVecVol[sv])[iv2]
									)
									{
										//	create a new Pyramid
										expVol = *m_grid.create<Tetrahedron>(
														TetrahedronDescriptor(sv->vertex(iv0),
															(m_aaVrtVecVol[sv])[iv0],
															(m_aaVrtVecVol[sv])[iv2],
															sv->vertex(iv1)));
									}
//									else if( ( m_aaVrtVecVol[sv])[iv0] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv0]));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv1] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv1]));
//									}
//									else if( ( m_aaVrtVecVol[sv])[iv2] )
//									{
//										//	create a new Tetrahedron
//										expVol = *m_grid.create<Tetrahedron>(
//														TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
//																			 (m_aaVrtVecVol[sv])[iv2]));
//									}
									else
									{
										//	Text from SR, still similar:
										//  this code-block should never be entered. If it is entered then
										//	we either selected the wrong faces (this shouldn't happen), or there
										//	are selected faces, which have fracture-boundary-vertices only.
										//	This is the same is if inner fracture edges exists, which are
										//	connected to two boundary vertices.
										//	Since we tried to remove those edges above, something went wrong.
										//	remove the temporary attachments and throw an error

										UG_LOG("Tetraeder ECC Fehlt eine Loesung " << std::endl);
										detachMarkers();
										throw(UGError("Error in ExpandFractures3d Arte Stasi. Implementation Error."));
										return false;
									}

								}
							}
							else
							{
								// test if face belongs to the closed faces

								std::pair<Face*,Face*> const & closedNeighbrs = ecfsi.spuckPairNeighbouredFractClosedManifEl();


								if( closedNeighbrs.first == tFace || closedNeighbrs.second == tFace )
								{

									// figure out the vertex that is not from the cut edge

									Vertex * freeVrtx = nullptr;

									// wasn das für ein Käse.....
									for( IndexType vrtxInd = 0; vrtxInd < triangVrtxNum; vrtxInd++ )
//									for( IndexType vrtxInd = 0; vrtxInd < sv->num_vertices(); vrtxInd++ )
									{
										if( ! EdgeContains(cutEdge, tFace->vertex(vrtxInd)))
											freeVrtx = tFace->vertex(vrtxInd);
									}

									if( freeVrtx == nullptr )
									{
										UG_LOG("Immer noch null " << std::endl);
										UG_THROW("Immer noch null " << std::endl);
									}

//									int freeVrtxInd = -1;
//									int secondCutEdgVrtxInd = -1;
//
//									for( IndexType vrtxInd = 0; vrtxInd < sv->num_vertices(); vrtxInd++ )
//									{
//										if( ! EdgeContains(cutEdge, freeVrtx))
//											freeVrtxInd = vrtxInd;
//
//										if( secondVrtxCutEdge == secondVrtxCutEdge)
//											secondCutEdgVrtxInd = vrtxInd;
//									}


									// for the "pyramid", exclude the base Vertex
									// for the other tetrahedron, exclude the shifted of the non-base of the cut edge
									// i.e. of secondVrtxCutEdge

//									Vertex * vrtxNotFromCutEdge = tFace->vertex(freeVrtxInd);

									// relate iv0, iv1, iv2 to vrtxNotFromCutEdge, base

									if( iv0 == indBasVrtx )
									{

										if(    ( m_aaVrtVecVol[sv] )[iv0]
											&& ( m_aaVrtVecVol[sv] )[iv1]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{
											expVol = *m_grid.create<Pyramid>(
															PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
																(m_aaVrtVecVol[sv])[iv2],
																(m_aaVrtVecVol[sv])[iv1],
																( m_aaVrtVecVol[sv] )[iv0]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv0]));

										}
										else if(    ( m_aaVrtVecVol[sv] )[iv0]
												 && ( m_aaVrtVecVol[sv] )[iv1]
										)
										{
											expVol = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), ( m_aaVrtVecVol[sv] )[iv0],
																				 (m_aaVrtVecVol[sv])[iv1]));


											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv0]));
										}
//										else if(    ( m_aaVrtVecVol[sv] )[iv1]
//											&& ( m_aaVrtVecVol[sv] )[iv2]
//										)
//										{
//
//										}
										else if(    ( m_aaVrtVecVol[sv] )[iv0]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{

											expVol = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1),  ( m_aaVrtVecVol[sv] )[iv0],
																				 (m_aaVrtVecVol[sv])[iv2]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv0]));
										} //
										else if(    ( m_aaVrtVecVol[sv] )[iv0]
										)
										{

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv0]));
										} //



									}

									if( iv1 == indBasVrtx )
									{

										if(    ( m_aaVrtVecVol[sv] )[iv0]
											&& ( m_aaVrtVecVol[sv] )[iv1]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{
											expVol = *m_grid.create<Pyramid>(
															PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
																(m_aaVrtVecVol[sv])[iv0],
																(m_aaVrtVecVol[sv])[iv2],
																( m_aaVrtVecVol[sv] )[iv1]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv1]));

										}
										else if(    ( m_aaVrtVecVol[sv] )[iv0]
											&& ( m_aaVrtVecVol[sv] )[iv1]
										)
										{
											expVol = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), ( m_aaVrtVecVol[sv] )[iv1], sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv0]));


											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv1]));

										}
										else if(    ( m_aaVrtVecVol[sv] )[iv1]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{
											expVol = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), ( m_aaVrtVecVol[sv] )[iv1], sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv2]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv1]));

										}
//										else if(    ( m_aaVrtVecVol[sv] )[iv0]
//											&& ( m_aaVrtVecVol[sv] )[iv2]
//										)
//										{
//
//										}
										else if(    ( m_aaVrtVecVol[sv] )[iv1]
										)
										{

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv1]));

										}


									}

									if( iv2 == indBasVrtx )
									{

										if(    ( m_aaVrtVecVol[sv] )[iv0]
											&& ( m_aaVrtVecVol[sv] )[iv1]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{

											expVol = *m_grid.create<Pyramid>(
															PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
																(m_aaVrtVecVol[sv])[iv1],
																(m_aaVrtVecVol[sv])[iv0],
																( m_aaVrtVecVol[sv] )[iv2]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv2]));
										}

	//									else if(    ( m_aaVrtVecVol[sv] )[iv0]
	//										&& ( m_aaVrtVecVol[sv] )[iv1]
	//									)
	//									{
	//
	//									}
										else if(    ( m_aaVrtVecVol[sv] )[iv1]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{

											expVol = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(( m_aaVrtVecVol[sv] )[iv2], sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv1]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor( sv->vertex(iv1) , sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv2]));

										}
										else if(    ( m_aaVrtVecVol[sv] )[iv0]
											&& ( m_aaVrtVecVol[sv] )[iv2]
										)
										{
											expVol = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), ( m_aaVrtVecVol[sv] )[iv0],
																				 (m_aaVrtVecVol[sv])[iv2]));

											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv2]));
										} //
										else if( ( m_aaVrtVecVol[sv] )[iv2]
										)
										{
											expVolTwo = *m_grid.create<Tetrahedron>(
															TetrahedronDescriptor( sv->vertex(iv1) , sv->vertex(iv1), sv->vertex(iv0),
																				 (m_aaVrtVecVol[sv])[iv2]));
										}



									}

								}


								std::vector<Face*> const & closedNotNeighbr = ecfsi.spuckVecClosedFracManifElNoNeighbr();

								for( Face * fac : closedNotNeighbr )
								{
									if( fac == tFace )
									{
										if( iv0 == indBasVrtx )
										{

											if(    ( m_aaVrtVecVol[sv] )[iv0]
												&& ( m_aaVrtVecVol[sv] )[iv1]
												&& ( m_aaVrtVecVol[sv] )[iv2]
											)
											{

												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv0]));

												expVolTwo = *m_grid.create<Pyramid>(
																PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
																	(m_aaVrtVecVol[sv])[iv2],
																	(m_aaVrtVecVol[sv])[iv1],
																	( m_aaVrtVecVol[sv] )[iv0]));

											}
											else if(    ( m_aaVrtVecVol[sv] )[iv0]
												&& ( m_aaVrtVecVol[sv] )[iv1]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv0]));

												expVolTwo = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), ( m_aaVrtVecVol[sv] )[iv0],
																					 (m_aaVrtVecVol[sv])[iv1]));

											}
//											else if(    ( m_aaVrtVecVol[sv] )[iv1]
//												&& ( m_aaVrtVecVol[sv] )[iv2]
//											)
//											{
//
//											}
//												&& ( m_aaVrtVecVol[sv] )[iv2]
//											)
//											{
//
//											}
											else if(    ( m_aaVrtVecVol[sv] )[iv0]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv0]));


											}




										}

										if( iv1 == indBasVrtx )
										{

											if(    ( m_aaVrtVecVol[sv] )[iv0]
												&& ( m_aaVrtVecVol[sv] )[iv1]
												&& ( m_aaVrtVecVol[sv] )[iv2]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv1]));

												expVolTwo = *m_grid.create<Pyramid>(
																PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
																	(m_aaVrtVecVol[sv])[iv0],
																	(m_aaVrtVecVol[sv])[iv2],
																	( m_aaVrtVecVol[sv] )[iv1]));

											}
											else if(    ( m_aaVrtVecVol[sv] )[iv0]
												&& ( m_aaVrtVecVol[sv] )[iv1]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv1]));

												expVolTwo = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), ( m_aaVrtVecVol[sv] )[iv1], sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv0]));

											}
											else if(    ( m_aaVrtVecVol[sv] )[iv1]
												&& ( m_aaVrtVecVol[sv] )[iv2]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv1]));

												expVolTwo = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), ( m_aaVrtVecVol[sv] )[iv1], sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv2]));

											}
//											else if(    ( m_aaVrtVecVol[sv] )[iv0]
//												&& ( m_aaVrtVecVol[sv] )[iv2]
//											)
//											{
//
//											}
											else if(    ( m_aaVrtVecVol[sv] )[iv1]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv1]));

											}


										}

										if( iv2 == indBasVrtx )
										{
											if(    ( m_aaVrtVecVol[sv] )[iv0]
												&& ( m_aaVrtVecVol[sv] )[iv1]
												&& ( m_aaVrtVecVol[sv] )[iv2]
											)
											{

												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv2]));

												expVolTwo = *m_grid.create<Pyramid>(
																PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
																	(m_aaVrtVecVol[sv])[iv1],
																	(m_aaVrtVecVol[sv])[iv0],
																	( m_aaVrtVecVol[sv] )[iv2]));

											}
	//										else if(    ( m_aaVrtVecVol[sv] )[iv0]
	//											&& ( m_aaVrtVecVol[sv] )[iv1]
	//										)
	//										{
	//
	//										}
											else if(    ( m_aaVrtVecVol[sv] )[iv1]
												&& ( m_aaVrtVecVol[sv] )[iv2]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv2]));

												expVolTwo = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor( ( m_aaVrtVecVol[sv] )[iv2], sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv1]));

											}
											else if(    ( m_aaVrtVecVol[sv] )[iv0]
												&& ( m_aaVrtVecVol[sv] )[iv2]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv2]));

												expVolTwo = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor( ( m_aaVrtVecVol[sv] )[iv2], sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv0]));

											}	//
											else if(   ( m_aaVrtVecVol[sv] )[iv2]
											)
											{
												expVol = *m_grid.create<Tetrahedron>(
																TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																					 (m_aaVrtVecVol[sv])[iv2]));


											}

										}

									}
								}

#if 0

								if(    ( m_aaVrtVecVol[sv] )[iv0]
									&& ( m_aaVrtVecVol[sv] )[iv1]
									&& ( m_aaVrtVecVol[sv] )[iv2]
								)
								{
									//	create a new prism
									expVol = *m_grid.create<Prism>(
													PrismDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																	(m_aaVrtVecVol[sv])[iv2],
																	(m_aaVrtVecVol[sv])[iv1],
																	(m_aaVrtVecVol[sv])[iv0]));
								}
								else if(    ( m_aaVrtVecVol[sv] )[iv0]
										 && ( m_aaVrtVecVol[sv] )[iv1]
								)
								{
									//	create a new Pyramid
									expVol = *m_grid.create<Pyramid>(
													PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
														(m_aaVrtVecVol[sv])[iv1],
														(m_aaVrtVecVol[sv])[iv0],
														sv->vertex(iv2)));
								}
								else if(    ( m_aaVrtVecVol[sv] )[iv1]
										 && ( m_aaVrtVecVol[sv] )[iv2]
								)
								{
									//	create a new Pyramid
									expVol = *m_grid.create<Pyramid>(
													PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
														(m_aaVrtVecVol[sv])[iv2],
														(m_aaVrtVecVol[sv])[iv1],
														sv->vertex(iv0)));
								}
								else if(    (m_aaVrtVecVol[sv])[iv0]
										 && (m_aaVrtVecVol[sv])[iv2]
								)
								{
									//	create a new Pyramid
									expVol = *m_grid.create<Pyramid>(
													PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
														(m_aaVrtVecVol[sv])[iv0],
														(m_aaVrtVecVol[sv])[iv2],
														sv->vertex(iv1)));
								}
								else if( ( m_aaVrtVecVol[sv])[iv0] )
								{
									//	create a new Tetrahedron
									expVol = *m_grid.create<Tetrahedron>(
													TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																		 (m_aaVrtVecVol[sv])[iv0]));
								}
								else if( ( m_aaVrtVecVol[sv])[iv1] )
								{
									//	create a new Tetrahedron
									expVol = *m_grid.create<Tetrahedron>(
													TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																		 (m_aaVrtVecVol[sv])[iv1]));
								}
								else if( ( m_aaVrtVecVol[sv])[iv2] )
								{
									//	create a new Tetrahedron
									expVol = *m_grid.create<Tetrahedron>(
													TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
																		 (m_aaVrtVecVol[sv])[iv2]));
								}
								else
								{
									//	Text from SR, still similar:
									//  this code-block should never be entered. If it is entered then
									//	we either selected the wrong faces (this shouldn't happen), or there
									//	are selected faces, which have fracture-boundary-vertices only.
									//	This is the same is if inner fracture edges exists, which are
									//	connected to two boundary vertices.
									//	Since we tried to remove those edges above, something went wrong.
									//	remove the temporary attachments and throw an error

									UG_LOG("Tetraeder Fehlt eine Loesung " << std::endl);
									detachMarkers();
									throw(UGError("Error in ExpandFractures3d Arte Stasi. Implementation Error."));
									return false;
								}


#endif

							}


						}
						else
						{
							//	traditionally only tetrahedrons are supported. This section thus raises an error
							// Markus tries to implement also Hexahedra
		//									grid.detach_from_vertices(aVrtVec);
		//								grid.detach_from_volumes(aVrtVec);
		//							grid.detach_from_vertices(aAdjMarker);
		//							grid.detach_from_edges(aAdjMarker);
							detachMarkers();
							throw(UGError("Incomplete implementation error in ExpandFractures3d Arte ending crossing clefts."));
							return false;
						}

						if(expVol)
						{

							IndexType newSubs = m_fracInfosBySubset.at(m_sh.get_subset_index(tFace)).newSubsetIndex;

							subsOfNewVolumes.push_back( newSubs );

							m_sh.assign_subset(expVol, newSubs);

							newFractureVolumes.push_back(expVol);
						}

						if(expVolTwo)
						{

							IndexType newSubs = m_fracInfosBySubset.at(m_sh.get_subset_index(tFace)).newSubsetIndex;

							subsOfNewVolumes.push_back( newSubs );

							m_sh.assign_subset(expVolTwo, newSubs);

							newFractureVolumes.push_back(expVolTwo);
						}

					}
				}
			}
		}
	}


	return true;
}


} /* namespace ug */
