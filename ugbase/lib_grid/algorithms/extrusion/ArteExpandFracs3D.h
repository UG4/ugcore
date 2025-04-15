/*
 * ArteExpandFracs3D.h
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
#include <stack>
#include <vector>
#include "lib_grid/lg_base.h"
#include "expand_layers.h"
	//#include "expand_layers_arte.h"
	//#include "expand_layers_arte3D.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"
//#include "lib_grid/util/simple_algebra/least_squares_solver.h"

#include <vector>

#include "support.h"
#include "support3D.h"



#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_

namespace ug {

class ArteExpandFracs3D
{

public:

	ArteExpandFracs3D( Grid & grid, SubsetHandler & sh,
		    std::vector<FractureInfo> const & fracInfos,
			bool useTrianglesInDiamonds, bool establishDiamonds );

	virtual ~ArteExpandFracs3D();


public:

	bool run( bool & needToRestart );

	using IndexType = unsigned short;


private:

	Grid & m_grid;
	SubsetHandler & m_sh;
	std::vector<FractureInfo> m_fracInfos;
	bool m_useTrianglesInDiamonds, m_establishDiamonds;

	Grid::VertexAttachmentAccessor<APosition> m_aaPos;

	//	objects for temporary results
	FaceDescriptor m_facDescr;
	VolumeDescriptor m_volDescr;

//	std::vector<Edge*> m_tmpEdges; // used for temporary results.
//	std::vector<Face*> m_tmpFaces; // used for temporary results.
//	std::vector<Volume*> m_tmpVols; // used for temporary results.

	std::vector<FractureInfo> m_fracInfosBySubset;

	Selector m_sel;

	bool initialize();

	bool setSelector();

	bool attachMarkers();

	bool detachMarkers();

	using NormalVectorFacIntoVol = vector3;

	using AttachedFractFaceEdgeSudo = support::AttachedFractElem<Face*,Edge*,IndexType,NormalVectorFacIntoVol>;

	using VecAttachedFractFaceEdgeSudo = std::vector<AttachedFractFaceEdgeSudo>;

	using AttachedGenerFaceEdgeSudo = support::AttachedGeneralElem<Face*,Edge*,IndexType>;

	using VecAttachedGenerFaceEdgeSudo = std::vector<AttachedGenerFaceEdgeSudo>;

	using AttachedBndryFaceEdgeSudo = support::AttachedBoundryElem<Face*,Edge*,IndexType,NormalVectorFacIntoVol>;

	using VecAttachedBndryFaceEdgeSudo = std::vector<AttachedBndryFaceEdgeSudo>;

	using EdgePair = std::pair<Edge*,Edge*>;

	using VertxFracPropts = support::VertexFracturePropertiesVol<IndexType, AttachedFractFaceEdgeSudo>;

	using AttVertFracProp = Attachment<VertxFracPropts>;

	AttVertFracProp m_aAdjMarkerVFP;

	// TODO FIXME verfehltes Konzept im 3D Fall!!! ERROR
	Grid::VertexAttachmentAccessor<AttVertFracProp> m_aaMarkVrtVFP;
	// TODO FIXME anstatt zu zählen, wieviele fractures angrenzen, was man ja lassen kann,
	// vielleicht irgendwo auch gebraucht, muss man vor allem zählen, wieviele subdomains
	// von fractures an dem Vertex zusammentreffen!!!!!
	// und ob an dem Vertex, wenn er innen ist, eine fracture aufhört, oder weiter geht
	// ob also der Vertex "rundum" von fracture faces umgeben ist, oder nur teilweise
	// das vertex fracture properties Konzept vom 2D Fall ist also nicht ausreichend

	//	AttVertFracProp m_aAdjMarkerVFP; // used to know if an edge is frac edge, suffix vfp misleading....

	Grid::EdgeAttachmentAccessor<AttVertFracProp> m_aaMarkEdgeVFP;
	// used to know if an edge is frac edge, suffix vfp misleading....

	ABool m_aAdjMarkerFaceIsFracB; // used to know if an face is frac face

	Grid::FaceAttachmentAccessor<ABool> m_aaMarkFaceIsFracB;

	// TODO FIXME die hier können alle entfernt werden im Prinzip
#if 0
	ABool m_aAdjMarkerFaceHasUnclosedFracSideB;

	Grid::FaceAttachmentAccessor<ABool> m_aaMarkFaceHasUnclosedFracSideB;

	ABool m_aAdjMarkerVrtxHasUnclosedFracB;

	Grid::VertexAttachmentAccessor<ABool> m_aaMarkVrtxHasUnclosedFracB;

	ABool m_aAdjMarkerVrtx2AtInnerEndOfEndingCrossingFract;

	Grid::VertexAttachmentAccessor<ABool> m_aaMarkVrtx2AtInnerEndOfEndingCrossingFract;

#endif

	ABool m_aAdjMarkerFaceWithEndingCrossingCleft;

	Grid::FaceAttachmentAccessor<ABool> m_aaMarkFaceWithEndingCrossingCleft;

	ABool m_aAdjMarkerVrtxAtEndingCrossingCleft;

	Grid::VertexAttachmentAccessor<ABool> m_aaMarkVrtxAtEndingCrossingCleft;



	bool countAndSelectFracBaseNums();

	bool distinguishSegments();

//	bool checkUnclosedFracFaces();

	bool shiftUnclosedFracFaces();

//	bool detectEndingCrossingClefts();

	bool detectEndingCrossingCleftsSegmBased();

	bool m_needToSplitEdgesConnectingNeighbrdEndingCrossCleftVrtx;

	bool seletForSegmented();

	std::vector<Face*> m_originalFractureFaces;

	bool assignOrigFracInfos();


//	using AttVrtVec = Attachment<std::vector<Vertex*> >;
//	AttVrtVec m_attVrtVec;
//	Grid::VolumeAttachmentAccessor<AttVrtVec> m_aaVrtVecVol;

	using AttVecEdge = Attachment<std::vector<Edge*>>;
	using AttVecFace = Attachment<std::vector<Face*>>;
//	using AttVecVol = Attachment<std::vector<Volume*>>;

	AttVecEdge m_aAdjInfoEdges;
	AttVecFace m_aAdjInfoFaces;
//	AttVecVol m_aAdjInfoVols;

	Grid::VertexAttachmentAccessor<AttVecEdge> m_aaVrtInfoAssoEdges;
	Grid::VertexAttachmentAccessor<AttVecFace> m_aaVrtInfoAssoFaces;
//	Grid::VertexAttachmentAccessor<AttVecVol> m_aaVrtInfoAssoVols;

	bool establishNewVrtBase();

//	bool generateVertexInfos();

	bool loop2EstablishNewVertices();

//	ug::support::VertexFractureQuadrupel<Face*,vector3,Volume*,Edge*,IndexType> bla();

//	using VrtxFractrQuadrplArte3D = support::VertexFractureQuadrupel<Face*,vector3,Volume*,Edge*,IndexType>;
//
//	using VrtxFractrQuadrplArte3DVec = std::vector<VrtxFractrQuadrplArte3D>;
//
//	VrtxFractrQuadrplArte3DVec m_vrtxFractrQuadrplVec;

//	using VertFracTrip = support::VertexFractureTripleMF<Face*, IndexType, Volume*, vector3, Edge*>;
//
//	using VecVertFracTrip = std::vector<VertFracTrip>;
//
//	using AttVecVertFracTrip = Attachment<VecVertFracTrip>;

//	AttVecVertFracTrip m_aAdjInfoAVVFT;
//
//	Grid::VertexAttachmentAccessor<AttVecVertFracTrip> m_aaVrtInfoFraTri;

	bool checkIfFacesVerticesCoincide( Face * const & facOne, Face * const & facTwo );

	bool collectFaceVertices( std::vector<Vertex*> & facVrt, Face * const & fac );

	bool createConditionForNewVrtcs();

	using AttVrtVec = Attachment<std::vector<Vertex*> >;

	AttVrtVec m_attVrtVec;

	Grid::VolumeAttachmentAccessor<AttVrtVec> m_aaVrtVecVol;

	using CrossVertInf = support::CrossingVertexInfoVol<Vertex*, IndexType >;

	std::vector<CrossVertInf> m_vecCrossVrtInf;

	using PairSudoBool = std::pair<IndexType,bool>;
	using VecPairSudoBool = std::vector<PairSudoBool>;

	// deprecated due to Stasi algo
//	bool isVrtxSurroundedByFracFaces( Vertex * const & vrt, VertxFracPropts & vrtxFracPrps );
//									  VecPairSudoBool & sudoSurrounded );

	// deprecated due to Stasi Algorithm
//	bool sortElemCircleIsClosed( VecAttachedFractFaceEdgeSudo const & vecAttFac,
//								 VecAttachedFractFaceEdgeSudo & vecSortedFac,
//								 int startFaceIndexUser = -1,
////								 int endFaceIndexUser = -1,
////								 IndexType startEdgeIndexUser = -1,
////								 IndexType endEdgeIndexUser = -1
////								 Face * const & startFacUser = nullptr,
////								 Face * const & endFacUser = nullptr,
//								 Edge * const & startEdgUser = nullptr,
//								 Edge * const & endEdgUser = nullptr
//								);

public:
	using VrtxFracProptsStatus = VertxFracPropts::VrtxFracStatus;

private:
//	static_assert< std::is_same< VrtxFracProptsStatus,support::VertexFracturePropertiesVol::VrtxFracStatus>::value );
//	static_assert( std::is_same<VrtxFracProptsStatus,support::VertexFracturePropertiesVol::VrtxFracStatus>::value );
//	static_assert( std::is_same<VrtxFracProptsStatus,support::VertexFracturePropertiesVol<IndexType, AttachedFractFaceEdgeSudo>::VrtxFracStatus>::value );

//	template<typename VOLUME_TYPE, VrtxFracProptsStatus vfps>
////	template<support::VertexFracturePropertiesVol::VrtxFracStatus vfp>
////	template<int I>
//	bool establishNewVertices( Vertex * const & oldVrt )
//	{
//		UG_THROW("too general, should not be called presently " << std::endl);
//		return false;
//	};

//	template< bool APPLY_GENERAL_SEGMENT_ORDERING,
//			  ArteExpandFracs3D::VrtxFracProptsStatus vfp,
//			  std::enable_if<APPLY_GENERAL_SEGMENT_ORDERING, bool>::type = true
//			>
//	bool establishNewVertices( Vertex * const & oldVrt );

//	template< bool APPLY_GENERAL_SEGMENT_ORDERING,
//			  ArteExpandFracs3D::VrtxFracProptsStatus vfp
//	>
//	bool establishNewVertices( Vertex * const & oldVrt );


//	template< bool APPLY_GENERAL_SEGMENT_ORDERING,
//			  ArteExpandFracs3D::VrtxFracProptsStatus vfp,
//			  std::enable_if<APPLY_GENERAL_SEGMENT_ORDERING, bool>::type = true
//			>
//	bool establishNewVertices( Vertex * const & oldVrt );
//
//	template< bool APPLY_GENERAL_SEGMENT_ORDERING,
//			  ArteExpandFracs3D::VrtxFracProptsStatus vfp,
//			  typename std::enable_if< std::integral_constant<bool,!APPLY_GENERAL_SEGMENT_ORDERING>>
//			>
//	bool establishNewVertices( Vertex * const & oldVrt );

	bool createNewElements();

	using AttachedVolumeElemInfo = support::AttachedFullDimElemInfo<Volume*, Face *, Edge *, IndexType, NormalVectorFacIntoVol>;

	using VecAttachedVolumeElemInfo = std::vector<AttachedVolumeElemInfo>;
	using AttVecAttachedVolumeElemInfo = Attachment<VecAttachedVolumeElemInfo>;

//	AttVecAttachedVolumeElemInfo m_aAdjVolElmInfo;
//	Grid::VertexAttachmentAccessor<AttVecAttachedVolumeElemInfo> m_aaVolElmInfo;

	using SegmentVolElmInfo = VecAttachedVolumeElemInfo;
	using VecSegmentVolElmInfo = std::vector<SegmentVolElmInfo>;

	bool stasiAlgo( Vertex * const & oldVrt );
	int prepareStasi( Vertex * const & vrt, AttachedVolumeElemInfo & attVolElmInfo );

	using AttVecSegmentVolElmInfo = Attachment<VecSegmentVolElmInfo>;

	AttVecSegmentVolElmInfo m_attAdjVecSegVolElmInfo;
	Grid::VertexAttachmentAccessor<AttVecSegmentVolElmInfo> m_accsAttVecSegVolElmInfo;

	bool computeNormalKuhVolProcedure( Volume * const & kuhVol, Face * const & fac, NormalVectorFacIntoVol & normalIntoVol );

	bool enableVolOptAutoGenFac();

	bool establishNewVertizesStasiBased( Vertex * const & oldVrt );

//	template< IndexType NUM_SURR_FRACS, bool isBndryVrtx >
//	bool expandWithinTheSegment( Vertex * const & oldVrt, SegmentVolElmInfo const & segmVolElmInfo );

	IndexType shiftUnclosedFracFacesToUnclosedFractFaces( Vertex * const & vrt );

//	bool extracFractSudosOfSegment(SegmentVolElmInfo const & segmVolElmInfo, std::vector<IndexType> & sudosInSegment );

public:
	using SegmentLimitingSides = support::SegmentSides<Volume*,Face*,Edge*,IndexType,vector3,Vertex*>;

//	using SegmentLimitSidesPairSudoNorml = SegmentLimitingSides::PairSudoNormlV;
//
//	using VecSegmentLimitSidesPairSudoNorml = SegmentLimitingSides::VecPairSudoNormlV;

	using SegmentVrtxFracStatus = SegmentLimitingSides::VrtxFracStatus;

	using VecSegmentLimitingSides = std::vector<SegmentLimitingSides>;

	using AttVecSegmLimSid = Attachment<VecSegmentLimitingSides>;

	using SegLimSidesFractFace = SegmentLimitingSides::AttFractElm;
	using VecSegLimSidesFractFace = SegmentLimitingSides::VecAttFractElm;

private:

	AttVecSegmLimSid m_attVecSegmLimSid;

	Grid::VertexAttachmentAccessor<AttVecSegmLimSid> m_vrtxAttAccsVecSegmLimSid;

	// establish an attechment of type SegmentLimitingSides for all vertices
	// transfer the functionality of establishNewVertizesStasiBased here
	// incorportate also the free ending clefts here and if it is a free ending cleft and so on
	bool establishSegmentLimitingSidesInfo();
	// unclear if useful to implement this
	// Ziel eigentlich: die frei endenden Faces, die als offene Faces da sind, jedem Segment zu zu ordnen
	// irgendwie sollen für jedes Segment die segment limiting sides bestimmt werden, als attachment vector
	// die Segment limiting sides sollen auch die in ihnen eingeschlossenen unvollendeten Faces kennen#
	// vielleicht auch gar nicht als attachment nötig, da die Segmente ihren Vertex kennen

	// the artificial normals are for the case of two crossing fractures inside, and the case
	// of boundary vertices, i.e. one fracture at one boundary sudo, two fracture at one boundary sudo, one fracture at two boundary sudos
//	template< SegmentVrtxFracStatus seVrtFracStat >
	bool expandWithinTheSegment( SegmentLimitingSides const & segmLimSides );

	using PlaneDescriptor = support::ManifoldDescriptor<vector3>;

	using VecPlaneDescriptor = std::vector<PlaneDescriptor>;

	using PlaneDescriptorType = PlaneDescriptor::ManifoldType;

	bool computeCrossingPointOf3Planes( VecPlaneDescriptor const & vecPlaneDescr, vector3 & crossingPoint );

//	template<SegmentVrtxFracStatus svfs>
	// For the case of one inner fracture
//	bool computeShiftVector( VecPlaneDescriptor const & vecPl );

	// For the case of two and three inner fractures, and
	// the case of one or two fractures at one outer boundary subdomain
	// for the case of one fracture at one outer boundary subdomain
//	bool computeShiftVector( VecSegmentLimitSidesPairSudoNorml const & vecSegmLimSidPrSudoNrml );

	int splitInnerFreeFracEdgs();

	template<typename ELEMTYP>
	bool addElem( std::vector<ELEMTYP> & knownElems, ELEMTYP elemToAdd );

	template< bool FACES_HAVE_SAME_SUDO >
	bool fractFacesAreNeighboured( SegLimSidesFractFace const & fractFaceOne,
								   SegLimSidesFractFace const & fractFaceTwo,
								   Edge * & commonEdge
								  );

	bool assignNotCuttingEdgeUnclosedCrossingCleft( SegLimSidesFractFace const & slsffUncl, Edge * const & commonEdge, Edge * & shiftEdge );

	bool findShiftEdgeUncuttingCrossingCleft( VecSegLimSidesFractFace const & vecSegmLimSiFFUnclosed,
											  Face * const & endingFractManifCutting,
											  Face * const & endingFractManifNotCutting,
											  Edge * const & uncuttingLowDimEl,
											  Edge * & shiftLowDimEl
											);

//	std::vector<Face*> m_d_endingCrossingCleftFaces;
//	std::vector<Vertex*> m_d_endingCrossingCleftVrtcs;
//	std::vector<Edge*> m_d_cuttingEdges;
//	std::vector<Face*> m_d_crossingNeighboredNotEndingFaces;
//	std::vector<Face*> m_d_crossingNeighboredNotEndingFacesCommEdg;
////	std::vector<Edge*> otherEdgeOfCrossingNotEndingFace;
////	std::vector<Face*> nextFaceOfCrossingNotEndingFaces;
//	std::vector<Face*> m_d_notEndingCrossingFacesNotNeighbour;
//	std::vector<Edge*> m_d_allContributingEdges;

	// edges that connect two ending crossing cleft vertices directly, need to be splitted
	std::vector<Edge *> m_vecEdgeDirectConnectingEndingCrossCleftVrtcs;

	void assignDebugSubsets( bool intermediate );

	bool splitEdgesOfNeighboredEndingCrossingFracVrtcs();


	using EndingCrossingFractureSegmentInfo = support::EndingCrossingFractSegmentInfo<Volume*, Face*, Edge*, Vertex*, IndexType >;

	using VecEndingCrossingFractureSegmentInfo = std::vector<EndingCrossingFractureSegmentInfo>;

	// TODO FIXME all need to be initialized in constructor and attached and detached

	VecEndingCrossingFractureSegmentInfo m_vecEndCrossFractSegmInfo;

	using AttVecEndingCrossingFractureSegmentInfo = Attachment<VecEndingCrossingFractureSegmentInfo>;

	AttVecEndingCrossingFractureSegmentInfo m_attAtVrtxVecEndingCrossFractSegmInfo;

	Grid::VertexAttachmentAccessor<AttVecEndingCrossingFractureSegmentInfo> m_vrtxAttAccsVecEndingCrossFractSegmInfo;

//	ABool m_attAtVrtxIfVrtxIsEndingCrossingCleftVrtx;
//
//	Grid::VertexAttachmentAccessor<ABool> m_vrtxAttAccsVrtxIsEndingCrossingCleftVrtx;

	ABool m_attAtFaceIfFaceIsSegmLimFaceEndingCrossingCleft;

	Grid::FaceAttachmentAccessor<ABool> m_facAttAccsIfFaceIsSegmLimFaceEndingCrossingCleft;

	ABool m_attAtVolIfVolTouchesEndingCrossingCleft;

	Grid::VolumeAttachmentAccessor<ABool> m_volAttAccsVolTouchesEndingCrossingCleft;

	template<typename ELMTYP>
	bool checkIfContentUnique( std::vector<ELMTYP> const & vecTest, std::vector<ELMTYP> & content, IndexType mandatoryDifferentElems );

	bool computeCrossPointOfPlaneWithLine( PlaneDescriptor const & shiftedPlane, Edge * const & shiftDirectionEdg, Vertex * const & oldVrt, vector3 & posCrossingPt );

	ABool m_attAtVrtxIfVrtxArisesFromExpandedEndingCrossingCleft;
	Grid::VertexAttachmentAccessor<ABool> m_vrtxAttAccsVrtxArisesFromExpandedEndingCrossingCleft;

	std::vector<Vertex*> m_vrtxArisesFromExpandedEndingCrossingCleft;

	bool etablishVolumesAtEndingCrossingClefts( std::vector<Volume*> & newFractureVolumes, std::vector<IndexType> & subsOfNewVolumes );

	IndexType deleteEndingCrossingCleftOrigFacs();

};

// specification has to be declared outside central class context, else compilation error

//template<ArteExpandFracs3D::SegmentVrtxFracStatus::oneFracSuDoAtt>
//bool computeShiftVector( );


// for only one surrounding subdom around the segment
//template<>
//bool ArteExpandFracs3D::expandWithinTheSegment<1,false>( Vertex * const & oldVrt, SegmentVolElmInfo const & segmVolElmInfo );


//using ArteOneFractCrossSegment = ArteExpandFracs3D::SegmentVrtxFracStatus::oneFracSuDoAtt;

//template<>
//bool ArteExpandFracs3D::expandWithinTheSegment<ArteOneFractCrossSegment>( SegmentLimitingSides const & segmLimSides );
//template<>
//bool ArteExpandFracs3D::expandWithinTheSegment<ArteExpandFracs3D::SegmentVrtxFracStatus::oneFracSuDoAtt>( SegmentLimitingSides const & segmLimSides );
//
//
//template<bool artificialNormalTwo, bool artificialNormalThree>
//bool ArteExpandFracs3D::expandWithinTheSegment<ArteExpandFracs3D::SegmentVrtxFracStatus::threeFracSuDoAtt>( SegmentLimitingSides const & segmLimSides, std::vector<Volume*> const & vecVolsOfSegment );


//template <>
//bool ArteExpandFracs3D::establishNewVertices< true,
//											  ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt
//											>( Vertex * const & oldVrt );
//
//template <>
//bool ArteExpandFracs3D::establishNewVertices< false,
//											  ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt
//											>( Vertex * const & oldVrt );


} /* namespace ug */

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_ */
