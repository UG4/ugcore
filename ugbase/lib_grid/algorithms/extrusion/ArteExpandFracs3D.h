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

	bool run();

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

	using IndexType = unsigned short;

	using AttachedFaceEdgeSudo = support::AttachedElem<Face*,Edge*,IndexType>;

	using VecAttachedFaceEdgeSudo = std::vector<AttachedFaceEdgeSudo>;

	using EdgePair = std::pair<Edge*,Edge*>;

	using VertxFracPropts = support::VertexFracturePropertiesVol<IndexType, AttachedFaceEdgeSudo>;

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

	ABool m_aAdjMarkerB; // used to know if an face is frac face

	Grid::FaceAttachmentAccessor<ABool> m_aaMarkFaceB;

	bool countAndSelectFracBaseNums();

	std::vector<Face*> m_originalFractureFaces;

	bool assignOrigFracInfos();


//	using AttVrtVec = Attachment<std::vector<Vertex*> >;
//	AttVrtVec m_attVrtVec;
//	Grid::VolumeAttachmentAccessor<AttVrtVec> m_aaVrtVecVol;

	using AttVecEdge = Attachment<std::vector<Edge*>>;
	using AttVecFace = Attachment<std::vector<Face*>>;
	using AttVecVol = Attachment<std::vector<Volume*>>;

	AttVecEdge m_aAdjInfoEdges;
	AttVecFace m_aAdjInfoFaces;
	AttVecVol m_aAdjInfoVols;

	Grid::VertexAttachmentAccessor<AttVecEdge> m_aaVrtInfoAssoEdges;
	Grid::VertexAttachmentAccessor<AttVecFace> m_aaVrtInfoAssoFaces;
	Grid::VertexAttachmentAccessor<AttVecVol> m_aaVrtInfoAssoVols;

	bool establishNewVrtBase();

	bool generateVertexInfos();

	bool loop2EstablishNewVertices();

//	ug::support::VertexFractureQuadrupel<Face*,vector3,Volume*,Edge*,IndexType> bla();

//	using VrtxFractrQuadrplArte3D = support::VertexFractureQuadrupel<Face*,vector3,Volume*,Edge*,IndexType>;
//
//	using VrtxFractrQuadrplArte3DVec = std::vector<VrtxFractrQuadrplArte3D>;
//
//	VrtxFractrQuadrplArte3DVec m_vrtxFractrQuadrplVec;

	using VertFracTrip = support::VertexFractureTripleMF<Face*, IndexType, Volume*, vector3, Edge*>;

	using VecVertFracTrip = std::vector<VertFracTrip>;

	using AttVecVertFracTrip = Attachment<VecVertFracTrip>;

	AttVecVertFracTrip m_aAdjInfoAVVFT;

	Grid::VertexAttachmentAccessor<AttVecVertFracTrip> m_aaVrtInfoFraTri;

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

	bool isVrtxSurroundedByFracFaces( Vertex * const & vrt, VertxFracPropts & vrtxFracPrps );
//									  VecPairSudoBool & sudoSurrounded );

	// transform to template soon
	bool sortElemCircleIsClosed( VecAttachedFaceEdgeSudo const & vecAttFac,
								 VecAttachedFaceEdgeSudo & vecSortedFac,
								 int startFaceIndexUser = -1,
//								 int endFaceIndexUser = -1,
//								 IndexType startEdgeIndexUser = -1,
//								 IndexType endEdgeIndexUser = -1
//								 Face * const & startFacUser = nullptr,
//								 Face * const & endFacUser = nullptr,
								 Edge * const & startEdgUser = nullptr,
								 Edge * const & endEdgUser = nullptr
								);

public:
	using VrtxFracProptsStatus = VertxFracPropts::VrtxFracStatus;

private:
//	static_assert< std::is_same< VrtxFracProptsStatus,support::VertexFracturePropertiesVol::VrtxFracStatus>::value );
//	static_assert( std::is_same<VrtxFracProptsStatus,support::VertexFracturePropertiesVol::VrtxFracStatus>::value );
//	static_assert( std::is_same<VrtxFracProptsStatus,support::VertexFracturePropertiesVol<IndexType, AttachedFaceEdgeSudo>::VrtxFracStatus>::value );

	template<VrtxFracProptsStatus vfps>
//	template<support::VertexFracturePropertiesVol::VrtxFracStatus vfp>
//	template<int I>
	bool establishNewVertices( Vertex * const & oldVrt )
	{
		UG_THROW("too general, should not be called presently " << std::endl);
		return false;
	};

	bool createNewElements();

	using AttachedVolumeElemInfo = support::AttachedFullDimElemInfo<Volume*, Face *, Edge *, IndexType>;
	using VecAttachedVolumeElemInfo = std::vector<AttachedVolumeElemInfo>;
	using AttVecAttachedVolumeElemInfo = Attachment<VecAttachedVolumeElemInfo>;

	AttVecAttachedVolumeElemInfo m_aAdjVolElmInfo;
	Grid::VertexAttachmentAccessor<AttVecAttachedVolumeElemInfo> m_aaVolElmInfo;


};

// specification has to be declared outside central class context, else compilation error
template <>
bool ArteExpandFracs3D::establishNewVertices<ArteExpandFracs3D::VrtxFracProptsStatus::oneFracSuDoAtt>( Vertex * const & oldVrt );

} /* namespace ug */

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_ */
