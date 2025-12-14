/*
 * DiamondsEstablish3D.h
 *
 *  Created on: 08.12.2025
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDSESTABLISH3D_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDSESTABLISH3D_H_

#include <boost/function.hpp>
#include <stack>
#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"

#include "DiamondInfo.h"

namespace ug
{

namespace arte
{

namespace diamonds
{


class DiamondsEstablish3D
{
public:

	using IndexType = unsigned short;

	using VolManifVrtxCombi = diamonds::VolManifVrtxCombi<Volume*,Face*,Edge*, Vertex*, IndexType>;

	using VecVolManifVrtxCombi = std::vector<VolManifVrtxCombi>;

	using VrtxPair = std::pair<Vertex*,Vertex*>;


	DiamondsEstablish3D( Grid & grid, SubsetHandler & sh, VecVolManifVrtxCombi const & vecVolManifVrtxC );

	virtual ~DiamondsEstablish3D();

	bool createTheDiamonds();


private:

	Grid & m_grid;
	SubsetHandler & m_sh;
	Grid::VertexAttachmentAccessor<APosition> m_aaPos;
	VecVolManifVrtxCombi m_vecVolManifVrtxCombiToShrink4Diams;

	bool initialize();
	bool figureOutTheEdges();

	bool findRegions2BShrinked();

	using VolumeElementTwin = FulldimLowdimTwin<Volume*, Edge*, IndexType>;

	using VolumeElementFaceQuintuplet = FullLowDimManifQuintuplet<Volume*, Face*, Edge*, Vertex*, IndexType>;

	using VecVolumeElementFaceQuintuplet = std::vector<VolumeElementFaceQuintuplet>;


	using Elems2BQuenched = ElemsToBeQuenched4DiamSpace<Volume*, Face*, Edge*, Vertex*, IndexType>;
	using VecElems2BQuenched = std::vector<Elems2BQuenched>;

	using PairVolFacVrtxCmb = std::pair<VolManifVrtxCombi,VolManifVrtxCombi>;

	bool trafoVolFacVrtxCombiPair2FullLowDimManifQuintuplet( PairVolFacVrtxCmb & prVolFacVrtxC,
															 VolumeElementFaceQuintuplet & vef5
															);

	bool establishElems2BeQuenched();

	VecVolumeElementFaceQuintuplet m_vecVolElmFac5;

	VecElems2BQuenched m_vecElems2BQuenched;

	std::vector<Volume*> m_disappearingVols;
	std::vector<Face*> m_disappearingFacs;
//	std::vector<Edge*> m_disappearingEdgs;

	using EdgePair = std::pair<Edge*,Edge*>;

	bool sortElems2BQuenched();

	Selector m_sel;

	bool setSelector();

	using ElemGroupVrtx2BQuenched4Diams = ElemGroupVrtxToBeQuenched4DiamSpace<Volume*, Face*, Edge*, Vertex*, IndexType>;

	using VecElemGroupVrtx2BQnchd4D = std::vector<ElemGroupVrtx2BQuenched4Diams>;

	VecElemGroupVrtx2BQnchd4D m_vecElemGroupVrtx2BQuenched;

	bool attachMarkers();
	bool detachMarkers();
	bool assignBasicAtts();

	using AttElemGrpVrtx2BQuenchd = Attachment<ElemGroupVrtx2BQuenched4Diams>;

	AttElemGrpVrtx2BQuenchd m_attElmGrpVrtx2BQnchd;

	Grid::VertexAttachmentAccessor<AttElemGrpVrtx2BQuenchd> m_attAccsElmGrpVrtx2BQnchd;

	bool trafoCollectedInfo2Attachments();

	void debugE2bQ(Elems2BQuenched & e2bq);

	ABool m_attMarkTwinFace; // used to know if an face is twin face

	Grid::FaceAttachmentAccessor<ABool> m_attAccsFacIsTwinFac;

	ABool m_attMarkShiftFace; // used to know if an face is shift face

	Grid::FaceAttachmentAccessor<ABool> m_attAccsFacIsShiftFac;

	ABool m_attMarkVolGetsShrinked;
	Grid::VolumeAttachmentAccessor<ABool> m_attAccsVolGetsShrinked;

	using AttVrtVec = Attachment<std::vector<Vertex*> >;

	AttVrtVec m_attVrtVec;

	Grid::VolumeAttachmentAccessor<AttVrtVec> m_attAccsVrtVecVol;

	ABool m_attMarkVrtxIsCenterVrtx;

	Grid::VertexAttachmentAccessor<ABool> m_attAccsVrtxIsCenterVrtx;

	ABool m_attMarkVrtxIsShiftVrtx;

	Grid::VertexAttachmentAccessor<ABool> m_attAccsVrtxIsShiftVrtx;

	ABool m_attMarkEdgeIsShiftEdge;

	Grid::EdgeAttachmentAccessor<ABool> m_attAccsEdgeIsShiftEdge;

	bool trafoQuintupleInfo2Attachments(VolumeElementFaceQuintuplet & vef5);

	using VolumeEdgeTwin = FulldimLowdimTwin<Volume*,Edge*,IndexType>;
	using PairVolumeEdgeTwin = std::pair<VolumeEdgeTwin, VolumeEdgeTwin>;

	using AttVecInt = Attachment<std::vector<IndexType>>;

	AttVecInt m_attInfoVecSudosTouchingVrtx;

	Grid::VertexAttachmentAccessor<AttVecInt> m_attAccsInfoVecSudosTouchingVrtx;

	bool assignMidPointOfShiftVrtcs( Elems2BQuenched & e2bq );

	ABool m_attVrtxIsMidPtOfShiftVrtx;

	Grid::VertexAttachmentAccessor<ABool> m_attAccsVrtxIsMidPtOfShiftVrtcs;

	bool distributeInfosForShrinkingVols();

	bool createConditionForNewVrtcs();
};

} /* namespace diamonds */

} /* namespace arte */

} /* namespace ug */

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDSESTABLISH3D_H_ */
