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
	VecVolManifVrtxCombi m_vecVolManifVrtxCombiToShrink4Diams;

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
};

} /* namespace diamonds */

} /* namespace arte */

} /* namespace ug */

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_DIAMONDSESTABLISH3D_H_ */
