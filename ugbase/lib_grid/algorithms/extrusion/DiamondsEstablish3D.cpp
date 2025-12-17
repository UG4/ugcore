/*
 * DiamondsEstablish3D.cpp
 *
 *  Created on: 08.12.2025
 *      Author: Markus Knodel
 */

#include <lib_grid/algorithms/extrusion/DiamondsEstablish3D.h>
#include <typeinfo>

namespace ug
{

namespace arte
{

namespace diamonds
{

DiamondsEstablish3D::DiamondsEstablish3D( Grid & grid,
										  SubsetHandler & sh,
										  DiamondsEstablish3D::VecVolManifVrtxCombi const & vecVolManifVrtxC
										 )
		:
			m_grid(grid),
			m_sh(sh),
			m_aaPos(Grid::VertexAttachmentAccessor<APosition>()),
			m_vecVolManifVrtxCombiToShrink4Diams(vecVolManifVrtxC),
			m_vecVolElmFac5(VecVolumeElementFaceQuintuplet()),
			m_vecElems2BQuenched(VecElems2BQuenched()),
			m_disappearingVols(std::vector<Volume*>()),
			m_disappearingFacs(std::vector<Face*>()),
			m_vecElemGroupVrtx2BQuenched(VecElemGroupVrtx2BQnchd4D()),
			m_attElmGrpVrtx2BQnchd(AttElemGrpVrtx2BQuenchd()),
			m_attAccsElmGrpVrtx2BQnchd(Grid::VertexAttachmentAccessor<AttElemGrpVrtx2BQuenchd>()),
			m_attMarkTwinFace(ABool()),
			m_attAccsFacIsTwinFac(Grid::FaceAttachmentAccessor<ABool>()),
			m_attMarkShiftFace(ABool()),
			m_attAccsFacIsShiftFac(Grid::FaceAttachmentAccessor<ABool>()),
			m_attMarkShiftTriangleFace(ABool()),
			m_attAccsFacIsShiftTriangleFac(Grid::FaceAttachmentAccessor<ABool>()),
			m_attMarkShiftQuadriliteralFace(ABool()),
			m_attAccsFacIsShiftQuadriliteralFac(Grid::FaceAttachmentAccessor<ABool>()),
			m_attMarkVolGetsShrinked(ABool()),
			m_attAccsVolGetsShrinked(Grid::VolumeAttachmentAccessor<ABool>()),
			m_attVrtVec(AttVrtVec()),
			m_attAccsVrtVecVol(Grid::VolumeAttachmentAccessor<AttVrtVec>()),
//			m_attVrtVecFace(AttVrtVec()),
//			m_attAccsVrtVecFace(Grid::VolumeAttachmentAccessor<AttVrtVec>()),
			m_attMarkVrtxIsCenterVrtx(ABool()),
			m_attAccsVrtxIsCenterVrtx(Grid::VertexAttachmentAccessor<ABool>()),
			m_attMarkVrtxIsShiftVrtx(ABool()),
			m_attAccsVrtxIsShiftVrtx(Grid::VertexAttachmentAccessor<ABool>()),
			m_attMarkEdgeIsShiftEdge(ABool()),
			m_attAccsEdgeIsShiftEdge(Grid::EdgeAttachmentAccessor<ABool>()),
			m_attInfoVecSudosTouchingVrtx(AttVecInt()),
			m_attAccsInfoVecSudosTouchingVrtx(Grid::VertexAttachmentAccessor<AttVecInt>()),
			m_attVrtxIsMidPtOfShiftVrtx(ABool()),
			m_attAccsVrtxIsMidPtOfShiftVrtcs(Grid::VertexAttachmentAccessor<ABool>()),
			m_vecCombiShiftVrtxMidVrtx(VecCombiShiftVrtxMidVrtx()),
			m_attMidPtVrtxOfShiftVrtx(AVertex()),
			m_attAccsMidPtVrtxOfShiftVrtx(Grid::VertexAttachmentAccessor<AVertex>()),
			m_attCenterVrtxOfShiftVrtx(AVertex()),
			m_attAccsCenterVrtxOfShiftVrtx(Grid::VertexAttachmentAccessor<AVertex>()),
			m_vecCombiNewVolsTwoCross(VecCombiNewVolsProps()),
			m_vecCombiNewVolsThreeCross(VecCombiNewVolsProps()),
			m_sudosTable(VecIndxVec()),
			m_vecCombiCntrVrtxSudo(VecCombiCntrVrtxSudo()),
			m_attEdgeCanBeRemoved(ABool()),
			m_attAccsEdgeCanBeRemoved(Grid::EdgeAttachmentAccessor<ABool>()),
			m_attNewSudoOfVrtx(AInt()),
			m_attAccsNewSudoOfVrtx(Grid::VertexAttachmentAccessor<AInt>()),
			m_faces2BDeletedAtLastStep(std::vector<Face*>()),
			m_edges2BDeletedAtLastStep(std::vector<Edge*>())
{
	//	// Notloesung, nicht in die erste Initialisierung vor geschweifter Klammer, da copy constructor privat
		m_sel = Selector();
}

//bool DiamondsEstablish3D::createTheDiamonds()
//{
//	IndexType sudosVols = m_sh.num_subsets();
//
//	IndexType sudosEdges = sudosVols + 1;
//
//	IndexType sudosFaces = sudosVols + 2;
//
//	for( auto & vmvcd : m_vecVolManifVrtxCombiToShrink4Diams )
//	{
//		Volume* vol;
//		vmvcd.spuckVol(vol);
//
//		m_sh.assign_subset(vol, sudosVols);
//
//		IndexType numLowdimElmsFnd = vmvcd.computeTheLowdimElm(m_grid);
//
//		if( numLowdimElmsFnd != 1 )
//		{
//			UG_LOG("number of lowdim elems found strange " << numLowdimElmsFnd << std::endl);
//			UG_THROW("number of lowdim elems found strange " << numLowdimElmsFnd << std::endl);
//		}
//
//		Edge* edge;
//		vmvcd.spuckLowdimElem( edge );
//
//		if( edge == nullptr )
//		{
//			UG_LOG("Edge nicht gefunden " << std::endl);
//			UG_THROW("Edge nicht gefunden " << std::endl);
//		}
//		m_sh.assign_subset( edge, sudosEdges);
//
//		Face * fac;
//		vmvcd.spuckManif(fac);
//		m_sh.assign_subset( fac, sudosFaces );
//	}
//
//	UG_LOG("Established diamonds" << std::endl);
//
//
//	return true;
//}

DiamondsEstablish3D::~DiamondsEstablish3D()
{
	// Auto-generated destructor stub
}

bool DiamondsEstablish3D::initialize()
{
	if(!m_grid.has_vertex_attachment(aPosition) )
	{
		UG_LOG("Error in ExpandFractures Arte 3D: Missing position attachment");
		return false;
	}

	m_aaPos = Grid::VertexAttachmentAccessor<APosition>(m_grid, aPosition);

	return true;
}

/////////////////////////////////////////////////////////////////


bool DiamondsEstablish3D::createTheDiamonds()
{
	if( ! initialize())
	{
		UG_LOG("initialization diamonds did not work " << std::endl);
		return false;
	}

	if( ! setSelector() )
		return false;

	UG_LOG("selektiert" << std::endl);

	if( ! figureOutTheEdges() )
	{
		UG_LOG("Edges not found " << std::endl);
		return false;
	}

	UG_LOG("edges figured out " << std::endl);

	if( ! findRegions2BShrinked())
	{
		UG_LOG("Regions to be shrinked not found " << std::endl);
		return false;
	}

	UG_LOG("regions to be shrinked found " << std::endl);

	if( ! establishElems2BeQuenched() )
	{
		UG_LOG("quenching elements not establishable " << std::endl);
		return false;
	}

	UG_LOG("quenching elements established " << std::endl);

	if( ! sortElems2BQuenched())
	{
		UG_LOG("sorting quench impossible" << std::endl);
		return false;
	}

	UG_LOG("sorted elements finished " << std::endl);

	if( ! attachMarkers())
	{
		UG_LOG("Markers do not want D" << std::endl);
		return false;
	}

	UG_LOG("markers attached" << std::endl);

	if( ! assignBasicAtts())
	{
		UG_LOG("unassignale basic att" << std::endl);
		return false;
	}

	UG_LOG("basic attcs assigned" << std::endl);

	if( ! trafoCollectedInfo2Attachments())
	{
		UG_LOG("untrafoable infos" << std::endl);
		return false;
	}

	UG_LOG("trafo collected infos " << std::endl);

	//	disable selection inheritance to avoid infinite recursion.
	m_sel.enable_selection_inheritance(false);

	if( ! createConditionForNewVrtcs())
	{
		UG_LOG("conditions for new vertices do not work " << std::endl);
		return false;
	}

	UG_LOG("conditions for new vertices created" << std::endl);

	if( ! distributeInfosForShrinkingVols())
	{
		UG_LOG("info distribution did not work " << std::endl);
		return false;
	}

	UG_LOG("info distributed" << std::endl);

	if( ! determineShiftFaces())
	{
		UG_LOG("shift faces not determinable " << std::endl);
		return false;
	}

	UG_LOG("detect removable edges " << std::endl);

	if( ! detectRemovableEdges())
	{
		UG_LOG("removable edges not detected " << std::endl);
		return false;
	}

	UG_LOG("shrink volumes " << std::endl);

	if( ! shrinkVolumes())
	{
		UG_LOG("shrinken schief gegangen " << std::endl);
		return false;
	}

	UG_LOG("volumes shrinked " << std::endl);

	if( ! postprocessNewDiamVols())
	{
		UG_LOG("postprocessing diam vols not working" << std::endl);
		return false;
	}

	UG_LOG("diam vols postprocessed" << std::endl);

	if( ! detachMarkers())
	{
		UG_LOG("Markers do not get detouched D" << std::endl);
		return false;
	}

	UG_LOG("markers detached" << std::endl);

	UG_LOG("noch nicht soweit: Established diamonds" << std::endl);



	return true;
}

/////////////////////////////////////////////////////////////////


bool DiamondsEstablish3D::setSelector()
{

	m_sel.assign_grid(m_grid);

	m_sel.enable_autoselection(false);
	m_sel.enable_selection_inheritance(true);	//required for select and mark, disabled later
	m_sel.enable_strict_inheritance(false);


	return true;
}

/////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::figureOutTheEdges()
{

	// compute the edges connecting the old and shift vertex

	for( auto & vmvcd : m_vecVolManifVrtxCombiToShrink4Diams )
	{
//		IndexType numLowdimElmsFnd = vmvcd.checkIntegrity(m_grid);
//		if( numLowdimElmsFnd != 1 )
		if( ! vmvcd.checkIntegrity(m_grid) )
		{
			UG_LOG("number of lowdim elems found strange " << std::endl);
			return false;
		}
	}

	UG_LOG("Edges found" << std::endl);

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::findRegions2BShrinked()
{
	UG_LOG("want to find regions to be shrinked" << std::endl);


	VecVolManifVrtxCombi vecVolManifVrtxCopy = m_vecVolManifVrtxCombiToShrink4Diams;

	int d_out = 0;

	for( typename VecVolManifVrtxCombi::iterator itVMVOuter = vecVolManifVrtxCopy.begin();
				  itVMVOuter != vecVolManifVrtxCopy.end(); )
	{

//		UG_LOG("out round " << d_out << std::endl);

		bool partnerFound = false;

		VolManifVrtxCombi outer = *itVMVOuter;

		VrtxPair oldAndShiftVrtxOuter;
		outer.spuckOldAndShiftVrtx( oldAndShiftVrtxOuter );

		Vertex * oldVrtxOuter = oldAndShiftVrtxOuter.first;
		Vertex * shiftVrtxOuter = oldAndShiftVrtxOuter.second;

		Face * faceOuter;
		outer.spuckManif(faceOuter);

		IndexType sudoOuter = outer.spuckSudo();

		int d_in = 0;

		for( typename VecVolManifVrtxCombi::iterator itVMVInner = itVMVOuter + 1;
					  itVMVInner != vecVolManifVrtxCopy.end(); )
		{
//			UG_LOG("start in " << d_in << std::endl);

			VolManifVrtxCombi inner = *itVMVInner;

			VrtxPair oldAndShiftVrtxInner;
			inner.spuckOldAndShiftVrtx( oldAndShiftVrtxInner );

			Vertex * oldVrtxInner = oldAndShiftVrtxInner.first;

			if( oldVrtxInner == oldVrtxOuter )
			{
				Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;

				// if connected by a face, twin to outer found, partner may exist or not
				// check if the face is identical, else is from the partner twin

				Face * faceInner;
				inner.spuckManif(faceInner);

				if( faceInner == faceOuter )
				{
					// twin found

					PairVolFacVrtxCmb prVolFacVrtxC( outer, inner );
					VolumeElementFaceQuintuplet vef5;

					trafoVolFacVrtxCombiPair2FullLowDimManifQuintuplet( prVolFacVrtxC, vef5 );

					if( ! vef5.checkIntegrity() )
					{
						UG_LOG("strange twin produced" << std::endl);
						return false;
					}

					m_vecVolElmFac5.push_back(vef5);

					itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);

					partnerFound = true;

				}
				else
				{
					itVMVInner++;
				}

				if( partnerFound )
				{
					break;
				}
			}
			else
			{
				itVMVInner++;
			}
			if( partnerFound )
			{
				break;
			}


//			UG_LOG("end in " << d_in << std::endl);
			d_in++;
		}

//		UG_LOG("first round finished " << std::endl);
		if( ! partnerFound )
		{
			UG_LOG("no partner found " << std::endl);
			return false;
		}
		else
		{
			itVMVOuter = vecVolManifVrtxCopy.erase(itVMVOuter);
		}

//		UG_LOG("end out " << d_out << std::endl);
		d_out++;
	}

	UG_LOG("end of search reached diams " << std::endl);

//	for( auto & v: m_vecVolElmFac5 )
//	{
//		IndexType sudoNum = m_sh.num_subsets();
//
//		IndexType sudoVols = sudoNum;
//		IndexType sudoFacs = sudoNum+1;
//		IndexType sudoEdgs = sudoNum+2;
//		IndexType sudoVrtx = sudoNum+3;
//
//		Volume * vol1;
//		Volume * vol2;
//		Face * fac;
//		Edge * edg1;
//		Edge * edg2;
//		Vertex * vrtC;
//		Vertex * vrtE1;
//		Vertex * vrtE2;
//
//		std::pair<VolumeElementTwin,VolumeElementTwin> pvv;
//
//		v.spuckPairFullLowDimTwin(pvv);
//
//		pvv.first.spuckFullDimElem(vol1);
//		pvv.second.spuckFullDimElem(vol2);
//
//		v.spuckManifElem( fac );
//
//		pvv.first.spuckLowDimElem(edg1);
//		pvv.second.spuckLowDimElem(edg2);
//
//		v.spuckCenterVertex(vrtC);
//
//		VrtxPair vp;
//		v.spuckShiftVrtcs(vp);
//
//		vrtE1 = vp.first;
//		vrtE2 = vp.second;
//
//		m_sh.assign_subset(vol1, sudoVols);
//		m_sh.assign_subset(vol2, sudoVols);
//		m_sh.assign_subset(fac, sudoFacs);
//		m_sh.assign_subset(edg1, sudoEdgs);
//		m_sh.assign_subset(edg2, sudoEdgs);
//		m_sh.assign_subset(vrtC, sudoVrtx);
//		m_sh.assign_subset(vrtE1, sudoVrtx);
//		m_sh.assign_subset(vrtE2, sudoVrtx);
//
//	}


	UG_LOG("found regions to be shrinked" << std::endl);


	return true;
}


////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::trafoVolFacVrtxCombiPair2FullLowDimManifQuintuplet(
		PairVolFacVrtxCmb & prVolFacVrtxC, VolumeElementFaceQuintuplet & vef5 )
{

	VolManifVrtxCombi & mvcOne = prVolFacVrtxC.first;
	VolManifVrtxCombi & mvcTwo = prVolFacVrtxC.second;

	IndexType sudoOne = mvcOne.spuckSudo();
	IndexType sudoTwo = mvcTwo.spuckSudo();

	if( sudoOne != sudoTwo )
	{
		UG_LOG("sudos differ " << std::endl);
		return false;
	}

	IndexType sudo = sudoOne;

	Face * faceOne;
	Face * faceTwo;

	mvcOne.spuckManif(faceOne);
	mvcTwo.spuckManif(faceTwo);

	if( faceOne != faceTwo )
	{
		UG_LOG("faces differ " << std::endl);
		return false;
	}

	Face * connectingFace = faceOne;

	VrtxPair oldAndShiftVrtxOne;
	mvcOne.spuckOldAndShiftVrtx( oldAndShiftVrtxOne );

	VrtxPair oldAndShiftVrtxTwo;
	mvcTwo.spuckOldAndShiftVrtx( oldAndShiftVrtxTwo );

	Vertex * oldVrtxOne = oldAndShiftVrtxOne.first;
	Vertex * oldVrtxTwo = oldAndShiftVrtxTwo.first;

	if( oldVrtxOne != oldVrtxTwo )
	{
		UG_LOG("center vertices not identical " << std::endl);
		return false;
	}

	Vertex * oldVrtx = oldVrtxOne;

	Vertex * shiftVrtxOne = oldAndShiftVrtxOne.second;
	Vertex * shiftVrtxTwo = oldAndShiftVrtxTwo.second;

	if( shiftVrtxOne == shiftVrtxTwo )
	{
		UG_LOG("shift vertices coincide but should not " << std::endl);
		return false;
	}

	Volume * volOne;
	Volume * volTwo;

	mvcOne.spuckFulldimElem( volOne );
	mvcTwo.spuckFulldimElem( volTwo );

	if( volOne == volTwo )
	{
		UG_LOG("volumes coincide but should not " << std::endl);
		return false;
	}

	Edge * edgeOne;
	Edge * edgeTwo;

	mvcOne.spuckLowDimElem( edgeOne );
	mvcTwo.spuckLowDimElem( edgeTwo );
	if( edgeOne == edgeTwo )
	{
		UG_LOG("edges coincide but should not " << std::endl);
		return false;
	}

	VolumeElementTwin volElTwinOne( volOne, edgeOne, sudo );
	VolumeElementTwin volElTwinTwo( volTwo, edgeTwo, sudo );

	if( ! volElTwinOne.checkIntegrity() || ! volElTwinTwo.checkIntegrity() )
	{
		UG_LOG("twins of vol edge not integer " << std::endl);
		return false;
	}

	std::pair<VolumeElementTwin,VolumeElementTwin> volElTwinPair( volElTwinOne, volElTwinTwo );

	vef5 = VolumeElementFaceQuintuplet( volElTwinPair, connectingFace );

	if( ! vef5.checkIntegrity() )
	{
		UG_LOG( "quitent not integer " << std::endl );
		return false;
	}

	return true;
}


////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::establishElems2BeQuenched()
{
	VecVolumeElementFaceQuintuplet vvef5 = m_vecVolElmFac5;

	for( typename VecVolumeElementFaceQuintuplet::iterator itOuter = vvef5.begin();
			      itOuter != vvef5.end();
	   )
	{
		VolumeElementFaceQuintuplet vef5Outer = *itOuter;

		VecVolumeElementFaceQuintuplet vfld5ThisEdgePr;

		vfld5ThisEdgePr.push_back(vef5Outer);

		// now collect all quintuplets belonging to the same edge pair

		EdgePair edgPrOuter;

		vef5Outer.spuckPairLowDimElem(edgPrOuter);

		for( typename VecVolumeElementFaceQuintuplet::iterator itInner = vvef5.begin() + 1;
				      itInner != vvef5.end();
		   )
		{
			VolumeElementFaceQuintuplet vef5Inner = *itInner;

			EdgePair edgPrInner;
			vef5Inner.spuckPairLowDimElem(edgPrInner);

			if( edgPrInner == edgPrOuter )
			{
				vfld5ThisEdgePr.push_back(vef5Inner);

				itInner = vvef5.erase(itInner);
			}
			else // if( edgPrInner != edgPrOuter )
			{
				EdgePair testEdges = edgPrInner;
				std::swap( testEdges.first, testEdges.second );

				if( testEdges == edgPrOuter )
				{
					if( ! vef5Inner.swapEntries() )
					{
						UG_LOG("swapping not worked " << std::endl);
						return false;
					}

					if( ! vef5Inner.checkIntegrity())
					{
						UG_LOG("not integer any more after sapping test" << std::endl);
						return false;
					}

					EdgePair testAgainEdgPr;
					vef5Inner.spuckPairLowDimElem(testAgainEdgPr);

					if( testAgainEdgPr == edgPrOuter )
					{
						vfld5ThisEdgePr.push_back(vef5Inner);

						itInner = vvef5.erase(itInner);
					}
					else
					{
						UG_LOG("should be correct but is not " << std::endl);
						return false;
					}
				}
				else
				{
					itInner++;
				}
			}
		}

		Elems2BQuenched elem2BQuenched( vfld5ThisEdgePr );

		if( ! elem2BQuenched.checkIntegrity())
		{
			UG_LOG("an elem to be quenched not integer" << std::endl);
			return false;
		}

		if( ! assignMidPointOfShiftVrtcs(elem2BQuenched))
		{
			UG_LOG("mid point not assignable" << std::endl);
			return false;
		}

		m_vecElems2BQuenched.push_back(elem2BQuenched);

		itOuter = vvef5.erase(itOuter);
	}

//	int d_q = 10;

	return true;
}

//////////////////////////////////////////////////////////////////////////////7

void DiamondsEstablish3D::debugE2bQ(Elems2BQuenched & e2bq)
{
//	for( Elems2BQuenched & e2bq : m_vecElems2BQuenched )
	{
		IndexType sudoNum = m_sh.num_subsets();

		IndexType sudoVols = sudoNum;
		IndexType sudoFacs = sudoNum+1;
		IndexType sudoEdgs = sudoNum+2;
		IndexType sudoVrtx = sudoNum+3;

		VecVolumeElementFaceQuintuplet vve5;

		EdgePair edgP;
		e2bq.spuckPairLowDimElem(edgP);

		e2bq.spuckVecFullLowDimManifQuintuplet(vve5);

		Vertex * centerVrtx;
		e2bq.spuckCenterVertex(centerVrtx);

		vector3 vertexLocation;

		if( centerVrtx != nullptr )
			vertexLocation = m_aaPos[centerVrtx];
		else
		{
			UG_LOG("center null " << std::endl);
		}

		number x = vertexLocation[0];
		number y = vertexLocation[1];
		number z = vertexLocation[2];

		number distSq = ( x - 0.5 )*( x - 0.5 ) + ( y - 0.5 )*( y - 0.5 ) + ( z - 0.5 )*( z - 0.5 );

		if( distSq < 0.02 )
		{

			for( VolumeElementFaceQuintuplet & v: vve5 )
			{
				Volume * vol1;
				Volume * vol2;
				Face * fac;
				Edge * edg1;
				Edge * edg2;
				Vertex * vrtC;
				Vertex * vrtE1;
				Vertex * vrtE2;

				std::pair<VolumeElementTwin,VolumeElementTwin> pvv;

				v.spuckPairFullLowDimTwin(pvv);

				pvv.first.spuckFullDimElem(vol1);
				pvv.second.spuckFullDimElem(vol2);

				v.spuckManifElem( fac );

				pvv.first.spuckLowDimElem(edg1);
				pvv.second.spuckLowDimElem(edg2);

				v.spuckCenterVertex(vrtC);

				VrtxPair vp;
				v.spuckShiftVrtcs(vp);

				vrtE1 = vp.first;
				vrtE2 = vp.second;


				Volume * newVol1;
				Volume * newVol2;

				newVol1 = *m_grid.create<Prism>(
								PrismDescriptor( vol1->vertex(0),
										vol1->vertex(1),
										vol1->vertex(2),
										vol1->vertex(3),
										vol1->vertex(4),
										vol1->vertex(5)
									)
									);

				newVol2 = *m_grid.create<Prism>(
								PrismDescriptor( vol2->vertex(0),
										vol2->vertex(1),
										vol2->vertex(2),
										vol2->vertex(3),
										vol2->vertex(4),
										vol2->vertex(5)
									)
									);


				m_sh.assign_subset(newVol1, sudoVols);
				m_sh.assign_subset(newVol2, sudoVols);

				Face * newFac;

				newFac = *m_grid.create<Triangle>(TriangleDescriptor( fac->vertex(0), fac->vertex(1), fac->vertex(2) ));

				m_sh.assign_subset(newFac, sudoFacs);


	//			m_sh.assign_subset(vol1, sudoVols);
	//			m_sh.assign_subset(vol2, sudoVols);

//				m_sh.assign_subset(fac, sudoFacs);


	//			m_sh.assign_subset(fac, sudoFacs);
	//			m_sh.assign_subset(edg1, sudoEdgs);
	//			m_sh.assign_subset(edg2, sudoEdgs);
	//			m_sh.assign_subset(vrtC, sudoVrtx);
	//			m_sh.assign_subset(vrtE1, sudoVrtx);
	//			m_sh.assign_subset(vrtE2, sudoVrtx);

				Edge * edgOne;
				Edge * edgTwo;

				edgOne = *m_grid.create<RegularEdge>(EdgeDescriptor( edgP.first->vertex(0), edgP.first->vertex(1) ));
				edgTwo = *m_grid.create<RegularEdge>(EdgeDescriptor( edgP.second->vertex(0), edgP.second->vertex(1) ));

				m_sh.assign_subset(edgOne,sudoEdgs);
				m_sh.assign_subset(edgTwo,sudoEdgs);

				Vertex * centerV;

				e2bq.spuckCenterVertex(centerV);

				Vertex * newCenter = *m_grid.create<RegularVertex>();

				m_aaPos[newCenter] = m_aaPos[centerV];

				m_sh.assign_subset(newCenter, sudoVrtx);

			}




//		if( d_q == 20 )
//			return true;
//
//		d_q++;
//
		}
	}
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::sortElems2BQuenched()
{
	VecElems2BQuenched ve2bq = m_vecElems2BQuenched;

	for( typename VecElems2BQuenched::iterator itOuter = ve2bq.begin();
			      itOuter != ve2bq.end();
	   )
	{
		Elems2BQuenched e2bqOuter = *itOuter;

		Vertex * centerVrtxOuter;
		VecElems2BQuenched vecElems2Q4ThisVrtx;
		e2bqOuter.spuckOrigCenterVertex( centerVrtxOuter );

		vecElems2Q4ThisVrtx.push_back(e2bqOuter);

		for( typename VecElems2BQuenched::iterator itInner = ve2bq.begin() + 1;
				      itInner != ve2bq.end();
		   )
		{
			Elems2BQuenched e2bqInner = *itInner;

			Vertex * centerVrtxInner;
			e2bqInner.spuckOrigCenterVertex( centerVrtxInner );

			if( centerVrtxOuter == centerVrtxInner )
			{
				vecElems2Q4ThisVrtx.push_back(e2bqInner);
				itInner = ve2bq.erase(itInner);
			}
			else
			{
				itInner++;
			}
		}

		ElemGroupVrtx2BQuenched4Diams egv2b( vecElems2Q4ThisVrtx );

		if( ! egv2b.checkIntegrity() )
		{
			UG_LOG("elem group not integer for vertex" << std::endl);
			return false;
		}

		m_vecElemGroupVrtx2BQuenched.push_back( egv2b );

		itOuter = ve2bq.erase(itOuter);

	}

	return true;
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::attachMarkers()
{
//	m_aAdjMarkerVFP = AttVertFracProp();
//
////	support::VertexFracturePropertiesVol<IndexType> vfp0; // false, 0 );
//	VertxFracPropts vfp0; // false, 0 );
//	// default value: no boundary fracture, no fractures crossing
//
//	m_grid.attach_to_vertices_dv( m_aAdjMarkerVFP, vfp0 );
//	m_aaMarkVrtVFP = Grid::VertexAttachmentAccessor<AttVertFracProp> ( m_grid, m_aAdjMarkerVFP );

	m_attElmGrpVrtx2BQnchd = AttElemGrpVrtx2BQuenchd();

	ElemGroupVrtx2BQuenched4Diams egv2bQEmpty;

	m_grid.attach_to_vertices_dv( m_attElmGrpVrtx2BQnchd, egv2bQEmpty );

	m_attAccsElmGrpVrtx2BQnchd = Grid::VertexAttachmentAccessor<AttElemGrpVrtx2BQuenchd>( m_grid, m_attElmGrpVrtx2BQnchd);

	m_attMarkTwinFace = ABool(); // used to know if an face is twin face

	m_grid.attach_to_faces_dv( m_attMarkTwinFace, false );

	m_attAccsFacIsTwinFac = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_attMarkTwinFace );

	m_attMarkShiftFace = ABool(); // used to know if an face is shift face

	m_grid.attach_to_faces_dv( m_attMarkShiftFace, false );

	m_attAccsFacIsShiftFac = Grid::FaceAttachmentAccessor<ABool>( m_grid, m_attMarkShiftFace );

	m_attMarkShiftTriangleFace = ABool(); // used to know if an face is shift face

	m_grid.attach_to_faces_dv(m_attMarkShiftTriangleFace, false);

	m_attAccsFacIsShiftTriangleFac = Grid::FaceAttachmentAccessor<ABool>(m_grid, m_attMarkShiftTriangleFace);

	m_attMarkShiftQuadriliteralFace = ABool(); // used to know if an face is shift face

	m_grid.attach_to_faces_dv(m_attMarkShiftQuadriliteralFace, false);

	m_attAccsFacIsShiftQuadriliteralFac = Grid::FaceAttachmentAccessor<ABool>(m_grid, m_attMarkShiftQuadriliteralFace);

	m_attMarkVolGetsShrinked = ABool();

	m_grid.attach_to_volumes_dv( m_attMarkVolGetsShrinked, false );

	m_attAccsVolGetsShrinked = Grid::VolumeAttachmentAccessor<ABool>( m_grid, m_attMarkVolGetsShrinked );

	m_attVrtVec = AttVrtVec();

	m_grid.attach_to_volumes_dv(m_attVrtVec, std::vector<Vertex*>());

//	m_attVrtVecFace = AttVrtVec();
//
//	m_grid.attach_to_faces(m_attVrtVecFace, std::vector<Vertex*>());
//
//	m_attAccsVrtVecFace = Grid::VolumeAttachmentAccessor<AttVrtVec>(m_grid, m_attVrtVecFace);

	m_attAccsVrtVecVol = Grid::VolumeAttachmentAccessor<AttVrtVec>(m_grid, m_attVrtVec);

	m_attMarkVrtxIsCenterVrtx = ABool();

	m_grid.attach_to_vertices_dv(m_attMarkVrtxIsCenterVrtx,false);

	m_attAccsVrtxIsCenterVrtx = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_attMarkVrtxIsCenterVrtx);

	m_attMarkVrtxIsShiftVrtx = ABool();

	m_grid.attach_to_vertices_dv(m_attMarkVrtxIsShiftVrtx,false);

	m_attAccsVrtxIsShiftVrtx = Grid::VertexAttachmentAccessor<ABool>( m_grid, m_attMarkVrtxIsShiftVrtx);

	m_attMarkEdgeIsShiftEdge = ABool();

	m_grid.attach_to_edges_dv(m_attMarkEdgeIsShiftEdge,false);

	m_attAccsEdgeIsShiftEdge = Grid::EdgeAttachmentAccessor<ABool>( m_grid, m_attMarkEdgeIsShiftEdge );

	m_attInfoVecSudosTouchingVrtx = AttVecInt();

	m_grid.attach_to_vertices_dv(m_attInfoVecSudosTouchingVrtx, std::vector<IndexType>());

	m_attAccsInfoVecSudosTouchingVrtx = Grid::VertexAttachmentAccessor<AttVecInt>( m_grid, m_attInfoVecSudosTouchingVrtx);

	m_attVrtxIsMidPtOfShiftVrtx = ABool();

	m_grid.attach_to_vertices_dv(m_attVrtxIsMidPtOfShiftVrtx,false);

	m_attAccsVrtxIsMidPtOfShiftVrtcs = Grid::VertexAttachmentAccessor<ABool>(m_grid, m_attVrtxIsMidPtOfShiftVrtx);

	m_attMidPtVrtxOfShiftVrtx = AVertex();

	m_grid.attach_to_vertices_dv(m_attMidPtVrtxOfShiftVrtx, nullptr);

	m_attAccsMidPtVrtxOfShiftVrtx = Grid::VertexAttachmentAccessor<AVertex>( m_grid, m_attMidPtVrtxOfShiftVrtx );

	m_attCenterVrtxOfShiftVrtx = AVertex();

	m_grid.attach_to_vertices_dv( m_attCenterVrtxOfShiftVrtx, nullptr);

	m_attAccsCenterVrtxOfShiftVrtx = Grid::VertexAttachmentAccessor<AVertex>( m_grid, m_attCenterVrtxOfShiftVrtx );

	m_attEdgeCanBeRemoved = ABool();

	m_grid.attach_to_edges_dv(m_attEdgeCanBeRemoved, false);

	m_attAccsEdgeCanBeRemoved = Grid::EdgeAttachmentAccessor<ABool>( m_grid, m_attEdgeCanBeRemoved);

	m_attNewSudoOfVrtx = AInt();

	m_grid.attach_to_vertices_dv( m_attNewSudoOfVrtx, -1 );

	m_attAccsNewSudoOfVrtx = Grid::VertexAttachmentAccessor<AInt>( m_grid, m_attNewSudoOfVrtx );

	return true;
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::detachMarkers()
{
	m_grid.detach_from_vertices(m_attElmGrpVrtx2BQnchd);

	m_grid.detach_from_faces(m_attMarkTwinFace);
	m_grid.detach_from_faces(m_attMarkShiftFace);

	m_grid.detach_from_faces(m_attMarkShiftTriangleFace);
	m_grid.detach_from_faces(m_attMarkShiftQuadriliteralFace);

	m_grid.detach_from_volumes(m_attMarkVolGetsShrinked);

	m_grid.detach_from_volumes( m_attVrtVec );

//	m_grid.detach_from_volumes(m_attVrtVecFace);

	m_grid.detach_from_vertices(m_attMarkVrtxIsCenterVrtx);

	m_grid.detach_from_vertices(m_attMarkVrtxIsShiftVrtx);

	m_grid.detach_from_edges(m_attMarkEdgeIsShiftEdge);

	m_grid.detach_from_vertices(m_attInfoVecSudosTouchingVrtx);

	m_grid.detach_from_vertices(m_attVrtxIsMidPtOfShiftVrtx);

	m_grid.detach_from_vertices(m_attMidPtVrtxOfShiftVrtx);

	m_grid.detach_from_vertices(m_attCenterVrtxOfShiftVrtx);

	m_grid.detach_from_edges(m_attEdgeCanBeRemoved);

	m_grid.detach_from_vertices(m_attNewSudoOfVrtx);

	return true;
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::assignBasicAtts()
{

	for( ElemGroupVrtx2BQuenched4Diams egv2bq : m_vecElemGroupVrtx2BQuenched )
	{
		Vertex * centerVrtx;
		egv2bq.spuckOrigCenterVertex(centerVrtx);

		m_sel.select(centerVrtx);

		m_attAccsElmGrpVrtx2BQnchd[centerVrtx] = egv2bq;
		m_attAccsVrtxIsCenterVrtx[centerVrtx] = true;
		IndxVec sudoList = egv2bq.spuckSudoList();
		m_attAccsInfoVecSudosTouchingVrtx[centerVrtx] = sudoList;

		m_sel.select( m_grid.associated_edges_begin(centerVrtx), m_grid.associated_edges_end(centerVrtx) );
		m_sel.select( m_grid.associated_faces_begin(centerVrtx), m_grid.associated_faces_end(centerVrtx) );
		m_sel.select( m_grid.associated_volumes_begin(centerVrtx), m_grid.associated_volumes_end(centerVrtx) );

		bool sudoAlreadyKnown = generateNewDiamSudos(centerVrtx, sudoList);

	}

	return true;
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::trafoCollectedInfo2Attachments()
{
	//for(VertexIterator iterV = m_sel.vertices_begin(); iterV != m_sel.vertices_end(); iterV++)
	for( ElemGroupVrtx2BQuenched4Diams & egv2bQ : m_vecElemGroupVrtx2BQuenched )
	{
		Vertex * centerV;

		egv2bQ.spuckOrigCenterVertex(centerV);

		VecElems2BQuenched ve2bq;
		egv2bQ.spuckVecElems2BQuenched4Diams(ve2bq);

		for( Elems2BQuenched & e2bq : ve2bq )
		{
			VecVolumeElementFaceQuintuplet vvef5;

			e2bq.spuckVecFullLowDimManifQuintuplet(vvef5);

			for( VolumeElementFaceQuintuplet & vef5 : vvef5 )
			{
				if( ! trafoQuintupleInfo2Attachments(vef5))
				{
					UG_LOG("Quintuples do not get to attachments" << std::endl);
					return false;
				}
			}
//			debugE2bQ(e2bq);

			Vertex * midPtVrtx;
			e2bq.spuckMidPointOfShiftVrtcs(midPtVrtx);
			m_attAccsVrtxIsMidPtOfShiftVrtcs[midPtVrtx] = true;

			m_sel.select(midPtVrtx);

			EdgePair ep;
			e2bq.spuckPairLowDimElem(ep);
			m_sel.select(ep.first);
			m_sel.select(ep.second);

			VrtxPair vp;
			e2bq.spuckShiftVrtcs(vp);
			Vertex * shiftVrtxOne = vp.first;
			Vertex * shiftVrtxTwo = vp.second;
			m_sel.select(shiftVrtxOne);
			m_sel.select(shiftVrtxTwo);

			m_attAccsMidPtVrtxOfShiftVrtx[shiftVrtxOne] = midPtVrtx;
			m_attAccsMidPtVrtxOfShiftVrtx[shiftVrtxTwo] = midPtVrtx;

			m_attAccsCenterVrtxOfShiftVrtx[shiftVrtxOne] = centerV;
			m_attAccsCenterVrtxOfShiftVrtx[shiftVrtxTwo] = centerV;

			// TODO FIXME hier die sudo Liste auch noch an den Center Vertex hängen als attachment

		}
	}



	return true;
}


//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::trafoQuintupleInfo2Attachments(VolumeElementFaceQuintuplet & vef5 )
{
	VrtxPair vp;
	vef5.spuckShiftVrtcs(vp);

	m_attAccsVrtxIsShiftVrtx[vp.first] = true;
	m_attAccsVrtxIsShiftVrtx[vp.second] = true;

	PairVolumeEdgeTwin pvet;
	vef5.spuckPairFullLowDimTwin(pvet);

	VolumeEdgeTwin vetOne = pvet.first;
	VolumeEdgeTwin vetTwo = pvet.second;

	Volume * volOne;
	Volume * volTwo;
	Edge * edgOne;
	Edge * edgTwo;

	vetOne.spuckFullDimElem(volOne);
	vetOne.spuckLowDimElem(edgOne);

	vetTwo.spuckFullDimElem(volTwo);
	vetTwo.spuckLowDimElem(edgTwo);

	m_attAccsVolGetsShrinked[volOne] = true;
	m_attAccsVolGetsShrinked[volTwo] = true;

	m_attAccsEdgeIsShiftEdge[edgOne] = true;
	m_attAccsEdgeIsShiftEdge[edgTwo] = true;

	m_sel.select(volOne);
	m_sel.select(volTwo);

	m_sel.select(edgOne);
	m_sel.select(edgTwo);

	Face * fac;

	vef5.spuckManifElem(fac);

	m_attAccsFacIsTwinFac[fac] = true;

	m_sel.select(fac);

	return true;
}

//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::assignMidPointOfShiftVrtcs(Elems2BQuenched & e2bq)
{
	VrtxPair shiftVrtcs;

	e2bq.spuckShiftVrtcs(shiftVrtcs);

	Vertex * midPtVrtx;

	if( !createMidVrtx(shiftVrtcs, midPtVrtx))
	{
		UG_LOG("mid vertex already created " << std::endl);
//		return false;
	}

	e2bq.assignMidPointOfShiftVrtcs(midPtVrtx);

	return true;
}

//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::createMidVrtx( VrtxPair const & shiftVrtcs, Vertex * & midPtVrtx )
{

	bool pairAlreadyCombined = false;

	for( CombiShiftVrtxMidVrtx & csvmv : m_vecCombiShiftVrtxMidVrtx )
	{
		VrtxPair prNoSwap;

		csvmv.spuckShiftVrtxPair(prNoSwap);

		if( prNoSwap == shiftVrtcs )
		{
			pairAlreadyCombined = true;
			csvmv.spuckSinglVrtx(midPtVrtx);

			break;
		}
		else
		{
			std::swap( prNoSwap.first, prNoSwap.second );

			if( prNoSwap == shiftVrtcs )
			{
				pairAlreadyCombined = true;
				csvmv.spuckSinglVrtx(midPtVrtx);

				break;
			}

			if( pairAlreadyCombined )
				break;

		}
	}

	if( ! pairAlreadyCombined )
	{

		Vertex * vrtxOne;
		Vertex * vrtxTwo;

		vrtxOne = shiftVrtcs.first;
		vrtxTwo = shiftVrtcs.second;


		vector3 posOne = m_aaPos[vrtxOne];
		vector3 posTwo = m_aaPos[vrtxTwo];

		vector3 sum;

		VecAdd(sum,posOne,posTwo);

		vector3 scal;

		VecScale(scal,sum,0.5);

		midPtVrtx = *m_grid.create<RegularVertex>();

		m_aaPos[midPtVrtx] = scal;

		CombiShiftVrtxMidVrtx csvmvNew( shiftVrtcs, midPtVrtx );

		m_vecCombiShiftVrtxMidVrtx.push_back(csvmvNew);
	}

	return (! pairAlreadyCombined);

}

//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::createConditionForNewVrtcs()
{

	//	iterate over all surrounding volumes to enable volume changes, this loop taken from SR but shortened
	for(VolumeIterator iterSurrVol = m_sel.volumes_begin(); iterSurrVol != m_sel.volumes_end(); iterSurrVol++ )
	{
		Volume * sv = *iterSurrVol;

		std::vector<Vertex*> & newVrts = m_attAccsVrtVecVol[sv];
		newVrts.resize(sv->num_vertices());

		for(size_t iVrt = 0; iVrt < sv->num_vertices(); iVrt++ )
		{
			newVrts[iVrt] = nullptr;
		}
			// erstmal so tun, als ob keine neuen Vertizes erzeugt werden an den alten Vertizes
	}

//	for( FaceIterator iterFac = m_sel.faces_begin(); iterFac != m_sel.faces_end(); iterFac++)
//	{
//		Face * fac = *iterFac;
//
//		std::vector<Vertex*> & newVrts = m_attAccsVrtVecFace[fac];
//
//		newVrts.resize(fac->num_vertices());
//
//		for(size_t iVrt = 0; iVrt < sv->num_vertices(); iVrt++ )
//		{
//			newVrts[iVrt] = nullptr;
//		}
//
//	}

	return true;
}


//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::distributeInfosForShrinkingVols()
{
	for( ElemGroupVrtx2BQuenched4Diams & egv2bq :  m_vecElemGroupVrtx2BQuenched )
	{
		VecElems2BQuenched ve2bq;

		egv2bq.spuckVecElems2BQuenched4Diams(ve2bq);

		for( Elems2BQuenched & e2bq : ve2bq )
		{
			VecVolumeElementFaceQuintuplet vvef5;

			e2bq.spuckVecFullLowDimManifQuintuplet(vvef5);

			Vertex * centerVrtx;

			egv2bq.spuckOrigCenterVertex(centerVrtx);

			Vertex * newMidVrtx;

			e2bq.spuckMidPointOfShiftVrtcs(newMidVrtx);

			for( VolumeElementFaceQuintuplet & vef5 : vvef5 )
			{
				PairVolumeEdgeTwin pvet;

				vef5.spuckPairFullLowDimTwin(pvet);

				VolumeEdgeTwin vetOne, vetTwo;

				vetOne = pvet.first;
				vetTwo = pvet.second;

				Volume * volOne;
				Volume * volTwo;

				vetOne.spuckFullDimElem(volOne);
				vetTwo.spuckFullDimElem(volTwo);

				if( !teachMidVrtx2Vol(volOne, centerVrtx, newMidVrtx) || ! teachMidVrtx2Vol(volTwo,centerVrtx,newMidVrtx) )
				{
					UG_LOG("not taughtable" << std::endl);
					return false;
				}
			}
		}
	}

	// TODO FIXME FACES detektieren, die verschoben werden müssen!!!

	return true;
}


//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::teachMidVrtx2Vol( Volume * const & vol, Vertex * const & origVrtx, Vertex * const & midVrtx )
{
	IndexType taught = 0;

	std::vector<Vertex*> & newVrts4Fac = m_attAccsVrtVecVol[ vol ];

	for( IndexType indVrt = 0; indVrt < (vol)->num_vertices();  indVrt++ )
	{
		Vertex* volVrt = (vol)->vertex(indVrt);

		if(  origVrtx == volVrt )
		{
			newVrts4Fac[ indVrt ] = midVrtx;
			taught++;
		}
	}

	if( taught != 1 )
	{
		UG_LOG("strange problem with the volume " << taught << std::endl);
		m_sh.assign_subset(vol, m_sh.num_subsets());
		return false;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::shrinkVolumes()
{
	//	holds local side vertex indices
	std::vector<size_t>	locVrtInds;

	//	first we create new edges from selected ones which are connected to
	//	inner vertices. This allows to preserve old subsets.
	//	Since we have to make sure that we use the right vertices,
	//	we have to iterate over the selected volumes and perform all actions on the edges
	//	of those volumes.

	int volNum = 0;

	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end(); ++iter_sv)
	{
		volNum++;

		Volume* sv = *iter_sv;

		UG_LOG("Diamond entering volume to create new elems " << CalculateCenter(sv, m_aaPos) << std::endl);

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
				Vertex* nv0 = (m_attAccsVrtVecVol[sv])[ind0];
				Vertex* nv1 = (m_attAccsVrtVecVol[sv])[ind1];

				if(nv0 || nv1)
				{
					//	if one vertex has no associated new one, then we use the vertex itself
					if(!nv0)
						nv0 = sv->vertex(ind0);
					if(!nv1)
						nv1 = sv->vertex(ind1);

					//	create the new edge if it not already exists.
					if( ! m_grid.get_edge(nv0, nv1))
						m_grid.create_by_cloning(e, EdgeDescriptor(nv0, nv1), e);
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

	int volNum2 = 0;

	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end(); ++iter_sv)
	{
		volNum2++;

		Volume* sv = *iter_sv;
		//	check for each face whether it has to be copied.

		UG_LOG("Diamond Face descriptor for vol " << CalculateCenter(sv, m_aaPos) << std::endl);

		for(size_t i_face = 0; i_face < sv->num_faces(); ++i_face)
		{
			Face* sf = m_grid.get_face(sv, i_face);

			if( m_sel.is_selected(sf))
			{
				//	check the associated vertices through the volumes aaVrtVecVol attachment.
				//	If no face between the new vertices already exists, we'll create the new face.
				sv->get_vertex_indices_of_face(locVrtInds, i_face);
				fd.set_num_vertices(sf->num_vertices());

				for(size_t i = 0; i < sf->num_vertices(); ++i)
				{
					Vertex* nVrt = (m_attAccsVrtVecVol[sv])[locVrtInds[i]];

					if(nVrt)
						fd.set_vertex(i, nVrt);
					else
						fd.set_vertex(i, sv->vertex(locVrtInds[i]));
				}

				//	if the new face does not already exist, we'll create it
				if(!m_grid.get_face(fd))
					m_grid.create_by_cloning(sf, fd, sf);
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

	VolumeDescriptor vd;

	// beim Volumen fängt das Abfangen an, es muss abgefragt werden, ob das Volumen aus
	// einem Segment ist, wo bisher falsch expandiert wird
	// erst beim Face ab zu fangen ist viel zu spät TODO FIXME
	// nicht zu vergessen, die Edges ordentlich sortiert zu sammeln, die zwei endende Vertizes enthalten,
	// und dann zu splitten, wenn alle Attachements entfernt sind, dann Neustart anfordern mit neuer Geometrie

	int sudoNewVols = m_sh.num_subsets();

//	std::vector<Volume*> newVolsTwoCut;
//	std::vector<Volume*> newVolsThreeCut;


	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end();)
	{

		Volume* sv = *iter_sv;
		++iter_sv;


		UG_LOG("Diamond Volume new creation try at " << CalculateCenter(sv, m_aaPos) << std::endl);


		//	now expand the fracture faces of sv to volumes.
		for(IndexType iSideVol = 0; iSideVol < sv->num_sides(); iSideVol++)
		{
			//	get the local vertex indices of the side of the volume
			sv->get_vertex_indices_of_face(locVrtInds, iSideVol);

			Face* sideFace = m_grid.get_side(sv, iSideVol);

			Volume * shiftVol = nullptr;

			if( sideFace )
			{

				if( m_attAccsFacIsShiftFac[sideFace] )
				{

					if( m_attAccsFacIsShiftQuadriliteralFac[sideFace] )
					{

						std::vector<Vertex*> centerVrtcs;
						std::vector<Vertex*> shiftVrtcs;
						std::vector<Vertex*> midPtVrtcs;


						for( IndexType vrtxInd = 0; vrtxInd < sideFace->num_vertices(); vrtxInd++ )
						{
							Vertex * sideVrtx = sideFace->vertex(vrtxInd);

							if( m_attAccsVrtxIsCenterVrtx[sideVrtx] )
							{
								centerVrtcs.push_back(sideVrtx);
							}
							else if( m_attAccsVrtxIsShiftVrtx[sideVrtx] )
							{
								shiftVrtcs.push_back(sideVrtx);
							}
							else
							{
								UG_LOG("not center not shift but shift face? "
										<< CalculateCenter( sideFace, m_aaPos ) << std::endl );
								return false;
							}
						}

						if( centerVrtcs.size() != 2 || shiftVrtcs.size() != 2 )
						{
							UG_LOG("strange face behaviour for shift face of fac and vol "
										<< CalculateCenter( sideFace, m_aaPos )
										<< " ---- "
										<< CalculateCenter( sv, m_aaPos )
										<< std::endl);
								return false;
						}

						std::swap( shiftVrtcs[0], shiftVrtcs[1] );

						if( ! findShiftFaceVertices( sv, centerVrtcs, midPtVrtcs ))
						{
							UG_LOG("vertices of shift face strange "
									<< CalculateCenter( sideFace, m_aaPos ) << std::endl);
							return false;
						}

						UG_LOG("Diamentenerzeugung fuer " << CalculateCenter( sideFace, m_aaPos ) << std::endl);

//						centers.push_back(centerVrtcs);
//						shifts.push_back(shiftVrtcs);
//						midPts.push_back(midPtVrtcs);

//						// sehr sonderbar, wieso folgende Reihenfolge notwendig, damit es klappt:
//
//						shiftVol = *m_grid.create<Prism>(
//								           PrismDescriptor( centerVrtcs[0], midPtVrtcs[0], shiftVrtcs[1],
//											            	centerVrtcs[1], midPtVrtcs[1], shiftVrtcs[0]
//										       	   	   	  )
//											   	   	   	 );
//						// sehr sonderbar, wieso folgende Reihenfolge notwendig, damit es klappt:

						shiftVol = *m_grid.create<Prism>(
								           PrismDescriptor( centerVrtcs[0], midPtVrtcs[0], shiftVrtcs[0],
											            	centerVrtcs[1], midPtVrtcs[1], shiftVrtcs[1]
										       	   	   	  )
											   	   	   	 );

						int isThreeCross = -1;
						IndexType foundThreeCross = 0;

						VrtxIndxCombi vrtxSudosCombiCenters;

						for( int i = 0; i < 2; i++ )
						{
							Vertex * cntrVrtx = centerVrtcs[i];

							IndxVec sudoVec = m_attAccsInfoVecSudosTouchingVrtx[cntrVrtx];

							IndexType sudoNums = sudoVec.size();

							VrtxIndxPair centerVrtxSudos( cntrVrtx, sudoVec );

							vrtxSudosCombiCenters.push_back(centerVrtxSudos);

							if( sudoNums < 2 )
							{
								UG_LOG("no cross but new prism " << std::endl);
								return false;
							}
							else if( sudoNums == 2 )
							{
								;
							}
							else if( sudoNums == 3 )
							{
								isThreeCross = i;
								foundThreeCross++;
							}
							else
							{
								UG_LOG("strange sudo number " << std::endl);
								return false;
							}
						}

						if( foundThreeCross > 0 )
						{
							if( foundThreeCross != 1 )
							{
								UG_LOG("three cross vertices too much " << std::endl);
								return false;
							}

							CombiNewVolsProps cnvp( shiftVol, sideFace, vrtxSudosCombiCenters, shiftVrtcs, midPtVrtcs, isThreeCross );

							m_vecCombiNewVolsThreeCross.push_back(cnvp);
						}
						else
						{
							CombiNewVolsProps cnvp( shiftVol, sideFace, vrtxSudosCombiCenters, shiftVrtcs, midPtVrtcs );

							m_vecCombiNewVolsTwoCross.push_back(cnvp);
						}

						//m_attAccsInfoVecSudosTouchingVrtx

//						newVols.push_back(shiftVol);
//						numFacs++;

//						if( numFacs == 1 )
//							return true;
//						Volume * shiftVol2 = *m_grid.create<Prism>(
//								           PrismDescriptor(
//											                shiftVrtcs[1], midPtVrtcs[1], centerVrtcs[1],
//															shiftVrtcs[0], midPtVrtcs[0], centerVrtcs[0]
//										       	   	   	  )
//											   	   	   	 );

//						m_sh.assign_subset(shiftVol,sudoNewVols);
//						m_sh.assign_subset(shiftVol2,m_sh.num_subsets()+2);

						UG_LOG("Diamentenerzeugung geschafft " << CalculateCenter( sideFace, m_aaPos ) << std::endl);

//						return true;

					}

//					if( shiftVol && shiftVol2 )
//					{
//						m_sh.assign_subset(shiftVol, sudoNewVols);
//						m_sh.assign_subset(shiftVol2, sudoNewVols+1);
//
//						m_sh.assign_subset( centerVrtcs[0], m_sh.num_subsets() );
//						m_sh.assign_subset( centerVrtcs[1], m_sh.num_subsets() );
//						m_sh.assign_subset( midPtVrtcs[0], m_sh.num_subsets() );
//						m_sh.assign_subset( midPtVrtcs[1], m_sh.num_subsets() );
//						m_sh.assign_subset( shiftVrtcs[0], m_sh.num_subsets() );
//						m_sh.assign_subset( shiftVrtcs[1], m_sh.num_subsets() );
//
//
//						return true;
//					}
				}
			}

//			if( shiftVol )
//			{
//				m_sh.assign_subset(shiftVol,m_sh.num_subsets());
//			}

		}

		//	now set up a new volume descriptor and replace the volume.
		if(vd.num_vertices() != sv->num_vertices())
			vd.set_num_vertices(sv->num_vertices());

		for(size_t i_vrt = 0; i_vrt < sv->num_vertices(); ++i_vrt)
		{
			if( (m_attAccsVrtVecVol[sv])[i_vrt] )
				vd.set_vertex(i_vrt, (m_attAccsVrtVecVol[sv])[i_vrt]);
			else
				vd.set_vertex(i_vrt, sv->vertex(i_vrt));
		}

		m_grid.create_by_cloning(sv, vd, sv);
		m_grid.erase(sv);

	}

//	for( Volume * v : newVols )
//	{
//		m_sh.assign_subset(v,m_sh.num_subsets());
//	}

	UG_LOG("Volumes erzeugt " << std::endl);

	for( EdgeIterator itEdg = m_sel.begin<Edge>(); itEdg != m_sel.end<Edge>(); )
	{
		Edge * edg = *itEdg;
		itEdg++;

		if( m_attAccsEdgeCanBeRemoved[edg] )
		{
			m_grid.erase(edg);
		}
	}

	UG_LOG("unnoetige Ecken entfernt " << std::endl);

	for(FaceIterator iter = m_sel.begin<Face>(); iter != m_sel.end<Face>();)
	{
		Face* fac = *iter;
		++iter;

		if( ! m_attAccsFacIsShiftFac[fac] )
			m_grid.erase(fac);
	}

	UG_LOG("Gesichter entfernt " << std::endl);


//	int suse = m_sh.num_subsets();
//
//	for( int i = 0; i < centers.size(); i++ )
//	{
//		std::vector<Vertex*> centerVrtcs = centers[i];
//		std::vector<Vertex*> shiftVrtcs = shifts[i];
//		std::vector<Vertex*> midPtVrtcs = midPts[i];
//
//		UG_LOG("versuche neues Volumen nr " << i << std::endl);
//
//		Volume * shiftVol = *m_grid.create<Prism>(
//				PrismDescriptor( centerVrtcs[0], shiftVrtcs[0], midPtVrtcs[0],
//						centerVrtcs[1], shiftVrtcs[1], midPtVrtcs[1]
//				)
//		);
//
//		m_sh.assign_subset(shiftVol, suse);
//
//		UG_LOG("neue Volumen " << i << CalculateCenter(shiftVol,m_aaPos) << std::endl);
//
//	}

	UG_LOG("neue Volumen erzeugt " << std::endl);

	return true;
}

//////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::determineShiftFaces()
{

	for(FaceIterator iterF = m_sel.begin<Face>(); iterF != m_sel.end<Face>();)
	{
		Face* fac = *iterF;

		IndexType numShiftVrcs = 0;
		IndexType numCentrVrtcs = 0;

		if( typeid(*fac) == typeid(Triangle) )
		{
			UG_LOG("we have a triangle " << CalculateCenter( fac, m_aaPos ) << " -> " << typeid(*fac).name() << std::endl);
		}
		else if( typeid(*fac) == typeid(Quadrilateral) )
		{
			UG_LOG("we have a quadriliteral " << CalculateCenter( fac, m_aaPos ) << " -> " << typeid(*fac).name() << std::endl);
		}
		else
		{
			UG_LOG("we have whatever at " << CalculateCenter( fac, m_aaPos ) << " -> " << typeid(*fac).name() << std::endl);
		}

		for( IndexType iVrt = 0; iVrt < fac->num_vertices(); iVrt++ )
		{
			Vertex * vrt = fac->vertex(iVrt);

			if( m_attAccsVrtxIsCenterVrtx[vrt] )
			{
				numCentrVrtcs++;
			}
			else if( m_attAccsVrtxIsShiftVrtx[vrt] )
			{
				numShiftVrcs++;
			}
		}

		if( numShiftVrcs == 2 && numCentrVrtcs == 2 )
		{
			m_attAccsFacIsShiftFac[fac] = true;
			m_attAccsFacIsShiftQuadriliteralFac[fac] = true;

			if( typeid(*fac) != typeid(Quadrilateral) )
			{
				UG_LOG("4 vertices but not Quadriliteral? " << CalculateCenter( fac, m_aaPos ) << " -> " << typeid(*fac).name() << std::endl);
			}

		}
		else if( numShiftVrcs == 2 && numCentrVrtcs == 1 )
		{
			m_attAccsFacIsShiftFac[fac] = true;
			m_attAccsFacIsShiftTriangleFac[fac] = true;
			UG_LOG("was ist das für ein Typ " << CalculateCenter( fac, m_aaPos ) << " -> " << typeid(*fac).name() << std::endl );
		}

		if( numCentrVrtcs > 2 || numShiftVrcs > 2 )
		{
			UG_LOG("too much shift or center vertices  in one face " << std::endl);
			UG_LOG("weia was ist das für ein Typ " << CalculateCenter( fac, m_aaPos ) << " -> " << typeid(*fac).name() << std::endl );
			return false;
		}

		iterF++;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////

//bool DiamondsEstablish3D::findShiftFaceVertices(
//		std::vector<Vertex*> & centerVrtcs,
//		std::vector<Vertex*> & shiftVrtcs,
//		std::vector<Vertex*> & midPtVrtcs
//)
//{
//
//	// figure out to which center vertex which shift vertex belongs
//
//	if( centerVrtcs.size() != 2 || shiftVrtcs.size() != 2 )
//	{
//		UG_LOG("strange face behaviour for shift face" << std::endl);
////				<< CalculateCenter( fac, m_aaPos ) << std::endl);
//		return false;
//	}
//
//	if( ! checkAttsOfShiftFaceVrtcs( centerVrtcs, shiftVrtcs) )
//	{
//		std::swap( shiftVrtcs[0], shiftVrtcs[1] );
//
//		if( ! checkAttsOfShiftFaceVrtcs( centerVrtcs, shiftVrtcs) )
//		{
//			UG_LOG("shift and swap vertices do not agree independent of ordering" << std::endl);
//			return false;
//		}
//	}
//
//	for( Vertex * & spv : shiftVrtcs )
//	{
//		midPtVrtcs.push_back( m_attAccsMidPtVrtxOfShiftVrtx[spv] );
//	}
//
//	return true;
//}

////////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::checkAttsOfShiftFaceVrtcs( std::vector<Vertex*> const & centerVrtcs, std::vector<Vertex*> const & shiftVrtcs )
{
	Vertex * const & centerVrtxOne = centerVrtcs[0];
	Vertex * const & centerVrtxTwo = centerVrtcs[1];

	Vertex * const & shiftVrtxOne = shiftVrtcs[0];
	Vertex * const & shiftVrtxTwo = shiftVrtcs[1];

	if( m_attAccsCenterVrtxOfShiftVrtx[shiftVrtxOne] == centerVrtxOne )
	{
		if( m_attAccsCenterVrtxOfShiftVrtx[shiftVrtxTwo] == centerVrtxTwo )
		{
			return true;
		}
		else
		{
//			UG_LOG("shift and center vertex do not fit together " << std::endl);
			return false;
		}
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::findShiftFaceVertices( Volume * & vol,
  	  	    std::vector<Vertex*> & centerVrtcs,
		std::vector<Vertex*> & midPtVrtcs
		 )
{
	midPtVrtcs = std::vector<Vertex*>(2);

	std::vector<Vertex*> & newVrts4Fac = m_attAccsVrtVecVol[ vol ];

	for( IndexType indVrt = 0; indVrt < (vol)->num_vertices();  indVrt++ )
	{
		Vertex* volVrt = (vol)->vertex(indVrt);

		for( int i = 0; i < 2; i++ )
		{
			Vertex * ceVe = centerVrtcs[i];

			if(  ceVe == volVrt )
			{
				midPtVrtcs[i] = newVrts4Fac[indVrt];
			}
		}
	}


	return true;
}

////////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::postprocessNewDiamVols()
{
	for( CombiNewVolsProps & cnvp : m_vecCombiNewVolsTwoCross )
	{
		VrtxIndxCombi viCombiCenter;
		cnvp.spuckCenterVrtcsSudos(viCombiCenter);
		Volume * vol;
		cnvp.spuckFulldimElem(vol);

		VrtxIndxPair vrtxIndxPr( viCombiCenter[0] );

		Vertex * center = vrtxIndxPr.first;
		IndxVec sudoVec = vrtxIndxPr.second;

		for( CombiCntrVrtxSudo & ccvs: m_vecCombiCntrVrtxSudo )
		{
			if( ccvs.spuckSudoVec() == sudoVec )
			{
//				m_sh.assign_subset(vol, ccvs.spuckNewSudo());
				if( ! assignSudoOfNewVols2VolAndSubElems(vol,ccvs.spuckNewSudo() ) )
				{
					UG_LOG("sudos not assignable" << std::endl);
					return false;
				}
				break;
			}
		}

	}

	for( CombiNewVolsProps & cnvp : m_vecCombiNewVolsThreeCross )
	{
		VrtxIndxCombi viCombiCenter;
		cnvp.spuckCenterVrtcsSudos(viCombiCenter);
		Volume * vol;
		cnvp.spuckFulldimElem(vol);

		IndexType i3c = cnvp.spuckThreeCrossIndex();

		VrtxIndxPair vrtxIndxPr( viCombiCenter[i3c] );

		Vertex * center = vrtxIndxPr.first;
		IndxVec sudoVec = vrtxIndxPr.second;

		for( CombiCntrVrtxSudo & ccvs: m_vecCombiCntrVrtxSudo )
		{
			if( ccvs.spuckSudoVec() == sudoVec )
			{
//				m_sh.assign_subset(vol, ccvs.spuckNewSudo());
				if( ! assignSudoOfNewVols2VolAndSubElems(vol,ccvs.spuckNewSudo() ) )
				{
					UG_LOG("sudos not assignable" << std::endl);
					return false;
				}
				break;
			}
		}

	}

	UG_LOG("try to split diam vols " << std::endl);

	for( CombiNewVolsProps & cnvp : m_vecCombiNewVolsThreeCross )
	{
		if( ! splitThreeCrossLargeDiams(cnvp) )
		{
			UG_LOG("splitting not possible");
			return false;
		}
	}

	for( Face * fac :  m_faces2BDeletedAtLastStep )
	{
		if( ! fac )
		{
			UG_LOG("want to delete null fac " << std::endl);
			return false;
		}

		m_grid.erase(fac);
	}

	for( Edge * edg : m_edges2BDeletedAtLastStep )
	{
		if( ! edg )
		{
			UG_LOG("want to delete edge unexisting at end " << std::endl);
			return false;
		}

		m_grid.erase(edg);
	}


	UG_LOG("all diam vols splitted " << std::endl);

	return true;
}

////////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::generateNewDiamSudos(Vertex * & centerV, IndxVec sudoList )
{
	bool sudoCombiUnKnown = addElem(m_sudosTable, sudoList);

	if( ! sudoCombiUnKnown )
	{
		for( CombiCntrVrtxSudo & ccvs : m_vecCombiCntrVrtxSudo )
		{
			if( ccvs.spuckSudoVec() == sudoList )
			{
				ccvs.schluckVertex(centerV);
				m_sh.assign_subset(centerV, ccvs.spuckNewSudo());
				m_attAccsNewSudoOfVrtx[centerV] = ccvs.spuckNewSudo();
				break;
			}
		}
	}
	else
	{
		IndexType newSudo = m_sh.num_subsets();
		m_sh.assign_subset(centerV,newSudo);

		std::string sudoName = std::string("diamond_");
		for( IndexType & sd : sudoList )
		{
			sudoName += std::string("_") + std::string( const_cast<char*>( m_sh.get_subset_name( sd ) ) );
		}

		m_sh.set_subset_name(sudoName.c_str(), newSudo);

		m_attAccsNewSudoOfVrtx[centerV] = newSudo;
		CombiCntrVrtxSudo ccvs( sudoList, newSudo );
		m_vecCombiCntrVrtxSudo.push_back(ccvs);
		UG_LOG("NNNNNNNNNNNNNNNNNNNN" << std::endl);
		UG_LOG("creating new sudo of " << std::endl);
		for( auto & i : sudoList )
		{
			UG_LOG("list part " << i << std::endl);
		}
		UG_LOG("DDDDDDDDDDDDDDDDDDDD" << std::endl);

	}

	return ( sudoCombiUnKnown );
}

////////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::assignSudoOfNewVols2VolAndSubElems(Volume * & vol, IndexType sudo)
{
	if( ! vol )
	{
		UG_LOG("vol to assign sudo null " << std::endl);
		return false;
	}

	m_sh.assign_subset(vol, sudo);

	for(IndexType iFace = 0; iFace < vol->num_faces(); ++iFace)
	{
		Face * fac = m_grid.get_face(vol, iFace);

		if( ! fac )
		{
			UG_LOG("face null sudo " << std::endl);
			return false;
		}

		m_sh.assign_subset( fac, sudo );

	}

	for(IndexType iEdge = 0; iEdge < vol->num_edges(); ++iEdge)
	{
		Edge* edg = m_grid.get_edge(vol, iEdge);

		if( ! edg )
		{
			UG_LOG("edge null sudo " << std::endl);
			return false;
		}

		m_sh.assign_subset( edg, sudo );

	}

	for( IndexType iVrt = 0; iVrt < vol->num_vertices(); iVrt++ )
	{
		Vertex * vrt = vol->vertex(iVrt);

		if( !vrt )
		{
			UG_LOG("vertex null sudo " << std::endl);
			return false;
		}

		m_sh.assign_subset( vrt, sudo );
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::detectRemovableEdges()
{

	for( FaceIterator itFac = m_sel.begin<Face>(); itFac != m_sel.end<Face>(); itFac++)
	{
		Face * fac = *itFac;

		if( ! fac )
		{
			UG_LOG("face does not exist " << std::endl);
			return false;
		}

		UG_LOG("looking into twin facs try " << std::endl);

		if( m_attAccsFacIsTwinFac[fac] )
		{
			UG_LOG("inside twin facs try " << std::endl);

			for( IndexType iEdg = 0; iEdg < fac->num_edges(); iEdg++)
			{
				UG_LOG("inside edges of twins try " << std::endl);

				Edge * edg = m_grid.get_edge(fac,iEdg);

				if( ! edg )
				{
					UG_LOG("edge does not exist " << std::endl);
					return false;
				}

				IndexType numAssoCenterVrcs = 0;

				for( IndexType iV = 0; iV < 2; iV++ )
				{
					Vertex * vrt = edg->vertex(iV);

					if( ! vrt )
					{
						UG_LOG("vertex does not exist " << std::endl);
						return false;
					}

					if( m_attAccsVrtxIsCenterVrtx[vrt] )
					{
						numAssoCenterVrcs++;
					}
				}

				UG_LOG("number of associatec center vertices is " << numAssoCenterVrcs << std::endl);

				if( numAssoCenterVrcs == 1 )
				{
					m_attAccsEdgeCanBeRemoved[edg] = true;
				}

			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::splitThreeCrossLargeDiams( CombiNewVolsProps & combiNewVolsProps )
{
	VrtxIndxCombi viCombiCenter;
	combiNewVolsProps.spuckCenterVrtcsSudos(viCombiCenter);
	Volume * vol;
	combiNewVolsProps.spuckFulldimElem(vol);

	IndexType ind3CrossSide = combiNewVolsProps.spuckThreeCrossIndex();

	IndexType indOtherSide = ( ind3CrossSide + 1 ) % 2;

	VrtxIndxPair diamCenterVrtxPr( viCombiCenter[ind3CrossSide] );

	Vertex * centerDiam = diamCenterVrtxPr.first;

	VrtxIndxPair diamOtherVrtxPr( viCombiCenter[ indOtherSide ] );

	Vertex * centerOther = diamOtherVrtxPr.first;

	std::vector<Vertex*> midPtVrtcs;
	combiNewVolsProps.spuckMidPtVrtcs(midPtVrtcs);

	Vertex * midPtDiam = midPtVrtcs[ind3CrossSide];
	Vertex * midPtOther = midPtVrtcs[indOtherSide];

	vector3 centerDiamVec3 = m_aaPos[centerDiam];
	vector3 midPtDiamVec3 = m_aaPos[midPtDiam];
	vector3 centerOtherVec3 = m_aaPos[centerOther];

	vector3 cutPtVec3;

	DropAPerpendicular(cutPtVec3, midPtDiamVec3, centerDiamVec3, centerOtherVec3);

	Vertex * cutVrtx = *m_grid.create<RegularVertex>();

	m_aaPos[cutVrtx] = cutPtVec3;

	combiNewVolsProps.schluckNewSplitVrtx(cutVrtx);

	std::vector<Vertex*> shiftVrtcs;
	combiNewVolsProps.spuckShiftVrtcs(shiftVrtcs);

	Vertex * shiftVrtxDiam = shiftVrtcs[ind3CrossSide];
	Vertex * shiftVrtxOther = shiftVrtcs[ indOtherSide];

	Volume * splitVolPrism = *m_grid.create<Prism>(
			           	   	   PrismDescriptor( cutVrtx, midPtDiam, shiftVrtxDiam,
						            			centerOther, midPtOther, shiftVrtxOther
					       	   	   	  )
						   	   	   	 );

	assignSudoOfNewVols2VolAndSubElems(splitVolPrism, m_attAccsNewSudoOfVrtx[centerOther]);
//	m_sh.assign_subset( splitVolPrism, m_attAccsNewSudoOfVrtx[centerOther] );

	Volume * splitVolTetra = *m_grid.create<Tetrahedron>(
							TetrahedronDescriptor( shiftVrtxDiam, midPtDiam, cutVrtx, centerDiam  ) );

//	m_sh.assign_subset( splitVolTetra,  m_attAccsNewSudoOfVrtx[centerDiam] );
	assignSudoOfNewVols2VolAndSubElems(splitVolTetra, m_attAccsNewSudoOfVrtx[centerDiam]);

	// figure out the faces and the edges which later should be erased

	IndexType numEdgesFound = 0;

	Edge* centersEdg = nullptr;

	for( IndexType iEdg = 0; iEdg < vol->num_edges(); iEdg++)
	{
		Edge * edg = m_grid.get_edge(vol, iEdg);

		if( EdgeContains(edg, centerDiam ) && EdgeContains(edg, centerOther))
		{
			centersEdg = edg;
			numEdgesFound++;
		}
	}

	if( numEdgesFound != 1 )
	{
		UG_LOG("edges found at end " << numEdgesFound << std::endl);
		return false;
	}

	addElem(m_edges2BDeletedAtLastStep,centersEdg);

	std::vector<Face*> faces2Del;

	for( IndexType iFac = 0; iFac < vol->num_faces(); iFac++)
	{
		Face * fac = m_grid.get_face(vol, iFac);

		if( FaceContains(fac, centersEdg) )
		{
			faces2Del.push_back(fac);
		}
	}

	if( faces2Del.size() != 2 )
	{
		UG_LOG("strange number of faces to del at fin " << faces2Del.size() << std::endl);
		return false;
	}

	for( Face * f : faces2Del )
	{
		addElem(m_faces2BDeletedAtLastStep, f);
	}

	m_grid.erase(vol);


	return true;
}


///////////////////////////////////////////////////////////////////////////////////


} /* namespace diamonds */

} /* namespace arte */

} /* namespace ug */
