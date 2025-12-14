/*
 * DiamondsEstablish3D.cpp
 *
 *  Created on: 08.12.2025
 *      Author: Markus Knodel
 */

#include <lib_grid/algorithms/extrusion/DiamondsEstablish3D.h>

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
			m_attMarkVolGetsShrinked(ABool()),
			m_attAccsVolGetsShrinked(Grid::VolumeAttachmentAccessor<ABool>()),
			m_attVrtVec(AttVrtVec()),
			m_attAccsVrtVecVol(Grid::VolumeAttachmentAccessor<AttVrtVec>()),
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
			m_vecCombiShiftVrtxMidVrtx(VecCombiShiftVrtxMidVrtx())
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

	if( ! shrinkVolumes())
	{
		UG_LOG("shrinken schief gegangen " << std::endl);
		return false;
	}

	UG_LOG("volumes shrinked " << std::endl);

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

	m_attMarkVolGetsShrinked = ABool();

	m_grid.attach_to_volumes_dv( m_attMarkVolGetsShrinked, false );

	m_attAccsVolGetsShrinked = Grid::VolumeAttachmentAccessor<ABool>( m_grid, m_attMarkVolGetsShrinked );

	m_attVrtVec = AttVrtVec();

	m_grid.attach_to_volumes_dv(m_attVrtVec, std::vector<Vertex*>());

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

	return true;
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::detachMarkers()
{
	m_grid.detach_from_vertices(m_attElmGrpVrtx2BQnchd);

	m_grid.detach_from_faces(m_attMarkTwinFace);
	m_grid.detach_from_faces(m_attMarkShiftFace);

	m_grid.detach_from_volumes(m_attMarkVolGetsShrinked);

	m_grid.detach_from_volumes( m_attVrtVec );

	m_grid.detach_from_vertices(m_attMarkVrtxIsCenterVrtx);

	m_grid.detach_from_vertices(m_attMarkVrtxIsShiftVrtx);

	m_grid.detach_from_edges(m_attMarkEdgeIsShiftEdge);

	m_grid.detach_from_vertices(m_attInfoVecSudosTouchingVrtx);

	m_grid.detach_from_vertices(m_attVrtxIsMidPtOfShiftVrtx);

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
		m_attAccsInfoVecSudosTouchingVrtx[centerVrtx] = egv2bq.spuckSudoList();

		m_sel.select( m_grid.associated_edges_begin(centerVrtx), m_grid.associated_edges_end(centerVrtx) );
		m_sel.select( m_grid.associated_faces_begin(centerVrtx), m_grid.associated_faces_end(centerVrtx) );
		m_sel.select( m_grid.associated_volumes_begin(centerVrtx), m_grid.associated_volumes_end(centerVrtx) );

	}

	return true;
}

////////////////////////////////////////////////////////////////////////////

bool DiamondsEstablish3D::trafoCollectedInfo2Attachments()
{
	//for(VertexIterator iterV = m_sel.vertices_begin(); iterV != m_sel.vertices_end(); iterV++)
	for( ElemGroupVrtx2BQuenched4Diams & egv2bQ : m_vecElemGroupVrtx2BQuenched )
	{
//		Vertex * centerV = *iterV;
//
//		ElemGroupVrtx2BQuenched4Diams egv2bQ = m_attAccsElmGrpVrtx2BQnchd[centerV];
//
//		Vertex * testVrtx;
//		egv2bQ.spuckOrigCenterVertex(testVrtx);
//
//		if( testVrtx != nullptr )
//		{
//			UG_LOG("found test vertex not null " << m_aaPos[testVrtx] << std::endl);
//
//			if( centerV != nullptr )
//				UG_LOG("also center not null " << m_aaPos[centerV] << std::endl);
//		}
//
//
//		if( testVrtx == nullptr )
//		{
//			UG_LOG("test vertex null " << std::endl);
//
//			if( centerV != nullptr )
//				UG_LOG("for center " << m_aaPos[centerV] << std::endl);
//
//			continue;
//		}
//
//		if( testVrtx != centerV )
//		{
//			UG_LOG("transfer did not work " << std::endl);
//
//			UG_LOG("different locations at " << m_aaPos[testVrtx] << " and " << m_aaPos[centerV] << std::endl);
//
//			return false;
//		}

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
			m_sel.select(vp.first);
			m_sel.select(vp.second);


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

	int volNum3 = 0;

	for(VolumeIterator iter_sv = m_sel.volumes_begin(); iter_sv != m_sel.volumes_end();)
	{
		volNum3++;

		Volume* sv = *iter_sv;
		++iter_sv;

		if( volNum3 == volNum2+1 || volNum3 == volNum+1 )
		{
			m_sh.assign_subset(sv,m_sh.num_subsets());
			UG_LOG("Abbrechen notwendig bei " << CalculateCenter(sv, m_aaPos) << std::endl);
//			break;
		}

		UG_LOG("Diamond Volume new creation try at " << CalculateCenter(sv, m_aaPos) << std::endl);

		//	now expand the fracture faces of sv to volumes.
//		for(size_t i_side = 0; i_side < sv->num_sides(); ++i_side)
//		{
//			//	get the local vertex indices of the side of the volume
//			sv->get_vertex_indices_of_face(locVrtInds, i_side);
//
//			Face* tFace = m_grid.get_side(sv, i_side);
//
//			if(tFace)
//			{}
//		}

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



	UG_LOG("Volumes erzeugt " << std::endl);


//	for(FaceIterator iter = m_sel.begin<Face>(); iter != m_sel.end<Face>();)
//	{
//		Face* f = *iter;
//		++iter;
//
//		m_grid.erase(f);
//	}

	UG_LOG("Gesichter entfernt " << std::endl);

	return true;
}

//////////////////////////////////////////////////////////////////////////////////




} /* namespace diamonds */

} /* namespace arte */

} /* namespace ug */
