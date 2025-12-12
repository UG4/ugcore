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
			m_vecVolManifVrtxCombiToShrink4Diams(vecVolManifVrtxC),
			m_vecVolElmFac5(VecVolumeElementFaceQuintuplet()),
			m_vecElems2BQuenched(VecElems2BQuenched()),
			m_disappearingVols(std::vector<Volume*>()),
			m_disappearingFacs(std::vector<Face*>())
//			m_disappearingEdgs(std::vector<Edge*>())
{
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


bool DiamondsEstablish3D::createTheDiamonds()
{
	if( ! figureOutTheEdges() )
	{
		UG_LOG("Edges not found " << std::endl);
		return false;
	}

	if( ! findRegions2BShrinked())
	{
		UG_LOG("Regions to be shrinked not found " << std::endl);
		return false;
	}

	if( ! establishElems2BeQuenched() )
	{
		UG_LOG("quenching elements not establishable " << std::endl);
		return false;
	}

	UG_LOG("Established diamonds" << std::endl);


	return true;
}

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

		UG_LOG("out round " << d_out << std::endl);

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
			UG_LOG("start in " << d_in << std::endl);

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


			UG_LOG("end in " << d_in << std::endl);
			d_in++;
		}

		UG_LOG("first round finished " << std::endl);
		if( ! partnerFound )
		{
			UG_LOG("no partner found " << std::endl);
			return false;
		}
		else
		{
			itVMVOuter = vecVolManifVrtxCopy.erase(itVMVOuter);
		}

		UG_LOG("end out " << d_out << std::endl);
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

		m_vecElems2BQuenched.push_back(elem2BQuenched);

		itOuter = vvef5.erase(itOuter);
	}

//	for( Elems2BQuenched & e2bq : m_vecElems2BQuenched )
//	{
//		IndexType sudoNum = m_sh.num_subsets();
//
//		IndexType sudoVols = sudoNum;
//		IndexType sudoFacs = sudoNum+1;
//		IndexType sudoEdgs = sudoNum+2;
//		IndexType sudoVrtx = sudoNum+3;
//
//
//
//	}


	return true;
}
////////////////////////////////////////////////////////////////////////////


} /* namespace diamonds */

} /* namespace arte */

} /* namespace ug */
