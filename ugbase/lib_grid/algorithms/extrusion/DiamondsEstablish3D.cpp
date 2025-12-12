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

//bool DiamondsEstablish3D::findRegions2BShrinked()
//{
//	UG_LOG("want to find regions to be shrinked" << std::endl);
//
//	IndexType sudoNum = m_sh.num_subsets();
//
//	IndexType sudoVols = sudoNum;
//	IndexType sudoFacs = sudoNum+1;
//	IndexType sudoEdgs = sudoNum+2;
//	IndexType sudoVrtx = sudoNum+3;
//
//	// collect the pairs of volumes connected by a face relevant for the shrinking of the volumes
//
////	for( auto & vmvcd : m_vecVolManifVrtxCombiToShrink4Diams )
////	{
//////		Elems2BQuenched e2bq;
////		VrtxPair oldAndShiftV;
////		vmvcd.spuckOldAndShiftVrtx( oldAndShiftV);
////		Vertex * oldCenter = oldAndShiftV.first;
////	}
//
//	VecVolManifVrtxCombi vecVolManifVrtxCopy = m_vecVolManifVrtxCombiToShrink4Diams;
//
//	int d_itE = 0;
//
//	for( typename VecVolManifVrtxCombi::iterator itVMVOuter = vecVolManifVrtxCopy.begin();
//				  itVMVOuter != vecVolManifVrtxCopy.end(); )
//	{
////		UG_LOG("entering E " << d_itE << std::endl);
//
//		VolManifVrtxCombi outer = *itVMVOuter;
//
//		VrtxPair oldAndShiftVrtxOuter;
//		outer.spuckOldAndShiftVrtx( oldAndShiftVrtxOuter );
//
//		Vertex * oldVrtxOuter = oldAndShiftVrtxOuter.first;
//		Vertex * shiftVrtxOuter = oldAndShiftVrtxOuter.second;
//
//		Face * faceOuter;
//		outer.spuckManif(faceOuter);
//
//		IndexType sudoOuter = outer.spuckSudo();
//
////		IndexType twin2OuterFound = 0;
////		IndexType partner2OuterFound = 0;
////		IndexType twin2InnerFound = 0;
//
//		// TODO FIXME von wegen Twin...... und sudo auch sehr wichtig!!!!
//		VecVolManifVrtxCombi twinCombi2Outer;
//		// should get length 2 in the procedure
//
//		twinCombi2Outer.push_back(outer); // in each case part a twin has to exist, sharing same face
//
//		// Quatsch: should get length 0 or 2 in the procedure: 0 in case of 3 cross, 2 in case of 2 cross
//		VecVolManifVrtxCombi partnerCombi2Outer;
//
////		bool partnerFound = false;
//
////		UG_LOG("entering inner " << d_itE << std::endl);
//
//		int d_itI = 0;
//
//		for( typename VecVolManifVrtxCombi::iterator itVMVInner = itVMVOuter + 1;
//					  itVMVInner != vecVolManifVrtxCopy.end(); )
//		{
////			UG_LOG("in the inner " << d_itI << std::endl);
//
//			VolManifVrtxCombi inner = *itVMVInner;
//
//			VrtxPair oldAndShiftVrtxInner;
//			inner.spuckOldAndShiftVrtx( oldAndShiftVrtxInner );
//
//			Vertex * oldVrtxInner = oldAndShiftVrtxInner.first;
////			Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;
//
//			if( oldVrtxInner == oldVrtxOuter )
//			{
//				Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;
//
//				if( shiftVrtxInner != shiftVrtxOuter )
//				{
//					// if connected by a face, twin to outer found, partner may exist or not
////					 twin2OuterFound++;
//					// check if the face is identical, else is from the partner twin
//
//					Face * faceInner;
//					inner.spuckManif(faceInner);
//
//					if( faceInner == faceOuter )
//					{
//						// twin found
//						twinCombi2Outer.push_back(inner);
//
////						twinFound = true;
//
//						itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);
//
//					}
//					// else can be any combination of a completely different side of the crossing point
//					else
//					{
//						itVMVInner++;
//					}
//				}
//				else if( shiftVrtxInner == shiftVrtxOuter )
//				{
//					// CHECK IF SUDO COINCIDE; ELSE NOT SUITABLE
//
//					IndexType sudoInner = inner.spuckSudo();
//
//					if( sudoInner == sudoOuter )
//					{
//						// side partner found, for this side partner, the twin is needed as well!
//	//					partner2OuterFound++;
//						partnerCombi2Outer.push_back(inner);
//						// will be able to find only one partner in this loop,
//						//its twin has to be found in an additional inner loop!
//						itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);
//
//					}
//					else
//					{
//						itVMVInner++;
//					}
//
//				}
//
//			}
//			else
//			{
//				// else: do nothing, as no relevant relation between the two objects
//				++itVMVInner;
//			}
//
////			UG_LOG("end one inner " << d_itI << std::endl);
//
//			d_itI++;
//		}
//
//	//		UG_LOG("left inner " << d_itE << std::endl);
//
//
//		if( partnerCombi2Outer.size() > 0 )
//		{
//			// additional inner loop to find the twin of the partner
//
//			for( VolManifVrtxCombi & partnerOuter : partnerCombi2Outer )
//			{
//	//			VolManifVrtxCombi partnerOuter = partnerCombi2Outer[0];
//
//				VrtxPair oldAndShiftVrtxPartnerOuter;
//				outer.spuckOldAndShiftVrtx( oldAndShiftVrtxPartnerOuter );
//
//				Vertex * oldVrtxPartnerOuter = oldAndShiftVrtxPartnerOuter.first;
//				Vertex * shiftVrtxPartnerOuter = oldAndShiftVrtxPartnerOuter.second;
//
//				Face * facePartnerOuter;
//				partnerOuter.spuckManif(facePartnerOuter);
//
//				IndexType sudoPartnerOuter = partnerOuter.spuckSudo();
//
//				for( typename VecVolManifVrtxCombi::iterator itVMVInner = itVMVOuter + 1;
//							  itVMVInner != vecVolManifVrtxCopy.end(); )
//				{
//					VolManifVrtxCombi inner = *itVMVInner;
//
//					VrtxPair oldAndShiftVrtxInner;
//					inner.spuckOldAndShiftVrtx( oldAndShiftVrtxInner );
//
//					Vertex * oldVrtxInner = oldAndShiftVrtxInner.first;
//	//				Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;
//
//					if( oldVrtxInner == oldVrtxPartnerOuter )
//					{
//						Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;
//
//						if( shiftVrtxInner != shiftVrtxPartnerOuter )
//						{
//							// if connected by a face, twin to outer found, partner may exist or not
//		//					 twin2OuterFound++;
//							// check if the face is identical, else is from the partner twin
//
//							Face * faceInner;
//							inner.spuckManif(faceInner);
//
//							if( faceInner == facePartnerOuter )
//							{
//								// twin found
//								partnerCombi2Outer.push_back(inner);
//
//		//						twinFound = true;
//
//								itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);
//
//							}
//							else
//							{
//								itVMVInner++;
//							}
//							// else can be any combination of a completely different side of the crossing point
//						}
//						else if( shiftVrtxInner == shiftVrtxOuter )
//						{
////							UG_LOG("too much partners or twins " << std::endl);
//							// darf nicht mehr passieren
//	//						UG_THROW("");
//							// side partner found, for this side partner, the twin is needed as well!
//		//					partner2OuterFound++;
//	//						partnerCombi2Outer.push_back(inner);
//							// will be able to find only one partner in this loop,
//							//its twin has to be found in an additional inner loop!
//
//							IndexType sudoInner = inner.spuckSudo();
//
//							if( sudoInner == sudoPartnerOuter )
//							{
//								partnerCombi2Outer.push_back(inner);
//
//								itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);
//
//							}
//							else
//							{
//								itVMVInner++;
//							}
//						}
//
//					}
//					else
//					{
//						// else: do nothing, as no relevant relation between the two objects
//						++itVMVInner;
//					}
//
//
//				}
//			}
//
//
//		}
////		else if( partnerCombi2Outer.size() != 0 )
////		{
////			UG_LOG("strange partner size " << partnerCombi2Outer.size() << std::endl);
////			return false;
////		}
//
//		if( twinCombi2Outer.size() %2 != 0 )
//		{
//			UG_LOG("strange twin size " << std::endl);
//			return false;
//		}
//
//		if( partnerCombi2Outer.size() % 2 != 0 )
//		{
//			UG_LOG("no partner found " << std::endl);
//			return false;
//		}
//
//		itVMVOuter = vecVolManifVrtxCopy.erase(itVMVOuter);
//
//		// TODO FIXME
//		// FullLowDimManifQntpl needs to be created from each of the two or four contributions
//		// something like transform VecVolManifVrtxCombi to VecVolManifVrtxCombi
//		// VecFullLowDimManifQuintuplet needs to be created
//		// to check from outside if the central vertex is kept!!!
//		// create the new combination and attach it to the vector
//		//		Elems2BQuenched
////		m_vecElems2BQuenched.push_back(...);
//
//		VolumeElementFaceQuintuplet vef5Twin;
//
//		if( ! trafoVolFacVrtxCombiPair2FullLowDimManifQuintuplet( twinCombi2Outer, vef5Twin, false ) )
//		{
//			UG_LOG("trafo failed twin " << std::endl);
//			return false;
//		}
//
//
//		if( ! vef5Twin.checkIntegrity())
//		{
//			UG_LOG("twins not integer " << std::endl);
//			return false;
//		}
//
//		VecVolumeElementFaceQuintuplet vecvef5;
//
//		vecvef5.push_back(vef5Twin);
//
//		if( partnerCombi2Outer.size() == 0 )
//		{
//			;
//		}
//		else if( partnerCombi2Outer.size() == 2 )
//		{
//			VolumeElementFaceQuintuplet vef5Partner;
//
//			if( ! trafoVolFacVrtxCombiPair2FullLowDimManifQuintuplet( partnerCombi2Outer, vef5Partner, true ) )
//			{
//				UG_LOG("trafo failed partner " << std::endl);
//				return false;
//			}
//
//			if( ! vef5Partner.checkIntegrity() )
//			{
//				UG_LOG("partner not integer " << std::endl);
//				return false;
//			}
//
//			vecvef5.push_back(vef5Partner);
//		}
//		else
//		{
//			UG_LOG("strange parter size " << std::endl);
//			return false;
//		}
//
//		for( auto & v: vecvef5 )
//		{
//			Volume * vol1;
//			Volume * vol2;
//			Face * fac;
//			Edge * edg1;
//			Edge * edg2;
//			Vertex * vrtC;
//			Vertex * vrtE1;
//			Vertex * vrtE2;
//
//			std::pair<VolumeElementTwin,VolumeElementTwin> pvv;
//
//			v.spuckPairFullLowDimTwin(pvv);
//
//			pvv.first.spuckFullDimElem(vol1);
//			pvv.second.spuckFullDimElem(vol2);
//
//			v.spuckManifElem( fac );
//
//			pvv.first.spuckLowDimElem(edg1);
//			pvv.second.spuckLowDimElem(edg2);
//
//			v.spuckCenterVertex(vrtC);
//
//			VrtxPair vp;
//			v.spuckShiftVrtcs(vp);
//
//			vrtE1 = vp.first;
//			vrtE2 = vp.second;
//
//			m_sh.assign_subset(vol1, sudoVols);
//			m_sh.assign_subset(vol2, sudoVols);
//			m_sh.assign_subset(fac, sudoFacs);
//			m_sh.assign_subset(edg1, sudoEdgs);
//			m_sh.assign_subset(edg2, sudoEdgs);
//			m_sh.assign_subset(vrtC, sudoVrtx);
//			m_sh.assign_subset(vrtE1, sudoVrtx);
//			m_sh.assign_subset(vrtE2, sudoVrtx);
//
//		}
//
//		Elems2BQuenched elem2BQuenched;
//
//		if( ! establishElems2BeQuenched( vecvef5, elem2BQuenched ) )
//		{
//			UG_LOG("establish quench failed " << std::endl);
//			return false;
//		}
//
//#if 0
//		if( ! elem2BQuenched.checkIntegrity())
//		{
//			UG_LOG("elems to be quenched not integer " << std::endl);
//			return false;
//		}
//
//		m_vecElems2BQuenched.push_back(elem2BQuenched);
//#endif
//
//		d_itE++;
//	}
//
//	UG_LOG("found regions to be shrinked" << std::endl);
//
//
//	return true;
//}

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

	mvcOne.spuckLowdimElem( edgeOne );
	mvcTwo.spuckLowdimElem( edgeTwo );

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

bool DiamondsEstablish3D::establishElems2BeQuenched( VecVolumeElementFaceQuintuplet & vvef5,
													 Elems2BQuenched & elem2BQuenched )
{



	return true;
}
////////////////////////////////////////////////////////////////////////////


} /* namespace diamonds */

} /* namespace arte */

} /* namespace ug */
