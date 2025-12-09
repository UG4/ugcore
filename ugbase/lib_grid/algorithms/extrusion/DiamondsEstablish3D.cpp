/*
 * DiamondsEstablish3D.cpp
 *
 *  Created on: 08.12.2025
 *      Author: Markus Knodel
 */

#include <lib_grid/algorithms/extrusion/DiamondsEstablish3D.h>

namespace ug
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
			m_vecElems2BQuenched(VecElems2BQuenched()),
			m_disappearingVols(std::vector<Volume*>()),
			m_disappearingFacs(std::vector<Face*>()),
			m_disappearingEdgs(std::vector<Edge*>())
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

bool DiamondsEstablish3D::findRegions2BShrinked()
{
	// collect the pairs of volumes connected by a face relevant for the shrinking of the volumes

//	for( auto & vmvcd : m_vecVolManifVrtxCombiToShrink4Diams )
//	{
////		Elems2BQuenched e2bq;
//		VrtxPair oldAndShiftV;
//		vmvcd.spuckOldAndShiftVrtx( oldAndShiftV);
//		Vertex * oldCenter = oldAndShiftV.first;
//	}

	VecVolManifVrtxCombi vecVolManifVrtxCopy = m_vecVolManifVrtxCombiToShrink4Diams;

	for( typename VecVolManifVrtxCombi::iterator itVMVOuter = vecVolManifVrtxCopy.begin();
				  itVMVOuter != vecVolManifVrtxCopy.end(); )
	{
		VolManifVrtxCombi outer = *itVMVOuter;

		VrtxPair oldAndShiftVrtxOuter;
		outer.spuckOldAndShiftVrtx( oldAndShiftVrtxOuter );

		Vertex * oldVrtxOuter = oldAndShiftVrtxOuter.first;
		Vertex * shiftVrtxOuter = oldAndShiftVrtxOuter.second;

		Face * faceOuter;
		outer.spuckManif(faceOuter);

//		IndexType twin2OuterFound = 0;
//		IndexType partner2OuterFound = 0;
//		IndexType twin2InnerFound = 0;

		VecVolManifVrtxCombi twinCombi2Outer;
		// should get length 2 in the procedure

		twinCombi2Outer.push_back(outer); // in each case part a twin has to exist, sharing same face

		// should get length 0 or 2 in the procedure: 0 in case of 3 cross, 2 in case of 2 cross
		VecVolManifVrtxCombi partnerCombi2Outer;

//		bool partnerFound = false;

		for( typename VecVolManifVrtxCombi::iterator itVMVInner = itVMVOuter + 1;
					  itVMVInner != vecVolManifVrtxCopy.end(); )
		{
			VolManifVrtxCombi inner = *itVMVInner;

			VrtxPair oldAndShiftVrtxInner;
			inner.spuckOldAndShiftVrtx( oldAndShiftVrtxInner );

			Vertex * oldVrtxInner = oldAndShiftVrtxInner.first;
//			Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;

			if( oldVrtxInner == oldVrtxOuter )
			{
				Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;

				if( shiftVrtxInner != shiftVrtxOuter )
				{
					// if connected by a face, twin to outer found, partner may exist or not
//					 twin2OuterFound++;
					// check if the face is identical, else is from the partner twin

					Face * faceInner;
					inner.spuckManif(faceInner);

					if( faceInner == faceOuter )
					{
						// twin found
						twinCombi2Outer.push_back(inner);

//						twinFound = true;

						itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);

					}
					// else can be any combination of a completely different side of the crossing point
				}
				else if( shiftVrtxInner == shiftVrtxOuter )
				{
					// side partner found, for this side partner, the twin is needed as well!
//					partner2OuterFound++;
					partnerCombi2Outer.push_back(inner);
					// will be able to find only one partner in this loop,
					//its twin has to be found in an additional inner loop!

					itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);

				}

			}
			else
			{
				// else: do nothing, as no relevant relation between the two objects
				++itVMVInner;
			}
		}

		if( partnerCombi2Outer.size() == 1 )
		{
			// additional inner loop to find the twin of the partner

			VolManifVrtxCombi partnerOuter = partnerCombi2Outer[0];

			VrtxPair oldAndShiftVrtxPartnerOuter;
			outer.spuckOldAndShiftVrtx( oldAndShiftVrtxPartnerOuter );

			Vertex * oldVrtxPartnerOuter = oldAndShiftVrtxPartnerOuter.first;
			Vertex * shiftVrtxPartnerOuter = oldAndShiftVrtxPartnerOuter.second;

			Face * facePartnerOuter;
			partnerOuter.spuckManif(facePartnerOuter);

			for( typename VecVolManifVrtxCombi::iterator itVMVInner = itVMVOuter + 1;
						  itVMVInner != vecVolManifVrtxCopy.end(); )
			{
				VolManifVrtxCombi inner = *itVMVInner;

				VrtxPair oldAndShiftVrtxInner;
				inner.spuckOldAndShiftVrtx( oldAndShiftVrtxInner );

				Vertex * oldVrtxInner = oldAndShiftVrtxInner.first;
//				Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;

				if( oldVrtxInner == oldVrtxPartnerOuter )
				{
					Vertex * shiftVrtxInner = oldAndShiftVrtxInner.second;

					if( shiftVrtxInner != shiftVrtxPartnerOuter )
					{
						// if connected by a face, twin to outer found, partner may exist or not
	//					 twin2OuterFound++;
						// check if the face is identical, else is from the partner twin

						Face * faceInner;
						inner.spuckManif(faceInner);

						if( faceInner == facePartnerOuter )
						{
							// twin found
							partnerCombi2Outer.push_back(inner);

	//						twinFound = true;

							itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);

						}
						// else can be any combination of a completely different side of the crossing point
					}
					else if( shiftVrtxInner == shiftVrtxOuter )
					{
						UG_LOG("too much partners or twins " << std::endl);
						// darf nicht mehr passieren
//						UG_THROW("");
						// side partner found, for this side partner, the twin is needed as well!
	//					partner2OuterFound++;
//						partnerCombi2Outer.push_back(inner);
						// will be able to find only one partner in this loop,
						//its twin has to be found in an additional inner loop!

						itVMVInner = vecVolManifVrtxCopy.erase(itVMVInner);

					}

				}
				else
				{
					// else: do nothing, as no relevant relation between the two objects
					++itVMVInner;
				}


			}

			if( partnerCombi2Outer.size() != 2 )
			{
				UG_LOG("no partner found " << std::endl);
				return false;
			}

		}
		else if( partnerCombi2Outer.size() != 0 )
		{
			UG_LOG("strange partner size " << std::endl);
			return false;
		}

		if( twinCombi2Outer.size() != 2 )
		{
			UG_LOG("strange twin size " << std::endl);
			return false;
		}

		itVMVOuter = vecVolManifVrtxCopy.erase(itVMVOuter);

		// TODO FIXME
		// FullLowDimManifQntpl needs to be created from each of the two or four contributions
		// something like transform VecVolManifVrtxCombi to VecVolManifVrtxCombi
		// VecFullLowDimManifQuintuplet needs to be created
		// to check from outside if the central vertex is kept!!!
		// create the new combination and attach it to the vector
		//		Elems2BQuenched
//		m_vecElems2BQuenched.push_back(...);
	}




	return true;
}



} /* namespace diamonds */


} /* namespace ug */
