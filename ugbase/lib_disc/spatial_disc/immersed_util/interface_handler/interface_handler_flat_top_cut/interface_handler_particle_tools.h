/*
 * interface_handler_local_tools.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_FLAT_TOP_TOOLS_H_
#define INTERFACE_HANDLER_FLAT_TOP_TOOLS_H_

namespace ug{
 
////////////////////////////////////////////////////////////////////////////////
//	FlatTopBase methods - ROID/CollectCorners...
////////////////////////////////////////////////////////////////////////////////


template <int TWorldDim>
bool InterfaceHandlerLocalParticle<TWorldDim>::
is_boundary_face_for2(const size_t sideID)
{
	if ( dim == 3 )
		UG_THROW("in 'is_boundary_face_for2()': only implemented for 2d!!\n");

// get data
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(this->m_roid);

	//	number of corners of side (special case bottom side pyramid)
		const int coOfSide = (this->m_roid != ROID_PYRAMID || sideID != 0)
							? rRefElem.num(dim-1, sideID, 0) : rRefElem.num(dim-1, sideID, 0) + 2;

	if ( coOfSide > 2 )
		UG_THROW("in 'is_boundary_face_for2()': hmm...coOfSide > 2 ??? \n coOfSide = " << coOfSide << "\n");

	std::vector<int> prtIndex(2);

	for(int co = 0; co < coOfSide; ++co)
	{
		size_t cornerID;
		if (this->m_roid != ROID_PYRAMID || sideID != 0)
			cornerID = rRefElem.id(dim-1, sideID, 0, co);
		else
			cornerID = rRefElem.id(dim-1, sideID, 0, (co % 3) + (co>3 ? 1 : 0));

		if ( !this->lies_onInterface(cornerID) )
			return false;
		else
			prtIndex[co] = getPrtIndex(co);

	}

	UG_LOG("prtIndex[0]: " << prtIndex[0] << "\n");
	UG_LOG("prtIndex[1]: " << prtIndex[1] << "\n");

	if (  prtIndex[0] ==  prtIndex[1] )
	{	UG_THROW("========= is_boundary_face = TRUE\n");

		return true;
	}
	UG_LOG("========= is_boundary_face = FALSE\n");

	return false;
}



/*
// compare implementation of 'DimFV1Geometry::update_boundary_faces()'
template <int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_faces()
{

	/////////////////////////////////////////////////////////////////////////////
	// get general data
	/////////////////////////////////////////////////////////////////////////////
 		const DimReferenceElement<dim>& rRefElem
			= ReferenceElementProvider::get<dim>(this->m_roid);

		DimReferenceMapping<dim, dim>& rMapping = ReferenceMappingProvider::get<dim, dim>(this->m_roid);
		rMapping.update(this->m_vCornerCoords);

		const LocalShapeFunctionSet<dim>& TrialSpace =
			LocalFiniteElementProvider::get<dim>(this->m_roid, LFEID(LFEID::LAGRANGE, dim, 1));

	/////////////////////////////////////////////////////////////////////////////
	//	compute local and global geom object midpoints for each dimension
	/////////////////////////////////////////////////////////////////////////////

		MathVector<dim> vvLocMid[dim+1][maxMid];
		MathVector<dim> vvGloMid[dim+1][maxMid];

 	// 	set corners of element as local centers of nodes
		for(size_t i = 0; i < rRefElem.num(0); ++i)
			vvLocMid[0][i] = rRefElem.corner(i);

	//	compute local midpoints
		interfaceComputeMidPoints<dim, DimReferenceElement<dim>, maxMid>(rRefElem, vvLocMid[0], vvLocMid);

	// 	remember global position of nodes
		for(size_t i = 0; i < rRefElem.num(0); ++i)
			vvGloMid[0][i] = this->m_vCornerCoords[i];
	 //	compute local midpoints
		interfaceComputeMidPoints<dim, DimReferenceElement<dim>, maxMid>(rRefElem, vvGloMid[0], vvGloMid);

	/////////////////////////////////////////////////////////////////////////////
	// collect boudary faces
	/////////////////////////////////////////////////////////////////////////////

	// get number of sides of element
 		size_t numSides = 0;
		numSides = rRefElem.num(dim-1);

	//	current number of bf
		size_t curr_bf = 0;

		m_vBF.clear();
	//	loop sides of element
		for(size_t side = 0; side < numSides; ++side)
		{
		// side is no boundary face => continue
			if ( !is_boundary_face(side) )
				 continue;

		//	number of corners of side (special case bottom side pyramid)
			const int coOfSide = (this->m_roid != ROID_PYRAMID || side != 0)
								? rRefElem.num(dim-1, side, 0) : rRefElem.num(dim-1, side, 0) + 2;

		//	resize vector
			m_vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				interfaceBF& bf = m_vBF[curr_bf];

			//	set node id == scv this bf belongs to
				if (this->m_roid != ROID_PYRAMID || side != 0)
					bf.nodeId = rRefElem.id(dim-1, side, 0, co);
				else
				{
					// map according to order defined in ComputeBFMidID
					bf.nodeId = rRefElem.id(dim-1, side, 0, (co % 3) + (co>3 ? 1 : 0));
				}

			//	Compute MidID for BF
				interfaceComputeBFMidID(rRefElem, side, bf.vMidID, co);

			// 	copy corners of bf
				interfaceCopyCornerByMidID<dim, maxMid>(bf.vLocPos, bf.vMidID, vvLocMid, interfaceBF::numCo);
				interfaceCopyCornerByMidID<dim, maxMid>(bf.vGloPos, bf.vMidID, vvGloMid, interfaceBF::numCo);

			// 	integration point
				AveragePositions(bf.localIP, bf.vLocPos, interfaceBF::numCo);
				AveragePositions(bf.globalIP, bf.vGloPos, interfaceBF::numCo);

			// 	normal on scvf
				traits::NormalOnSCVF(bf.Normal, bf.vGloPos, vvGloMid[0]);

 			//	compute volume
				bf.Vol = VecTwoNorm(bf.Normal);

			//	compute shapes and grads
				bf.numSH = TrialSpace.num_sh();
				TrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
				TrialSpace.grads(&(bf.vLocalGrad[0]), bf.localIP);

			//	get reference mapping
				rMapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
				bf.detj = rMapping.sqrt_gram_det(bf.localIP);

			//	compute global gradients
				for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
					MatVecMult(bf.vGlobalGrad[sh], bf.JtInv, bf.vLocalGrad[sh]);

			//	increase curr_bf
				++curr_bf;

			} // end loop of corners of side

		} // end loop sides of element


}
*/





/*
template <int TWorldDim>
int InterfaceHandlerLocalParticle<TWorldDim>::
CollectCorners_StdFV_for2(GridObject* elem)
{

	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	this->m_vCornerCoords.clear();
	this->m_vInterfaceID.clear();
	this->m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

	// get prtIndex:
	bool isOutsideVertex = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
 		if ( is_OutsideVertex(vVertex[i]) || is_FTVertex(vVertex[i]) )
 		{
 			isOutsideVertex = true;
 			break;
 		}
	}

	if ( !isOutsideVertex )
		UG_THROW("Error in 'CollectCorners_FlatTop_2d': no vertex is FTVertex: should be true for at least 1 vertex!\n");

// loop vertices
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// get element
		Vertex* vrtRoot = vVertex[i];
		int prtIndex;

		//////////////////////////////////////////////
		// case 1:
		// vertex insideFluid => Position und index puschen
		if ( !is_FTVertex(vrtRoot, i) )

		//if ( !is_OutsideVertex(vrtRoot) && !is_FTVertex(vVertex[i]))
		{
			if ( is_nearInterfaceVertex(vrtRoot) )
			{
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");
			}
			else
			{
				this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
				this->m_vOriginalCornerID.push_back(i);

				vInsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
			}

		}
		else
		{
 			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
			this->m_vOriginalCornerID.push_back(i);
			this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'this->m_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));

		}

	}

	// skip whole element, since only FT points are included
 	if ( dim == 2 && this->m_vCornerCoords.size() != 3 )
 		UG_THROW("in 'CollectCorners_StdFV()': this->m_vCornerCoords.size() "
					"= " << this->m_vCornerCoords.size() << "not possible!\n");
 	if ( dim == 3 && this->m_vCornerCoords.size() != 4 )
 		UG_THROW("in 'CollectCorners_StdFV()': this->m_vCornerCoords.size() "
					"= " << this->m_vCornerCoords.size() << "not possible!\n");

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	const char* filename;
 	filename = "output_StdFV.";

	std::string name(filename);
	char ext[50]; sprintf(ext, "txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");

	fprintf(outputFile,"\n_________________________ begin 'CollectIBCorners_StdFV' ________________________\n");

	for ( size_t i = 0; i < vVertex.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", this->m_aaPos[vVertex[i]][0] , this->m_aaPos[vVertex[i]][1], this->m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: this->m_vCornerCoords.size() %lu \n", this->m_vCornerCoords.size());

	for ( size_t i = 0; i < this->m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", this->m_vCornerCoords[i][0], this->m_vCornerCoords[i][1], this->m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < this->m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: this->m_vInterfaceID = %lu\n", this->m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < this->m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: this->m_vOriginalCornerID = %lu\n", this->m_vOriginalCornerID[i]);


	fprintf(outputFile,"\n_________________________ end 'CollectIBCorners_StdFV()' ______________________________\n\n");

	fclose(outputFile);



	return this->m_vCornerCoords.size();

}


template <int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
ResortQuadrilateral_for2(std::vector<MathVector<dim> > vQuadriCorners)
{

	UG_LOG("start ResortQuadrilateral_for2()...\n");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check for CCW
	// VecDot(vNormOut1,zDir) < 0 => corners of quadrilateral are not CCW:
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	MathVector<3> vNormOut1, vNormOut2;
	std::vector<MathVector<3> > vCornerCoordsBlow; vCornerCoordsBlow.resize(4);
	for ( size_t i = 0; i < vQuadriCorners.size(); ++i )
	{
		for ( size_t d = 0; d < dim; ++d )
			vCornerCoordsBlow[i][d] = vQuadriCorners[i][d];
		if ( dim == 2 ) vCornerCoordsBlow[i][2] = 0.0;
	}
	CalculateTriangleNormalNoNormalize(vNormOut1, vCornerCoordsBlow[0], vCornerCoordsBlow[1], vCornerCoordsBlow[2]);
	CalculateTriangleNormalNoNormalize(vNormOut2, vCornerCoordsBlow[2], vCornerCoordsBlow[3], vCornerCoordsBlow[0]);

	UG_LOG("vNormOut1: " << vNormOut1 << "\n");
	UG_LOG("vNormOut2: " << vNormOut2 << "\n");

	if ( VecDot(vNormOut1, vNormOut2) < 0 ) UG_THROW("----> check 1 failed after second check!\n");


	MathVector<3> normalDir;
	for ( size_t i = 0; i < 2; ++i )
		normalDir[i] = 0.0;
	normalDir[2] = 1.0;

	if ( VecDot(vNormOut1, normalDir) < 0 && VecDot(vNormOut2, normalDir) > 0)
	{
		UG_THROW("case1: inconsistent orientation: VecDot(vNormOut1, normalDir): " << VecDot(vNormOut1, normalDir) <<
				", VecDot(vNormOut2, normalDir): " << VecDot(vNormOut2, normalDir) <<"\n");
	}
	else if ( VecDot(vNormOut1, normalDir) > 0 && VecDot(vNormOut2, normalDir) < 0)
	{
		UG_THROW("case1: inconsistent orientation: VecDot(vNormOut1, normalDir): " << VecDot(vNormOut1, normalDir) <<
				", VecDot(vNormOut2, normalDir): " << VecDot(vNormOut2, normalDir) <<"\n");
	}
	else
	{
		std::vector<MathVector<dim> > vCornerCoordsBuffer;
		vCornerCoordsBuffer.clear();
		m_vQuadriCorners_for2.clear();
		std::vector<size_t> vQuadriIDBuffer;
		vQuadriIDBuffer.clear();
		if ( vQuadriCorners.size() != 4 )
			UG_THROW("wrong size of vQuadriCorners: " << vQuadriCorners.size() << "\n");

	// A: change orientation OR simply copy:
		if ( VecDot(vNormOut1, normalDir) < 0 && VecDot(vNormOut2, normalDir) < 0)
		{
			for ( size_t i = 0; i < 4; ++i )
			{
				vCornerCoordsBuffer.push_back(vQuadriCorners[3-i]);
				vQuadriIDBuffer.push_back(this->m_vQuadriOrigID[3-i]);
			}
		}
		else
		{
			for ( size_t i = 0; i < 4; ++i )
			{
				vCornerCoordsBuffer.push_back(vQuadriCorners[i]);
				vQuadriIDBuffer.push_back(this->m_vQuadriOrigID[i]);
			}
		}
	// B: re-order, s.t. 1.,2. node with prtIndex = 0; 3.,4. node with prtIndex = 1:
		this->m_vQuadriOrigID.clear();
		for ( size_t i = 0; i < 3; ++i )
		{
			m_vQuadriCorners_for2.push_back(vCornerCoordsBuffer[i+1]);
			this->m_vQuadriOrigID.push_back(vQuadriIDBuffer[i+1]);
		}
		m_vQuadriCorners_for2.push_back(vCornerCoordsBuffer[0]);
		this->m_vQuadriOrigID.push_back(vQuadriIDBuffer[0]);

	}


}

template <int TWorldDim>
int InterfaceHandlerLocalParticle<TWorldDim>::
CollectTriangle_and_Quadri_for2(GridObject* elem, Vertex* vrtInside)
{
	UG_LOG("start CollectTriangle_and_Quadri_for2()...\n");

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

	for ( size_t i = 0; i < vVertex.size(); ++i )
		UG_LOG("pos[" << i << "]: " << this->m_aaPos[vVertex[i]] << "\n");

//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrtRoot = vVertex[i];

		if ( !is_FTVertex(vrtRoot) )
		{
			if ( vrtRoot != vrtInside )
				UG_THROW("in 'CollectTriangle_and_Quadri_for2()': error in collecting inside vertex\n");
			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
			this->m_vOriginalCornerID.push_back(i);
			this->m_vNOInterfaceID.push_back(i);
		}
		else
		{

		// compute intersection point of edge, build by vrtRoot and vrtInside:
			MathVector<dim> intersectionPnt;
 			get_intersection_point(intersectionPnt, vrtInside, vrtRoot);

			UG_LOG("posRoot: " << this->m_aaPos[vrtRoot] << "posInside: " << this->m_aaPos[vrtInside] << "\n");
			UG_LOG("intersectionPnt: " << intersectionPnt << "\n");

			// check for correct inersectionPnt
			if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
				UG_THROW("in 'CollectTriangle_and_Quadri_for2()': Error in computation of 'intersectionPnt':\n "
						" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");

			this->m_vCornerCoords.push_back(intersectionPnt);
			this->m_vOriginalCornerID.push_back(i);
			this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);
			m_vQuadriCorners_for2.push_back(intersectionPnt);
			this->m_vQuadriOrigID.push_back(i);

 		}

	} // end vrt-loop

	CollectQuadrilateral_besideTri_for2(elem);

}



template <int TWorldDim>
int InterfaceHandlerLocalParticle<TWorldDim>::
CollectQuadrilateral_besideTri_for2(GridObject* elem)
{
	UG_LOG("start CollectQuadrilateral_besideTri_for2()...\n");

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	//////////////////////////////////////////////////////////////////////////////////////////
	//	A: start with ordering 1 0 for the first corners in 'm_vQuardi_for2'
	//		=> if ordering is 0 1 ---> switch
	//	B: add the two further corners in the order 0 1
	//		=> we end up with the setting 1 0 0 1
	//  C: reorder as done in 'ResortQuadrilateral_for2():
	// 		=> final order as it ought to be: 0 0 1 1 :-)
	//////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////
	//	A: check for ordering, weather it is 0 1:
	if ( m_vQuadriCorners_for2.size() != 2 )
		UG_THROW("worng start in 'CollectQuadrilateral_besideTri_for2'!\n");

	size_t index_of_first_quadri_corner = this->m_vInterfaceID[0];
	Vertex* vrt = vVertex[index_of_first_quadri_corner];
	int prtIndex_of_first_quadri_corner;
// get prtIndex of first quadri corner:
	if ( !is_FTVertex(prtIndex_of_first_quadri_corner, vrt) )
		UG_THROW("error in 'CollectQuadrilateral_besideTri_for2()'...\n");

	UG_LOG("prtIndex_of_first_quadri_corner: " << prtIndex_of_first_quadri_corner << "\n");

	if ( prtIndex_of_first_quadri_corner == 0 )
		this->switch_order();

	UG_LOG("============= after A: CHECK order ===============\n");
	for ( size_t i = 0; i < m_vQuadriCorners_for2.size(); ++i  )
		UG_LOG("m_vQuadriCorners_for2[" << i << "]: " << m_vQuadriCorners_for2[i] << "\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i  )
		UG_LOG("this->m_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

	//////////////////////////////////////////////////////////////////////////////////////////
	//	B: add the two further corners in the order 0 1
	//		=> we end up with the setting 1 0 0 1

	// loop edges and get edge, cut by both particles:
	bool pushOrderSwitch = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		size_t j = i+1;
		if ( j == vVertex.size() ) j = 0;

		Vertex* vrt1 = vVertex[i];
		Vertex* vrt2 = vVertex[j];
		int prtIndex1, prtIndex2;

		MathVector<dim> intersectionPnt1, intersectionPnt2;

		if ( !is_FTVertex(prtIndex1, vrt1) || !is_FTVertex(prtIndex2, vrt2) )
			continue;

		UG_LOG("this->m_aaPos[vrt1] = " << this->m_aaPos[vrt1] << "\n");
		UG_LOG("this->m_aaPos[vrt2] = " << this->m_aaPos[vrt2] << "\n");
		UG_LOG("prtIndex1 = " << prtIndex1 << ", prtIndex2 = " << prtIndex2 << "\n");

		if ( prtIndex1 == prtIndex2 )
			UG_THROW("prtIndex1 == prtIndex2 NOT POSSIBLE on this element!\n");

	/// intersectionPnt1 := intersection with prtIndex1 and vrt1 = inside circle of particle with prtIndex1
		get_intersection_point_for2(intersectionPnt1, vrt2, vrt1, prtIndex1);
		// check for correct inersectionPnt
		if ( fabs(get_LSvalue_byPosition_for2(intersectionPnt1, prtIndex1)) > 1e-6  )
			UG_THROW("in 'CollectQuadrilateral_besideTri_for2()': Error in computation of 'intersectionPnt1':\n "
					" intersectionPnt1 = " << intersectionPnt1 << "\n distance from interace = " << fabs(get_LSvalue_byPosition_for2(intersectionPnt1, prtIndex1)) << "\n");

	/// intersectionPnt2 := intersection with prtIndex2 and vrt2 = inside circle of particle with prtIndex2
		get_intersection_point_for2(intersectionPnt2, vrt1, vrt2, prtIndex2);
		// check for correct inersectionPnt
		if ( fabs(get_LSvalue_byPosition_for2(intersectionPnt2, prtIndex2)) > 1e-6  )
			UG_THROW("in 'CollectQuadrilateral_besideTri_for2()': Error in computation of 'intersectionPnt2':\n "
					" intersectionPnt2 = " << intersectionPnt2 << "\n distance from interace = " << fabs(get_LSvalue_byPosition_for2(intersectionPnt2, prtIndex2)) << "\n");

		UG_LOG("intersectionPnt1 = " << intersectionPnt1 << "\n");
		UG_LOG("intersectionPnt2 = " << intersectionPnt2 << "\n");

	// FIRST: push intersectionpoint of particle 0, SECOND: push intersectionpoint of particle 1
		if ( prtIndex1 == 0 )
		{
			m_vQuadriCorners_for2.push_back(intersectionPnt1);
			m_vQuadriCorners_for2.push_back(intersectionPnt2);

			this->m_vQuadriOrigID.push_back(i);
			this->m_vQuadriOrigID.push_back(j);

		}
		else
		{
			m_vQuadriCorners_for2.push_back(intersectionPnt2);
			m_vQuadriCorners_for2.push_back(intersectionPnt1);

			this->m_vQuadriOrigID.push_back(j);
			this->m_vQuadriOrigID.push_back(i);
		}
		// => resulting order according to prtIndex: 1 0 0 1 :-)

	} // end edge-loop

	UG_LOG("============= after B: CHECK order ===============\n");
	for ( size_t i = 0; i < m_vQuadriCorners_for2.size(); ++i  )
		UG_LOG("m_vQuadriCorners_for2[" << i << "]: " << m_vQuadriCorners_for2[i] << "\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i  )
		UG_LOG("this->m_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

	//////////////////////////////////////////////////////////////////////////////////////////
	//  C: reorder as done in 'ResortQuadrilateral_for2():
	// 		=> final order as it ought to be: 0 0 1 1 :-)

	ResortQuadrilateral_for2(m_vQuadriCorners_for2);
	//  => nach reordering: 0 0 1 1 ;-)


	UG_LOG("============= after C: CHECK order ===============\n");
	for ( size_t i = 0; i < m_vQuadriCorners_for2.size(); ++i )
		UG_LOG("m_vQuadriCorners_for2[" << i << "]: " << m_vQuadriCorners_for2[i] << "\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i  )
		UG_LOG("this->m_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

	UG_LOG("end CollectQuadrilateral_besideTri_for2()...\n");

}


template <int TWorldDim>
int InterfaceHandlerLocalParticle<TWorldDim>::
CollectQuadrilateral_for2(GridObject* elem)
{
	UG_LOG("start CollectQuadrilateral_for2()...\n");

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vQuadriCorners;
	std::vector<MathVector<dim> > buffer_vQuadriCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

	for ( size_t i = 0; i < vVertex.size(); ++i )
		UG_LOG("pos[" << i << "]: " << this->m_aaPos[vVertex[i]] << "\n");

//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	bool pushOrderSwitch = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		size_t j = i+1;
		if ( j == vVertex.size() ) j = 0;


		Vertex* vrt1 = vVertex[i];
		Vertex* vrt2 = vVertex[j];
		int prtIndex1, prtIndex2;

		MathVector<dim> intersectionPnt1, intersectionPnt2;

		if ( !is_FTVertex(prtIndex1, vrt1) )
			UG_THROW("wrong method! vrt1 = " << this->m_aaPos[vrt1] << " lies NOT outside!\n");
		if ( !is_FTVertex(prtIndex2, vrt2) )
			UG_THROW("wrong method! vrt2 = " << this->m_aaPos[vrt2] << " lies NOT outside!\n");

		UG_LOG("this->m_aaPos[vrt1] = " << this->m_aaPos[vrt1] << "\n");
		UG_LOG("this->m_aaPos[vrt2] = " << this->m_aaPos[vrt2] << "\n");

		UG_LOG("prtIndex1 = " << prtIndex1 << ", prtIndex2 = " << prtIndex2 << "\n");

	// no intersection by interface for prtIndex1 == prtIndex2
		if ( prtIndex1 == prtIndex2 )
		{
			UG_LOG("continue...\n");

			continue;
		}
	/// intersectionPnt1 := intersection with prtIndex1 and vrt1 = inside circle of particle with prtIndex1
		get_intersection_point_for2(intersectionPnt1, vrt2, vrt1, prtIndex1);
		// check for correct inersectionPnt
		if ( fabs(get_LSvalue_byPosition_for2(intersectionPnt1, prtIndex1)) > 1e-6  )
			UG_THROW("in 'CollectQuadrilateral_for2()': Error in computation of 'intersectionPnt1':\n "
					" intersectionPnt1 = " << intersectionPnt1 << "\n distance from interace = " << fabs(get_LSvalue_byPosition_for2(intersectionPnt1, prtIndex1)) << "\n");

	/// intersectionPnt2 := intersection with prtIndex2 and vrt2 = inside circle of particle with prtIndex2
		get_intersection_point_for2(intersectionPnt2, vrt1, vrt2, prtIndex2);
		// check for correct inersectionPnt
		if ( fabs(get_LSvalue_byPosition_for2(intersectionPnt2, prtIndex2)) > 1e-6  )
			UG_THROW("in 'CollectQuadrilateral_for2()': Error in computation of 'intersectionPnt2':\n "
					" intersectionPnt2 = " << intersectionPnt2 << "\n distance from interace = " << fabs(get_LSvalue_byPosition_for2(intersectionPnt2, prtIndex2)) << "\n");

		UG_LOG("intersectionPnt1 = " << intersectionPnt1 << "\n");
		UG_LOG("intersectionPnt2 = " << intersectionPnt2 << "\n");

		if ( !pushOrderSwitch )
		{
			pushOrderSwitch = true;
		// FIRST: push intersectionpoint of particle 1, SECOND: push intersectionpoint of particle 0
			if ( prtIndex1 == 1 )
			{
				vQuadriCorners.push_back(std::make_pair(intersectionPnt1, prtIndex1));
				vQuadriCorners.push_back(std::make_pair(intersectionPnt2, prtIndex2));

				this->m_vQuadriOrigID.push_back(i);
				this->m_vQuadriOrigID.push_back(j);
			}
			else
			{
				vQuadriCorners.push_back(std::make_pair(intersectionPnt2, prtIndex2));
				vQuadriCorners.push_back(std::make_pair(intersectionPnt1, prtIndex1));

				this->m_vQuadriOrigID.push_back(j);
				this->m_vQuadriOrigID.push_back(i);
			}
		}
		else
		{
		// FIRST: push intersectionpoint of particle 0, SECOND: push intersectionpoint of particle 1
			if ( prtIndex1 == 1 )
			{
				vQuadriCorners.push_back(std::make_pair(intersectionPnt2, prtIndex2));
				vQuadriCorners.push_back(std::make_pair(intersectionPnt1, prtIndex1));

				this->m_vQuadriOrigID.push_back(j);
				this->m_vQuadriOrigID.push_back(i);
			}
			else
			{
				vQuadriCorners.push_back(std::make_pair(intersectionPnt1, prtIndex1));
				vQuadriCorners.push_back(std::make_pair(intersectionPnt2, prtIndex2));

				this->m_vQuadriOrigID.push_back(i);
				this->m_vQuadriOrigID.push_back(j);
			}
		}

		// => resulting order according to prtIndex: 0 1 1 0 :-)

	} // end edge-loop

	// copy data for call of 'ResortQuadrilateral_for2()':
	buffer_vQuadriCorners.clear();
	for ( size_t i = 0; i < vQuadriCorners.size(); ++i )
	{
		buffer_vQuadriCorners.push_back(vQuadriCorners[i].first);
		UG_LOG("vQuadriCorners[" << i << "]: " << vQuadriCorners[i].first << ", index: " << vQuadriCorners[i].second << "\n");
	}
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i )
		UG_LOG("this->m_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

	ResortQuadrilateral_for2(buffer_vQuadriCorners);


	UG_LOG("\n\n------------> NACH ordering:\n");
	for ( size_t i = 0; i < vQuadriCorners.size(); ++i )
		UG_LOG("m_vQuadriCorners_for2[" << i << "]: " << m_vQuadriCorners_for2[i] << "\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i )
		UG_LOG("this->m_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

	UG_LOG("end ResortQuadri...\n");

}


template <int TWorldDim>
int InterfaceHandlerLocalParticle<TWorldDim>::
CollectCorners_FlatTop_2d_for2(GridObject* elem)
{
	UG_LOG("start CollectCorners_FlatTop_2d_for2()...\n");


	m_vQuadriCorners_for2.clear();
	this->m_vQuadriOrigID.clear();
	this->m_vCornerCoords.clear();
	this->m_vInterfaceID.clear();
	this->m_vNOInterfaceID.clear();
	this->m_vOriginalCornerID.clear();


//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

	bool inside = false;
	Vertex* vrtInside;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		if ( ! is_FTVertex(vVertex[i]) )
		{
 			inside = true;
 			vrtInside = vVertex[i];
 			break;
		}
	}

	if ( inside )
	{
		UG_LOG("vrtInside: " << this->m_aaPos[vrtInside] << "\n");
		this->m_roid = ROID_TRIANGLE;
		CollectTriangle_and_Quadri_for2(elem, vrtInside);
	}
	else
	{
		this->m_roid = ROID_QUADRILATERAL;
		CollectQuadrilateral_for2(elem);
	/////////////////////////////////////////////////////////////////////////////
	// write data seperately, since 'CollectQuadrilateral_for2()'
	//	is also called during 'CollectTriangle_and_Quadri_for2()':
	/////////////////////////////////////////////////////////////////////////////
	// A: copy corners into this->m_vCornerCoords:
		for ( size_t i = 0; i < m_vQuadriCorners_for2.size(); ++i )
			this->m_vCornerCoords.push_back(m_vQuadriCorners_for2[i]);

	// B: set all corners to InterfaceIDs:
		for ( size_t i = 0; i < 4; ++i )
		{
			this->m_vInterfaceID.push_back(i);
			this->m_vOriginalCornerID.push_back(this->m_vQuadriOrigID[i]);
		}
	}

	if ( inside )
	{
		for ( size_t i = 0; i < this->m_vNOInterfaceID.size(); ++i )
		{
			UG_LOG("this->m_vNOInterfaceID[" << i << "]: " << this->m_vNOInterfaceID[i] << "\n");
		}
	UG_LOG("--------------------> m_vQuadriCorners_for2.size(): " << m_vQuadriCorners_for2.size() << "\n");
	for ( size_t i = 0; i < m_vQuadriCorners_for2.size(); ++i )
	{
		UG_LOG("m_vQuadriCorners_for2[" << i << "]: " << m_vQuadriCorners_for2[i] << "\n");
	}
	for ( size_t i = 0; i < this->m_vCornerCoords.size(); ++i )
	{
		UG_LOG("this->m_vCornerCoords[" << i << "]: " << this->m_vCornerCoords[i] << "\n");
	}
	if ( this->m_vCornerCoords.size() == 3 )
			UG_LOG("++++++++++++++++ end CollectCorners_FlatTop_2d_for2()...\n");
	}


}

*/

template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_old()
{
// VORSICHT: update_inner_boundary_radial() wird in fv1FT_geom_impl.h aufgerufen, BEVOR remap_BF() stattfindet!
//	=> bf.node_id() immer mit interfaceID, nicht OriginalCornerID!

	if (m_prtIndex == -1)
		UG_THROW(
				"InterfaceHandlerLocalParticle::update_inner_boundary_radial_old(): m_prtIndex not set!\n");
	const MathVector<dim> center = m_spCutElementHandler->get_center(
			m_prtIndex);

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

//  Loop the boundary faces to assemble impulse equations
	for (size_t i = 0; i < this->m_vBF.size(); ++i) {
		interfaceBF bf = this->m_vBF[i];
		VecSubtract(m_vRadialAtCo[bf.node_id()], this->corner(bf.node_id()),
				center);

		// original version for 'm_vRadialAtIP'-computation;
		//
		//			VecSubtract(m_vRadialAtIP[bf.nodeId], bf.globalIP, center);
		//
		// BUT: after multiplication with normal of particle boundary (:= n = bf.normal) the
		// integral '\int (w x r)*n dS' must be 0.0!!
		//  	=> r := m_vRadialAtIP != bf.normal
		// 		=> (wxr)*n = 0.0, since r = n
		//
		// Remark: m_massDefect = 0.0 INDEPENDENT of choice of r!

		if (1) //m_bRadial_forMassEq_equals_Normal ) // = default case!
		{
			VecScale(m_vRadialAtIP[bf.node_id()], bf.normal(), -1.0);
			//	UG_LOG("m_vRadialAtIP[" << bf.node_id() << "] = " << m_vRadialAtIP[bf.node_id()] << "\n");
			//	UG_LOG("bf.normal() = " << bf.normal() << "\n");
		} else if (1) {
			VecSubtract(m_vRadialAtIP[bf.node_id()], bf.global_ip(), center);
			//	VecNormalize(m_vRadialAtIP[bf.node_id()], m_vRadialAtIP[bf.node_id()]);
			//	VecScale(m_vRadialAtIP[bf.node_id()], m_vRadialAtIP[bf.node_id()], VecLength(bf.normal()));
			//		UG_LOG("bf.global_ip() " << bf.global_ip() << "\n");
			//		UG_LOG("* m_vRadialAtIP[" << bf.node_id() << "] = " << m_vRadialAtIP[bf.node_id()] << "\n");
		} else {
            typedef DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> > TFVGeom;
			TFVGeom& geo = GeomProvider<TFVGeom>::get(
					LFEID(LFEID::LAGRANGE, dim, 1), 1);

			int ip = (bf.node_id() - 1);
			if (ip < 0)
				ip = ip + 3;

			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
			const MathVector<dim>& SCVFip_Pos = scvf.global_ip();
			VecSubtract(m_vRadialAtIP[ip], SCVFip_Pos, center);
			VecNormalize(m_vRadialAtIP[ip], m_vRadialAtIP[ip]);
			VecScale(m_vRadialAtIP[ip], m_vRadialAtIP[ip],
					VecLength(bf.normal()));

		}

	}

	//	UG_LOG("_____________\n");
	//////////////////////////////////////////////////////////////////////////////////////////
	// REMARK:
	// 	bf.node_id() is still the node id of the flat top element, i.e. the value
	//		can potentially be 2 and 3!
	//			=> rewriting m_vRadialAtIP/AtCo to the original node id necessary
	// 		=> for Quadrilaterals: 2 m_vRadialAtIP/AtCo-values, BUT only 1 original corner
	//					---> averaging over both values and write this average to single corner
	//////////////////////////////////////////////////////////////////////////////////////////

	// 0. initialize buffer array
	std::vector<MathVector<dim> > vRadialAtCo;
	std::vector<MathVector<dim> > vRadialAtIP;
	vRadialAtCo.clear();
	vRadialAtCo.resize(this->numCo(), 0.0);
	vRadialAtIP.clear();
	vRadialAtIP.resize(this->numCo(), 0.0);
	std::vector<int> vCounter;
	vCounter.clear();
	vCounter.resize(this->numCo(), 0);

	// ToDo for 3d?

// 1. add and copy
	for (size_t i = 0; i < this->interface_id_all().size(); ++i) {
		size_t interfaceID = this->interface_id(i);
		size_t originalID = this->corner_orig(interfaceID);

		vRadialAtCo[originalID] += m_vRadialAtCo[interfaceID];
		vRadialAtIP[originalID] += m_vRadialAtIP[interfaceID];

		vCounter[originalID]++;
	}
// 2. rescale and copy
	/*	for(size_t i = 0; i < this->interface_id_all().size(); ++i)
	 {
	 size_t interfaceID = this->interface_id(i);
	 size_t originalID = this->corner_orig(interfaceID);
	 if ( vCounter[originalID] > 1 )
	 {
	 VecScale(m_vRadialAtCo[originalID], vRadialAtCo[originalID], 1.0/vCounter[originalID]);
	 VecScale(m_vRadialAtIP[originalID], vRadialAtIP[originalID], 1.0); // /vCounter[originalID]);
	 //UG_LOG("m_vRadialAtIP[" << originalID << "] = " << m_vRadialAtIP[originalID] << "\n");
	 }
	 }*/
	// 2. rescale and copy
	for (size_t i = 0; i < vRadialAtCo.size(); ++i)
		if (vCounter[i] > 0) {
			VecScale(m_vRadialAtCo[i], vRadialAtCo[i], 1.0 / vCounter[i]);
			VecScale(m_vRadialAtIP[i], vRadialAtIP[i], 1.0 / vCounter[i]);

		}

}

// VORSICHT: update_inner_boundary_radial() wird in fv1FT_geom_impl.h aufgerufen, BEVOR remap_BF() stattfindet!
//	=> bf.node_id() immer mit interfaceID, nicht OriginalCornerID!
template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial(const MathVector<TWorldDim>* vCornerCoords)
{
	if (dim == 2) {
		if (this->m_vBF.size() != 2 && this->m_vBF.size() != 0)
			UG_THROW(
					"oha, this->m_vBF.size() != 2:" << this->m_vBF.size() << "\n");
	}

	if (0)
		UG_LOG("start update_inner_boundary_radial()...\n");
	const int modusIP = 0;
	const bool useResized = true;

	if (m_prtIndex == -1)
		UG_THROW(
				"InterfaceHandlerLocalParticle::update_inner_boundary_radial(): m_prtIndex not set!\n");

	const MathVector<dim> center = m_spCutElementHandler->get_center(m_prtIndex);

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// REMARK: for elements with only ONE outside corner ON the interface, bf.size() = 0!
// 			=> for assembling connections u_fluid -> u_prt, although 'm_vRadialAtCo' is needed!
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// A. set 'm_vRadialAtCo' --> independent of boundary faces!
	for (size_t i = 0; i < this->numCo(); ++i) {
		if (useResized) {
			VecSubtract(m_vRadialAtCo[i], this->corner(i), center);
        } else { UG_THROW("in 'update_inner_boundary_radial()': case not implemented anymore!\n");}
/*
        {
			size_t nodeID = this->corner_orig(i);
			VecSubtract(m_vRadialAtCo[nodeID], vCornerCoords[nodeID], center);
			// project in radial direction onto interface
			VecNormalize(m_vRadialAtCo[nodeID], m_vRadialAtCo[nodeID]);
			VecScale(m_vRadialAtCo[nodeID], m_vRadialAtCo[nodeID], radius);
		}
*/
		if (0) {
			UG_LOG("i = " << i << "\n");
			UG_LOG("this->corner(" << i << ") = " << this->corner(i) << "\n");
			UG_LOG(
					"---> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
		}
	}

	// B. set 'm_vRadialAtIP' --> for boundary faces!
	for (size_t i = 0; i < this->m_vBF.size(); ++i) {
		interfaceBF bf = this->m_vBF[i];

		size_t nodeID = bf.node_id();
		if (!useResized)
			nodeID = this->corner_orig(nodeID);

        typedef DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> > TFVGeom;
		TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

		int ip = (bf.node_id() - 1);
		if (ip < 0)
			ip = ip + 3;

		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		const MathVector<dim>& SCVFip_Pos = scvf.global_ip();

		switch (modusIP) {
		case 0: 	// AtIP == AtCo
			m_vRadialAtIP[nodeID] = m_vRadialAtCo[nodeID];
			break;
		case 1: 	// AtIP == bf.global_ip()
			VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center);
			//		VecNormalize(m_vRadialAtIP[nodeID], m_vRadialAtIP[nodeID]);
			//		VecScale(m_vRadialAtIP[nodeID], m_vRadialAtIP[nodeID], radius);
			break;
		case 2: 	// AtIP == bf.normal()
            UG_THROW("in 'update_inner_boundary_radial()': case not implemented anymore!\n");
            break;
/*
			VecScale(m_vRadialAtIP[nodeID], bf.normal(), -1.0);
			// rescale
			VecNormalize(m_vRadialAtIP[nodeID], m_vRadialAtIP[nodeID]);
			VecScale(m_vRadialAtIP[nodeID], m_vRadialAtIP[nodeID], radius);
			break;
 */
		case 3:
			VecSubtract(m_vRadialAtIP[nodeID], SCVFip_Pos, center);
			//	VecNormalize(m_vRadialAtIP[ip], m_vRadialAtIP[ip]);
			//	VecScale(m_vRadialAtIP[ip], m_vRadialAtIP[ip], VecLength(bf.normal()));
			break;
		default:
			throw(UGError(
					"Error in IInterfaceMapper::update_inner_boundary_radial()!"));
		}

		if (0) {

			UG_LOG("i = " << i << "\n");
			UG_LOG("bf.node_id() = " << bf.node_id() << "\n");
			UG_LOG(
					"this->corner(" << nodeID << ") = " << this->corner(nodeID) << "\n");
			UG_LOG(
					"vCornerCoords[" << nodeID << "] = " << vCornerCoords[nodeID] << "\n");
			UG_LOG(
					"* m_vRadialAtCo[" << nodeID << "] = " << m_vRadialAtCo[nodeID] << "\n");
			UG_LOG(
					"* m_vRadialAtIP[" << nodeID << "] = " << m_vRadialAtIP[nodeID] << "\n");
			//VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center);
			UG_LOG(
					"** m_vRadialAtIP[" << nodeID << "] = " << m_vRadialAtIP[nodeID] << "\n");

			UG_LOG("end update_inner_boundary_radial()...\n");

		}
	}

}

template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_for2_StdFV()
{
	UG_LOG("\n start update_inner_boundary_radial_for2_StdFV()...\n");

	if (m_prtIndex == -1)
		UG_THROW(
				"InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2_StdFV(): m_prtIndex not set!\n");

	const MathVector<dim> center1 = m_spCutElementHandler->get_center(0);
 	const MathVector<dim> center2 = m_spCutElementHandler->get_center(1);
 
	UG_LOG("center1: " << center1 << "\n");
	UG_LOG("center2: " << center2 << "\n");

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// REMARK: for elements with only ONE outside corner ON the interface, bf.size() = 0!
	// 			=> for assembling connections u_fluid -> u_prt, although 'm_vRadialAtCo' is needed!
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// A. set 'm_vRadialAtCo' --> independent of boundary faces!
	// fill m_vRadialAtCO for center1:

	MathVector<dim> RadialIPRef1;
	MathVector<dim> RadialIPRef2;
	RadialIPRef1 = 0.0;
	RadialIPRef2 = 0.0;
	size_t numOutCo1 = 0;
	size_t numOutCo2 = 0;

	if (dim == 2 && this->numCo() != 3)
		UG_THROW(
				"' in update_inner_boundary_radial_for2_StdFV(): wrong number of corners! -> numCo = " << this->numCo() << "\n");

	for (size_t i = 0; i < this->numCo(); ++i) {
		const int prtIndex = getPrtIndex(i);

		UG_LOG("prtIndex = " << prtIndex << "\n");

		// set 'm_vRadialAtCo'
		if (prtIndex == 0) {
			VecSubtract(m_vRadialAtCo[i], this->corner(i), center1);
			RadialIPRef1 += this->corner(i);
			numOutCo1++;

			if (!this->lies_onInterface(i))
				UG_THROW(
						"inconsistent in 'update_inner_boundary_radial_for2_StdFV()'!\n");
		} else if (prtIndex == 1) {
			VecSubtract(m_vRadialAtCo[i], this->corner(i), center2);
			RadialIPRef2 += this->corner(i);
			numOutCo2++;

			if (!this->lies_onInterface(i))
				UG_THROW(
						"inconsistent in 'update_inner_boundary_radial_for2_StdFV()'!\n");
		}

		if (1) {
			UG_LOG("i = " << i << "\n");

			UG_LOG("this->corner(" << i << ") = " << this->corner(i) << "\n");
			UG_LOG(
					"---> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
		}
	}

	UG_LOG(
			"____________________________________________________________________________\n");
	UG_LOG(
			"____________________________________________________________________________\n");
	UG_LOG("numOutCo1 = " << numOutCo1 << "\n");
	UG_LOG("numOutCo2 = " << numOutCo2 << "\n");

	// set 'm_vRadialAtIP'
	VecScale(RadialIPRef1, RadialIPRef1, 1.0 / numOutCo1);
	VecScale(RadialIPRef2, RadialIPRef2, 1.0 / numOutCo2);
	UG_LOG("RadialIPRef1 = " << RadialIPRef1 << "\n");
	UG_LOG("RadialIPRef2 = " << RadialIPRef2 << "\n");

	for (size_t i = 0; i < this->numCo(); ++i) {
		const int prtIndex = getPrtIndex(i);

		UG_LOG("IP: prtIndex = " << prtIndex << "\n");

		if (prtIndex == 0) {
			VecSubtract(m_vRadialAtIP[i], RadialIPRef1, center1);
			UG_LOG("this->corner(" << i << ") = " << this->corner(i) << "\n");
			UG_LOG(
					"---> m_vRadialAtIP[" << i << "] = " << m_vRadialAtIP[i] << "\n");
		} else if (prtIndex == 1) {
			VecSubtract(m_vRadialAtIP[i], RadialIPRef2, center2);
			UG_LOG("this->corner(" << i << ") = " << this->corner(i) << "\n");
			UG_LOG(
					"---> m_vRadialAtIP[" << i << "] = " << m_vRadialAtIP[i] << "\n");
		}

	}

	UG_LOG("END update_inner_boundary_radial_for2_StdFV()...\n\n");

}

// VORSICHT: update_inner_boundary_radial() wird in fv1FT_geom_impl.h aufgerufen, BEVOR remap_BF() stattfindet!
//	=> bf.node_id() immer mit interfaceID, nicht OriginalCornerID!
template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_for2()
{
	UG_LOG("start update_inner_boundary_radial_for2()...\n");

// get data:
	std::vector<MathVector<dim> > vCornerCoords;
	for (size_t i = 0; i < 4; ++i)
		vCornerCoords.push_back(this->m_vQuadriCorners_for2[i]);

	for (size_t i = 0; i < 4; ++i)
		UG_LOG("vCornerCoords[" << i << "] = " << vCornerCoords[i] << "\n");

	const int modusIP = 1;
	const bool useResized = true;

	if (m_prtIndex == -1)
		UG_THROW(
				"InterfaceHandlerLocalParticle::update_inner_boundary_radial(): m_prtIndex not set!\n");

	if (this->m_vBF.size() != 4)
		UG_THROW(
				"in 'update_inner_boundary_radial_for2()': this->m_vBF.size() should be 4, but is " << this->m_vBF.size() << "\n");

	const MathVector<dim> center1 = m_spCutElementHandler->get_center(0);
 	const MathVector<dim> center2 = m_spCutElementHandler->get_center(1);
 
	UG_LOG("center1: " << center1 << "\n");
	UG_LOG("center2: " << center2 << "\n");

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(4, 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(4, 0.0);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// REMARK: for elements with only ONE outside corner ON the interface, bf.size() = 0!
	// 			=> for assembling connections u_fluid -> u_prt, although 'm_vRadialAtCo' is needed!
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// A. set 'm_vRadialAtCo' --> independent of boundary faces!
	// fill m_vRadialAtCO for center1:
	for (size_t i = 0; i < 2; ++i) {
		UG_LOG(
				"setting for center1 = " << center1 << " and corner[" << i << "]: " << vCornerCoords[i] << "\n");

		if (useResized)
			VecSubtract(m_vRadialAtCo[i], vCornerCoords[i], center1);
		else
			UG_THROW(
					"in 'update_inner_boundary_radial_for2()': not implemented in the for2-case!\n");

		UG_LOG("-> vCornerCoords[" << i << "] = " << vCornerCoords[i] << "\n");
		UG_LOG("-> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
	}
	// fill m_vRadial for center2:
	for (size_t i = 2; i < 4; ++i) {
		UG_LOG(
				"setting for center2 = " << center2 << " and corner[" << i << "]: " << vCornerCoords[i] << "\n");
		if (useResized)
			VecSubtract(m_vRadialAtCo[i], vCornerCoords[i], center2);
		else
			UG_THROW(
					"in 'update_inner_boundary_radial_for2()': not implemented in the for2-case!\n");

		UG_LOG("-> vCornerCoords[" << i << "] = " << vCornerCoords[i] << "\n");
		UG_LOG("-> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// B. set 'm_vRadialAtIP' --> for boundary faces!
	// fill m_vRadialAtCO for center1:
	for (size_t i = 0; i < 2; ++i) {
		interfaceBF bf = this->m_vBF[i];

		size_t nodeID = bf.node_id(); // ok, since bf.nodeID = 0,1,2,3
		if (!useResized)
			UG_THROW(
					"in 'update_inner_boundary_radial_for2()': case '!useResized' not implemented!\n");

		int ip = (bf.node_id() - 1);
		if (ip < 0)
			ip = ip + 3;

		switch (modusIP) {
		case 0: 	// AtIP == AtCo
			m_vRadialAtIP[nodeID] = m_vRadialAtCo[nodeID];
			break;
		case 1: 	// AtIP == bf.global_ip()
			VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center1);
			break;
		default:
			throw(UGError(
					"Error in InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2()!"));
		}

		if (1) {
			UG_LOG("i = " << i << "\n");
			UG_LOG("bf.node_id() = " << bf.node_id() << "\n");
			UG_LOG("bf.normal() = " << bf.normal() << "\n");
			UG_LOG("bf.global_ip()" << bf.global_ip() << "\n");
			UG_LOG(
					"vCornerCoords[" << nodeID << "] = " << vCornerCoords[nodeID] << "\n");
			UG_LOG(
					"* m_vRadialAtCo[" << nodeID << "] = " << m_vRadialAtCo[nodeID] << "\n");
			UG_LOG(
					"** m_vRadialAtIP[" << nodeID << "] = " << m_vRadialAtIP[nodeID] << "\n");
		}
	}

	// fill m_vRadialAtCO for center2:
	for (size_t i = 2; i < 4; ++i) {
		interfaceBF bf = this->m_vBF[i];

		size_t nodeID = bf.node_id(); // ok, since bf.nodeID = 0,1,2,3
		if (!useResized)
			UG_THROW(
					"in 'update_inner_boundary_radial_for2()': case '!useResized' not implemented!\n");

		int ip = (bf.node_id() - 1);
		if (ip < 0)
			ip = ip + 3;

		switch (modusIP) {
		case 0: 	// AtIP == AtCo
			m_vRadialAtIP[nodeID] = m_vRadialAtCo[nodeID];
			break;
		case 1: 	// AtIP == bf.global_ip()
			VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center2);
			break;
		default:
			throw(UGError(
					"Error in InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2()!"));
		}

		if (1) {
			UG_LOG("i = " << i << "\n");
			UG_LOG("bf.node_id() = " << bf.node_id() << "\n");
			UG_LOG("bf.normal() = " << bf.normal() << "\n");
			UG_LOG("bf.global_ip()" << bf.global_ip() << "\n");
			UG_LOG(
					"vCornerCoords[" << nodeID << "] = " << vCornerCoords[nodeID] << "\n");
			UG_LOG(
					"* m_vRadialAtCo[" << nodeID << "] = " << m_vRadialAtCo[nodeID] << "\n");
			UG_LOG(
					"** m_vRadialAtIP[" << nodeID << "] = " << m_vRadialAtIP[nodeID] << "\n");
		}
	}
	UG_LOG("end update_inner_boundary_radial()...\n");

	for (size_t i = 0; i < 4; ++i) {
		UG_LOG("i = " << i << "\n");
		UG_LOG("--> vCornerCoords[" << i << "] = " << vCornerCoords[i] << "\n");
		UG_LOG("--> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
		UG_LOG("--> m_vRadialAtIP[" << i << "] = " << m_vRadialAtIP[i] << "\n");

	}
}

template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_StdFV(const MathVector<TWorldDim>* vCornerCoords)
{

	if (m_prtIndex == -1)
		UG_THROW(
				"InterfaceHandlerLocalParticle::update_inner_boundary_radial(): m_prtIndex not set!\n");

	const MathVector<dim> center = m_spCutElementHandler->get_center(
			m_prtIndex);

	/*
	 typedef DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> > TFVGeom;
	 TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);
	 */
// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

//  Loop the boundary faces to assemble impulse equations
	MathVector<dim> RadialIPRef;
	RadialIPRef = 0.0;
	size_t numOutCo = 0;
	for (size_t i = 0; i < this->numCo(); ++i) {
		// set 'm_vRadialAtCo'
		VecSubtract(m_vRadialAtCo[i], this->corner(i), center);

		// set 'm_vRadialAtIP'
		/*		const typename TFVGeom::SCVF& scvf = geo.scvf(i);
		 const MathVector<dim>& SCVFip_Pos = scvf.global_ip();
		 VecSubtract(m_vRadialAtIP[i], SCVFip_Pos, center);
		 */
		//VecSubtract(m_vRadialAtIP[i], this->corner(i), center);
		// compute radialAtIP for current element:
		if (this->lies_onInterface(i)) {
			RadialIPRef += this->corner(i);
			numOutCo++;
		}

		if (0) {
			UG_LOG("i = " << i << "\n");

			UG_LOG("this->corner(" << i << ") = " << this->corner(i) << "\n");
			UG_LOG(
					"---> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
		}
	}

// set 'm_vRadialAtIP'
	VecScale(RadialIPRef, RadialIPRef, 1.0 / numOutCo);
	for (size_t i = 0; i < this->numCo(); ++i)
		VecSubtract(m_vRadialAtIP[i], RadialIPRef, center);



}

} // end namespace ug



#endif /* INTERFACE_HANDLER_FLAT_TOP_TOOLS_H_ */
