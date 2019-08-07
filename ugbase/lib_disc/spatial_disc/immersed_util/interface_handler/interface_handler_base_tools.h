/*
 * interface_handler_local_tools.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_BASE_TOOLS_H_
#define INTERFACE_HANDLER_LOCAL_BASE_TOOLS_H_

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	FlatTopBase methods - ROID/CollectCorners...
////////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
number InterfaceHandlerLocalBase<TWorldDim>::
get_edge_intersection(Vertex* vrt1, Vertex* vrt2)
{

	number LSvalue1 = get_LSvalue(vrt1, 0);
	number LSvalue2 = get_LSvalue(vrt2, 0);

 	number alpha = LSvalue1 / (LSvalue1 - LSvalue2);

	return alpha;
}

// used for boundary implementation
template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
lies_onInterface(const size_t newID)
{
	for ( size_t i = 0; i < m_vInterfaceID.size(); ++i)
		if ( m_vInterfaceID[i] == newID )
			return true;

	return false;
}


template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
remapped_fromInterface(const size_t origID, size_t& get_interfaceID)
{
 	for ( size_t i = 0; i < m_vQuadriOrigID.size(); ++i)
	{
 		size_t origID_of_interface_corners = m_vQuadriOrigID[i];
		if ( origID_of_interface_corners == origID )
		{
			get_interfaceID = i;
			return true;
		}
  	}

	return false;
}

template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
remapped_fromInterface(const size_t origID)
{
 	for ( size_t i = 0; i < m_vInterfaceID.size(); ++i)
	{
		size_t interfaceID = m_vInterfaceID[i];
		size_t origID_of_interface_corners = m_vOriginalCornerID[interfaceID];
		if ( origID_of_interface_corners == origID )
			return true;
  	}

	return false;
}


template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
is_boundary_face(const size_t sideID)
{
// get data
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

//	number of corners of side (special case bottom side pyramid)
	const int coOfSide = (m_roid != ROID_PYRAMID || sideID != 0)
						? rRefElem.num(dim-1, sideID, 0) : rRefElem.num(dim-1, sideID, 0) + 2;

	for(int co = 0; co < coOfSide; ++co)
	{
		size_t cornerID;
		if (m_roid != ROID_PYRAMID || sideID != 0)
			cornerID = rRefElem.id(dim-1, sideID, 0, co);
		else
			cornerID = rRefElem.id(dim-1, sideID, 0, (co % 3) + (co>3 ? 1 : 0));

		if ( !lies_onInterface(cornerID) )
			return false;
	}

 	return true;

}

/*
// compare implementation of 'DimFV1Geometry::update_boundary_faces()'
template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
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
template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
set_roid_2d()
{
	size_t numCo = m_vCornerCoords.size();

	switch(numCo)
	{
		case 0: m_roid = ROID_UNKNOWN; 			break;
		case 3: m_roid = ROID_TRIANGLE; 		break;
		case 4: m_roid = ROID_QUADRILATERAL;	break;

		default: throw(UGError("during switch in 'set_roid_2d()': "
				"m_numCo invalid => default case chosen"));
	}

}

template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
set_roid_3d()
{
	size_t numCo = m_vCornerCoords.size();

	switch(numCo)
	{
		case 0: m_roid = ROID_UNKNOWN; 		break;
		case 4: m_roid = ROID_TETRAHEDRON; 	break;
		case 5: m_roid = ROID_PYRAMID; 		break;
		case 6: m_roid = ROID_PRISM; 		break;

		default: throw(UGError("during switch in 'set_roid_3d()': "
				"m_numCo invalid => default case chosen"));
	}


}

// parameter 'normal' is the direction from the triangle into the inner of the prism
template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
isCCW(std::vector<MathVector<dim> > vCornerCoords, MathVector<dim> normal)
{

	if ( dim != 3 )
		UG_THROW("stop: in 2d the method VecCross() may not be used!\n");

	MathVector<dim> e1, e2;
	VecSubtract(e1, vCornerCoords[1], vCornerCoords[0]);
	VecSubtract(e2, vCornerCoords[2], vCornerCoords[0]);

	MathVector<dim> n;
	VecCross(n, e1, e2);
	VecNormalize(n,n);
	VecNormalize(normal,normal);

	if ( fabs(VecDot(n, normal)) < 1e-6 )
		UG_THROW("in 'isCCW()': vectors are perpendicular.\n");

	if ( VecDot(n, normal) > 0 ) 	return true; // CCW
	else							return false;

}

template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
isIncluded(std::vector<MathVector<dim> > vCheckList, MathVector<dim> checkPoint)
{
	for(size_t i = 0; i < vCheckList.size(); ++i)
	{
		bool isEqual = true;

		if ( VecDistance(vCheckList[i],checkPoint) > 1e-10)
			isEqual = false;

		if ( isEqual )
			return true;

	}

	return false;
}

template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
ResortQuadrilateral(std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners,
					std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners,
					MathVector<dim> normalDir)
{
// some checks:
	if ( vOutsideCorners.size() != 2 ) UG_THROW("vOutsideCorners.size() has to be 2 in 2d case, but is " << vOutsideCorners.size() << ".\n");
	if ( vInsideCorners.size()  != 2 ) UG_THROW("vInsideCorners.size() has to be 2 in 2d case, but is " << vInsideCorners.size() << ".\n");
	if ( vOutsideCorners[0].second != vOutsideCorners[1].second ) UG_THROW("Outside indices must be identical!\n");

////////////////////////////////////////////////////////////////////////////////////////////
//	1) Copy 'vOutsideCorners'/'vOutsideCorners' to 'm_vCornerCoords' in special order:
// 		 -> ordering: 1. insideCorner, 2./3. outsideCorner, 4. insideCorner
////////////////////////////////////////////////////////////////////////////////////////////

	m_vCornerCoords.clear(); m_vCornerCoords.resize(4);
	m_vOriginalCornerID.clear(); m_vOriginalCornerID.resize(4);
// fill with insideCorners+indices:
	m_vCornerCoords[0] 		= vInsideCorners[0].first;
	m_vOriginalCornerID[0]  = vInsideCorners[0].second;

	m_vCornerCoords[3] 		= vInsideCorners[1].first;
	m_vOriginalCornerID[3]  = vInsideCorners[1].second;

// fill with outsideCorners+indices:
	m_vCornerCoords[1] 		= vOutsideCorners[0].first;
	m_vOriginalCornerID[1]  = vOutsideCorners[0].second;

	m_vCornerCoords[2] 		= vOutsideCorners[1].first;
	m_vOriginalCornerID[2]  = vOutsideCorners[1].second;

// remember outside Index in order to fill 'vm_vInterfaceID' after re-ordering.
	size_t outInd = vOutsideCorners[0].second;

////////////////////////////////////////////////////////////////////////////////////////////
//	2) Check for potential cross-ordering or not-CCW
////////////////////////////////////////////////////////////////////////////////////////////
	MathVector<3> vNormOut1, vNormOut2;
	std::vector<MathVector<3> > vCornerCoordsBlow; vCornerCoordsBlow.resize(4);
	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
	{
		for ( size_t d = 0; d < dim; ++d )
			vCornerCoordsBlow[i][d] = m_vCornerCoords[i][d];
		vCornerCoordsBlow[i][2] = 0.0;
	}

	CalculateTriangleNormalNoNormalize(vNormOut1, vCornerCoordsBlow[0], vCornerCoordsBlow[1], vCornerCoordsBlow[2]);
	CalculateTriangleNormalNoNormalize(vNormOut2, vCornerCoordsBlow[2], vCornerCoordsBlow[3], vCornerCoordsBlow[0]);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Check 1 - cross-ordering
// VecDot(vNormOut1,vNormOut2) < 0 => corners of quadrilateral are not ordered along the
//		endges, i.e. across the diagonal => switch order of the two outideCorners, i.e. index [1],[2]
/////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( VecDot(vNormOut1, vNormOut2) < 0 )
	{
 		MathVector<dim> bufferCoord = m_vCornerCoords[1];
		size_t bufferInd = m_vOriginalCornerID[1];

		m_vCornerCoords[1] 		= m_vCornerCoords[2];
		m_vOriginalCornerID[1]  = m_vOriginalCornerID[2];

		m_vCornerCoords[2] 		= bufferCoord;
		m_vOriginalCornerID[2]  = bufferInd;
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Check 2 - CCW?
// VecDot(vNormOut1,zDir) < 0 => corners of quadrilateral are not CCW:
/////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
	{
		for ( size_t d = 0; d < dim; ++d )
			vCornerCoordsBlow[i][d] = m_vCornerCoords[i][d];
		if ( dim == 2 ) vCornerCoordsBlow[i][2] = 0.0;


	}
	CalculateTriangleNormalNoNormalize(vNormOut1, vCornerCoordsBlow[0], vCornerCoordsBlow[1], vCornerCoordsBlow[2]);
	CalculateTriangleNormalNoNormalize(vNormOut2, vCornerCoordsBlow[2], vCornerCoordsBlow[3], vCornerCoordsBlow[0]);

	if ( VecDot(vNormOut1, vNormOut2) < 0 ) UG_THROW("check 1 failed after second check!\n");


// pure copy of 'normalDir' data because of data type inconsistence (MathVector<3> vs. MathVector<dim>)
	MathVector<3> _normalDir;
	for ( size_t i = 0; i < 2; ++i )
		_normalDir[i] = normalDir[i];
	if ( dim == 2 ) _normalDir[2] = 1.0;
	if ( dim == 3 ) _normalDir[2] = normalDir[2];

	if ( VecDot(vNormOut1, _normalDir) < 0 )
	{

 		MathVector<dim> bufferCoord = m_vCornerCoords[1];
		size_t bufferInd = m_vOriginalCornerID[1];

		m_vCornerCoords[1] 		= m_vCornerCoords[3];
		m_vOriginalCornerID[1]  = m_vOriginalCornerID[3];

		m_vCornerCoords[3] 	= bufferCoord;
		m_vOriginalCornerID[3]  = bufferInd;
	}

// fill 'vInterfaceID':
	m_vInterfaceID.clear();
	for ( size_t i = 0; i < m_vOriginalCornerID.size(); ++i )
		if ( m_vOriginalCornerID[i] == outInd )
			m_vInterfaceID.push_back(i);

}



template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_StdFV(GridObject* elem)
{
	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);


	bool isOutsideVertex = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
 		if ( is_OutsideVertex(vVertex[i], i) || is_FTVertex(vVertex[i], i) )
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

		//////////////////////////////////////////////
		// case 1:
		// vertex insideFluid => Position und index puschen
		if ( !is_FTVertex(vVertex[i], i) )
		{
			if ( is_nearInterfaceVertex(vrtRoot, i) )
			{
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");
			}
			else
			{
				m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);

				vInsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));
			}

		}
		else
		{
 			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));

		}

	}

	// skip whole element, since only FT points are included
 	if ( dim == 2 && m_vCornerCoords.size() != 3 )
 		UG_THROW("in 'CollectCorners_StdFV()': m_vCornerCoords.size() "
					"= " << m_vCornerCoords.size() << "not possible!\n");
 	if ( dim == 3 && m_vCornerCoords.size() != 4 )
 		UG_THROW("in 'CollectCorners_StdFV()': m_vCornerCoords.size() "
					"= " << m_vCornerCoords.size() << "not possible!\n");

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
		fprintf(outputFile,"%e, %e, %e \n", m_aaPos[vVertex[i]][0] , m_aaPos[vVertex[i]][1], m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: m_vCornerCoords.size() %lu \n", m_vCornerCoords.size());

	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1], m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vInterfaceID = %lu\n", m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vOriginalCornerID = %lu\n", m_vOriginalCornerID[i]);


	fprintf(outputFile,"\n_________________________ end 'CollectIBCorners_StdFV()' ______________________________\n\n");

	fclose(outputFile);



	return m_vCornerCoords.size();
}



template <int TWorldDim>
size_t InterfaceHandlerLocalBase<TWorldDim>::
get_vertex_index(Vertex* vrt, GridObject* elem)
{
  	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

 	for(size_t i = 0; i < vVertex.size(); ++i)
		if ( vrt == vVertex[i])
			return i;

	UG_THROW("in CutElementHandlerImmersed::get_vertex_index: no index found!\n");
}

/*
template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_2d(GridObject* elem)
{
	bool bOutput = false;

	if ( bOutput ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_2d()' ______________________________\n\n");

	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

 	bool isFTVertex = false;
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
		isFTVertex = is_FTVertex(vVertex[i], i);
		if ( isFTVertex )
			break;
	}

	if ( !isFTVertex )
		UG_THROW("Error in 'CollectCorners_FlatTop_2d': no vertex is FTVertex: should be true for at least 1 vertex!\n");

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *m_spMG, elem);

	// loop vertices
	//////////////////////////////////////////////
	// REMARK:
	// order is the same as in 'vCornerCoords', therefore we can be sure, that the
	// order of the new 'vCornerIBCoords' will be consistent with the grid standard
	//////////////////////////////////////////////

	bool bNearInterface = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{

		// get element
		Vertex* vrtRoot = vVertex[i];

		//////////////////////////////////////////////
		// case 1:
		// vertex insideFluid
		// 		=> Position und index puschen
		if ( ! is_FTVertex(vrtRoot, i) )
		{

			if ( is_nearInterfaceVertex(vrtRoot, i) )
			{
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");
			}
			else
			{
				m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);

				vInsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));
			}

		}
		//////////////////////////////////////////////
  		// case 2:
		// vertex = FT + ON interface
		// 		=> KEINE Berechnung von 'intersectionPoint' notwendig! -> pushen und alten index pushen

		// REMARK: is_nearInterfaceVerx = false per default, if m_vThresholdOnLevel = 0.0
		else if ( is_nearInterfaceVertex(vrtRoot, i) )
		{
 			bNearInterface = true;
			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));
			vNearIntCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));

		}
		//////////////////////////////////////////////
  		// case 3:
		// vertex 'outsideFluid'
		// 		=> NEUE Position berechen+pushen und alten index pushen
		else
		{
 			//////////////////////////////////////////////////////////////////////////////////////////
			// loop alle edges, die interface schneiden und damit einen neuen intersectionPnt
			// beitragen zum damit assoziierten alten index
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

				MathVector<dim> intersectionPnt;

	 		///////////////////////////////////////////////////////////////////
 			// lies vrtRoot on a cutted edge?
		 	///////////////////////////////////////////////////////////////////
			// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
				if ( is_nearInterfaceVertex(vrt2, vrtInd2) || is_nearInterfaceVertex(vrt1, vrtInd1) )
				{ bNearInterface = true; continue; }
			 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt1 && !is_FTVertex(vrt2, vrtInd2) ){
					get_intersection_point(intersectionPnt, vrt2, vrt1);
 				}
			// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt2 && !is_FTVertex(vrt1, vrtInd1) )
					get_intersection_point(intersectionPnt, vrt1, vrt2);
				else
 					continue;


			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");


	 		///////////////////////////////////////////////////////////////////
	 		// only push_back, if not included yet!
			// 	-> can be ONLY the case, if the intersectionPoint is a node
	 			if ( ! isIncluded(m_vCornerCoords, intersectionPnt) )
	 			{

	 				m_vCornerCoords.push_back(intersectionPnt);
	 				m_vOriginalCornerID.push_back(i);
 	 				m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

	 				vOutsideCorners.push_back(std::make_pair(intersectionPnt, i));
   	 			}


 			} // end edge-loop

		} // end else-case

 	} // end vrt-loop

////////////////////////////////////////////////////////////////////////////////////////////
// Postprecessing for quadrilaterals ( <=>  vOutsideCorners == 2 )
// (vInsideCorners.size() == 2) && (bNearInterface)	 => ALL nodes insideFluid, BUT one ON surface
//		=> no Quadrilateral, but Triangle!!
////////////////////////////////////////////////////////////////////////////////////////////
	MathVector<dim> normalDir(0.0);
	if ( (m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
 		ResortQuadrilateral(vInsideCorners, vOutsideCorners, normalDir);
	else if ( bNearInterface )
	{
	// Quadrilateral -> Triangle
		if ( vInsideCorners.size() == 1 ) // case 1
		{
			// do nothing, since re-sorting not necessary...???
		}
	// skip whole element, since only FT points are included
		else if ( vInsideCorners.size() == 0 )
			UG_THROW("in 'CollectCorners_FlatTop_2d()': vInsideCorners.size() "
					"= " << vInsideCorners.size() << "not possible!\n");
	}


	if ( bOutput)
	{
////////////////////////////////////////////////////////////////////////////////////////////
// Wirte to files:
////////////////////////////////////////////////////////////////////////////////////////////
	if ( (m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
		{
 			const char* filename = "File_Quadrilateral.";
				std::string name(filename);
		 	char ext[50]; sprintf(ext, "txt");
			name.append(ext);
			FILE* printFile = fopen(name.c_str(), "a");
			fprintf(printFile, "------------------------------------------------------\n\n new triangle: \n");
			for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
				fprintf(printFile," vCornerCoords: %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1]);
				for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
				fprintf(printFile,"OriginalCornerID = %lu\n", m_vOriginalCornerID[i]);
			for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
				fprintf(printFile,"m_vInterfaceID = %lu\n", m_vInterfaceID[i]);
			for ( size_t i = 0; i < vInsideCorners.size(); ++i )
				fprintf(printFile," vInsideCorners: %e, %e \n", vInsideCorners[i].first[0], vInsideCorners[i].first[1]);
			for ( size_t i = 0; i < vNearIntCorners.size(); ++i )
				fprintf(printFile," vNearIntCorners: %e, %e \n", vNearIntCorners[i].first[0], vNearIntCorners[i].first[1]);
			for ( size_t i = 0; i < vOutsideCorners.size(); ++i )
				fprintf(printFile," vOutsideCorners: %e, %e \n", vOutsideCorners[i].first[0], vOutsideCorners[i].first[1]);
			fclose(printFile);
		}
		else if ( bNearInterface )
		{
		// Quadrilateral -> Triangle
			if ( vInsideCorners.size() == 1 ) // case 1
			{
				// do nothing, since re-sorting not necessary...???

				const char* filename = "File_newTirangle.";
	  			std::string name(filename);
			 	char ext[50]; sprintf(ext, "txt");
				name.append(ext);
				FILE* printFile = fopen(name.c_str(), "a");
				fprintf(printFile, "------------------------------------------------------\n\n new triangle: \n");
				for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
					fprintf(printFile," vCornerCoords: %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1]);
	 			for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
					fprintf(printFile,"OriginalCornerID = %lu\n", m_vOriginalCornerID[i]);
				for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
					fprintf(printFile,"m_vInterfaceID = %lu\n", m_vInterfaceID[i]);
				for ( size_t i = 0; i < vInsideCorners.size(); ++i )
					fprintf(printFile," vInsideCorners: %e, %e \n", vInsideCorners[i].first[0], vInsideCorners[i].first[1]);
				for ( size_t i = 0; i < vNearIntCorners.size(); ++i )
					fprintf(printFile," vNearIntCorners: %e, %e \n", vNearIntCorners[i].first[0], vNearIntCorners[i].first[1]);
				for ( size_t i = 0; i < vOutsideCorners.size(); ++i )
					fprintf(printFile," vOutsideCorners: %e, %e \n", vOutsideCorners[i].first[0], vOutsideCorners[i].first[1]);
				fclose(printFile);
			}

		}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	if (bOutput)
	{
	const char* filename = "output.";
	std::string name(filename);
 	char ext[50];
	if ( dim == 2 )
		sprintf(ext, "2d.txt");
	if ( dim == 3 )
		sprintf(ext, "3d.txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");

	fprintf(outputFile,"\n_________________________ begin 'CollectIBCorners2d' ________________________\n");

	for ( size_t i = 0; i < vVertex.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_aaPos[vVertex[i]][0] , m_aaPos[vVertex[i]][1], m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: m_vCornerCoords.size() %lu \n", m_vCornerCoords.size());

	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1], m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vInterfaceID = %lu\n", m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vOriginalCornerID = %lu\n", m_vOriginalCornerID[i]);


	fprintf(outputFile,"\n_________________________ end 'CollectIBCorners2d()' ______________________________\n\n");

	fclose(outputFile);

	}
	}
	if ( bOutput ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_2d()' ______________________________\n\n");

	return m_vCornerCoords.size();

}
*/




template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_2d(GridObject* elem)
{
	bool bOutput = false;
    bool bOutputToFile = false;
	if ( bOutput ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_2d()' ______________________________\n\n");

	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

 	bool isFTVertex = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
		isFTVertex = is_FTVertex(vVertex[i], i);
		if ( isFTVertex )
			break;
	}
    

	if ( !isFTVertex )
		UG_THROW("Error in 'CollectCorners_FlatTop_2d': no vertex is FTVertex: should be true for at least 1 vertex!\n");

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *m_spMG, elem);

	// loop vertices
	//////////////////////////////////////////////
	// REMARK:
	// order is the same as in 'vCornerCoords', therefore we can be sure, that the
	// order of the new 'vCornerIBCoords' will be consistent with the grid standard
	//////////////////////////////////////////////

	bool bNearInterface = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{

		// get element
		Vertex* vrtRoot = vVertex[i];

        if ( (fabs(m_aaPos[vVertex[i]][0] - 1.08941) < 1e-4) && (fabs(m_aaPos[vVertex[i]][1] + 0.140733) < 1e-4) )
            UG_LOG("stop here..\n");

		//////////////////////////////////////////////
		// case 1:
		// vertex insideFluid
		// 		=> Position und index puschen
		if ( ! is_FTVertex(vrtRoot, i) )
		{

			if ( is_nearInterfaceVertex(vrtRoot, i) )
			{
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");
			}
			else
			{
				m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);

				vInsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));
			}

		}
		//////////////////////////////////////////////
  		// case 2:
		// vertex = FT + ON interface
		// 		=> KEINE Berechnung von 'intersectionPoint' notwendig! -> pushen und alten index pushen

		// REMARK: is_nearInterfaceVerx = false per default, if m_vThresholdOnLevel = 0.0
		else if ( is_nearInterfaceVertex(vrtRoot, i) )
		{
 			bNearInterface = true;
			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));
			vNearIntCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));

		}
		//////////////////////////////////////////////
  		// case 3:
		// vertex 'outsideFluid'
		// 		=> NEUE Position berechen+pushen und alten index pushen
		else
		{
 			//////////////////////////////////////////////////////////////////////////////////////////
			// loop alle edges, die interface schneiden und damit einen neuen intersectionPnt
			// beitragen zum damit assoziierten alten index
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

				MathVector<dim> intersectionPnt;

	 		///////////////////////////////////////////////////////////////////
 			// lies vrtRoot on a cutted edge?
		 	///////////////////////////////////////////////////////////////////
			// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
				if ( is_nearInterfaceVertex(vrt2, vrtInd2) || is_nearInterfaceVertex(vrt1, vrtInd1) )
				{ bNearInterface = true; continue; }
			 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt1 && !is_FTVertex(vrt2, vrtInd2) ){
					get_intersection_point(intersectionPnt, vrt2, vrt1);
 				}
			// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt2 && !is_FTVertex(vrt1, vrtInd1) )
					get_intersection_point(intersectionPnt, vrt1, vrt2);
				else
 					continue;


			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");


	 		///////////////////////////////////////////////////////////////////
	 		// only push_back, if not included yet!
			// 	-> can be ONLY the case, if the intersectionPoint is a node
	 			if ( ! isIncluded(m_vCornerCoords, intersectionPnt) )
	 			{

	 				m_vCornerCoords.push_back(intersectionPnt);
	 				m_vOriginalCornerID.push_back(i);
 	 				m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

	 				vOutsideCorners.push_back(std::make_pair(intersectionPnt, i));
   	 			}


 			} // end edge-loop

		} // end else-case

 	} // end vrt-loop

////////////////////////////////////////////////////////////////////////////////////////////
// Postprecessing for quadrilaterals ( <=>  vOutsideCorners == 2 )
// (vInsideCorners.size() == 2) && (bNearInterface)	 => ALL nodes insideFluid, BUT one ON surface
//		=> no Quadrilateral, but Triangle!!
////////////////////////////////////////////////////////////////////////////////////////////
	MathVector<dim> normalDir(0.0);
	if ( (m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
 		ResortQuadrilateral(vInsideCorners, vOutsideCorners, normalDir);
	else if ( bNearInterface )
	{
	// Quadrilateral -> Triangle
		if ( vInsideCorners.size() == 1 ) // case 1
		{
			// do nothing, since re-sorting not necessary...???
		}
	// skip whole element, since only FT points are included
		else if ( vInsideCorners.size() == 0 )
			UG_THROW("in 'CollectCorners_FlatTop_2d()': vInsideCorners.size() "
					"= " << vInsideCorners.size() << "not possible!\n");
	}


////////////////////////////////////////////////////////////////////////////////////////////
// Wirte to files:
////////////////////////////////////////////////////////////////////////////////////////////
if ( bOutputToFile )
{
	if ( (m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
		{
 			const char* filename = "File_Quadrilateral.";
				std::string name(filename);
		 	char ext[50]; sprintf(ext, "txt");
			name.append(ext);
			FILE* printFile = fopen(name.c_str(), "a");
			fprintf(printFile, "------------------------------------------------------\n\n new triangle: \n");
			for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
				fprintf(printFile," vCornerCoords: %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1]);
				for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
				fprintf(printFile,"OriginalCornerID = %lu\n", m_vOriginalCornerID[i]);
			for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
				fprintf(printFile,"m_vInterfaceID = %lu\n", m_vInterfaceID[i]);
			for ( size_t i = 0; i < vInsideCorners.size(); ++i )
				fprintf(printFile," vInsideCorners: %e, %e \n", vInsideCorners[i].first[0], vInsideCorners[i].first[1]);
			for ( size_t i = 0; i < vNearIntCorners.size(); ++i )
				fprintf(printFile," vNearIntCorners: %e, %e \n", vNearIntCorners[i].first[0], vNearIntCorners[i].first[1]);
			for ( size_t i = 0; i < vOutsideCorners.size(); ++i )
				fprintf(printFile," vOutsideCorners: %e, %e \n", vOutsideCorners[i].first[0], vOutsideCorners[i].first[1]);
			fclose(printFile);
		}
		else if ( bNearInterface )
		{
		// Quadrilateral -> Triangle
			if ( vInsideCorners.size() == 1 ) // case 1
			{
				// do nothing, since re-sorting not necessary...???

				const char* filename = "File_newTirangle.";
	  			std::string name(filename);
			 	char ext[50]; sprintf(ext, "txt");
				name.append(ext);
				FILE* printFile = fopen(name.c_str(), "a");
				fprintf(printFile, "------------------------------------------------------\n\n new triangle: \n");
				for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
					fprintf(printFile," vCornerCoords: %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1]);
	 			for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
					fprintf(printFile,"OriginalCornerID = %lu\n", m_vOriginalCornerID[i]);
				for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
					fprintf(printFile,"m_vInterfaceID = %lu\n", m_vInterfaceID[i]);
				for ( size_t i = 0; i < vInsideCorners.size(); ++i )
					fprintf(printFile," vInsideCorners: %e, %e \n", vInsideCorners[i].first[0], vInsideCorners[i].first[1]);
				for ( size_t i = 0; i < vNearIntCorners.size(); ++i )
					fprintf(printFile," vNearIntCorners: %e, %e \n", vNearIntCorners[i].first[0], vNearIntCorners[i].first[1]);
				for ( size_t i = 0; i < vOutsideCorners.size(); ++i )
					fprintf(printFile," vOutsideCorners: %e, %e \n", vOutsideCorners[i].first[0], vOutsideCorners[i].first[1]);
				fclose(printFile);
			}

		}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	const char* filename = "output.";
	std::string name(filename);
 	char ext[50];
	if ( dim == 2 )
		sprintf(ext, "2d.txt");
	if ( dim == 3 )
		sprintf(ext, "3d.txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");

	fprintf(outputFile,"\n_________________________ begin 'CollectIBCorners2d' ________________________\n");

	for ( size_t i = 0; i < vVertex.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_aaPos[vVertex[i]][0] , m_aaPos[vVertex[i]][1], m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: m_vCornerCoords.size() %lu \n", m_vCornerCoords.size());

	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1], m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vInterfaceID = %lu\n", m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vOriginalCornerID = %lu\n", m_vOriginalCornerID[i]);


	fprintf(outputFile,"\n_________________________ end 'CollectIBCorners2d()' ______________________________\n\n");

	fclose(outputFile);


} // end if ( bOutputToFile )

	if ( bOutput ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_2d()' ______________________________\n\n");

	return m_vCornerCoords.size();

}

template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_Prism3(GridObject* elem)
{
///////////////////////////////////////////////////////////////////
// 	1) push_back into 'vInsideCorners'
// 	2) collect according outsideCorner (lying on the commom edge)
// 	3) push_back into 'vOutsideCorners'
///////////////////////////////////////////////////////////////////

	bool output = false;

	if ( output ) UG_LOG("_________________________ begin 'CollectIBCorners3d_Prism3()' ______________________________\n\n");

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

 	std::vector<MathVector<dim> > vOutsideCorners;
	std::vector<MathVector<dim> > vInsideCorners;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *m_spMG, elem);

	// loop vertices
 	size_t outsideInd;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		// get element
		Vertex* vrtRoot = vVertex[i];

		if ( ! is_FTVertex(vrtRoot, i) )
		{
	///////////////////////////////////////////////////////////////////
	// 	1) push_back into 'vInsideCorners'
			vInsideCorners.push_back(m_aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);

	///////////////////////////////////////////////////////////////////
	// 	2) collect according outsideCorner (lying on the commom edge)
			bool outsideCornerFound = false;
			MathVector<dim> intersectionPnt;
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("in 'NEWcall_Extrapolation': error in collecting vertices associated to an edge!....EXIT!...\n");
				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

			// lies vrtRoot on a cutted edge?
			// case1: vrt1 = outsideParticle && vrt2 = insideParticle:
				if ( vrtRoot == vrt1 && is_FTVertex(vrt2, vrtInd2) )
  					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt1, vrt2);}
			// case2: vrt2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt2 && is_FTVertex(vrt1, vrtInd1) )
					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt2, vrt1);}
			// NO nearInterface handling necessary, since this is Pyramid or Tetrahedron case!
				else if ( is_nearInterfaceVertex(vrt2, vrtInd2) || is_nearInterfaceVertex(vrt1, vrtInd1) )
					{UG_THROW("in 'CollectCorners_Prism3()': no nearInterface case possible!\n");}
				else
					continue;

 			// check for correct intersectionPnt
				if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners3d_Prism3()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");

			} // end edge-loop

			if ( !outsideCornerFound )
				UG_THROW("During edge-loop no outsideCorner found!\n");

	///////////////////////////////////////////////////////////////////
	// 	3) push_back into 'vOutsideCorners'
	//		check for 'isIndluded(), BUT: ALL 3 intersectionPonits need to be included for PRISM!!
			if ( isIncluded(vOutsideCorners, intersectionPnt) )
				UG_THROW("intersectionPnt allready included? ERROR\n");
			vOutsideCorners.push_back(intersectionPnt);

		} // end inside-case
		else
		// recall ousideIndex for 'm_vOriginalCornerID'
			outsideInd = i;

	} // end vertex-loop

	for(size_t i = 0; i < 3; ++i)
	{
		m_vOriginalCornerID.push_back(outsideInd);
		m_vInterfaceID.push_back(3+i);
	}

// WRITE 'vCornerCoords': FIRST insideCorners, SECOND outsideCorners
	m_vCornerCoords.clear();
	for(size_t i = 0; i < 3; ++i)
		m_vCornerCoords.push_back(vInsideCorners[i]);
	for(size_t i = 0; i < 3; ++i)
		m_vCornerCoords.push_back(vOutsideCorners[i]);

// some checks:
	if ( m_vInterfaceID.size() 	!= 3 ) UG_THROW("wrong size for m_vInterfaceID! Should be 3 and is " << m_vInterfaceID.size() << "\n");
	if ( m_vOriginalCornerID.size()!= 6 ) UG_THROW("wrong size for m_vOriginalCornerID! Should be 6 and is " << m_vOriginalCornerID.size() << "\n");
	if ( vInsideCorners.size() 	!= 3 ) UG_THROW("wrong size for vInsideCorners! Should be 3 and is " << vInsideCorners.size() << "\n");
	if ( vOutsideCorners.size() != 3 ) UG_THROW("wrong size for vOutsideCorners! Should be 3 and is " << vOutsideCorners.size() << "\n");
	if ( m_vCornerCoords.size()!= 6 ) UG_THROW("wrong size for m_vCornerCoords! Should be 6 and is " << m_vCornerCoords.size() << "\n");

	if ( output )
	{
		UG_LOG("---------- vor SWITCH ---------- \n");
		for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
			UG_LOG("m_vCornerCoords = " << m_vCornerCoords[i] << "\n");
		for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
			UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	}
// check 'vOutsideCorners'-array for CCW:
// get normal to triangle, built by 'vOutsideCorners'-array:
	MathVector<dim> normal = vInsideCorners[0];
 	VecSubtract(normal, vOutsideCorners[0], normal);

	if ( !isCCW(vOutsideCorners, normal) )
	{
	// switch '	m_vCornerCoords'
 		m_vCornerCoords[4] = vOutsideCorners[2];
		m_vCornerCoords[5] = vOutsideCorners[1];
	// switch '	m_vOriginalCornerID'
		size_t buffer = m_vOriginalCornerID[4];
		m_vOriginalCornerID[4] = m_vOriginalCornerID[5];
		m_vOriginalCornerID[5] = buffer;
		if ( output ) UG_LOG("SWITCH vOutsideCorneers\n");
	}

// check 'vInsideCorners'-array for CCW:
	if ( !isCCW(vInsideCorners, normal) )
	{
	// switch '	m_vCornerCoords'
 		m_vCornerCoords[1] = vInsideCorners[2];
		m_vCornerCoords[2] = vInsideCorners[1];
	// switch '	m_vOriginalCornerID'
		size_t buffer = m_vOriginalCornerID[1];
		m_vOriginalCornerID[1] = m_vOriginalCornerID[2];
		m_vOriginalCornerID[2] = buffer;
		if ( output ) UG_LOG("SWITCH vInsideCorneers\n");
 	}

	if ( output )
	{
		UG_LOG("---------- nach SWITCH ---------- \n");
		for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
			UG_LOG("m_vCornerCoords = " << m_vCornerCoords[i] << "\n");
		for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
			UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	}

/*
 	for(size_t i = 0; i < 3; ++i)
		UG_LOG("vOutsideCorners[i] = " << vOutsideCorners[i] << "\n");
	for(size_t i = 0; i < 3; ++i)
		UG_LOG("vInsideCorners[i] = " << vInsideCorners[i] << "\n");
	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		UG_LOG("m_vInterfaceID[i] = " << m_vInterfaceID[i] << "\n");
	for(size_t i = 0; i < 6; ++i)
		UG_LOG("m_vCornerCoords[i] = " << m_vCornerCoords[i] << "\n");
*/

	if ( output )
	{
	const char* filename = "output3d.";
	std::string name(filename);
 	char ext[50]; sprintf(ext, "txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");


	fprintf(outputFile,"\n_________________________ begin 'Prism3' ________________________\n");

	for ( size_t i = 0; i < vVertex.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_aaPos[vVertex[i]][0] , m_aaPos[vVertex[i]][1], m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: m_vCornerCoords.size() %lu \n", m_vCornerCoords.size());

	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1], m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vInterfaceID = %lu\n", m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vOriginalCornerID = %lu\n", m_vOriginalCornerID[i]);


	fprintf(outputFile,"\n_________________________ end 'Prism3' ______________________________\n\n");

	fclose(outputFile);
	//UG_THROW("first check: Prism3\n");
	}
	return m_vCornerCoords.size();

}

template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_Prism4(GridObject* elem)
{
///////////////////////////////////////////////////////////////////
// 	1) push_back 'insideCorner' into 'm_vCornerCoords'
// 	2) collect according outsideCorner (lying on the commom edge)
// 	3) push_back 'intersectionPoints' into 'm_vCornerCoords'
//  4) check for CCW-property
///////////////////////////////////////////////////////////////////

	bool output = false;
	if ( output ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_Prism4()' ______________________________\n\n");

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *m_spMG, elem);

	// collect outside indices:
	std::vector<size_t> outsideInd;
	for(size_t i = 0; i < vVertex.size(); ++i)
		if ( is_FTVertex(vVertex[i], i) )
			outsideInd.push_back(i);

	if ( outsideInd.size() != 2 ) UG_THROW("Wrong number of outsideInd!\n");
	for(size_t i = 0; i < outsideInd.size(); ++i)
		if ( output ) UG_LOG("outsideInd = " << outsideInd[i] << "\n");

	// loop vertices
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrtRoot = vVertex[i];
		if ( ! is_FTVertex(vrtRoot, i) )
		{
	///////////////////////////////////////////////////////////////////
	// 	1) push_back 'insideCorner' into 'm_vCornerCoords'
			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);

	///////////////////////////////////////////////////////////////////
	// 	2) collect according outsideCorner (lying on the commom edge)
 			MathVector<dim> intersectionPnt;
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("error in collecting vertices associated to an edge!\n");
				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

			// lies vrtRoot on a cutted edge?
			// case1: vrt1 = outsideParticle && vrt2 = insideParticle:
				if ( vrtRoot == vrt1 && is_FTVertex(vrt2, vrtInd2) )
				{
	///////////////////////////////////////////////////////////////////
	// 	3) push_back 'intersectionPnt' into 'm_vCornerCoords'
					get_intersection_point(intersectionPnt, vrt1, vrt2);
					m_vCornerCoords.push_back(intersectionPnt);

					if ( vrt2 == vVertex[outsideInd[0]] )
						m_vOriginalCornerID.push_back(outsideInd[0]);
					else if ( vrt2 == vVertex[outsideInd[1]] )
						m_vOriginalCornerID.push_back(outsideInd[1]);
					else UG_THROW("No outsideInd found.\n");
				}
			// case2: vrt1 = insideParticle && vrt2 = outsideParticle:
				else if ( vrtRoot == vrt2 && is_FTVertex(vrt1, vrtInd1) )
				{
	///////////////////////////////////////////////////////////////////
	// 	3) push_back 'intersectionPnt' into 'm_vCornerCoords'
					get_intersection_point(intersectionPnt, vrt2, vrt1);
					m_vCornerCoords.push_back(intersectionPnt);

					if ( vrt1 == vVertex[outsideInd[0]] )
						m_vOriginalCornerID.push_back(outsideInd[0]);
					else if ( vrt1 == vVertex[outsideInd[1]] )
						m_vOriginalCornerID.push_back(outsideInd[1]);
					else UG_THROW("No outsideInd found.\n");
				}
				else
					continue;

			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners3d_Prism2()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = "
							<< fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");
			} // end edge-loop
		} // end inside-case
 	} // end vertex-loop

	m_vInterfaceID.push_back(1);
	m_vInterfaceID.push_back(2);
	m_vInterfaceID.push_back(4);
	m_vInterfaceID.push_back(5);


// some checks:
	if ( m_vInterfaceID.size() 	!= 4 ) UG_THROW("wrong size for m_vInterfaceID! Should be 4 and is " << m_vInterfaceID.size() << "\n");
	if ( m_vOriginalCornerID.size()!= 6 ) UG_THROW("wrong size for m_vOriginalCornerID! Should be 6 and is " << m_vOriginalCornerID.size() << "\n");
 	if ( m_vCornerCoords.size()!= 6 ) UG_THROW("wrong size for m_vCornerCoords! Should be 6 and is " << m_vCornerCoords.size() << "\n");

///////////////////////////////////////////////////////////////////
//  4) check for CCW-property

	if ( output )
	{
		UG_LOG("---------- vor SWITCH ---------- \n");
		for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
			UG_LOG("m_vCornerCoords = " << m_vCornerCoords[i] << "\n");
		for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
			UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	}

	std::vector<MathVector<dim> > vTriangleCorners1;
	std::vector<MathVector<dim> > vTriangleCorners2;
	for(size_t i = 0; i < 3; ++i){
		vTriangleCorners1.push_back(m_vCornerCoords[i]);
		vTriangleCorners2.push_back(m_vCornerCoords[3+i]);
	}
// get normal to triangle, built by 'vOutsideCorners'-array:
	MathVector<dim> normal = m_vCornerCoords[0];
	VecSubtract(normal, m_vCornerCoords[3], normal);

	if ( !isCCW(vTriangleCorners1, normal) )
	{
	// switch '	m_vCornerCoords'
		m_vCornerCoords[1] = vTriangleCorners1[2];
		m_vCornerCoords[2] = vTriangleCorners1[1];
	// switch '	m_vOriginalCornerID'
		size_t buffer = m_vOriginalCornerID[1];
		m_vOriginalCornerID[1] = m_vOriginalCornerID[2];
		m_vOriginalCornerID[2] = buffer;
		if ( output ) UG_LOG("SWITCH vTriangleCorners1\n");

	}

// check 'vInsideCorners'-array for CCW:
	if ( !isCCW(vTriangleCorners2, normal) )
	{
	// switch '	m_vCornerCoords'
		m_vCornerCoords[4] = vTriangleCorners2[2];
		m_vCornerCoords[5] = vTriangleCorners2[1];
	// switch '	m_vOriginalCornerID'
		size_t buffer = m_vOriginalCornerID[4];
		m_vOriginalCornerID[4] = m_vOriginalCornerID[5];
		m_vOriginalCornerID[5] = buffer;
		if ( output ) UG_LOG("SWITCH vTriangleCorners2\n");
	}

	if ( output )
	{
		UG_LOG("---------- nach SWITCH ---------- \n");
		for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
			UG_LOG("m_vCornerCoords = " << m_vCornerCoords[i] << "\n");
		for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
			UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	}

/*
	for(size_t i = 0; i < 3; ++i)
		UG_LOG("vTriangleCorners1[i] = " << vTriangleCorners1[i] << "\n");
	for(size_t i = 0; i < 3; ++i)
		UG_LOG("vTriangleCorners2[i] = " << vTriangleCorners2[i] << "\n");
	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		UG_LOG("m_vInterfaceID[i] = " << m_vInterfaceID[i] << "\n");
	for(size_t i = 0; i < 6; ++i)
		UG_LOG("m_vCornerCoords[i] = " << m_vCornerCoords[i] << "\n");
*/

	if ( output )
	{
		UG_LOG("_________________________ end 'CollectCorners_FlatTop_Prism4()' ______________________________\n\n");


	const char* filename = "output3d.";
	std::string name(filename);
 	char ext[50]; sprintf(ext, "txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");

	fprintf(outputFile,"\n_________________________ begin 'Prism4' ________________________\n");

	for ( size_t i = 0; i < vVertex.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_aaPos[vVertex[i]][0] , m_aaPos[vVertex[i]][1], m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: m_vCornerCoords.size() %lu \n", m_vCornerCoords.size());

	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1], m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vInterfaceID = %lu\n", m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vOriginalCornerID = %lu\n", m_vOriginalCornerID[i]);


	fprintf(outputFile,"\n_________________________ end 'Prism4' ______________________________\n\n");

	fclose(outputFile);
	//UG_THROW("first check: Prism2\n");
	}
	return m_vCornerCoords.size();

}

template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_Pyramid(GridObject* elem)
{
///////////////////////////////////////////////////////////////////
// 	1) push_back into 'vInsideCorners'
// 	2) push_back into 'vOutsideCorners' the 2 outsideCorners on a cut edge
//  3) write nearInterfaceCorner on m_vCornerCoords[4], since in this setting it
//		is ALWAYS the top of the Pyramid! --> see grid_object_3d.h for ordering!
///////////////////////////////////////////////////////////////////

	bool output = false;

	if ( output )
        UG_LOG("_________________________ begin 'CollectCorners_FlatTop_Pyramid()' ______________________________\n\n");

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

// buffer vectors for (cornerCoords, cornerIndex) - for application of 'ResortQuadrilateral'!
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::pair<MathVector<dim>, size_t > nearInterfaceData;

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *m_spMG, elem);

	// loop vertices
	size_t outsideInd;
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
 		Vertex* vrtRoot = vVertex[i];

	// get nearInterfaceData and directly write cornerCoords to data:
		if ( is_nearInterfaceVertex(vrtRoot, i) )
			nearInterfaceData = std::make_pair(m_aaPos[vrtRoot], i);
		else if ( ! is_FTVertex(vrtRoot, i) )
		{
	///////////////////////////////////////////////////////////////////
	// 	1) push_back into 'vInsideCorners'
			vInsideCorners.push_back(std::make_pair(m_aaPos[vrtRoot], i));
			m_vOriginalCornerID.push_back(i);

	///////////////////////////////////////////////////////////////////
	// 	2) collect according outsideCorner (lying on the commom edge)
			bool outsideCornerFound = false;
			MathVector<dim> intersectionPnt;
			for(size_t e = 0; e < vEdges.size(); ++e)
			{
				Edge* edge = vEdges[e];
				std::vector<Vertex*> vVertexEdge;
				CollectVertices(vVertexEdge, *m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("in 'NEWcall_Extrapolation': error in collecting vertices associated to an edge!....EXIT!...\n");
				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = get_vertex_index(vrt1, elem);
				size_t vrtInd2 = get_vertex_index(vrt2, elem);

			// lies vrtRoot on a cutted edge?
			// do nothing for edges with 1 insideCorner and 1 nearInterface corner:
				if ( is_nearInterfaceVertex(vrt1,vrtInd1) || is_nearInterfaceVertex(vrt2,vrtInd2) )
					continue;
  			// case2: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt1 && is_FTVertex(vrt2, vrtInd2) )
					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt1, vrt2);}
			// case3: vrt2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt2 && is_FTVertex(vrt1, vrtInd1) )
					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt2, vrt1);}
				else
					continue;

			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners3d_Prism3()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition(intersectionPnt)) << "\n");

			} // end edge-loop

			if ( !outsideCornerFound )
				UG_THROW("During edge-loop no outsideCorner found!\n");

	///////////////////////////////////////////////////////////////////
	// 	3) push_back into 'vOutsideCorners'
			vOutsideCorners.push_back(std::make_pair(intersectionPnt, 0));

		} // end inside-case
		else
		// recall ousideIndex for 'm_vOriginalCornerID'
			outsideInd = i;

	} // end vertex-loop

// some checks:
	if ( vInsideCorners.size() 	!= 2 ) UG_THROW("wrong size for vInsideCorners! Should be 2 and is " << vInsideCorners.size() << "\n");
	if ( vOutsideCorners.size() != 2 ) UG_THROW("wrong size for vOutsideCorners! Should be 2 and is " << vOutsideCorners.size() << "\n");

	for(size_t i = 0; i < vOutsideCorners.size(); ++i)
		vOutsideCorners[i].second = outsideInd;

// during ResortQuadrilateral: m_vCornerCoords.resize(4) and m_vOriginalCornerID.resize(4)!
	MathVector<dim> normalDir;
	VecSubtract(normalDir, nearInterfaceData.first, vOutsideCorners[0].first);
	ResortQuadrilateral(vInsideCorners, vOutsideCorners, normalDir);

	m_vCornerCoords.push_back(nearInterfaceData.first);
	m_vOriginalCornerID.push_back(nearInterfaceData.second);
	m_vInterfaceID.push_back(4);

// some checks:
	if ( m_vInterfaceID.size() != 3 ) UG_THROW("wrong size for m_vInterfaceID! Should be 3 and is " << m_vInterfaceID.size() << "\n");
	if ( m_vOriginalCornerID.size()!= 5 ) UG_THROW("wrong size for m_vOriginalCornerID! Should be 5 and is " << m_vOriginalCornerID.size() << "\n");
	if ( m_vCornerCoords.size()!= 5 ) UG_THROW("wrong size for m_vCornerCoords! Should be 5 and is " << m_vCornerCoords.size() << "\n");

/*
	for(size_t i = 0; i < 2; ++i)
		UG_LOG("vOutsideCorners[i] = " << vOutsideCorners[i] << "\n");
	for(size_t i = 0; i < 2; ++i)
		UG_LOG("vInsideCorners[i] = " << vInsideCorners[i] << "\n");
	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		UG_LOG("m_vOriginalCornerID[i] = " << m_vOriginalCornerID[i] << "\n");
	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		UG_LOG("m_vInterfaceID[i] = " << m_vInterfaceID[i] << "\n");
	for(size_t i = 0; i < 5; ++i)
		UG_LOG("m_vCornerCoords[i] = " << m_vCornerCoords[i] << "\n");
*/

	const char* filename = "output3d_Prism.";
	std::string name(filename);
	char ext[50]; sprintf(ext, "txt");
	name.append(ext);
	FILE* outputFile = fopen(name.c_str(), "a");


	fprintf(outputFile,"\n_________________________ begin 'Pyramid' ________________________\n");

	for ( size_t i = 0; i < vVertex.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_aaPos[vVertex[i]][0] , m_aaPos[vVertex[i]][1], m_aaPos[vVertex[i]][2]);

	fprintf(outputFile,"\n****Nachher: m_vCornerCoords.size() %lu \n", m_vCornerCoords.size());
/*
	for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
		fprintf(outputFile,"%e, %e, %e \n", m_vCornerCoords[i][0], m_vCornerCoords[i][1], m_vCornerCoords[i][2]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vInterfaceID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vInterfaceID = %lu\n", m_vInterfaceID[i]);

	fprintf(outputFile,"\n");

	for(size_t i = 0; i < m_vOriginalCornerID.size(); ++i)
		fprintf(outputFile,"Check in CollectIBCorner: m_vOriginalCornerID = %lu\n", m_vOriginalCornerID[i]);

*/
	fprintf(outputFile,"\n_________________________ end 'Pyramid' ______________________________\n\n");

	fclose(outputFile);
	//UG_THROW("first check: Prism3\n");

	return m_vCornerCoords.size();

}

template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_originalTet(GridObject* elem)
{
	if ( 0 ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_originalTet()' ______________________________\n\n");

	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

	for(size_t i = 0; i < vVertex.size(); ++i)
  		is_FTVertex(vVertex[i], i);

	// loop vertices
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		// get element
		Vertex* vrtRoot = vVertex[i];

		if ( ! is_FTVertex(vrtRoot, i) )
		{
			if ( is_nearInterfaceVertex(vrtRoot, i) )
			{	UG_THROW("NearInterface BUT !is_FTVertex() => error!\n"); }
			else
			{
				m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);
  			}

		}
		else if ( is_nearInterfaceVertex(vrtRoot, i) )
		{
 			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);
 		}
		else
			UG_THROW("in 'CollectCorners_FlatTop_originalTet()': outside fluid and NOT nearInterface case not possible!\n");
	}

	if ( 0 ) UG_LOG("_________________________ end 'CollectCorners_FlatTop_originalTet()' ______________________________\n\n");

	return m_vCornerCoords.size();

}

template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
get_cutMode(std::vector<Vertex*> vVertex)
{
    size_t numOutside = 0;
    size_t numNearInterface = 0;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
 		if ( is_FTVertex(vVertex[i], i) )
		{
			numOutside++;
			if ( is_nearInterfaceVertex(vVertex[i], i) )
				numNearInterface++;
		}
	}

	bool logMode = false;

// if all outside corners are nearInterface, no computation of intersectionPoints
//  => purely collecting of corners necessary:
	if ( numOutside == numNearInterface )
	{
		if ( logMode ) UG_LOG("numOutside = " << numOutside << ", numNearInterface = " << numNearInterface << "\n");
		if ( logMode ) UG_LOG("cutMode = 0\n");
		return 0; // original
	}
	else if ( numOutside == 1 )
	{
		if ( logMode ) UG_LOG("numOutside = " << numOutside << ", numNearInterface = " << numNearInterface << "\n");
		if ( logMode ) UG_LOG("cutMode = 1\n");
		return 1; // Prism3
	}
	else if ( numOutside == 2 && numNearInterface == 0 )
	{
		if ( logMode ) if ( logMode ) UG_LOG("numOutside = " << numOutside << ", numNearInterface = " << numNearInterface << "\n");
		if ( logMode ) UG_LOG("cutMode = 2\n");
		return 2; // Prism4
	}
	else if ( numOutside == 2 && numNearInterface == 1 )
	{
		if ( logMode ) UG_LOG("numOutside = " << numOutside << ", numNearInterface = " << numNearInterface << "\n");
		if ( logMode ) UG_LOG("cutMode = 3\n");
		return 3; // Pyramid
	}
	else if ( numOutside == 3 )
	{
		if ( logMode ) UG_LOG("numOutside = " << numOutside << ", numNearInterface = " << numNearInterface << "\n");
		if ( logMode ) UG_LOG("cutMode = 4\n");
		return 4; // Tetrahedron
	}
	return -1;

}

template <int TWorldDim>
int InterfaceHandlerLocalBase<TWorldDim>::
CollectCorners_FlatTop_3d(GridObject* elem)
{
	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

	size_t cutMode = get_cutMode(vVertex);

 	switch(cutMode)
	{
		case 0: return CollectCorners_FlatTop_originalTet(elem); 	// Tetrahedron with 1 corner nearInterface
		case 1: return CollectCorners_FlatTop_Prism3(elem); 		// Prism3 with 3 corners outside
		case 2: return CollectCorners_FlatTop_Prism4(elem); 	 	// Prism4 with 4 corners outside
		case 3: return CollectCorners_FlatTop_Pyramid(elem); 	 	// Pyramid with 2 corners outside and 1 corner nearInterface
		case 4: return CollectCorners_FlatTop_2d(elem); 			// Tetrahedron with 3 corners outside
		default: { UG_LOG("cutMode = " << cutMode << "\n");
			throw(UGError("Error in calling case dependent CollectIBCorners3d()!"));
		}
	}


}

template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
print_InterfaceDdata()
{
	UG_LOG(" ------------ print_InterfaceIDdata -----------\n");

		if ( m_roid == ROID_QUADRILATERAL )
			UG_LOG("ín add_jac_A_elem ---------ROID_QUADRILATERAL--------\n");
	if ( m_roid == ROID_TRIANGLE )
			UG_LOG("ín add_jac_A_elem ---------ROID_TRIANGLE--------\n");

		for ( size_t i = 0; i < m_vCornerCoords.size(); ++i )
			UG_LOG("m_vCornerCoords = " << m_vCornerCoords[i] << "\n");
		for ( size_t i = 0; i < m_vOriginalCornerID.size(); ++i )
			UG_LOG("Original: id = " << m_vOriginalCornerID[i] << "\n");
		UG_LOG("\n");
		for ( size_t i = 0; i < m_vInterfaceID.size(); ++i )
			UG_LOG("Interface: id = " << m_vInterfaceID[i] << "\n");
		UG_LOG("\n");

}

} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_BASE_TOOLS_H_ */
