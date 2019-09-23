/*
 * integrate_tools.h
 *
 *  Created on: 03.07.2015
 *      Author: suze
 */

#ifndef INTEGRATE_TOOLS_H_
#define INTEGRATE_TOOLS_H_

namespace ug{

template<int dim, typename TGridFunction>
number get_value(SmartPtr<TGridFunction> spGridFct, const int level, const size_t fct, const MathVector<dim> Corner)
{
 	std::vector<size_t> vTransInd(6);
	vTransInd[0] = 0;
	vTransInd[1] = 143;
	vTransInd[2] = 549;
	vTransInd[3] = 3029;
	vTransInd[4] = 8529;
	vTransInd[5] = 62201;

 	std::vector<size_t> vRotInd(6);
 	vRotInd[0] = 0;
 	vRotInd[1] = 12;
 	vRotInd[2] = 12;
 	vRotInd[3] = 12;
 	vRotInd[4] = 2153;
 	vRotInd[5] = 2153;

	MathVector<dim> RotVec;
	RotVec[0] = -Corner[1];
	RotVec[1] =  Corner[0];

	const size_t transInd = vTransInd[level];
	const size_t rotInd = vRotInd[level];

	const number transVel = DoFRef(*spGridFct, DoFIndex(transInd, fct));
	const number rotVel = DoFRef(*spGridFct, DoFIndex(rotInd, fct));

	//UG_LOG("transVel = " << transVel << ", rotVel = " << rotVel << "\n");

	const number val = transVel + rotVel*RotVec[fct];

	return val;
}

template<int dim>
bool set_nearInterface(MathVector<dim> vrtPos, int prtIndex)
{
	const number radius = 0.05;
	const MathVector<dim> center(0.2);
	const number threshold = 1e-9;
 	// ToDo threshold!

  // default value of threshold = 0.0 => do nothing in this case
//	if (threshold == 0.0)
//		return false;

  	const number dist = VecDistance(vrtPos, center);

	if (fabs(dist - radius) < threshold)
		return true;

	return false;
}


template<int dim>
number get_insideParticle(MathVector<dim> vrtPos, int prtIndex)
{
	const number radius = 0.05;
	const MathVector<dim> center(0.2);

	bool outsideFluid = false;
	if ( (VecDistance(vrtPos, center)-radius) < 1e-10)
	{
 		outsideFluid = true;
	}
	else
	{
 		outsideFluid = set_nearInterface<dim>(vrtPos, prtIndex);
	}

	return (VecDistance(vrtPos, center)-radius);
}

template<int dim>
bool is_insideParticle(MathVector<dim> vrtPos, int prtIndex)
{
	const number radius = 0.05;
	const MathVector<dim> center(0.2);

	bool outsideFluid = false;
	if ( (VecDistance(vrtPos, center)-radius) < 1e-10)
		outsideFluid = true;
	else
		outsideFluid = set_nearInterface<dim>(vrtPos, prtIndex);

	return outsideFluid;
}

template<int dim>
number get_LSvalue_byPosition(MathVector<dim> vrtPos, int PrtIndex)
{
	const number radius = 0.05;
	const MathVector<dim> center(0.2);

	MathVector<dim> localCoords;
	VecSubtract(localCoords, vrtPos, center);

	number dist = VecDot(localCoords, localCoords);
	dist = sqrt(dist);

	return radius - dist;
}

template<int dim>
bool get_LineCircle_Intersection(MathVector<dim>& Intersect, Vertex* vrtOutsideCirc,
		Vertex* vrtInsideCirc, int PrtIndex,
		typename domain_traits<dim>::position_accessor_type& aaPos)
{
	const number radius = 0.05;
	const MathVector<dim> center(0.2);

	//if ( dim == 3 ) UG_THROW("in 'get_LineCircle_Intersection()': not implemented for 3d case!\n");

	///////////////////////////////////////////////////////////////////////////////////////
	//
	// 'vrtPosOut':= the starting point of the ray,
	// 'vrtPosIn' := the end point of the ray,
	// 'center'   := the center of sphere you're testing against
	// 'radius'   := the radius of that sphere
	// 'lineDir'  := direction vector of ray from start to end
	// 'rayDir'   := direction vector of ray from center to start
	//
	// 	Ansatz:
	//  (1) Intersect = vrtPosOut + alpha*lineDir
	//	(2) Intersect - center = radius
	//
	// 		=> (1)  Intersect[0] = vrtPosOut[0] + alpha*lineDir[0]
	// 			    Intersect[1] = vrtPosOut[1] + alpha*lineDir[1]
	// 		=> (2)  (Intersect[0] - center[0])^2 + (Intersect[1] - center[1])^2 = radius^2
	//
	// 	Plug (1) into (2) => ... => quadratic equation for alpha:
	//
	//			alpha^2 * <lineDir,lineDir> + alpha * 2*<lineDir, rayDir> + ( <rayDir, rayDir>-radius^2 ) = 0
	//
	// 			-> a =  <lineDir,lineDir>, b =  <lineDir,linrayDireDir>, c =  <rayDir,rayDir> - radius^2
	//
	//
	// 	=> from (1): Intersect = vrtPosOut + alpha*(vrtPosIn - vrtPosOut)
	///////////////////////////////////////////////////////////////////////////////////////

	number alpha;

 	const MathVector<dim>& vrtPosOut = aaPos[vrtOutsideCirc];
	const MathVector<dim>& vrtPosIn = aaPos[vrtInsideCirc];

	MathVector<dim> lineDir;
	MathVector<dim> rayDir;

// lineDir = vrtPosIn - vrtPosOut;
	VecSubtract(lineDir, vrtPosIn, vrtPosOut);
// rayDir = vrtPosOut - center;
	VecSubtract(rayDir, vrtPosOut, center);

	const number a = VecDot(lineDir, lineDir);
	const number b = 2.0 * VecDot(lineDir, rayDir);
	const number c = VecDot(vrtPosOut, vrtPosOut) + VecDot(center, center)
			- 2 * VecDot(vrtPosOut, center) - radius * radius;

	const number discriminant = b * b - 4 * a * c;

// check that 'vrtPosOut' and 'vrtPosIn' really lie on different sides of the circle:
	if (discriminant < -1e-8)
 		UG_THROW("Value of discriminant = " << discriminant << "\n");


// discriminant = 0!
	const number alpha1 = (-b - sqrt(discriminant)) / (2.0 * a);
	const number alpha2 = (-b + sqrt(discriminant)) / (2.0 * a);

	if (alpha1 <= alpha2)
		alpha = alpha1;
	else
		alpha = alpha2;

	if (alpha < 0 || (alpha - 1.0) > 1e-8)
		UG_THROW(
				"Error in 'get_LineCircle_Intersection()': alpha not valid; should lie between 0 and 1: " << alpha << "\n");

	for (size_t d = 0; d < dim; ++d)
		Intersect[d] = vrtPosOut[d] + alpha * lineDir[d];

	return true;
}
/*
bool isIncluded2(const size_t newID, std::vector<size_t> vInterfaceID)
{
	for ( size_t i = 0; i < vInterfaceID.size(); ++i)
		if ( vInterfaceID[i] == newID )
			return true;

	return false;
}
*/
/*
void reset_sh_if_on_interface(size_t& sh,
							  std::vector<size_t> vInterfaceID,
							  std::vector<size_t>& vOriginalCornerID)
{
	for ( size_t i = 0; i < vInterfaceID.size(); ++i)
		if ( vInterfaceID[i] == sh )
			sh = vOriginalCornerID[i];

}
*/
template <int dim>
bool isIncluded(std::vector<MathVector<dim> > vCheckList, MathVector<dim> checkPoint)
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

template <int dim>
ReferenceObjectID set_roid_2d(const int numCo)
{
	ReferenceObjectID roid;

 	switch(numCo)
	{
		case 0: roid = ROID_UNKNOWN; 			break;
		case 3: roid = ROID_TRIANGLE; 			break;
		case 4: roid = ROID_QUADRILATERAL;		break;

		default: throw(UGError("during switch in 'set_roid_2d()': "
				"m_numCo invalid => default case chosen"));
	}
	return roid;

}

template <int dim>
ReferenceObjectID set_roid_3d(const int numCo)
{
	ReferenceObjectID roid;

	switch(numCo)
	{
		case 0: roid = ROID_UNKNOWN; 			break;
		case 4: roid = ROID_TETRAHEDRON;  		break;
		case 5: roid = ROID_PYRAMID; 			break;
		case 6: roid = ROID_PRISM; 				break;

		default: throw(UGError("during switch in 'set_roid_3d()': "
				"m_numCo invalid => default case chosen"));
	}
	return roid;

}

// parameter 'normal' is the direction from the triangle into the inner of the prism
template <int dim>
bool isCCW(std::vector<MathVector<dim> > vCornerCoords, MathVector<dim> normal)
{

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

template <int dim>
ReferenceObjectID reset_sh_if_on_interface(size_t& sh,
							  std::vector<size_t> vInterfaceID,
							  std::vector<size_t>& vOriginalCornerID)
{

	for ( size_t i = 0; i < vInterfaceID.size(); ++i)
		if ( vInterfaceID[i] == sh )
			sh = vOriginalCornerID[i];

	ReferenceObjectID roid;
	int numCo = 4;
	switch(numCo)
	{
		case 0: roid = ROID_UNKNOWN; 			break;
		case 4: roid = ROID_TETRAHEDRON;  		break;
		case 5: roid = ROID_PYRAMID; 			break;
		case 6: roid = ROID_PRISM; 				break;

		default: throw(UGError("during switch in 'set_roid_3d()': "
				"m_numCo invalid => default case chosen"));
	}

	return roid;

}

template <int dim>
void ResortQuadrilateral(std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners,
					std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners,
					MathVector<dim> normalDir,
					std::vector<MathVector<dim> >& m_vCornerCoords,
					std::vector<size_t>& m_vOriginalCornerID,
					std::vector<size_t>& m_vInterfaceID)
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

template <int dim>
bool CollectCorners_FlatTop_2d(GridObject* elem,
		SmartPtr<MultiGrid> m_spMG,
		typename domain_traits<dim>::position_accessor_type& aaPos,
		std::vector<MathVector<dim> >& m_vCornerCoords,
		std::vector<size_t>& m_vOriginalCornerID,
		std::vector<size_t>& m_vInterfaceID,
		ReferenceObjectID& coarseROID
)
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


	// get prtIndex:
	bool isFTVertex = false;
    size_t numCoOut = 0;
 
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
		isFTVertex = is_insideParticle<dim>(aaPos[vVertex[i]], 0);
		if ( isFTVertex )
            numCoOut++;
	}


	if ( !isFTVertex || numCoOut == vVertex.size() )
        return false;
    

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *m_spMG, elem);

	// loop vertices
	//////////////////////////////////////////////
	// REMARK:
	// order is the same as in 'vCornerCoords', therefore we can be sure, that the
	// order of the new 'vCornerIBCoords' will be consistent with the grid standard
	///////////////////////// /////////////////////

	bool bNearInterface = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{

		// get element
		Vertex* vrtRoot = vVertex[i];
		//////////////////////////////////////////////
		// case 1:
		// vertex insideFluid
		// 		=> Position und index puschen
		if ( ! is_insideParticle<dim>(aaPos[vrtRoot], 0) )
		{
 
			if ( set_nearInterface<dim>(aaPos[vrtRoot], 0) )
			{
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");
			}
			else
			{
				m_vCornerCoords.push_back(aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);

				vInsideCorners.push_back(std::make_pair(aaPos[vrtRoot], i));
 			}

		}
		//////////////////////////////////////////////
  		// case 2:
		// vertex = FT + ON interface
		// 		=> KEINE Berechnung von 'intersectionPoint' notwendig! -> pushen und alten index pushen

		// REMARK: is_nearInterfaceVerx = false per default, if m_vThresholdOnLevel = 0.0
		else if ( set_nearInterface<dim>(aaPos[vrtRoot], 0) )
		{
 
 			bNearInterface = true;
			m_vCornerCoords.push_back(aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(aaPos[vrtRoot], i));
			vNearIntCorners.push_back(std::make_pair(aaPos[vrtRoot], i));

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

				MathVector<dim> intersectionPnt;

	 		///////////////////////////////////////////////////////////////////
 			// lies vrtRoot on a cutted edge?
		 	///////////////////////////////////////////////////////////////////
			// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
				if ( set_nearInterface<dim>(aaPos[vrt2], 0) || set_nearInterface<dim>(aaPos[vrt1], 0) )
                {	bNearInterface = true; continue; }
			 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt1 && !is_insideParticle<dim>(aaPos[vrt2], 0) ){
 					get_LineCircle_Intersection<dim>(intersectionPnt, vrt2, vrt1, 0, aaPos);
 				}
			// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt2 && !is_insideParticle<dim>(aaPos[vrt1], 0) )
				{ get_LineCircle_Intersection<dim>(intersectionPnt, vrt1, vrt2, 0, aaPos);}
				else
				{continue;}

			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) << "\n");


	 		///////////////////////////////////////////////////////////////////
	 		// only push_back, if not included yet!
			// 	-> can be ONLY the case, if the intersectionPoint is a node
	 			if ( ! isIncluded<dim>(m_vCornerCoords, intersectionPnt) )
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
	{
 		ResortQuadrilateral<dim>(vInsideCorners, vOutsideCorners, normalDir, m_vCornerCoords, m_vOriginalCornerID, m_vInterfaceID);
	}
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

	const int numCo = m_vCornerCoords.size();

	if ( dim == 2 )
		coarseROID = set_roid_2d<dim>(numCo);

    if ( numCo == 0 )
        UG_THROW("hmmm: coarseROID = " << coarseROID << "\n");
    
	return true;

}


template <int dim>
bool CollectCorners_FlatTop_Prism3(GridObject* elem,
		SmartPtr<MultiGrid> m_spMG,
		typename domain_traits<dim>::position_accessor_type& aaPos,
		std::vector<MathVector<dim> >& m_vCornerCoords,
		std::vector<size_t>& m_vOriginalCornerID,
		std::vector<size_t>& m_vInterfaceID,
		ReferenceObjectID& coarseROID)
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

		if ( ! is_insideParticle<dim>(aaPos[vrtRoot], 0) )
		{
	///////////////////////////////////////////////////////////////////
	// 	1) push_back into 'vInsideCorners'
			vInsideCorners.push_back(aaPos[vrtRoot]);
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

			// lies vrtRoot on a cutted edge?
			// case1: vrt1 = outsideParticle && vrt2 = insideParticle:
				if ( vrtRoot == vrt1 && is_insideParticle<dim>(aaPos[vrt2], 0) )
  					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt1, vrt2);}
			// case2: vrt2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt2 && is_insideParticle<dim>(aaPos[vrt1], 0) )
					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt2, vrt1);}
			// NO nearInterface handling necessary, since this is Pyramid or Tetrahedron case!
				else if ( set_nearInterface<dim>(aaPos[vrt2], 0) || set_nearInterface<dim>(aaPos[vrt1], 0) )
					{UG_THROW("in 'CollectCorners_Prism3()': no nearInterface case possible!\n");}
				else
					continue;

			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners3d_Prism3()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) << "\n");

			} // end edge-loop

			if ( !outsideCornerFound )
				UG_THROW("During edge-loop no outsideCorner found!\n");

	///////////////////////////////////////////////////////////////////
	// 	3) push_back into 'vOutsideCorners'
			if ( ! isIncluded<dim>(vOutsideCorners, intersectionPnt) )
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

// check 'vOutsideCorners'-array for CCW:
// get normal to triangle, built by 'vOutsideCorners'-array:
	MathVector<dim> normal = vInsideCorners[0];
 	VecSubtract(normal, vOutsideCorners[0], normal);

	if ( !isCCW(vOutsideCorners, normal) )
	{
 		m_vCornerCoords[4] = vOutsideCorners[2];
		m_vCornerCoords[5] = vOutsideCorners[1];
	}

// check 'vInsideCorners'-array for CCW:
	if ( !isCCW(vInsideCorners, normal) )
	{
 		m_vCornerCoords[1] = vInsideCorners[2];
		m_vCornerCoords[2] = vInsideCorners[1];
 	}


	return true;

}

template <int dim>
bool CollectCorners_FlatTop_Prism4(GridObject* elem,
				SmartPtr<MultiGrid> m_spMG,
				typename domain_traits<dim>::position_accessor_type& aaPos,
				std::vector<MathVector<dim> >& m_vCornerCoords,
				std::vector<size_t>& m_vOriginalCornerID,
				std::vector<size_t>& m_vInterfaceID,
				ReferenceObjectID& coarseROID)
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
		if ( is_insideParticle<dim>(aaPos[vVertex[i]], 0) )
			outsideInd.push_back(i);

	if ( outsideInd.size() != 2 ) UG_THROW("Wrong number of outsideInd!\n");
	for(size_t i = 0; i < outsideInd.size(); ++i)
		if ( output ) UG_LOG("outsideInd = " << outsideInd[i] << "\n");

	// loop vertices
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		// get element
		Vertex* vrtRoot = vVertex[i];

		if ( ! is_insideParticle<dim>(aaPos[vrtRoot], 0) )
		{
	///////////////////////////////////////////////////////////////////
	// 	1) push_back 'insideCorner' into 'm_vCornerCoords'
			m_vCornerCoords.push_back(aaPos[vrtRoot]);
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

			// lies vrtRoot on a cutted edge?
			// case1: vrt1 = outsideParticle && vrt2 = insideParticle:
				if ( vrtRoot == vrt1 && is_insideParticle<dim>(aaPos[vrt2], 0) )
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
				else if ( vrtRoot == vrt2 && is_insideParticle<dim>(aaPos[vrt1], 0) )
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
				if ( fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners3d_Prism2()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = "
							<< fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) << "\n");
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
	}


	return true;

}

template <int dim>
bool CollectCorners_FlatTop_Pyramid(GridObject* elem,
		SmartPtr<MultiGrid> m_spMG,
		typename domain_traits<dim>::position_accessor_type& aaPos,
		std::vector<MathVector<dim> >& m_vCornerCoords,
		std::vector<size_t>& m_vOriginalCornerID,
		std::vector<size_t>& m_vInterfaceID,
		ReferenceObjectID& coarseROID)
{
///////////////////////////////////////////////////////////////////
// 	1) push_back into 'vInsideCorners'
// 	2) push_back into 'vOutsideCorners' the 2 outsideCorners on a cut edge
//  3) write nearInterfaceCorner on m_vCornerCoords[4], since in this setting it
//		is ALWAYS the top of the Pyramid! --> see grid_object_3d.h for ordering!
///////////////////////////////////////////////////////////////////

	bool output = false;

	if ( output ) UG_LOG("_________________________ begin 'CollectCorners_FlatTop_Pyramid()' ______________________________\n\n");

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
	// get element
		Vertex* vrtRoot = vVertex[i];

	// get nearInterfaceData and directly write cornerCoords to data:
		if ( set_nearInterface<dim>(aaPos[vrtRoot], 0) )
			nearInterfaceData = std::make_pair(aaPos[vrtRoot], i);
		else if ( ! is_insideParticle<dim>(aaPos[vrtRoot], 0) )
		{
	///////////////////////////////////////////////////////////////////
	// 	1) push_back into 'vInsideCorners'
			vInsideCorners.push_back(std::make_pair(aaPos[vrtRoot], i));
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

			// lies vrtRoot on a cutted edge?
			// do nothing for edges with 1 insideCorner and 1 nearInterface corner:
				if ( set_nearInterface<dim>(aaPos[vrt1], 0) || set_nearInterface<dim>(aaPos[vrt2], 0) )
					continue;
  			// case2: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt1 && is_insideParticle<dim>(aaPos[vrt2], 0) )
					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt1, vrt2);}
			// case3: vrt2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt2 && is_insideParticle<dim>(aaPos[vrt1], 0) )
					{outsideCornerFound = get_intersection_point(intersectionPnt, vrt2, vrt1);}
				else
					continue;

			// check for correct inersectionPnt
				if ( fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners3d_Prism3()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(get_LSvalue_byPosition<dim>(intersectionPnt, 0)) << "\n");

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

	return true;

}

template <int dim>
bool CollectCorners_FlatTop_originalTet(GridObject* elem,
		SmartPtr<MultiGrid> m_spMG,
		typename domain_traits<dim>::position_accessor_type& aaPos,
		std::vector<MathVector<dim> >& m_vCornerCoords,
		std::vector<size_t>& m_vOriginalCornerID,
		std::vector<size_t>& m_vInterfaceID,
		ReferenceObjectID& coarseROID)
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

	// loop vertices
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		// get element
		Vertex* vrtRoot = vVertex[i];

		if ( ! is_insideParticle<dim>(aaPos[vrtRoot], 0) )
		{
			if ( set_nearInterface<dim>(aaPos[vrtRoot], 0) )
			{	UG_THROW("NearInterface BUT !is_onInterfaceVertex() => error!\n"); }
			else
			{
				m_vCornerCoords.push_back(aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);
  			}

		}
		else if ( set_nearInterface<dim>(aaPos[vrtRoot], 0) )
		{
 			m_vCornerCoords.push_back(aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);
 		}
		else
			UG_THROW("in 'CollectCorners_FlatTop_originalTet()': outside fluid and NOT nearInterface case not possible!\n");
	}

	if ( 0 ) UG_LOG("_________________________ end 'CollectCorners_FlatTop_originalTet()' ______________________________\n\n");

	return true;

}

template <int dim>
int get_cutMode(std::vector<Vertex*> vVertex, typename domain_traits<dim>::position_accessor_type& aaPos)
{
    size_t numOutside = 0;
    size_t numNearInterface = 0;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		if ( is_insideParticle<dim>(aaPos[vVertex[i]], 0) )
		{
			numOutside++;
			if ( set_nearInterface<dim>(aaPos[vVertex[i]], 0) )
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

template <int dim>
bool CollectCorners_FlatTop_3d(GridObject* elem,
		SmartPtr<MultiGrid> m_spMG,
		typename domain_traits<dim>::position_accessor_type& aaPos,
		std::vector<MathVector<dim> >& m_vCornerCoords,
		std::vector<size_t>& m_vOriginalCornerID,
		std::vector<size_t>& m_vInterfaceID,
		ReferenceObjectID& coarseROID)
{
	bool output = false;
	if ( output ) UG_LOG("_________________________ begin 'get_LSvalue_byPosition()' ______________________________\n\n");

	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

// remember boolian for check, weather there existe at least one vertex, which is FT!
	bool isFTVertex = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		isFTVertex = is_insideParticle<dim>(aaPos[vVertex[i]], 0);
		if ( isFTVertex )
			break;
	}
	if ( !isFTVertex )
		return false;


	size_t cutMode = get_cutMode<dim>(vVertex, aaPos);

 	switch(cutMode)
	{
		case 0: return CollectCorners_FlatTop_originalTet<dim>(elem,m_spMG,aaPos,m_vCornerCoords,m_vOriginalCornerID,m_vInterfaceID,coarseROID);
		case 1: return CollectCorners_FlatTop_Prism3<dim>(elem,m_spMG,aaPos,m_vCornerCoords,m_vOriginalCornerID,m_vInterfaceID,coarseROID);
		case 2: return CollectCorners_FlatTop_Prism4<dim>(elem,m_spMG,aaPos,m_vCornerCoords,m_vOriginalCornerID,m_vInterfaceID,coarseROID);
		case 3: return CollectCorners_FlatTop_Pyramid<dim>(elem,m_spMG,aaPos,m_vCornerCoords,m_vOriginalCornerID,m_vInterfaceID,coarseROID);
		case 4: return CollectCorners_FlatTop_2d<dim>(elem,m_spMG,aaPos,m_vCornerCoords,m_vOriginalCornerID,m_vInterfaceID,coarseROID);
		default: { UG_LOG("cutMode = " << cutMode << "\n");
			throw(UGError("Error in calling case dependent CollectIBCorners3d()!"));
		}
	}


}


template <int dim>
bool CollectCorners_StdFV(GridObject* elem,
		SmartPtr<MultiGrid> m_spMG,
		typename domain_traits<dim>::position_accessor_type& aaPos,
		std::vector<MathVector<dim> >& m_vCornerCoords,
		std::vector<size_t>& m_vOriginalCornerID,
		std::vector<size_t>& m_vInterfaceID,
		ReferenceObjectID& coarseROID)
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

	// get prtIndex:
	bool isOutsideVertex = false;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
 		if ( is_insideParticle<dim>(vVertex[i], 0) )
 		{
 			isOutsideVertex = true;
 			break;
 		}
	}

	if ( !isOutsideVertex )
		return false;

// loop vertices
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// get element
		Vertex* vrtRoot = vVertex[i];

		//////////////////////////////////////////////
		// case 1:
		// vertex insideFluid => Position und index puschen
		if ( !is_insideParticle<dim>(vVertex[i]) )
		{
			if ( set_nearInterface<dim>(aaPos[vrtRoot], 0) )
			{
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");
			}
			else
			{
				m_vCornerCoords.push_back(aaPos[vrtRoot]);
				m_vOriginalCornerID.push_back(i);

				vInsideCorners.push_back(std::make_pair(aaPos[vrtRoot], i));
			}

		}
		else
		{
 			m_vCornerCoords.push_back(aaPos[vrtRoot]);
			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(aaPos[vrtRoot], i));

		}

	}

	// skip whole element, since only FT points are included
 	if ( dim == 2 && m_vCornerCoords.size() != 3 )
 		UG_THROW("in 'CollectCorners_StdFV()': m_vCornerCoords.size() "
					"= " << m_vCornerCoords.size() << "not possible!\n");
 	if ( dim == 3 && m_vCornerCoords.size() != 4 )
 		UG_THROW("in 'CollectCorners_StdFV()': m_vCornerCoords.size() "
					"= " << m_vCornerCoords.size() << "not possible!\n");


	return true;
}

} // end namespace ug


#endif /* INTEGRATE_TOOLS_H_ */
