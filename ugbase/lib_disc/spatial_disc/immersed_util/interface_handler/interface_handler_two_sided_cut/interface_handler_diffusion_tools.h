/*
 * interface_handler_local_tools.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_DIFFUSION_TOOLS_H_
#define INTERFACE_HANDLER_LOCAL_DIFFUSION_TOOLS_H_

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	FlatTopBase methods - ROID/CollectCorners...
////////////////////////////////////////////////////////////////////////////////

// used for local_to_global mapping

template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
lies_onInterface_tri(const size_t newID)
{
	for ( size_t i = 0; i < m_vInterfaceID_tri.size(); ++i)
		if ( m_vInterfaceID_tri[i] == newID )
			return true;

	return false;
}

template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
lies_onInterface_quad(const size_t newID)
{
	for ( size_t i = 0; i < m_vInterfaceID_quad.size(); ++i)
		if ( m_vInterfaceID_quad[i] == newID )
			return true;

	return false;
}

template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
lies_onInterface_size(const size_t newID, size_t size)
{
	if ( size == 3 )
		return lies_onInterface_tri(newID);
	else if ( size == 4 )
		return lies_onInterface_quad(newID);
	else UG_THROW("in 'lies_onInterface_sized()': invalid size (shold be 3 or 4!):" << size << "\n");

}


// used during 'diffusion_interface/diffusion_interface.h:initialize_interface_Nitsche()':
template <int TWorldDim>
size_t InterfaceHandlerLocalDiffusion<TWorldDim>::
get_or_insert_indexPair_Nitsche(const size_t global_index)
{

 	std::pair<typename std::map<size_t,size_t>::iterator,bool> ret;
	ret = m_MapInserted_Nitsche.insert ( std::pair<size_t,size_t>(global_index, m_MapInserted_Nitsche.size()) );

	if (ret.second==false) {
		//UG_LOG("element already existed with a value of " << ret.first->second << "\n");
	}
	else{
		m_verticesGlobalIndex.push_back(global_index);
        UG_LOG("global_index: " << global_index << "\n");

        UG_LOG("ret.first->first: " << ret.first->first << "\n");
        UG_LOG("ret.first->second: " << ret.first->second << "\n");
 	}

	return ret.first->second;

}

template <int TWorldDim>
size_t InterfaceHandlerLocalDiffusion<TWorldDim>::
get_index_for_global_index_Nitsche(const size_t global_index)
{
	for ( size_t i = 0; i < m_verticesGlobalIndex.size(); ++i)
		if ( m_verticesGlobalIndex[i] == global_index )
			return i;

	UG_THROW("in 'global_index_of_interface_DoF_Nitsche': no index in m_verticesGlobalIndex-vector found for global_index = " << global_index << "\n");
}

template <int TWorldDim>
size_t InterfaceHandlerLocalDiffusion<TWorldDim>::
get_or_insert_vertex(Vertex* vrt)
{
    MathVector<dim> vrtPos = this->m_aaPos[vrt];
        
    // near interface vertices are remapped onto grid points => NO new DoF! => No addint into list!
    if (is_nearInterface(vrt) )
        return -1;
        
    get_or_insert_vertex(vrtPos);
        
}
    
template <int TWorldDim>
bool InterfaceHandlerLocalDiffusion<TWorldDim>::
find_vrtPos(const MathVector<dim>& vrtPos)
{
    typename std::map<MathVector<dim>,size_t>::iterator it;
    it = m_MapInserted.find(vrtPos);
    
    // if ret.second == false, then vrtPos allready exists in n_MapInserted
    if (it != m_MapInserted.end())
        return true;
        
    return false;
        
}

template <int TWorldDim>
size_t InterfaceHandlerLocalDiffusion<TWorldDim>::
get_or_insert_vertex(const MathVector<dim>& vrtPos)
{
 	std::pair<typename std::map<MathVector<dim>,size_t>::iterator,bool> ret;
	ret = m_MapInserted.insert ( std::pair<MathVector<dim>,size_t>(vrtPos,m_MapInserted.size()) );

	if (ret.second==false) {
		//UG_LOG("element already existed with a value of " << ret.first->second << "\n");
	}
	else{
		m_verticesPos.push_back(vrtPos);
 	}

	return ret.first->second;
}

template <int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
Resort_RealCornerID()
{
	m_vRealCornerID_quad.clear();

	for ( size_t i = 0; i < this->m_vOriginalCornerID.size(); ++i )
	{
		size_t origID = this->m_vOriginalCornerID[i];
		if ( !this->lies_onInterface(i) )
			m_vRealCornerID_quad.push_back(origID);
        else if ( !find_vrtPos(this->m_vCornerCoords[i]) )
            m_vRealCornerID_quad.push_back(origID);
		else
			m_vRealCornerID_quad.push_back(get_or_insert_vertex(this->m_vCornerCoords[i]));
	}
}

template <int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
Collect_Data_Nitsche(GridObject* elem)
{
	this->m_vInterfaceID.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

 	this->m_vCornerCoords.clear();
	for(size_t i = 0; i < vVertex.size(); ++i)
		this->m_vCornerCoords.push_back(this->m_aaPos[vVertex[i]]);

	int numFTVertex = 0;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{

		// remember boolian for check, weather there existe at least one vertex, which is FT!
 		if ( this->is_onInterfaceVertex(vVertex[i], i) )
			numFTVertex++;
	}

	m_vAlpha.clear();
	m_vIntersectionPnts.clear();

	MathVector<dim> insidePnt_check;

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	for(size_t e = 0; e < vEdges.size(); ++e)
	{
		Edge* edge = vEdges[e];
		std::vector<Vertex*> vVertexEdge;
		CollectVertices(vVertexEdge, *this->m_spMG, edge);
		if ( vVertexEdge.size() != 2 )
			UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");

		Vertex* vrt1 = vVertexEdge[0];
		Vertex* vrt2 = vVertexEdge[1];
		size_t vrtInd1 = this->get_vertex_index(vrt1, elem);
		size_t vrtInd2 = this->get_vertex_index(vrt2, elem);

		UG_LOG("m_aaPos[vrt1] = " <<  this->m_aaPos[vrt1][0] << ", m_aaPos[vrt1][1] = " <<  this->m_aaPos[vrt1][1] << "\n");
		UG_LOG("m_aaPos[vrt2] = " <<  this->m_aaPos[vrt2][0] << ", m_aaPos[vrt2][1] = " <<  this->m_aaPos[vrt2][1] << "\n");

		std::vector<number> alphaOut;
		MathVector<dim> intersectionPnt;
        bool bNearInterface = false;


	///////////////////////////////////////////////////////////////////
	// lies vrtRoot on a cutted edge?
	///////////////////////////////////////////////////////////////////
	// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
		if ( this->is_nearInterfaceVertex(vrt2, vrtInd2) || this->is_nearInterfaceVertex(vrt1, vrtInd1) )
		{ bNearInterface = true; continue; }
	 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
		else if ( this->is_onInterfaceVertex(vrt1, vrtInd1) && !this->is_onInterfaceVertex(vrt2, vrtInd2) )
		{
			get_intersection_point(intersectionPnt, vrt2, vrt1, alphaOut);
			UG_LOG("* alphaOut[0] = " << alphaOut[0] << ", alphaOut[1] = " << alphaOut[1] << "\n");
			UG_LOG("* intersectionPnt[0] = " << intersectionPnt[0] << ", intersectionPnt[1] = " << intersectionPnt[1] << "\n");
			m_vAlpha.push_back(alphaOut);
			m_vIntersectionPnts.push_back(intersectionPnt);
			m_vShapeValues[vrtInd2][vrtInd1] = alphaOut[0];
			m_vShapeValues[vrtInd1][vrtInd2] = alphaOut[1];
			if ( numFTVertex > 1 )
			{
				m_insidePnt[0] = this->m_aaPos[vrt2][0];
				m_insidePnt[1] = this->m_aaPos[vrt2][1];
			}
			else
			{
				m_insidePnt[0] = this->m_aaPos[vrt1][0];
				m_insidePnt[1] = this->m_aaPos[vrt1][1];
			}
 			this->m_vInterfaceID.push_back(vrtInd2);
		}
	// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
		else if ( this->is_onInterfaceVertex(vrt2, vrtInd2) && !this->is_onInterfaceVertex(vrt1, vrtInd1) )
		{
			get_intersection_point(intersectionPnt, vrt1, vrt2, alphaOut);
			UG_LOG("# alphaOut[0] = " << alphaOut[0] << ", alphaOut[1] = " << alphaOut[1] << "\n");
			UG_LOG("# intersectionPnt[0] = " << intersectionPnt[0] << ", intersectionPnt[1] = " << intersectionPnt[1] << "\n");
			m_vAlpha.push_back(alphaOut);
			m_vIntersectionPnts.push_back(intersectionPnt);
			m_vShapeValues[vrtInd1][vrtInd2] = alphaOut[0];
			m_vShapeValues[vrtInd2][vrtInd1] = alphaOut[1];
			if ( numFTVertex > 1 )
			{
				insidePnt_check[0] = this->m_aaPos[vrt1][0];
				insidePnt_check[1] = this->m_aaPos[vrt1][1];
			}
			else
			{
				insidePnt_check[0] = this->m_aaPos[vrt2][0];
				insidePnt_check[1] = this->m_aaPos[vrt2][1];
			}
 			this->m_vInterfaceID.push_back(vrtInd1);
		}
		else
		{
			alphaOut.clear(); alphaOut.push_back(0.0); alphaOut.push_back(0.0);
 			m_vAlpha.push_back(alphaOut);
			m_vShapeValues[vrtInd1][vrtInd2] = 0.0;
			m_vShapeValues[vrtInd2][vrtInd1] = 0.0;
 			continue;
		}

	// check for correct inersectionPnt
		if ( fabs(this->get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
			UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
					" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(this->get_LSvalue_byPosition(intersectionPnt)) << "\n");

	} // end edge-loop

/*
// and copy data to m_vInterfaceID_tri:
	m_vInterfaceID_tri.clear();
	for ( size_t i = 0; i < this->m_vInterfaceID.size(); ++i)
		m_vInterfaceID_1.push_back(m_vInterfaceID[i]);
*/
// check
	if ( insidePnt_check[0] != m_insidePnt[0])
	{
		UG_LOG("insidePnt_check[0] = " << insidePnt_check[0] << "\n");
 		UG_LOG("m_insidePnt[0] = " << m_insidePnt[0] << "\n");
 	}
	if ( insidePnt_check[1] != m_insidePnt[1])
	{
 		UG_LOG("insidePnt_check[1] = " << insidePnt_check[1] << "\n");
 		UG_LOG("m_insidePnt[1] = " << m_insidePnt[1] << "\n");
	}

// set values to zero
	m_vShapeValues[0][0] = m_vShapeValues[1][1] = m_vShapeValues[2][2] = 0.0;

// compute area
	if ( dim != 2 )
		UG_THROW(" for dim = 3 the implementatino needs to be adapted here!\n");

	MathVector<dim> _x10, _x20;
	MathVector<dim+1> x10, x20, n;

	VecSubtract(_x10, m_vIntersectionPnts[0], m_insidePnt);
	VecSubtract(_x20, m_vIntersectionPnts[1], m_insidePnt);
	for(size_t i = 0; i < dim; ++i)
	{
		x10[i] = _x10[i];
		x20[i] = _x20[i];
	}
	x10[2] = x20[2] = 0.0;

	VecCross(n, x10, x20);

	m_Area = 0.5 * VecTwoNorm(n);

	VecSubtract(_x10, this->m_aaPos[vVertex[1]], this->m_aaPos[vVertex[0]]);
	VecSubtract(_x20, this->m_aaPos[vVertex[2]], this->m_aaPos[vVertex[0]]);
	for(size_t i = 0; i < dim; ++i)
	{
		x10[i] = _x10[i];
		x20[i] = _x20[i];
	}
	x10[2] = x20[2] = 0.0;

	VecCross(n, x10, x20);

	m_AreaOrig = 0.5 * VecTwoNorm(n);

// if numFTVertex = 1: m_Area needs to be the area of the remaining quadrilateral!
	if ( numFTVertex == 1 )
		m_Area = m_AreaOrig - m_Area;

	m_AreaScale = m_Area/m_AreaOrig;

// compute normal and Gamma
	MathVector<dim> normal;
 	VecSubtract(normal, m_vIntersectionPnts[1], m_vIntersectionPnts[0]);
 	m_NormalToFace[0] = -normal[1];
 	m_NormalToFace[1] = normal[0];

 	m_Gamma = sqrt(VecDot(m_NormalToFace,m_NormalToFace));

// compute integrals along Gamma
 	m_vIntegralGamma.clear(); m_vIntegralGamma.resize(3);
	for(size_t i = 0; i < 3; ++i)
	{
		number valueSum = 0.0;
		for(size_t j = 0; j < 3; ++j)
			valueSum += m_vShapeValues[i][j];
		m_vIntegralGamma[i] = 0.5 * valueSum; //* m_Gamma;
	}

    UG_LOG("m_insidePnt = " << m_insidePnt[0] << " , " << m_insidePnt[1] << "\n");
	UG_LOG("m_NormalToFace = " << m_NormalToFace[0] << " , " << m_NormalToFace[1] << "\n");
	UG_LOG("m_Gamma = " << m_Gamma << "\n");
	UG_LOG("m_Area = " << m_Area << "\n");
    UG_LOG("m_AreaOrig = " << m_AreaOrig << "\n");
    UG_LOG("m_AreaScale = " << m_AreaScale << "\n");

	for(size_t e = 0; e < vEdges.size(); ++e)
	{
		UG_LOG("alphaOut[" << e << "] = " << m_vAlpha[e][0] << " , " << m_vAlpha[e][1] << "\n");
		UG_LOG("vIntersectionPnt[" << e << "] = " << m_vIntersectionPnts[e][0] << " , " << m_vIntersectionPnts[e][1] << "\n");
	}

	for(size_t vrtInd1 = 0; vrtInd1 < 3; ++vrtInd1)
	{
		for(size_t vrtInd2 = 0; vrtInd2 < 3; ++vrtInd2)
			UG_LOG("m_vShapeValues[vrtInd1][vrtInd2] = " << m_vShapeValues[vrtInd1][vrtInd2] << "\n");
		UG_LOG("\n");
	}
	for(size_t i = 0; i < 3; ++i)
		UG_LOG("m_vIntegralGamma = " << m_vIntegralGamma[i] << "\n");

}

template <int TWorldDim>
int InterfaceHandlerLocalDiffusion<TWorldDim>::
CollectCorners_FlatTop_2d(GridObject* elem)
{
	//////////////////////////////////////////////
	// 1) fill vector with fluid corners:
	//////////////////////////////////////////////

	this->m_vCornerCoords.clear();
	this->m_vInterfaceID.clear();
	this->m_vOriginalCornerID.clear();
	m_vRealCornerID_tri.clear();


// buffer vectors for (cornerCoords, cornerIndex)
	std::vector<std::pair<MathVector<dim>, size_t > > vOutsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vInsideCorners;
	std::vector<std::pair<MathVector<dim>, size_t > > vNearIntCorners;
	std::vector<size_t> vRealCornerID_tri; vRealCornerID_tri.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

 	bool isFTVertex = false;
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// remember boolian for check, weather there existe at least one vertex, which is FT!
		isFTVertex = this->is_onInterfaceVertex(vVertex[i], i);
		if ( isFTVertex )
			break;
	}

	if ( !isFTVertex )
		UG_THROW("Error in 'CollectCorners_FlatTop_2d': no vertex is FTVertex: should be true for at least 1 vertex!\n");

	//	collect all edges of the element
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, *this->m_spMG, elem);

	// loop vertices
	//////////////////////////////////////////////
	// REMARK:
	// order is the same as in 'vCornerCoords', therefore we can be sure, that the
	// order of the new 'vCornerIBCoords' will be consistent with the grid standard
	//////////////////////////////////////////////

    bool bNearInterface = false; this->set_bNearInterface(false);
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrtRoot = vVertex[i];

		//////////////////////////////////////////////
		// case 1:
		// vertex insideDomain
 		if ( !this->is_onInterfaceVertex(vrtRoot, i) )
 		{
			if ( this->is_nearInterfaceVertex(vrtRoot, i) )
				UG_THROW("NearInterface BUT !is_FT => neuerdings Fehler!!....\n");

			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
			this->m_vOriginalCornerID.push_back(i);
			vRealCornerID_tri.push_back(i);

			vInsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
		}
		//////////////////////////////////////////////
  		// case 2:
		// vertex = FT + ON interface
		// 		=> KEINE Berechnung von 'intersectionPoint' notwendig! -> pushen und alten index pushen

		// REMARK: is_nearInterfaceVerx = false per default, if m_vThresholdOnLevel = 0.0
		else if ( this->is_nearInterfaceVertex(vrtRoot, i) )
		{
            bNearInterface = true; this->set_bNearInterface(true);

 			this->m_vCornerCoords.push_back(this->m_aaPos[vrtRoot]);
 			this->m_vOriginalCornerID.push_back(i);
			vRealCornerID_tri.push_back(i);

 			this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

			vOutsideCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));
			vNearIntCorners.push_back(std::make_pair(this->m_aaPos[vrtRoot], i));

			get_or_insert_vertex(this->m_aaPos[vrtRoot]);

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
				CollectVertices(vVertexEdge, *this->m_spMG, edge);
				if ( vVertexEdge.size() != 2 )
					UG_THROW("error in collecting vertices associated to an edge!....EXIT!...\n");

				Vertex* vrt1 = vVertexEdge[0];
				Vertex* vrt2 = vVertexEdge[1];
				size_t vrtInd1 = this->get_vertex_index(vrt1, elem);
				size_t vrtInd2 = this->get_vertex_index(vrt2, elem);

				MathVector<dim> intersectionPnt;

	 		///////////////////////////////////////////////////////////////////
 			// does vrtRoot lie on a cutted edge?
		 	///////////////////////////////////////////////////////////////////
			// case1: vrtRoot is intersectionPnt with insideCorner = near_interface_corner => remove!
				if ( this->is_nearInterfaceVertex(vrt2, vrtInd2) || this->is_nearInterfaceVertex(vrt1, vrtInd1) )
				{ bNearInterface = true;  this->set_bNearInterface(true); continue; }
			 // case2: vert2 = outsideParticle && vrt1 = insideParticle:
				else if ( vrtRoot == vrt1 && !this->is_onInterfaceVertex(vrt2, vrtInd2) ){
					get_intersection_point(intersectionPnt, vrt2, vrt1);
 				}
			// case3: vrt1 = outsideParticle && vrt2 = insideParticle:
				else if ( vrtRoot == vrt2 && !this->is_onInterfaceVertex(vrt1, vrtInd1) )
					get_intersection_point(intersectionPnt, vrt1, vrt2);
				else
 					continue;

			// check for correct inersectionPnt
				if ( fabs(this->get_LSvalue_byPosition(intersectionPnt)) > 1e-6  )
					UG_THROW("in 'CollectIBCorners2d()': Error in computation of 'intersectionPnt':\n "
							" intersectionPnt = " << intersectionPnt << "\n distance from interace = " << fabs(this->get_LSvalue_byPosition(intersectionPnt)) << "\n");

	 		///////////////////////////////////////////////////////////////////
	 		// only push_back, if not included yet!
			// 	-> can be ONLY the case, if the intersectionPoint is a node
	 			if ( ! this->isIncluded(this->m_vCornerCoords, intersectionPnt) )
	 			{

	 				this->m_vCornerCoords.push_back(intersectionPnt);
	 				this->m_vOriginalCornerID.push_back(i);
	 				this->m_vInterfaceID.push_back(this->m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!

	 				vOutsideCorners.push_back(std::make_pair(intersectionPnt, i));

	 				size_t index = get_or_insert_vertex(intersectionPnt);
	 				vRealCornerID_tri.push_back(index);

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
	if ( (this->m_vCornerCoords.size() == 4) && (!bNearInterface) && (dim == 2) )
	{
		this->ResortQuadrilateral(vInsideCorners, vOutsideCorners, normalDir);
	// and resort m_vRealCornerID_quad:
		Resort_RealCornerID();
	// and copy data to m_vInterfaceID_quad:
		m_vInterfaceID_quad.clear();
		for ( size_t i = 0; i < this->m_vInterfaceID.size(); ++i)
			m_vInterfaceID_quad.push_back(this->m_vInterfaceID[i]);
	}
    else if ( bNearInterface && this->get_orientation() == -1 )
    {
        // m_vRealCornerID_quad was not written and needs to be called:
        // and copy data to m_vInterfaceID_tri:
        Resort_RealCornerID();
        // and copy data to m_vInterfaceID_quad:
        m_vInterfaceID_quad.clear();
        for ( size_t i = 0; i < this->m_vInterfaceID.size(); ++i)
            m_vInterfaceID_quad.push_back(this->m_vInterfaceID[i]);
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

	if ( this->m_vCornerCoords.size() == 3 )
	{
	// and copy data to m_vInterfaceID_tri:
		m_vInterfaceID_tri.clear();
		for ( size_t i = 0; i < this->m_vInterfaceID.size(); ++i)
			m_vInterfaceID_tri.push_back(this->m_vInterfaceID[i]);
	// and copy data to m_vRealCornerID_tri:
		m_vRealCornerID_tri.clear();
		for ( size_t i = 0; i < vRealCornerID_tri.size(); ++i)
			m_vRealCornerID_tri.push_back(vRealCornerID_tri[i]);
	}

    if ( this->m_vCornerCoords.size() == 0 )
        UG_THROW("m_vCornerCoords.size() = " << this->m_vCornerCoords.size() << "not possible!\n");

	return this->m_vCornerCoords.size();

}

template <int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
print_CutElementData()
{

    const char* filename = "CutElementData.";
    std::string name(filename);
    char ext[50]; sprintf(ext, "txt");
    name.append(ext);
    FILE* printFile = fopen(name.c_str(), "a");
    if ( dim == 3 )
        fprintf(printFile, "--------- New cut element --------\n\n");
    else if ( this->m_roid == ROID_QUADRILATERAL )
        fprintf(printFile, "--------- ROID_QUADRILATERAL --------\n\n");
    else if ( this->m_roid == ROID_TRIANGLE )
       fprintf(printFile, "------------ ROID_TRIANGLE ----------\n\n");
    
    for ( size_t i = 0; i < this->m_vCornerCoords.size(); ++i )
        fprintf(printFile,"Cut element corner %lu: %e, %e \n", i, this->m_vCornerCoords[i][0], this->m_vCornerCoords[i][1]);
    fprintf(printFile,"\n");
    
    for(size_t i = 0; i < this->m_vOriginalCornerID.size(); ++i)
        fprintf(printFile,"Original corner ID: %lu\n", this->m_vOriginalCornerID[i]);
    fprintf(printFile,"\n");
    
    for(size_t i = 0; i < this->m_vInterfaceID.size(); ++i)
        fprintf(printFile,"Interface corner ID: %lu\n", this->m_vInterfaceID[i]);
    fprintf(printFile,"\n");
    
    
	if ( InterfaceHandlerLocalBase<dim>::m_roid == ROID_TRIANGLE )
	{
		for ( size_t i = 0; i < m_vInterfaceID_tri.size(); ++i )
            fprintf(printFile,"Interface corner ID - tri: %lu\n", m_vInterfaceID_tri[i]);

		for ( size_t i = 0; i < m_vRealCornerID_tri.size(); ++i )
            fprintf(printFile,"Real corner ID - tri: %lu\n", m_vRealCornerID_tri[i]);
    }

	if ( InterfaceHandlerLocalBase<dim>::m_roid == ROID_QUADRILATERAL )
	{
		for ( size_t i = 0; i < m_vInterfaceID_quad.size(); ++i )
            fprintf(printFile,"Interface corner ID - quad: %lu\n", m_vInterfaceID_quad[i]);
 
		for ( size_t i = 0; i < m_vRealCornerID_quad.size(); ++i )
            fprintf(printFile,"Real corner ID - quad: %lu\n", m_vRealCornerID_quad[i]);
 	}

    fclose(printFile);

}

template <int TWorldDim>
void InterfaceHandlerLocalDiffusion<TWorldDim>::
print_Nitsche_Data()
{
	for ( size_t i = 0; i < 3; ++i )
		UG_LOG("m_vAlpha = " << m_vAlpha[i][0] << "\t" << m_vAlpha[i][1] << "\n");
	UG_LOG("\n");

	UG_LOG("vIntersectionPnt 1:" << m_vIntersectionPnts[0][0] << "\t" << m_vIntersectionPnts[0][1] << "\n");
	UG_LOG("vIntersectionPnt 2:" << m_vIntersectionPnts[1][0] << "\t" << m_vIntersectionPnts[1][1] << "\n");
	for(size_t vrtInd1 = 0; vrtInd1 < 3; ++vrtInd1)
	{
		for(size_t vrtInd2 = 0; vrtInd2 < 3; ++vrtInd2)
			UG_LOG("m_vShapeValues[vrtInd1][vrtInd2] = " << m_vShapeValues[vrtInd1][vrtInd2] << "\n");
		UG_LOG("\n");
	}
	UG_LOG("m_NormalToFace = " << m_NormalToFace[0] << " , " << m_NormalToFace[1] << "\n");
	UG_LOG("m_Gamma = " << m_Gamma << "\n");
	UG_LOG("m_Area = " << m_Area << "\n");
	UG_LOG("m_AreaOrig = " << m_AreaOrig << "\n");
	UG_LOG("m_AreaScale = " << m_AreaScale << "\n");

}

} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_DIFFUSION_TOOLS_H_ */
