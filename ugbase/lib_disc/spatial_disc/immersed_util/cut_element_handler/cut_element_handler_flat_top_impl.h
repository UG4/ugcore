/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef CUT_ELEMENT_HANDLER_IMPL_H_
#define CUT_ELEMENT_HANDLER_IMPL_H_

#ifdef UG_PARALLEL
    #include "pcl/pcl_interface_communicator.h"
    #include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
    #include "lib_grid/multi_grid.h"
#endif

namespace ug{

template <int TWorldDim>
CutElementHandlerFlatTop<TWorldDim>::
CutElementHandlerFlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
                         SmartPtr<ParticleProviderSphere<dim> > interfaceProvider)
    : CutElementHandlerBase<dim>(mg, interfaceProvider),
      m_fctNames(fctNames),
 	  m_spInterfaceProvider(interfaceProvider)
{
        m_vPrtIndex.resize(3);
}
    
template <int TWorldDim>
CutElementHandlerFlatTop<TWorldDim>::
CutElementHandlerFlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
                         SmartPtr<ParticleProviderEllipse<dim> > interfaceProvider)
    : CutElementHandlerBase<dim>(mg, interfaceProvider),
      m_fctNames(fctNames),
  	  m_spInterfaceProvider(interfaceProvider)
{
	m_vPrtIndex.resize(3);
}
    
// checks weather node is transDoF OR rotDoF
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_extraDoF(DoFIndex dofIndex, int levIndex)
{
// pressure dof is always NO extraDoF:
    if ( dofIndex[1] == dim )
        return false;
        
    for (size_t p = 0; p < num_particles(); ++p)
    {
        
#ifdef UG_PARALLEL
        if ( m_vvvElemListCut[levIndex][p].size() == 0 ) {
            continue;
        }
#endif
    // get data
        DoFIndex transIndex = get_transInd_Comp(levIndex, p, 0);
        DoFIndex rotIndex = get_rotInd_Comp(levIndex, p, 0);
        
    /// extra index = translational velocity DoF
        if ( dofIndex[0] == transIndex[0] )
        {
            if ( get_DoF_modus_linear(p) )
            {
                return false;
            }
            else
                return true;
        }
    /// extra index = angular velocity DoF: 2d case => only ONE velocity DoF
        else if ( dofIndex[0] == rotIndex[0] )
        {
            if ( get_DoF_modus_angular(p) )
            {
                return false;
            }
            else
            {
                if ( dim == 3 ) return true;
                if ( dim == 2 && dofIndex[1] == 0 ) return true;
            }
        }
    }
        
    return false;
}
    
template<int TWorldDim>
MathMatrix<TWorldDim,TWorldDim> CutElementHandlerFlatTop<TWorldDim>::
get_rotationMat(MathVector<TWorldDim> radialVector)
{
    MathMatrix<TWorldDim,TWorldDim> rotationMat;
    if ( dim == 2 )
    {
        rotationMat[0][0] = -radialVector[1];
        rotationMat[0][1] = 0.0;
        
        rotationMat[1][0] = radialVector[0];
        rotationMat[1][1] = 0.0;
    }
    else if ( dim == 3 )
    {
        rotationMat[0][0] = 0.0;
        rotationMat[0][1] = radialVector[2];
        rotationMat[0][2] = -radialVector[1];
        
        rotationMat[1][0] = -radialVector[2];
        rotationMat[1][1] = 0.0;
        rotationMat[1][2] = radialVector[0];
        
        rotationMat[2][0] = radialVector[1];
        rotationMat[2][1] = -radialVector[0];
        rotationMat[2][2] = 0.0;
    }
        
    return rotationMat;

}
    
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_nearInterface(Vertex* vrt, const number threshold)
{
// loop all particles:
    for (size_t p = 0; p < m_spInterfaceProvider->num_particles(); ++p)
    {
    // compute the distance between the location of vrt and the interface:
        const number LS_value = get_LSvalue(vrt, p);
            
        if (fabs(LS_value) < threshold)
        {
            set_prtIndex(p);
            return true;
        }
            
    } // end particle loop
        
    return false;
}

    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid(Vertex* vrt, const number threshold)
{
    const int orientation = this->get_orientation();
        
// loop all particles
    for (size_t p = 0; p < m_spInterfaceProvider->num_particles(); ++p)
    {
    // compute the distance between the location of vrt and the interface:
    // level set value: LS_value := radius - distance
        const number LS_value = get_LSvalue(vrt, p);
            
    // compute
        if ( LS_value > threshold && orientation == 1 )
        {
            set_prtIndex(p);
            return true;
        }
        else if ( LS_value < -threshold && orientation == -1 )
        {
          // set particle index to data
            set_prtIndex(p);
            return true;
        }
        else if ( is_nearInterface(vrt) )
        {
        // set particle index to data
            set_prtIndex(p);
            return true;
        }
            
    } // end particle loop
    
    return false;
}
  
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid(int& PrtIndex, Vertex* vrt, const number threshold)
{
    bool outsideFluid = false;
        
//	loop over all centers and pick the index with minimal distance
    for (size_t p = 0; p < this->num_particles(); ++p)
    {
    // compute the distance between the location of vrt and the interface:
    // level set value: LS_value := radius - distance
        const number LS_value = get_LSvalue(vrt, p);
            
        if ( LS_value > threshold )
        {
        // set particle index to data
            PrtIndex = p;
            outsideFluid = true;
        }
        else if ( is_nearInterface(vrt) )
        {
        // set particle index to data
            PrtIndex = p;
            outsideFluid = true;
        }
    }
    
    return outsideFluid;
}
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid_prtIndex(const int prtIndex, Vertex* vrt, const number threshold)
{
        
    bool outsideFluid = false;
        
// compute the distance between the location of vrt and the interface:
// level set value: LS_value := radius - distance
    const number LS_value = get_LSvalue(vrt, prtIndex);
        
    if ( LS_value > threshold )
    {
        outsideFluid = true;
    }
    else if ( is_nearInterface(vrt) )
    {
        outsideFluid = true;
    }
        
    return outsideFluid;
}

   
template<int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
compute_and_set_prtIndex(GridObject* elem)
{
    m_vPrtIndices.clear();
        
//	collect all vertices of the element
    std::vector<Vertex*> vVertex;
    CollectVertices(vVertex, *this->m_spMG, elem);
        
    std::vector<int> vPrtIndex(vVertex.size(), -1);
    int isOutside = 0;
    int refIndex = -1;
    
 // loop all vertices to get prtIndex
    for(size_t i = 0; i < vVertex.size(); ++i)
    {
        if ( is_outsideFluid(vPrtIndex[i], vVertex[i]) )
        {
            refIndex = vPrtIndex[i];
            m_vPrtIndices[vVertex[i]] = vPrtIndex[i];
            isOutside++;
        }
    }
    
// SPECIAL case: 1 element is cut by 2 DIFFERENT particles:
    bool prtIndexSet = true;
    if ( isOutside > 0 )
    {
        for(size_t i = 0; i < vVertex.size(); ++i)
        {
            if ( vPrtIndex[i] != refIndex && vPrtIndex[i] != -1)
            {
                prtIndexSet = false;
                UG_LOG("SPECIAL case: 1 element is cut by 2 DIFFERENT particles! EXIT...\n");
            }
        }
    }
        
// set m_prtIndex to unique prtIndex found:
    if ( prtIndexSet )
        set_prtIndex(refIndex); // set_prtIndex() already called during is_outsideFluid()!
    
}

    
template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_prtCoords(const int topLevel, const number deltaT)
{    
    if ( deltaT ==  0.0 )
        UG_THROW("ParticleProvider:update_prtCoords: deltaT = " << deltaT << " => no update necessary!\n");
    
// update center
    for (size_t p = 0; p < num_particles(); ++p)
    {
#ifdef UG_PARALLEL
// use size of member 'CutElementHandlerFlatTop::m_vvvElemListCut' in order to
// indicate, whether a particle lies on a processor or not
    const int levIndex = this->get_Index(GridLevel(topLevel, GridLevel::LEVEL));
    std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
    UG_LOG("1 CutElementHandlerFlatTop::update_prtCoords() ElemList.size(): " << ElemList.size() << "\n");
    if ( ElemList.size() == 0 ) {
        UG_LOG("2 CutElementHandlerFlatTop::update_prtCoords() ElemList.size(): "
               << ElemList.size() << " => skip assembling! \n");
        return;
    }
#endif
    
        MathVector<dim> transSol = get_transSol(p, 0);
        MathVector<dim> rotSol = get_rotSol(p, 0);
    
     // finally update the center of the p-th particle
        m_spInterfaceProvider->update(deltaT, transSol, rotSol, p);
        
    } // end particle loop

}



template<int TWorldDim>
ElementModus CutElementHandlerFlatTop<TWorldDim>::
compute_element_modus(int prtIndex, GridObject* elem)
{
    bool insideFluid = false;
    bool outsideFluid = false;
    std::vector<Vertex*> vVertex;
    CollectVertices(vVertex, *this->m_spMG, elem);
    this->m_aaPos.access(*this->m_spMG, this->m_aPos);
    
//	loop vertices
    for(size_t i = 0; i < vVertex.size(); ++i)
    {
        Vertex* vrt = vVertex[i];
        if ( is_outsideFluid_prtIndex(prtIndex, vrt) )
        {
            outsideFluid = true;
            this->m_spOutsideMarker->mark(vrt);
        }
        else
        {
            insideFluid = true;
            if ( is_nearInterface(vrt) )
                UG_THROW("CutElementHandlerFlatTop::compute_element_modus(): case 'is_nearInterface(vrt) = true' not possible!\n");
        }
            
    } // vertex loop
        
    if (insideFluid && outsideFluid) 		return CUT_BY_INTERFACE;
    else if (outsideFluid) 					return OUTSIDE_DOM;
        
    return INSIDE_DOM;
        
}
    
template<int TWorldDim>
ElementModus CutElementHandlerFlatTop<TWorldDim>::
compute_element_modus(GridObject* elem)
{
 	bool insideFluid = false;
	bool outsideFluid = false;
  	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);
	this->m_aaPos.access(*this->m_spMG, this->m_aPos);


 //	loop vertices
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrt = vVertex[i];
 		if ( is_outsideFluid(vrt) )
		{
			outsideFluid = true;
 			this->m_spOutsideMarker->mark(vrt);
 		}
		else
		{
			insideFluid = true;
			if ( is_nearInterface(vrt) )
				UG_THROW("CutElementHandlerFlatTop::compute_element_modus(): case 'is_nearInterface(vrt) = true' not possible!\n");
		}

	} // vertex loop

	if (insideFluid && outsideFluid) 		return CUT_BY_INTERFACE;
	else if (outsideFluid) 					return OUTSIDE_DOM;

	return INSIDE_DOM;

}


template <int TWorldDim>
ElementModus CutElementHandlerFlatTop<TWorldDim>::
get_element_modus(GridObject* elem)
{
    if ( !this->m_bBoolMarkerInit )
       UG_THROW("global element modus not initialized! Call the method update_interface_data()'\n");

    if ( this->m_spCutMarker->is_marked(elem) )
	{
 		return CUT_BY_INTERFACE;
	}
	else if ( this->m_spOutsideMarker->is_marked(elem) )
	{
 		return OUTSIDE_DOM;
	}
	else
		return INSIDE_DOM;

}


template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_global_indices(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
{
// get data
	FunctionGroup fctGrp(dd->function_pattern(), TokenizeString(m_fctNames));

// initialize vectors
	m_vvvGlobalIndices_linearVel.resize(levIndex + 1);
	m_vvvGlobalIndices_angularVel.resize(levIndex + 1);

	m_vvvGlobalIndices_linearVel[levIndex].clear();
	m_vvvGlobalIndices_angularVel[levIndex].clear();

	m_vvvGlobalIndices_linearVel[levIndex].resize(num_particles());
	m_vvvGlobalIndices_angularVel[levIndex].resize(num_particles());

// initialize boolian for return value
	std::vector<bool> ExtraDoFs_found;
	for (size_t p = 0; p < this->num_particles(); ++p)
		ExtraDoFs_found.push_back(false);

 // initialize distance array
	std::vector < number > distance;
	distance.clear();
	distance.resize(num_particles(), 10000.0);
   	typedef typename std::vector<grid_base_object*>::iterator ListIter;


/*
 * hier keine #ifdef UG_PARALLEL-Abfrage notwendig, da 'm_vvvElemList' leer ist
 *  auf dem entsprechendem proc => der proc ohne particle ueberspringt schleife sowieso! :)
*/

	for (size_t p = 0; p < num_particles(); ++p)
	{
#ifdef UG_PARALLEL
    // use size of member 'CutElementHandlerFlatTop::m_vvvElemListCut' in order to
    // indicate, whether a particle lies on a processor or not
		std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
 		UG_LOG("1_ update_global_indices() ElemList.size(): " << ElemList.size() << "\n");
		if ( ElemList.size() == 0 ) {
 			UG_LOG("Process " << pcl::ProcRank() << ": 2_ update_global_indices() ElemList.size(): "
                   << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif


 		grid_base_object* transVelElem;
		Vertex* transVelVertex;

    // loop all elements covered by the particle in order to pick the DoF-nodes
		for (ListIter elemIter = m_vvvElemList[levIndex][p].begin();
				     elemIter != m_vvvElemList[levIndex][p].end(); ++elemIter)
		{
		//	get element
			grid_base_object* elem = *elemIter;

        // collect vertices
		 	std::vector<Vertex*> vVertex;
			CollectVertices(vVertex, *this->m_spMG, elem);

		//	loop vertices in order to:
        //      1) compute sum of distances
        //      2) check weather at least two vrt are 'outside_Fluid'
			number dist = 0.0;
			std::vector<Vertex*> vrtArray;
			for(size_t i = 0; i < elem->num_vertices(); ++i)
			{
            //	get vertex
				Vertex* vrt = elem->vertex(i);
				const MathVector<dim>& vrtPos = this->m_aaPos[vrt];
				const MathVector<dim>& center = get_center(p);

				dist += VecDistance(vrtPos, center);
				if (is_outsideFluid_prtIndex(p, vrt))
					vrtArray.push_back(vrt);
				if (is_outsideFluid_prtIndex(p, vrt) && !this->m_spOutsideMarker->is_marked(vrt))
					UG_THROW("Mist, immeroch falsche marker:-(\n");
			}

        //  2) check weather at least two vrt are 'outside_Fluid'
        //      and pick the node nearest to the center
			if (dist < distance[p] && vrtArray.size() >= 2)
			{
				ExtraDoFs_found[p] = true;

				transVelElem = elem;
				transVelVertex = vrtArray[0];

				distance[p] = dist;
				m_vvvGlobalIndices_linearVel[levIndex][p].clear();
				m_vvvGlobalIndices_angularVel[levIndex][p].clear();

			//	create multi index
				std::vector < DoFIndex > vInd1;
				std::vector < DoFIndex > vInd2;

			//	loop all velocity components
				for (int d = 0; d < dim; ++d)
				{
                //	get fct id for comopent
					const size_t fct = fctGrp[d];
                    
                //	get multi indices
					if (dd->inner_dof_indices(vrtArray[0], fct, vInd1) != 1)
						UG_THROW("Only one index expected.");
                    
                //	get multi indices
					if (dd->inner_dof_indices(vrtArray[1], fct, vInd2) != 1)
						UG_THROW("Only one index expected.");
                    
                // set multiIndex for ExtraDoFs
					m_vvvGlobalIndices_linearVel[levIndex][p].push_back(vInd2[0]);
					m_vvvGlobalIndices_angularVel[levIndex][p].push_back(vInd1[0]);

 				} // end dim-loop

			}


		} // end elem-loop

    // security check
		if ( !ExtraDoFs_found[p] )
			UG_THROW("ParticleHandlerGlobal:update_global_indices: no ExtraDoFs found for the "
                     << p << "-th particle!\n");

	} // end prt-loop


}

    
template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_interface_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
{
    this->m_bBoolMarkerInit = true;

 // initialize vectors
	m_vvvElemList.resize(levIndex + 1);
	m_vvvElemListCut.resize(levIndex + 1);
 	m_vvvElemListOutside.resize(levIndex + 1);

 	m_vvvElemList[levIndex].clear();
	m_vvvElemListCut[levIndex].clear();
 	m_vvvElemListOutside[levIndex].clear();

	m_vvvElemList[levIndex].resize(num_particles());
 	m_vvvElemListCut[levIndex].resize(num_particles());
	m_vvvElemListOutside[levIndex].resize(num_particles());

// get data
	typedef typename domain_traits<dim>::grid_base_object grid_base_object;
	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

//	loop elements in order to compute the volume and set rhs:
	for( ; iter != iterEnd; iter++)
	{
 	//	get element
		grid_base_object* elem = *iter;

    for (size_t prtIndex = 0; prtIndex < num_particles(); ++prtIndex)
    {
		ElementModus elemModus = compute_element_modus(prtIndex, elem);

		if ( elemModus == CUT_BY_INTERFACE )
		{
			m_vvvElemList[levIndex][prtIndex].push_back(elem);
			m_vvvElemListCut[levIndex][prtIndex].push_back(elem);

			this->m_spCutMarker->mark(elem);

		// for dim = 3: exclude Pyramids from assembling, i.e. set outside nodes to 'near_interface', so that element lies inside:
			bool isPyramid = false;
			if ( dim == 3 )
			{
				if ( element_is_pyramid(elem) )
			 		isPyramid = true;
			}

		// mark vrt in order to remove them from 'm_spOutsideMarker'-list!
			for(size_t i = 0; i < elem->num_vertices(); ++i)
			{
   				if ( is_outsideFluid_prtIndex(prtIndex, elem->vertex(i)) )
 				{
					this->m_spCutMarker->mark(elem->vertex(i));
 					this->m_spInterfaceVrtMarker->mark(elem->vertex(i));

 				// for pyramids, set ALL outside nodes to 'nearInterface' nodes
 				// => in 'get_cutMode()': numOutside == numNearInterface => original tetrahedron
 					if ( 0 ) //isPyramid )
 					{
 						this->m_spNearInterfaceVrtMarker->mark(elem->vertex(i));
 					}
 					else if ( is_nearInterface(elem->vertex(i)))
 					{
 						if ( !this->m_spNearInterfaceVrtMarker->is_marked(elem->vertex(i)))
 							UG_THROW("hmmm...muesste schon laengst markiert sein...oder noch nicht implementiert in 'is_nearInterface!!\n");
 					}
 				}
			}
		}
		else if ( elemModus == OUTSIDE_DOM )
		{
 			this->m_spOutsideMarker->mark(elem);
			m_vvvElemList[levIndex][prtIndex].push_back(elem);
			m_vvvElemListOutside[levIndex][prtIndex].push_back(elem);
   		}
		else // INSIDE_DOM
		{
			if ( elemModus != INSIDE_DOM )
				UG_THROW("in 'update_interface_data()': no case found for 'elemModus'!\n");
 		}
 
        } // end particle loop
	} // end elem-loop

}



// same fuction as in class 'FlatTopHandler', but since m_spParticleHandlerLocal is not a member of
// CutElementHandlerFlatTop, there is no option for using it
template <int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
element_is_pyramid(grid_base_object* elem)
{

    size_t numOutside = 0;
    size_t numNearInterface = 0;

	for(size_t i = 0; i < elem->num_vertices(); ++i)
	{
 		if ( is_outsideFluid(elem->vertex(i)) )
		{
			numOutside++;
			if ( is_nearInterfaceVertex(elem->vertex(i), i) )
				numNearInterface++;
		}
	}

// Pyramid
	if ( numOutside == 2 && numNearInterface == 1 )
 		return true;

	return false;

}



template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_multigrid_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
{
// 1. 'update_interface_data()':
//  --> all BoolMarker and the element lists are written
	update_interface_data(dd, levIndex);
    
// 2. set 'm_vvvGlobalIndices_linearVel/m_vvvGlobalIndices_angularVel'
	update_global_indices(dd, levIndex);

}


template<int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
print_elem_lists(ConstSmartPtr<DoFDistribution> dd)
{

// get data
  	typedef typename DoFDistribution::traits<grid_base_object>::const_iterator ElemIter;
	ElemIter iterBegin = dd->begin<grid_base_object>();
	ElemIter iterEnd = dd->end<grid_base_object>();

	const char* filename;
	std::string name;
 	char ext[50]; sprintf(ext, "txt");
 	FILE* printFile;

 //	loop elements in order to compute the volume and set rhs:
	for(ElemIter elemIter = iterBegin; elemIter != iterEnd; ++elemIter)
	{
	//	get element
		grid_base_object* elem = *elemIter;

		//	get all corner coordinates
		std::vector<MathVector<dim> > vCornerCoords;
		CollectCornerCoordinates(vCornerCoords, *elem, this->m_aaPos);

		if ( this->m_spOutsideMarker->is_marked(elem) )
		{
			filename = "File_OutsideMarker.";
			name = filename;
			name.append(ext);
			printFile = fopen(name.c_str(), "a");
			fprintf(printFile, "------------------------------------------------------\n\n new triangle: \n");
			for ( size_t i = 0; i < 3; ++i )
				fprintf(printFile," vCornerCoords: %e, %e \n", vCornerCoords[i][0], vCornerCoords[i][1]);
			fclose(printFile);
		}
		if ( this->m_spInterfaceVrtMarker->is_marked(elem) )
		{
			filename = "File_TrueCutElemMarker.";
			name = filename;
			name.append(ext);
			printFile = fopen(name.c_str(), "a");
			fprintf(printFile, "------------------------------------------------------\n\n new triangle: \n");
			for ( size_t i = 0; i < 3; ++i )
				fprintf(printFile," vCornerCoords: %e, %e \n", vCornerCoords[i][0], vCornerCoords[i][1]);
			fclose(printFile);
		}


	//	loop vertices
		for(size_t i = 0; i < elem->num_vertices(); ++i)
		{
			//	get vertex
			Vertex* vrt = elem->vertex(i);
			if ( this->m_spInterfaceVrtMarker->is_marked(vrt) )
			{
				filename = "File_m_pFlatTopMarker_Vrt.";
				name = filename;
				name.append(ext);
				printFile = fopen(name.c_str(), "a");
				fprintf(printFile," vCornerCoords: %e, %e \n", this->m_aaPos[vrt][0], this->m_aaPos[vrt][1]);
				fclose(printFile);
			}
			if ( this->m_spOutsideMarker->is_marked(vrt) )
			{
				filename = "File_m_pOutsideMarker_Vrt.";
				name = filename;
				name.append(ext);
				printFile = fopen(name.c_str(), "a");
				fprintf(printFile," vCornerCoords: %e, %e \n", this->m_aaPos[vrt][0], this->m_aaPos[vrt][1]);
				fclose(printFile);
			}
			if ( !this->m_spOutsideMarker->is_marked(vrt) && !this->m_spInterfaceVrtMarker->is_marked(vrt) )
			{
				filename = "File_m_pInsideMarker_Vrt.";
				name = filename;
				name.append(ext);
				printFile = fopen(name.c_str(), "a");
				fprintf(printFile," vCornerCoords: %e, %e \n", this->m_aaPos[vrt][0], this->m_aaPos[vrt][1]);
				fclose(printFile);
			}
			if ( this->m_spNearInterfaceVrtMarker->is_marked(vrt) )
				{
					filename = "File_m_pNearInterfaceVrtMarker_Vrt.";
					name = filename;
					name.append(ext);
					printFile = fopen(name.c_str(), "a");
					fprintf(printFile," vCornerCoords: %e, %e \n", this->m_aaPos[vrt][0], this->m_aaPos[vrt][1]);
					fclose(printFile);
				}

		}
	}


}


#ifdef UG_PARALLEL
template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
synchronize_particles(int levIndex) {
        
        bool verbose = false;
        
#ifdef UG_DEBUG
        verbose = true;
#endif
        
        pcl::ProcessCommunicator com;
        
        if (active_mpi_routine == 1) {
            ///////////////////////////////////////
            //	Synchronisation using Allreduce	 //
            ///////////////////////////////////////
            for (size_t p = 0; p < num_particles(); ++p)
            {
                if (verbose) {
                    UG_LOG("Synchronize particle " << p << ".\n");
                    UG_LOG("m_vvvElemListCut.size(): " << m_vvvElemListCut.size() << ".\n");
                    UG_LOG("m_vvvElemListCut[levIndex].size(): " << m_vvvElemListCut[levIndex].size() << ".\n");
                    UG_LOG("m_vvvElemListCut[levIndex][" << p << "].size(): " << m_vvvElemListCut[levIndex][p].size() << ".\n");
                }
                std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
                if (verbose)
                    UG_LOG("Process " << pcl::ProcRank() << ": 1_ synchronize_particles() ElemList.size(): " << ElemList.size() << "\n");
                if ( ElemList.size() == 0 ) {
                    // send zero vector as center to all other processes
                    std::vector<double> zero_values(3*dim,0.0);
                    std::vector<double> values(3*dim,0.0);
                    com.allreduce(&zero_values[0], &values[0], 3*dim, MPI_DOUBLE, PCL_RO_SUM);
                    MathVector<dim> center;
                    MathVector<dim> linearVelocity;
                    MathVector<dim> angularVelocity;
                    for (size_t i = 0; i < dim; ++i) {
                        center[i] = values[i];
                        linearVelocity[i] = values[dim+i];
                        angularVelocity[i] = values[2*dim+i];
                    }
                    if (verbose) {
                        UG_LOG("Process " << pcl::ProcRank() << ": setting center to " << center << "\n");
                        UG_LOG("Process " << pcl::ProcRank() << ": setting linear velocity to " << linearVelocity << "\n");
                        UG_LOG("Process " << pcl::ProcRank() << ": setting angular velocity to " << angularVelocity << ")\n");
                    }
                    m_spInterfaceProvider->set_center(center, p);
                    m_spInterfaceProvider->set_linear_velocity(linearVelocity, p, 0);
                    m_spInterfaceProvider->set_angular_velocity(angularVelocity, p, 0);
                } else {
                    // send correct center value to all other processes. Since all other processes should send zeros the allreduce with PCL_RO_SUM results in the correct center values on all processes.
                    MathVector<dim> center = m_spInterfaceProvider->get_center(p);
                    MathVector<dim> linearVelocity = m_spInterfaceProvider->get_linear_velocity(p, 0);
                    MathVector<dim> angularVelocity = m_spInterfaceProvider->get_angular_velocity(p, 0);
                    std::vector<double> values(3*dim,0.0);
                    for (size_t i = 0; i < dim; ++i) {
                        values[i] = center[i];
                        values[dim+i] = linearVelocity[i];
                        values[2*dim+i] = angularVelocity[i];
                    }
                    MathVector<3*dim> new_values;
                    com.allreduce(&values[0], &new_values[0], 3*dim, MPI_DOUBLE, PCL_RO_SUM);
                    MathVector<dim> new_center;
                    MathVector<dim> new_linearVelocity;
                    MathVector<dim> new_angularVelocity;
                    
                    for (size_t i = 0; i < dim; ++i) {
                        new_center[i] = new_values[i];
                        new_linearVelocity[i] = new_values[dim+i];
                        new_angularVelocity[i] = new_values[2*dim+i];
                    }
                    if (verbose)
                        UG_LOG("Process " << pcl::ProcRank() << ": setting center to " << new_center << "\n");
                    m_spInterfaceProvider->set_center(new_center, p);
                    if (verbose)
                        UG_LOG("Process " << pcl::ProcRank() << ": setting linear velocity to " << new_linearVelocity << "\n");
                    m_spInterfaceProvider->set_linear_velocity(new_linearVelocity, p, 0);
                    if (verbose)
                        UG_LOG("Process " << pcl::ProcRank() << ": setting angular velocity to " << new_angularVelocity << "\n");
                    m_spInterfaceProvider->set_angular_velocity(new_angularVelocity, p, 0);
                }
            }
        } else {
            
            ////////////////////////////////////////
            //	Synchronisation using Interfaces  //
            ////////////////////////////////////////
            // Create up to date particle data struct
            if (verbose)
                UG_LOG("Create data struct.\n");
            CombinedParticleData particleValues;
            ParticleData dataSet;
            for (size_t p = 0; p < num_particles(); ++p) {
                std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
                if ( ElemList.size() == 0 ) {
                    MathVector<dim> center;
                    if (verbose)
                        UG_LOG(pcl::ProcRank() << " sends zeros.\n");
                    for (std::size_t i = 0; i < dim; ++i){
                        center[i] = 0.0;
                        dataSet.transVel[i] = 0.0;
                        dataSet.rotVel[i] = 0.0;
                    }
                    dataSet.center = center;
                    dataSet.valid = false;
                } else {
                    dataSet.center = m_spInterfaceProvider->get_center(p);
                    dataSet.transVel = m_spInterfaceProvider->get_linear_velocity(p, 0);
                    dataSet.rotVel = m_spInterfaceProvider->get_angular_velocity(p, 0);
                    dataSet.valid = true;
                    if (verbose) {
                        UG_LOG(pcl::ProcRank() << " sends\n" );
                        UG_LOG("(" << dataSet.center[0] << ", " << dataSet.center[1] << ")\n");
                        UG_LOG("(" << dataSet.transVel[0] << ", " << dataSet.transVel[1] << ")\n");
                        UG_LOG("(" << dataSet.rotVel[0] << ", " << dataSet.rotVel[1] << ")\n");
                        UG_LOG(std::boolalpha << dataSet.valid);
                    }
                }
                particleValues.push_back(dataSet);
            }
            if (verbose)
                UG_LOG("Struct created.\n");
            com.barrier();
            
            ParticleData default_new_data;
            MathVector<dim> zero_vec;
            for (std::size_t i = 0; i < dim; ++i){
                zero_vec[i] = 0.0;
            }
            default_new_data.valid = false;
            default_new_data.center = zero_vec;
            default_new_data.transVel = zero_vec;
            default_new_data.rotVel = zero_vec;
            CombinedParticleData new_particleData(num_particles(),default_new_data);
            
            const GridLayoutMap& glm = this->m_spMG->distributed_grid_manager()->grid_layout_map();
            for (std::size_t proc = 0; proc < pcl::NumProcs(); ++proc) {
                if (verbose)
                    UG_LOG("proc " << proc << " sends to ");
                if (glm.has_layout<Vertex>(INT_H_MASTER))
                {
                    const typename GridLayoutMap::Types<Vertex>::Layout& vrt_hm_layout = glm.get_layout<Vertex>(INT_H_MASTER);
                    if (vrt_hm_layout.interface_exists(proc, levIndex))
                    {
                        for (GridLayoutMap::Types<Vertex>::Layout::const_iterator it = vrt_hm_layout.begin(levIndex); it != vrt_hm_layout.begin(levIndex); ++it) {
                            com.send_data(&particleValues, sizeof(particleValues), vrt_hm_layout.proc_id(it), proc);
                            //UG_LOG(vrt_hm_layout.proc_id(it) << " ");
                        }
                    }
                }
                if (glm.has_layout<Vertex>(INT_H_SLAVE))
                {
                    const typename GridLayoutMap::Types<Vertex>::Layout& vrt_hm_layout = glm.get_layout<Vertex>(INT_H_SLAVE);
                    if (vrt_hm_layout.interface_exists(proc, levIndex))
                    {
                        for (GridLayoutMap::Types<Vertex>::Layout::const_iterator it = vrt_hm_layout.begin(levIndex); it != vrt_hm_layout.begin(levIndex); ++it) {
                            com.receive_data(&new_particleData, sizeof(new_particleData), vrt_hm_layout.proc_id(it), proc);
                        }
                    }
                }
            }
            if (verbose)
                UG_LOG("\n")
                for (size_t p = 0; p < num_particles(); ++p) {
                    if (new_particleData[p].valid) {
                        m_spInterfaceProvider->set_center(new_particleData[p].center, p);
                        m_spInterfaceProvider->set_linear_velocity(new_particleData[p].transVel, p, 0);
                        m_spInterfaceProvider->set_angular_velocity(new_particleData[p].rotVel, p, 0);
                    }
                }
        }
    }
#endif

} // end ug namespace



#endif /* CUT_ELEMENT_HANDLER_IMPL_H_ */
