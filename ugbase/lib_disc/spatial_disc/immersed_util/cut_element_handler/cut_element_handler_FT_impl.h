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
    : m_spMG(mg.operator->()), m_fctNames(fctNames),
 	  m_spInterfaceProvider(interfaceProvider)
    {
        m_vPrtIndex.resize(3);
        
        
        m_vThresholdOnLevel.resize(mg->num_levels(), 0.0);
        
        //	get position attachment
        m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
        
        // 	let position accessor access Vertex Coordinates
        if(!m_spMG->has_attachment<Vertex>(m_aPos))
            m_spMG->attach_to<Vertex>(m_aPos);
        m_aaPos.access(*m_spMG, m_aPos);
        
        // initialize data
        m_spCutMarker = make_sp(new BoolMarker);
        m_spOutsideMarker = make_sp(new BoolMarker);
        m_spFlatTopVrtMarker = make_sp(new BoolMarker);
        m_spNearInterfaceVrtMarker = make_sp(new BoolMarker);
        
        
        m_spCutMarker->assign_grid(*m_spMG);
        m_spOutsideMarker->assign_grid(*m_spMG);
        m_spFlatTopVrtMarker->assign_grid(*m_spMG);
        m_spNearInterfaceVrtMarker->assign_grid(*m_spMG);
        
        m_spCutMarker->enable_mark_inheritance(false);
        m_spOutsideMarker->enable_mark_inheritance(false);
        m_spFlatTopVrtMarker->enable_mark_inheritance(false);
        m_spNearInterfaceVrtMarker->enable_mark_inheritance(false);
        
        UG_LOG("init bool_marker done!\n");
        
    }
    
template <int TWorldDim>
CutElementHandlerFlatTop<TWorldDim>::
CutElementHandlerFlatTop(SmartPtr<MultiGrid> mg, const char* fctNames,
		SmartPtr<ParticleProviderEllipse<dim> > interfaceProvider)
	: m_spMG(mg.operator->()), m_fctNames(fctNames),
 	  m_spInterfaceProvider(interfaceProvider)
{
	m_vPrtIndex.resize(3);


	m_vThresholdOnLevel.resize(mg->num_levels(), 0.0);

 //	get position attachment
	m_aPos = GetDefaultPositionAttachment<position_attachment_type>();

// 	let position accessor access Vertex Coordinates
	if(!m_spMG->has_attachment<Vertex>(m_aPos))
		m_spMG->attach_to<Vertex>(m_aPos);
	m_aaPos.access(*m_spMG, m_aPos);

// initialize data
	m_spCutMarker = make_sp(new BoolMarker);
 	m_spOutsideMarker = make_sp(new BoolMarker);
	m_spFlatTopVrtMarker = make_sp(new BoolMarker);
	m_spNearInterfaceVrtMarker = make_sp(new BoolMarker);


	m_spCutMarker->assign_grid(*m_spMG);
 	m_spOutsideMarker->assign_grid(*m_spMG);
	m_spFlatTopVrtMarker->assign_grid(*m_spMG);
	m_spNearInterfaceVrtMarker->assign_grid(*m_spMG);

	m_spCutMarker->enable_mark_inheritance(false);
 	m_spOutsideMarker->enable_mark_inheritance(false);
	m_spFlatTopVrtMarker->enable_mark_inheritance(false);
	m_spNearInterfaceVrtMarker->enable_mark_inheritance(false);

	UG_LOG("init bool_marker done!\n");

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
        // checks weather node is transDoF OR rotDoF
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
/// call during update_marker(); necessary for time-dependent case
template<int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
clear_bool_marker()
{
    m_spNearInterfaceVrtMarker->clear();
    m_spFlatTopVrtMarker->clear();
    m_spOutsideMarker->clear();
    m_spCutMarker->clear();
    
    UG_LOG("clear bool marker done...\n");
}
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
set_nearInterface(Vertex* vrt)
{
    const number threshold = this->get_threshold(vrt);
    if ( set_nearInterface(vrt, threshold) )
    {
        this->m_spNearInterfaceVrtMarker->mark(vrt);
        return true;
    }
    else return false;
}
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
set_nearInterface(Vertex* vrt, const number threshold)
{
    number localThres;
    // default value of threshold = 0.0 => set minimal thres for 'ON interface'-case!
    if (threshold == 0.0)
        localThres = 1e-9;
    else localThres = threshold;
        
    const number LS_value = get_LSvalue(vrt, 0);
        
    if (fabs(LS_value) < localThres)
    {
        const char* filename = "0_near_vertex";
        std::string name(filename);
        char ext[50]; sprintf(ext, ".txt");
        name.append(ext);
        FILE* outputFile = fopen(name.c_str(), "a");
        fprintf(outputFile,"nearPos: %e \t %e (%e, %e)\n", this->m_aaPos[vrt][0], this->m_aaPos[vrt][1], fabs(LS_value), threshold);
        fclose(outputFile);
        
        return true;
    }
    return false;
}
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid_inverse(Vertex* vrt)
{
    bool outsideFluid = false;
    
    //	loop over all centers and pick the index with minimal distance
    for (size_t p = 0; p < this->num_particles(); ++p)
    {
        const number LS_value = get_LSvalue(vrt, p);
        
        if ( LS_value < 1e-10)
        {
            outsideFluid = true;
        }
        else if ( set_nearInterface(vrt) )
        {
            outsideFluid = true;
        }
    }
    
    return outsideFluid;
}
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid(Vertex* vrt)
{
    bool outsideFluid = false;
    
    //	loop over all centers and pick the index with minimal distance
    for (size_t p = 0; p < this->num_particles(); ++p)
    {
        const number LS_value = get_LSvalue(vrt, p);
        
        if ( LS_value > 1e-10 )
        {
            outsideFluid = true;
        }
        else if ( set_nearInterface(vrt) )
        {
            outsideFluid = true;
        }
    }
        
    return outsideFluid;
}
    
template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid_prtIndex(const int prtIndex, Vertex* vrt)
{
        
    bool outsideFluid = false;
        
    //	loop over all centers and pick the index with minimal distance
        
    const number LS_value = get_LSvalue(vrt, prtIndex);
        
    if ( LS_value > 1e-10 )
    {
        outsideFluid = true;
    }
    else if ( set_nearInterface(vrt) )
    {
        outsideFluid = true;
    }
        
    return outsideFluid;
}

template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid(Vertex* vrt, const int interfaceOrientation)
{
    if (interfaceOrientation == 1 )
        return is_outsideFluid(vrt);
    else if (interfaceOrientation == -1 )
        return is_outsideFluid_inverse(vrt);
    else UG_THROW("in 'CutElementHandler::is_outsideFluid()': no valid orientation given!\n");
}

template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_outsideFluid(int& PrtIndex, Vertex* vrt)
{
    bool outsideFluid = false;
        
    //	loop over all centers and pick the index with minimal distance
    for (size_t p = 0; p < this->num_particles(); ++p)
    {
        const number LS_value = get_LSvalue(vrt, p);
            
        if ( LS_value > 1e-10 )
        {
            PrtIndex = p;
            outsideFluid = true;
        }
        else if ( set_nearInterface(vrt) )
        {
            PrtIndex = p;
            outsideFluid = true;
        }
    }
        
    return outsideFluid;
}
    
template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_prtCoords(const int topLevel, const number deltaT)
{    
    if ( deltaT ==  0.0 )
        UG_THROW("ParticleProvider:update_prtCoords: deltaT = " << deltaT << " => no update necessary!\n");
        
    // get level index
    const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
        
    // output data:
    UG_LOG("update_prtCoords() for levIndex = " << levIndex << "\n");
    UG_LOG("update_prtCoords() for deltaT = " << deltaT << "\n");
    
    // update center
    for (size_t p = 0; p < num_particles(); ++p)
    {
#ifdef UG_PARALLEL
    std::vector<grid_base_object*> ElemList = this->m_vvvElemListCut[levIndex][p];
    UG_LOG("1_ update_prtCoords() ElemList.size(): " << ElemList.size() << "\n");
    if ( ElemList.size() == 0 ) {
        UG_LOG("2_ update_prtCoords() ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
        return;
    }
#endif
    
        MathVector<dim> transSol = get_transSol(p, 0);
        MathVector<dim> rotSol = get_rotSol(p, 0);
    
        m_spInterfaceProvider->update(deltaT, transSol, rotSol, p);
        UG_LOG("END update_prtCoords for prtIndex = " << p << "\n\n");
        
    } // end particle loop

}
    



template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_FTVertex(Vertex* vrt)
{
// some checks: (should be ok, since all markers are cleaned during 'update()')
	if ( is_outsideFluid(vrt, 1) && !m_spOutsideMarker->is_marked(vrt) )
		UG_THROW("in 'is_FTVertex(vrt, prtIndex)': inconsistent boolians case 1!\n");
	if ( !is_outsideFluid(vrt, 1) && m_spOutsideMarker->is_marked(vrt) )
		UG_THROW("in 'is_FTVertex(vrt, prtIndex)': inconsistent boolians case 2!\n");

	return m_spFlatTopVrtMarker->is_marked(vrt);
}

template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_FTVertex(Vertex* vrt, size_t vrtIndex)
{
// some checks: (should be ok, since all markers are cleaned during 'update()')
	if ( is_outsideFluid(vrt, 1) && !m_spOutsideMarker->is_marked(vrt) )
		UG_THROW("in 'is_FTVertex(vrt, prtIndex)': inconsistent boolians case 1!\n");
	if ( !is_outsideFluid(vrt, 1) && m_spOutsideMarker->is_marked(vrt) )
		UG_THROW("in 'is_FTVertex(vrt, prtIndex)': inconsistent boolians case 2!\n");

	return m_spFlatTopVrtMarker->is_marked(vrt);
}

template<int TWorldDim>
bool CutElementHandlerFlatTop<TWorldDim>::
is_FTVertex(int& PrtIndex, Vertex* vrt)
{
// some checks: (should be ok, since all markers are cleaned during 'update()')
	if ( is_outsideFluid(PrtIndex, vrt) && !m_spOutsideMarker->is_marked(vrt) )
		UG_THROW("in 'is_FlatTopVrt(vrt, prtIndex)': inconsistent boolians case 1!\n");
	if ( !is_outsideFluid(PrtIndex, vrt) && m_spOutsideMarker->is_marked(vrt) )
        UG_THROW("in 'is_FlatTopVrt(vrt, prtIndex)': inconsistent boolians case 2!\n");
    
	return m_spFlatTopVrtMarker->is_marked(vrt);
}

template<int TWorldDim>
ElementModus CutElementHandlerFlatTop<TWorldDim>::
compute_element_modus(int prtIndex, GridObject* elem, const int interfaceOrientation)
{
    bool insideFluid = false;
    bool outsideFluid = false;
    std::vector<Vertex*> vVertex;
    CollectVertices(vVertex, *m_spMG, elem);
    m_aaPos.access(*m_spMG, m_aPos);
        
        
    //	loop vertices
    for(size_t i = 0; i < vVertex.size(); ++i)
    {
        Vertex* vrt = vVertex[i];
        if ( is_outsideFluid_prtIndex(prtIndex, vrt) )
        {
            outsideFluid = true;
            m_spOutsideMarker->mark(vrt);
        }
        else
        {
            insideFluid = true;
            if ( set_nearInterface(vrt) )
                UG_THROW("CutElementHandlerFlatTop::compute_element_modus(): case 'set_nearInterface(vrt) = true' not possible!\n");
        }
            
    } // vertex loop
        
    if (insideFluid && outsideFluid) 		return CUT_BY_INTERFACE;
    else if (outsideFluid) 					return OUTSIDE_DOM;
        
    return INSIDE_DOM;
        
}
    
template<int TWorldDim>
ElementModus CutElementHandlerFlatTop<TWorldDim>::
compute_element_modus(GridObject* elem, const int interfaceOrientation)
{
 	bool insideFluid = false;
	bool outsideFluid = false;
  	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);
	m_aaPos.access(*m_spMG, m_aPos);


 //	loop vertices
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
		Vertex* vrt = vVertex[i];
 		if ( is_outsideFluid(vrt, interfaceOrientation) )
		{
			outsideFluid = true;
 			m_spOutsideMarker->mark(vrt);
 		}
		else
		{
			insideFluid = true;
			if ( set_nearInterface(vrt) )
				UG_THROW("CutElementHandlerFlatTop::compute_element_modus(): case 'set_nearInterface(vrt) = true' not possible!\n");
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
	if ( m_spCutMarker->is_marked(elem) )
	{
 		return CUT_BY_INTERFACE;
	}
	else if ( m_spOutsideMarker->is_marked(elem) )
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
#ifdef UG_PARALLEL
    UG_LOG("Process " << pcl::ProcRank() << ": num_particles(): " << num_particles() << ".\n");
#endif
   	// WIRKLICH? --> provider->num_particles() = 1...
	for (size_t p = 0; p < num_particles(); ++p)
	{
#ifdef UG_PARALLEL
		std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
 		UG_LOG("1_ update_global_indices() ElemList.size(): " << ElemList.size() << "\n");
		if ( ElemList.size() == 0 ) {
 			UG_LOG("Process " << pcl::ProcRank() << ": 2_ update_global_indices() ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif


 		grid_base_object* transVelElem;
		Vertex* transVelVertex;

		for (ListIter elemIter = m_vvvElemList[levIndex][p].begin();
				elemIter != m_vvvElemList[levIndex][p].end(); ++elemIter)
		{
		//	get element
			grid_base_object* elem = *elemIter;

		 	std::vector<Vertex*> vVertex;
			CollectVertices(vVertex, *m_spMG, elem);
			//compute_threshold(vVertex);


		//	loop vertices in order to:
			//  1) compute sum of distances
			//  2) check weather at least two vrt are 'outside_Fluid'
			number dist = 0.0;
			std::vector<Vertex*> vrtArray;
			for(size_t i = 0; i < elem->num_vertices(); ++i)
			{
				//	get vertex
				Vertex* vrt = elem->vertex(i);
				const MathVector<dim>& vrtPos = m_aaPos[vrt];
				const MathVector<dim>& center = get_center(p);

				dist += VecDistance(vrtPos, center);
				if (is_outsideFluid_prtIndex(p, vrt))
					vrtArray.push_back(vrt);
				if (is_outsideFluid_prtIndex(p, vrt) && !m_spOutsideMarker->is_marked(vrt))
					UG_THROW("Mist, immeroch falsche marker:-(\n");
			}

			if (dist < distance[p] && vrtArray.size() >= 2)
			{
				ExtraDoFs_found[p] = true;

				transVelElem = elem;
				transVelVertex = vrtArray[0];

				distance[p] = dist;
				m_vvvGlobalIndices_linearVel[levIndex][p].clear();
				m_vvvGlobalIndices_angularVel[levIndex][p].clear();

			//	create multi index
				std::vector < DoFIndex > vInd1;  // ToDo: hier std::vector?
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


		if ( !ExtraDoFs_found[p] )
			UG_THROW("ParticleHandlerGlobal:update_global_indices: no ExtraDoFs found for the " << p << "-th particle!\n");

	} // end prt-loop


}


template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_marker(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
{

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
// 	typedef typename geometry_traits<grid_base_object>::const_iterator ElemIter;
//	ElemIter iterBegin = m_spMG->begin<grid_base_object>(gridLevel.level());
//	ElemIter iterEnd = m_spMG->end<grid_base_object>(gridLevel.level());
//	ElemIter iterBegin = m_spMG->begin<grid_base_object>(m_spMG->top_level());
//	ElemIter iterEnd = m_spMG->end<grid_base_object>(m_spMG->top_level());

	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

	//	loop elements in order to compute the volume and set rhs:
	for( ; iter != iterEnd; iter++)
	{
 	//	get element
		grid_base_object* elem = *iter;

    for (int prtIndex = 0; prtIndex < num_particles(); ++prtIndex)
    {
            //int prtIndex = 0;
            
		ElementModus elemModus = compute_element_modus(prtIndex, elem);

		if ( elemModus == CUT_BY_INTERFACE )
		{
			m_vvvElemList[levIndex][prtIndex].push_back(elem);
			m_vvvElemListCut[levIndex][prtIndex].push_back(elem);

			m_spCutMarker->mark(elem);

		// for dim = 3: exclude Pyramids from assembling, i.e. set outside nodes to 'near_interface', so that element lies inside:
			bool isPyramid = false;
			if ( dim == 3 )
			{
				if ( element_is_pyramid(elem) )
				{
			 		isPyramid = true;
			 		//UG_LOG("---------------------------------> is Pyramid\n");
				}
			}

		// mark vrt in order to remove them from 'm_spOutsideMarker'-list!
			for(size_t i = 0; i < elem->num_vertices(); ++i)
			{
   				if ( is_outsideFluid_prtIndex(prtIndex, elem->vertex(i)) )
 				{
					m_spCutMarker->mark(elem->vertex(i));
 					m_spFlatTopVrtMarker->mark(elem->vertex(i));

 				// for pyramids, set ALL outside nodes to 'nearInterface' nodes
 				// => in 'get_cutMode()': numOutside == numNearInterface => original tetrahedron
 					if ( 0 ) //isPyramid )
 					{
 						m_spNearInterfaceVrtMarker->mark(elem->vertex(i));
 					}
 					else if ( set_nearInterface(elem->vertex(i)))
 					{
 						if ( !m_spNearInterfaceVrtMarker->is_marked(elem->vertex(i)))
 							UG_THROW("hmmm...muesste schon laengst markiert sein...oder noch nicht implementiert in 'set_nearInterface!!\n");
 					}
 				}
			}
		}
		else if ( elemModus == OUTSIDE_DOM )
		{
 			m_spOutsideMarker->mark(elem);
			m_vvvElemList[levIndex][prtIndex].push_back(elem);
			m_vvvElemListOutside[levIndex][prtIndex].push_back(elem);
   		}
		else // INSIDE_DOM
		{
			if ( elemModus != INSIDE_DOM )
				UG_THROW("in 'update_marker()': no case found for 'elemModus'!\n");
 		}

        } // end particle loop
	} // end elem-loop

	UG_LOG("update_marker done...\n");


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
int CutElementHandlerFlatTop<TWorldDim>::
get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd)
{
	std::pair<std::map<GridLevel,size_t>::iterator,bool> ret;
	ret = m_Map.insert ( std::pair<GridLevel,size_t>(gridLevel,m_Map.size()) );
//	if ( gridLevel == TOPLEVEL )
//		m_myProvider->set_indexTopGrid(m_Map.size());

	if (ret.second==false) {
//		std::cout << "element already existed";
//		std::cout << " with a value of " << ret.first->second << '\n';
	}
	else{
		update_multigrid_data(dd, ret.first->second);
	}

 	return ret.first->second;
}

template <int TWorldDim>
int CutElementHandlerFlatTop<TWorldDim>::
get_Index(const GridLevel& gridLevel)
{
	std::map<GridLevel,size_t>::iterator it;
	it = m_Map.find(gridLevel);
// if NO element is found for the key 'gridLevel', the iterator points to the end:
	if (it == m_Map.end())
		UG_THROW("in CutElementHandlerFlatTop::get_Index(): no data available on gridLevel " << gridLevel << "!\n");

 	return it->second;
}


template <int TWorldDim>
int CutElementHandlerFlatTop<TWorldDim>::
get_Index_old(const GridLevel& gridLevel)
{
// get level for computation:
	// if gridlevel = TOP: gridLevel->level() = -1 => get regular topLevel
	GridLevel level;
	if ( gridLevel.level() == -1 )
	{
		size_t topLev = m_spMG->num_levels()-1;
		level = GridLevel(topLev, GridLevel::LEVEL);
 	}
	else
		level = gridLevel;

// get according map entry
	std::pair<std::map<GridLevel,size_t>::iterator,bool> ret;
	ret = m_Map.insert ( std::pair<GridLevel,size_t>(level,m_Map.size()) );

// Error if entry did not exist for given level
	if ( ret.second )
		UG_THROW("CutElementHandlerFlatTop::get_Index():"
				" no data given for grid level " << level.level() << "\n");

	return ret.first->second;

}


template <int TWorldDim>
void CutElementHandlerFlatTop<TWorldDim>::
update_multigrid_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
{
// 1. 'update_marker()', because also the element lists are written,
	update_marker(dd, levIndex);
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
		CollectCornerCoordinates(vCornerCoords, *elem, m_aaPos);

		if ( m_spOutsideMarker->is_marked(elem) )
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
		if ( m_spFlatTopVrtMarker->is_marked(elem) )
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
			if ( m_spFlatTopVrtMarker->is_marked(vrt) )
			{
				filename = "File_m_pFlatTopMarker_Vrt.";
				name = filename;
				name.append(ext);
				printFile = fopen(name.c_str(), "a");
				fprintf(printFile," vCornerCoords: %e, %e \n", m_aaPos[vrt][0], m_aaPos[vrt][1]);
				fclose(printFile);
			}
			if ( m_spOutsideMarker->is_marked(vrt) )
			{
				filename = "File_m_pOutsideMarker_Vrt.";
				name = filename;
				name.append(ext);
				printFile = fopen(name.c_str(), "a");
				fprintf(printFile," vCornerCoords: %e, %e \n", m_aaPos[vrt][0], m_aaPos[vrt][1]);
				fclose(printFile);
			}
			if ( !m_spOutsideMarker->is_marked(vrt) && !m_spFlatTopVrtMarker->is_marked(vrt) )
			{
				filename = "File_m_pInsideMarker_Vrt.";
				name = filename;
				name.append(ext);
				printFile = fopen(name.c_str(), "a");
				fprintf(printFile," vCornerCoords: %e, %e \n", m_aaPos[vrt][0], m_aaPos[vrt][1]);
				fclose(printFile);
			}
			if ( m_spNearInterfaceVrtMarker->is_marked(vrt) )
				{
					filename = "File_m_pNearInterfaceVrtMarker_Vrt.";
					name = filename;
					name.append(ext);
					printFile = fopen(name.c_str(), "a");
					fprintf(printFile," vCornerCoords: %e, %e \n", m_aaPos[vrt][0], m_aaPos[vrt][1]);
					fclose(printFile);
				}

		}
	}


}

template <int TWorldDim>
template <typename TDomain>
void CutElementHandlerFlatTop<TWorldDim>::
init(ConstSmartPtr<DoFDistribution> dd, const int baseLevel, const int topLevel)
{
// clear all marks in grid
	clear_bool_marker();

// clear data associated to gridlevel
	m_Map.clear();

	update_multigrid_data(dd, topLevel);
/*
// loop all levels: update data associated to gridlevel
	for (size_t lev = topLevel; (int)lev >= baseLevel; lev--)
		print_elem_lists(dd);

*/
}



} // end ug namespace



#endif /* CUT_ELEMENT_HANDLER_IMPL_H_ */
