/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef CUT_ELEMENT_HANDLER_IMMERSED_IMPL_H_
#define CUT_ELEMENT_HANDLER_IMMERSED_IMPL_H_



namespace ug{


template <int TWorldDim>
CutElementHandler_TwoSided<TWorldDim>::
CutElementHandler_TwoSided(SmartPtr<MultiGrid> mg, const char* fctNames,
		SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider)
    :     CutElementHandlerBase<dim>(mg, interfaceProvider),
          m_spInterfaceProvider(interfaceProvider),
          m_bElementNearInterface(false)
{
    m_vPrtIndex.resize(3);
        
    m_verticesNearPos.clear();
    m_MapNearVertices.clear();

}


template<int TWorldDim>
bool CutElementHandler_TwoSided<TWorldDim>::
is_nearInterface(Vertex* vrt, const number threshold)
{
// loop all particles:
    for (size_t p = 0; p < m_spInterfaceProvider->num_particles(); ++p)
    {
    // compute the distance between the location of vrt and the interface:
    // level set value: LS_value := radius - distance
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
bool CutElementHandler_TwoSided<TWorldDim>::
is_outsideDomain(Vertex* vrt, const number threshold)
{
// get data
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
            set_prtIndex(p);
            return true;
        }
        else if ( is_nearInterface(vrt) )
        {
            set_prtIndex(p);
            return true;
        }
    
    } // end particle loop

    return false;
}
    
template<int TWorldDim>
bool CutElementHandler_TwoSided<TWorldDim>::
is_outsideDomain(int& PrtIndex, Vertex* vrt, const number threshold)
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
            PrtIndex = p;
            outsideFluid = true;
        }
        else if ( is_nearInterface(vrt) )
        {
            PrtIndex = p;
            outsideFluid = true;
        }
    }
    
    return outsideFluid;
}

template<int TWorldDim>
void CutElementHandler_TwoSided<TWorldDim>::
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
        if ( is_outsideDomain(vPrtIndex[i], vVertex[i]) )
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
size_t CutElementHandler_TwoSided<TWorldDim>::
get_or_insert_vertex_near(const MathVector<dim>& vrtPos)
{
    std::pair<typename std::map<MathVector<dim>,size_t>::iterator,bool> ret;
    ret = m_MapNearVertices.insert ( std::pair<MathVector<dim>,size_t>(vrtPos,m_MapNearVertices.size()) );
        
    if (ret.second==false) {
        //UG_LOG("element already existed with a value of " << ret.first->second << "\n");
    }
    else{
        m_verticesNearPos.push_back(vrtPos);
    }
        
        
    return ret.first->second;
}
    



} // end ug namespace



#endif /* CUT_ELEMENT_HANDLER_IMMERSED_IMPL_H_ */
