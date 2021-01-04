/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef CUT_ELEMENT_HANDLER_BASE_IMPL_H_
#define CUT_ELEMENT_HANDLER_BASE_IMPL_H_


namespace ug{
    
template <int TWorldDim>
CutElementHandlerBase<TWorldDim>::
CutElementHandlerBase(SmartPtr<MultiGrid> mg,
                          SmartPtr<ParticleProviderEllipse<dim> > interfaceProvider)
    : m_spMG(mg),
    m_elementModus(INSIDE_DOM),
    m_bBoolMarkerInit(false),
    m_spInterfaceProvider((SmartPtr<ParticleProviderEllipse<dim> >) interfaceProvider)
{
        
        m_vvVertexMode.resize(2);
        m_vThresholdOnLevel.resize(mg->num_levels(), 1e-10);
        
        m_vvVertexMode.clear();
    
        //	get position attachment
        m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
        
        // 	let position accessor access Vertex Coordinates
        if(!m_spMG->has_attachment<Vertex>(m_aPos))
            m_spMG->attach_to<Vertex>(m_aPos);
        m_aaPos.access(*m_spMG, m_aPos);
        
        
        // initialize data
        m_spCutMarker = make_sp(new BoolMarker);
        m_spOutsideMarker = make_sp(new BoolMarker);
        m_spInsideMarker = make_sp(new BoolMarker);
        m_spInterfaceVrtMarker = make_sp(new BoolMarker);
        m_spNearInterfaceVrtMarker = make_sp(new BoolMarker);
        
        m_spCutMarker->assign_grid(*m_spMG);
        m_spOutsideMarker->assign_grid(*m_spMG);
        m_spInsideMarker->assign_grid(*m_spMG);
        m_spInterfaceVrtMarker->assign_grid(*m_spMG);
        m_spNearInterfaceVrtMarker->assign_grid(*m_spMG);
        
        m_spCutMarker->enable_mark_inheritance(false);
        m_spOutsideMarker->enable_mark_inheritance(false);
        m_spInsideMarker->assign_grid(*m_spMG);
        m_spInterfaceVrtMarker->enable_mark_inheritance(false);
        m_spNearInterfaceVrtMarker->enable_mark_inheritance(false);
        
        UG_LOG("init bool_marker done!\n");
        
}
    
template <int TWorldDim>
CutElementHandlerBase<TWorldDim>::
CutElementHandlerBase(SmartPtr<MultiGrid> mg,
                      SmartPtr<ParticleProviderSphere<dim> > interfaceProvider)
    : m_spMG(mg),
    m_elementModus(INSIDE_DOM),
    m_bBoolMarkerInit(false),
    m_spInterfaceProvider((SmartPtr<ParticleProviderSphere<dim> >) interfaceProvider)
{
        
        m_vvVertexMode.resize(2);
        m_vThresholdOnLevel.resize(mg->num_levels(), 1e-10);
        
        m_vvVertexMode.clear();
    
        //	get position attachment
        m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
        
        // 	let position accessor access Vertex Coordinates
        if(!m_spMG->has_attachment<Vertex>(m_aPos))
            m_spMG->attach_to<Vertex>(m_aPos);
        m_aaPos.access(*m_spMG, m_aPos);
        
        
        // initialize data
        m_spCutMarker = make_sp(new BoolMarker);
        m_spOutsideMarker = make_sp(new BoolMarker);
        m_spInsideMarker = make_sp(new BoolMarker);
        m_spInterfaceVrtMarker = make_sp(new BoolMarker);
        m_spNearInterfaceVrtMarker = make_sp(new BoolMarker);
        
        m_spCutMarker->assign_grid(*m_spMG);
        m_spOutsideMarker->assign_grid(*m_spMG);
        m_spInsideMarker->assign_grid(*m_spMG);
        m_spInterfaceVrtMarker->assign_grid(*m_spMG);
        m_spNearInterfaceVrtMarker->assign_grid(*m_spMG);
        
        m_spCutMarker->enable_mark_inheritance(false);
        m_spOutsideMarker->enable_mark_inheritance(false);
        m_spInsideMarker->assign_grid(*m_spMG);
        m_spInterfaceVrtMarker->enable_mark_inheritance(false);
        m_spNearInterfaceVrtMarker->enable_mark_inheritance(false);
        
        UG_LOG("init bool_marker done!\n");
        
    }

template <int TWorldDim>
CutElementHandlerBase<TWorldDim>::
CutElementHandlerBase(SmartPtr<MultiGrid> mg,
                      SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider)
    : m_spMG(mg),
      m_elementModus(INSIDE_DOM),
      m_bBoolMarkerInit(false),
      m_spInterfaceProvider((SmartPtr<DiffusionInterfaceProvider<dim> >) interfaceProvider)
{
    
    m_vvVertexMode.resize(2);    
    m_vThresholdOnLevel.resize(mg->num_levels(), 1e-10);
    
    m_vvVertexMode.clear();
 
    //	get position attachment
    m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
    
    // 	let position accessor access Vertex Coordinates
    if(!m_spMG->has_attachment<Vertex>(m_aPos))
        m_spMG->attach_to<Vertex>(m_aPos);
        m_aaPos.access(*m_spMG, m_aPos);
        
        
    // initialize data
    m_spCutMarker = make_sp(new BoolMarker);
    m_spOutsideMarker = make_sp(new BoolMarker);
    m_spInsideMarker = make_sp(new BoolMarker);
    m_spInterfaceVrtMarker = make_sp(new BoolMarker);
    m_spNearInterfaceVrtMarker = make_sp(new BoolMarker);
        
    m_spCutMarker->assign_grid(*m_spMG);
    m_spOutsideMarker->assign_grid(*m_spMG);
    m_spInsideMarker->assign_grid(*m_spMG);
    m_spInterfaceVrtMarker->assign_grid(*m_spMG);
    m_spNearInterfaceVrtMarker->assign_grid(*m_spMG);
        
    m_spCutMarker->enable_mark_inheritance(false);
    m_spOutsideMarker->enable_mark_inheritance(false);
    m_spInsideMarker->enable_mark_inheritance(false);
    m_spInterfaceVrtMarker->enable_mark_inheritance(false);
    m_spNearInterfaceVrtMarker->enable_mark_inheritance(false);
    
    clear_bool_marker();
    
    UG_LOG("CutElementHandlerBase constructor done!\n");
    
}
   
    
template <int TWorldDim>
template <typename TDomain>
void CutElementHandlerBase<TWorldDim>::
init(ConstSmartPtr<DoFDistribution> dd, const int baseLevel, const int topLevel)
{
// clear all marks in grid
    clear_bool_marker();
        
// clear data associated to gridlevel
    m_Map.clear();

// update cut element data
    const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL), dd);
    update_multigrid_data(dd, levIndex);

}
    
    
template <int TWorldDim>
int CutElementHandlerBase<TWorldDim>::
get_Index(const GridLevel& gridLevel, ConstSmartPtr<DoFDistribution> dd)
{
    std::pair<std::map<GridLevel,size_t>::iterator,bool> ret;
    ret = m_Map.insert ( std::pair<GridLevel,size_t>(gridLevel,m_Map.size()) );
   
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
int CutElementHandlerBase<TWorldDim>::
get_Index(const GridLevel& gridLevel)
{
    std::map<GridLevel,size_t>::iterator it;
    it = m_Map.find(gridLevel);
    // if NO element is found for the key 'gridLevel', the iterator points to the end:
    if (it == m_Map.end())
        UG_THROW("in CutElementHandler_FlatTop::get_Index(): no data available on gridLevel " << gridLevel << "!\n");
        
    return it->second;
}

    
    
template<int TWorldDim>
bool CutElementHandlerBase<TWorldDim>::
is_nearInterface(Vertex* vrt)
{
// get threshold (either default = 1e-10 or set by user via .lua-call:
    const number threshold = get_threshold(vrt);

// compute the distance between the location of vrt and the interface:
// level set value: LS_value := radius - distance
    const number LS_value = get_LSvalue(vrt);
    
    if (fabs(LS_value) < threshold)
        return true;
        
    return false;
}
    
template<int TWorldDim>
bool CutElementHandlerBase<TWorldDim>::
is_outsideDomain(Vertex* vrt)
{
// get data
    const int orientation = get_orientation();
   
// get threshold (either default = 1e-10 or set by user via .lua-call:
    const number threshold = get_threshold(vrt);

// compute the distance between the location of vrt and the interface:
// level set value: LS_value := radius - distance
    const number LS_value = get_LSvalue(vrt);
    
// a vertex is considered outside, if:
//      (A) it lies accross the interface OR
//      (B) it lies near to it according to either side
// check (A)
    if ( LS_value > threshold && orientation == 1 )
        return true;
// check (A)
    else if ( LS_value < -threshold && orientation == -1 )
        return true;
// check (B)
    else if ( is_nearInterface(vrt) )
        return true;

    
    return false;
}
    
    
template<int TWorldDim>
ElementModus CutElementHandlerBase<TWorldDim>::
get_element_modus(GridObject* elem) // --> ToDo: siehe FT impl!!
{
    if ( !m_bBoolMarkerInit )
        UG_THROW("global element modus not initialized! Call the method update_interface_data()'\n");
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
    
// This method can be used for call during 'compute_element_modes():
//      --> check 'VertexModus' instead of 'boolian' during  is_outside()/is_nearInterface()
template<int TWorldDim>
VertexModus CutElementHandlerBase<TWorldDim>::
compute_vertex_modus(Vertex* vrt)
{
    if ( is_nearInterface(vrt) )
        return ON_INTERFACE;
    else if ( is_outsideDomain(vrt) )
        return OUTSIDE;
    else
        return INSIDE;
}
    
    
template<int TWorldDim>
ElementModus CutElementHandlerBase<TWorldDim>::
compute_element_modus(GridObject* elem)
{
    
    size_t localInd = -0.5*get_orientation() + 0.5;

    m_vvVertexMode[localInd].clear();
    m_vVertexMode.clear();

    bool insideDomain = false;
    bool outsideDomain = false;
    std::vector<Vertex*> vVertex;
    CollectVertices(vVertex, *m_spMG, elem);
    m_aaPos.access(*m_spMG, m_aPos);
    
    size_t num_onInterface = 0;
    size_t num_outside = 0;
    size_t num_inside = 0;
        
//	loop vertices
    for(size_t i = 0; i < vVertex.size(); ++i)
    {
        Vertex* vrt = vVertex[i];
        if ( is_nearInterface(vrt) )
        {
        // A. update counter:
            num_onInterface++;
            
  //  get_or_insert_vertex_near(m_aaPos[vrt]); // ToDo: brauche ich diese Funktion und daten 'm_MapNearVertices' (in Immersed-Klasse) noch??
            
        // B. set local modus
            m_vVertexMode.push_back(ON_INTERFACE);
            m_vvVertexMode[localInd].push_back(ON_INTERFACE);

        // C. set global modus: BoolMarker
            m_spNearInterfaceVrtMarker->mark(vrt);
            m_spOutsideMarker->mark(vrt);

        // D. set boolian
            outsideDomain = true;

        }
        else if ( is_outsideDomain(vrt) )
        {
         // A. update counter
            num_outside++;
            
         // B. set local modus
            m_vVertexMode.push_back(OUTSIDE);
            m_vvVertexMode[localInd].push_back(OUTSIDE);

        // C. set global modus: BoolMarker
            m_spOutsideMarker->mark(vrt);
            m_spCutMarker->mark(vrt);

        // D. set boolian
            outsideDomain = true;
        }
        else
        {
        // A. update counter
            num_inside++;
            
        // B. set local modus
            m_vVertexMode.push_back(INSIDE);
            m_vvVertexMode[localInd].push_back(INSIDE);

        // C. set global modus: BoolMarker
            m_spInsideMarker->mark(vrt);

        // D. set boolian
            insideDomain = true;
                
            if ( is_nearInterface(vrt) )
                UG_THROW("CutElementHandlerBase::compute_element_modus(): case 'is_nearInterface(vrt) = true' not possible!\n");
        }
    } // vertex loop
    
    
    size_t checkOut = num_onInterface + num_outside;
    size_t checkIn = num_onInterface + num_inside;
        
    if ( checkOut == vVertex.size() )      m_elementModus = OUTSIDE_DOM;
    else if ( checkIn == vVertex.size() )  m_elementModus = INSIDE_DOM;
    else if (insideDomain && outsideDomain)
    {
        m_elementModus = CUT_BY_INTERFACE;
        m_spCutMarker->mark((grid_base_object*)elem);
    }
    else if (outsideDomain)
    {
        m_elementModus = OUTSIDE_DOM;
        m_spOutsideMarker->mark((grid_base_object*)elem);
    }
    else
    {
        m_elementModus = INSIDE_DOM;       
        m_spInsideMarker->mark((grid_base_object*)elem);
    }
    
    return m_elementModus;
        
}

    
template<int TWorldDim>
void CutElementHandlerBase<TWorldDim>::
clear_bool_marker()
{
    m_spNearInterfaceVrtMarker->clear();
    m_spInterfaceVrtMarker->clear();
    m_spOutsideMarker->clear();
    m_spInsideMarker->clear();
    m_spCutMarker->clear();
}
    
    
template <int TWorldDim>
void CutElementHandlerBase<TWorldDim>::
update_interface_data(ConstSmartPtr<DoFDistribution> dd, const int levIndex)
{
    m_bBoolMarkerInit = true;

// initialize vectors
    m_vvElemList.resize(levIndex + 1);
    m_vvElemListCut.resize(levIndex + 1);
    m_vvElemListOutside.resize(levIndex + 1);
    
    m_vvElemList[levIndex].clear();
    m_vvElemListCut[levIndex].clear();
    m_vvElemListOutside[levIndex].clear();
        
        
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
        
        ElementModus elemModus = compute_element_modus(elem);
            
        if ( elemModus == CUT_BY_INTERFACE )
        {
            m_vvElemList[levIndex].push_back(elem);
            m_vvElemListCut[levIndex].push_back(elem);
                
        // mark vrt in order to remove them from 'm_spOutsideMarker'-list!
            for(size_t i = 0; i < elem->num_vertices(); ++i)
            {
                if ( m_spOutsideMarker->is_marked(elem->vertex(i)) )
                    this->m_spInterfaceVrtMarker->mark(elem->vertex(i));
                
            // check: all nearInterface-vertices also are outside-vertices:
                if ( !m_spOutsideMarker->is_marked(elem->vertex(i)) &&
                      m_spNearInterfaceVrtMarker->is_marked(elem->vertex(i)) )
                    UG_THROW("Error: all nearInterface-vertices also are outside-vertices!!\n");
            }
        }
        else if ( elemModus == OUTSIDE_DOM )
        {
            m_vvElemList[levIndex].push_back(elem);
            m_vvElemListOutside[levIndex].push_back(elem);
        }
        else // INSIDE_DOM
        {
            if ( elemModus != INSIDE_DOM )
                UG_THROW("in 'update_interface_data()': no case found for 'elemModus'!\n");
        }
 
    } // end elem-loop
    

}
    
    
} // end ug namespace

#endif /* CUT_ELEMENT_HANDLER_BASE_IMPL_H_ */
