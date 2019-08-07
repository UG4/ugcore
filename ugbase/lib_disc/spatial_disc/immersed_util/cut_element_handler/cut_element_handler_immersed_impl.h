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
CutElementHandlerImmersed<TWorldDim>::
CutElementHandlerImmersed(SmartPtr<MultiGrid> mg, const char* fctNames,
		SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider)
    :  	  m_spInterfaceProvider(interfaceProvider),
          m_fctNames(fctNames), m_spMG(mg.operator->()),
          m_bElementNearInterface(false)
{
    m_vPrtIndex.resize(3);
    m_vvVertexMode.resize(2);
    
    m_vThresholdOnLevel.resize(mg->num_levels(), 0.0);
    
    m_verticesNearPos.clear();
    m_MapNearVertices.clear();

 //	get position attachment
	m_aPos = GetDefaultPositionAttachment<position_attachment_type>();

// 	let position accessor access Vertex Coordinates
	if(!m_spMG->has_attachment<Vertex>(m_aPos))
		m_spMG->attach_to<Vertex>(m_aPos);
	m_aaPos.access(*m_spMG, m_aPos);

}


template<int TWorldDim>
bool CutElementHandlerImmersed<TWorldDim>::
is_nearInterface(Vertex* vrt, const number threshold)
{
	number localThres;
 // default value of threshold = 0.0 => set minimal thres for 'ON interface'-case!
	if (threshold == 0.0)
		localThres = 1e-9;
	else localThres = threshold;

	const MathVector<dim>& center = get_center(0);
	const number radius = m_spInterfaceProvider->get_radius(0);
	const number dist = VecDistance(m_aaPos[vrt], center);

	if (fabs(dist - radius) < localThres)
	{
		const char* filename = "0_near_vertex";
		std::string name(filename);
		char ext[50]; sprintf(ext, ".txt");
		name.append(ext);
		FILE* outputFile = fopen(name.c_str(), "a");
		fprintf(outputFile,"nearPos: %e \t %e (%e, %e)\n", m_aaPos[vrt][0], m_aaPos[vrt][1], fabs(dist - radius), threshold);
		fclose(outputFile);

		return true;
	}
	return false;
}

template<int TWorldDim>
bool CutElementHandlerImmersed<TWorldDim>::
is_outsideDomain_inverse(Vertex* vrt)
{
// get data
	const MathVector<dim>& vrtPos = m_aaPos[vrt];
	const MathVector<dim>& center = get_center(0);
	const number radius			  = m_spInterfaceProvider->get_radius(0);

// compute
	if ( is_nearInterface(vrt) )
		return true;
	else if ( (VecDistance(vrtPos, center)-radius) > 1e-10)
		return true;

	return false;
}

template<int TWorldDim>
bool CutElementHandlerImmersed<TWorldDim>::
is_outsideDomain(Vertex* vrt)
{
// get data
 	const MathVector<dim>& vrtPos = m_aaPos[vrt];
	const MathVector<dim>& center = get_center(0);
	const number radius			  = m_spInterfaceProvider->get_radius(0);

// compute
    if ( is_nearInterface(vrt) )
    	return true;
    else if ( (VecDistance(vrtPos, center)-radius) < 1e-10)
		return true;

	return false;
}

// methode wegen komplizierten doppelt gemoppelt in 'compute_element_modes(): return bool
// fuer is_outside() und set_interface nicht optimal!
template<int TWorldDim>
VertexModus CutElementHandlerImmersed<TWorldDim>::
compute_vertex_modus(Vertex* vrt, const int interfaceOrientation)
{
	if ( is_nearInterface(vrt) )
		return ON_INTERFACE;
	else if ( is_outsideDomain(vrt, interfaceOrientation) )
		return OUTSIDE;
	else
		return INSIDE;
}
    
template <int TWorldDim>
size_t CutElementHandlerImmersed<TWorldDim>::
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
        
    if ( 0 )
    {
        UG_LOG("m_verticesNearPos.size(): " << m_verticesNearPos.size() << "\n");
        
        for ( size_t i = 0; i < m_verticesNearPos.size(); ++i)
            UG_LOG("m_verticesNearPos: " << m_verticesNearPos[i][0] << " , " << m_verticesNearPos[i][1] << "\n");
    }
        
    return ret.first->second;
}


template<int TWorldDim>
ElementModus CutElementHandlerImmersed<TWorldDim>::
compute_element_modus(GridObject* elem, const int interfaceOrientation)
{
    size_t localInd = -0.5*interfaceOrientation + 0.5;
    m_vvVertexMode[localInd].clear();
    
    bool insideFluid = false;
    bool outsideFluid = false;
    std::vector<Vertex*> vVertex;
    CollectVertices(vVertex, *m_spMG, elem);
    m_aaPos.access(*m_spMG, m_aPos);
    m_vVertexMode.clear();
    
    m_bElementNearInterface = false;
    
    size_t num_onInterface = 0;
    size_t num_outside = 0;
    size_t num_inside = 0;
        
    //	loop vertices
    for(size_t i = 0; i < vVertex.size(); ++i)
    {
        Vertex* vrt = vVertex[i];
        if ( is_nearInterface(vrt) )
        {
            num_onInterface++;
                
            get_or_insert_vertex_near(m_aaPos[vrt]);
            UG_LOG("is near Interface: " << m_aaPos[vrt][0] << "\t" << m_aaPos[vrt][1] << "\n");
            m_vVertexMode.push_back(ON_INTERFACE);
            m_vvVertexMode[localInd].push_back(ON_INTERFACE);
            outsideFluid = true;
            m_bElementNearInterface = true;
        }
        else if ( is_outsideDomain(vrt, interfaceOrientation) )
        {
            num_outside++;
            m_vVertexMode.push_back(OUTSIDE);
            m_vvVertexMode[localInd].push_back(OUTSIDE);
            outsideFluid = true;
        }
        else
        {
            num_inside++;
            m_vVertexMode.push_back(INSIDE);
            m_vvVertexMode[localInd].push_back(INSIDE);
            insideFluid = true;
                
            if ( is_nearInterface(vrt) )
                UG_THROW("CutElementHandlerImmersed::compute_element_modus(): case 'is_nearInterface(vrt) = true' not possible!\n");
        }
    } // vertex loop
        
    size_t checkOut = num_onInterface + num_outside;
    size_t checkIn = num_onInterface + num_inside;
        
    if ( checkOut == vVertex.size() )      m_elementModus = OUTSIDE_DOM;
    else if ( checkIn == vVertex.size() )  m_elementModus = INSIDE_DOM;
    else if (insideFluid && outsideFluid)
    {
        m_elementModus = CUT_BY_INTERFACE;
    }
    else if (outsideFluid)                  m_elementModus = OUTSIDE_DOM;
    else                                    m_elementModus = INSIDE_DOM;
        
    return m_elementModus;
        
}



template <int TWorldDim>
int CutElementHandlerImmersed<TWorldDim>::
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
		//update_multigrid_data(dd, ret.first->second);
	}

 	return ret.first->second;
}

template <int TWorldDim>
int CutElementHandlerImmersed<TWorldDim>::
get_Index(const GridLevel& gridLevel)
{
	std::map<GridLevel,size_t>::iterator it;
	it = m_Map.find(gridLevel);
// if NO element is found for the key 'gridLevel', the iterator points to the end:
	if (it == m_Map.end())
		UG_THROW("in CutElementHandlerImmersed::get_Index(): no data available on gridLevel " << gridLevel << "!\n");

 	return it->second;
}


template <int TWorldDim>
int CutElementHandlerImmersed<TWorldDim>::
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
		UG_THROW("CutElementHandlerImmersed::get_Index():"
				" no data given for grid level " << level.level() << "\n");

	return ret.first->second;

}



} // end ug namespace



#endif /* CUT_ELEMENT_HANDLER_IMMERSED_IMPL_H_ */
