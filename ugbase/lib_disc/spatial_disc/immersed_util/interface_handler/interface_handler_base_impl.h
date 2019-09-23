/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_BASE_IMPL_H_
#define INTERFACE_HANDLER_LOCAL_BASE_IMPL_H_



namespace ug{


///////////////////////////////////////////////////////////////
//	Constructor
///////////////////////////////////////////////////////////////
template <int TWorldDim>
InterfaceHandlerLocalBase<TWorldDim>::
InterfaceHandlerLocalBase(SmartPtr<CutElementHandler_FlatTop<dim> > cutElementHandler) :
		m_spMG(cutElementHandler->m_spMG),
		m_roid(ROID_UNKNOWN),
		m_bUseStdFVAssembling(false),
        m_bPrintCutElemDataToFile(false),
        m_spCutElementHandler((SmartPtr<CutElementHandler_FlatTop<dim> >)cutElementHandler)
{
//	get position attachment
	m_aPos = GetDefaultPositionAttachment<position_attachment_type>();

// 	let position accessor access Vertex Coordinates
	if(!m_spMG->has_attachment<Vertex>(m_aPos))
		m_spMG->attach_to<Vertex>(m_aPos);
	m_aaPos.access(*m_spMG, m_aPos);

// clear member
	m_vCornerCoords.clear();
	m_vInterfaceID.clear();
	m_vNOInterfaceID.clear();
	m_vOriginalCornerID.clear();
}

template <int TWorldDim>
InterfaceHandlerLocalBase<TWorldDim>::
InterfaceHandlerLocalBase(SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler) :
    m_spMG(cutElementHandler->m_spMG),
    m_roid(ROID_UNKNOWN),
    m_bUseStdFVAssembling(false),
    m_spCutElementHandler((SmartPtr<CutElementHandler_TwoSided<dim> >)cutElementHandler)
{
    //	get position attachment
    m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
        
    // 	let position accessor access Vertex Coordinates
    if(!m_spMG->has_attachment<Vertex>(m_aPos))
        m_spMG->attach_to<Vertex>(m_aPos);
    m_aaPos.access(*m_spMG, m_aPos);
        
    // clear member
    m_vCornerCoords.clear();
    m_vInterfaceID.clear();
    m_vNOInterfaceID.clear();
    m_vOriginalCornerID.clear();
}
    
    
template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
compute_cut_element_data(GridObject* elem)
{
    
// get new element corners and according element type
	if ( StdFV_assembling() ) // Version an Stelle von 'm_bUsualAss = true'
	{
		CollectCorners_StdFV(elem);

		if ( dim == 2 ) set_roid_2d();
		if ( dim == 3 ) set_roid_3d();
	}
	else
	{
		if ( dim == 2 ) { CollectCorners_FlatTop_2d(elem); set_roid_2d(); }
		if ( dim == 3 ) { CollectCorners_FlatTop_3d(elem); set_roid_3d(); }
	}

 // some checks:
	if ( m_vCornerCoords.size() == 0 )
		UG_THROW("m_vCornerCoords.size() = " << m_vCornerCoords.size() << "not possible!\n");

	if ( dim == 2 )
		if ( elem->reference_object_id() != ROID_TRIANGLE )
			UG_THROW("Discretisation only coded for triangular elements!\n");

	if ( dim == 3 )
		if ( elem->reference_object_id() != ROID_TETRAHEDRON )
			UG_THROW("Discretisation only coded for tetrahedral elements!\n");

// output computed data to file
    if ( print_cutElment_data() )
        print_CutElementData();

}

template <int TWorldDim>
void InterfaceHandlerLocalBase<TWorldDim>::
set_element_data(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords, ReferenceObjectID roid)
{
// for elements, lying inside the domain, set the number of corners
// to the standard case, which is dim+1, since code is only valid
// for Triangles and Tetrahedra
	size_t numCo = dim+1;

// set 'm_vCornerCoords'
	m_vCornerCoords.clear();
	for ( size_t i = 0; i < numCo; ++i )
		m_vCornerCoords.push_back(vCornerCoords[i]);

// set 'm_roid'
	m_roid = roid;

// for OUTSIDE_DOM, also fill 'vOriginalCornerID' and 'vInterfaceID':
	m_vInterfaceID.clear();
	m_vOriginalCornerID.clear();
	m_vCornerCoords.clear();

//	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *m_spMG, elem);

// loop vertices and fill cut element data with the usual nodes of the element
 	for(size_t i = 0; i < vVertex.size(); ++i)
	{
	// get element
		Vertex* vrtRoot = vVertex[i];

  		if ( !is_onInterfaceVertex(vrtRoot, i) )
 		{
			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
   			m_vOriginalCornerID.push_back(i);
 		}
		else
		{
 			m_vCornerCoords.push_back(m_aaPos[vrtRoot]);
 			m_vOriginalCornerID.push_back(i);
			m_vInterfaceID.push_back(m_vCornerCoords.size()-1);  // attention: push AFTER 'm_vCornerCoords.push_back()'!!
		}

	}

}



// called by geo.update()!!
template <int TWorldDim>
bool InterfaceHandlerLocalBase<TWorldDim>::
update_elem(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords)
{
	bool do_update_local = false;
    
// compute and set the element modus:
	m_elemModus = compute_element_modus(elem);

//   m_elemModus = get_element_modus(elem);
//       --> not possible to use, since local data will be computed during 'compute_element_modus()'

	switch(m_elemModus)
	{
 		case INSIDE_DOM:	   if ( dim == 2 ) set_element_data(elem, vCornerCoords, ROID_TRIANGLE);
 							   if ( dim == 3 ) set_element_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case OUTSIDE_DOM: 	   if ( dim == 2 ) set_element_data(elem, vCornerCoords, ROID_TRIANGLE);
							   if ( dim == 3 ) set_element_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case CUT_BY_INTERFACE: compute_cut_element_data(elem);
								//if ( m_roid == ROID_PYRAMID )UG_THROW("PYRAMID\n");
							   do_update_local = true;
						  	   break;  // cut element assembling
		default:
			throw(UGError("Error in InterfaceHandlerLocalBase::update(): switch(m_elemModus)!"));
	}

 	return do_update_local;

}


} // end ug namespace



#endif /* INTERFACE_HANDLER_LOCAL_BASE_IMPL_H_ */
