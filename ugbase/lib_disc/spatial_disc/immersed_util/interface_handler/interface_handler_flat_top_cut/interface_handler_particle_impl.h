/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_PARTICLE_IMPL_H_
#define INTERFACE_HANDLER_LOCAL_PARTICLE_IMPL_H_



namespace ug{


///////////////////////////////////////////////////////////////
//	Constructor
///////////////////////////////////////////////////////////////
template <int TWorldDim>
InterfaceHandlerLocalParticle<TWorldDim>::
InterfaceHandlerLocalParticle(SmartPtr<CutElementHandler_FlatTop<dim> > cutElementHandler,
			   number fluidDensity, number fluidKinVisc) :
		  InterfaceHandlerLocalBase<dim>(cutElementHandler),
          m_fluidDensity(fluidDensity),
		  m_fluidKinVisc(fluidKinVisc),
          m_numFct(0),
          m_numCo(0),
 		  m_bRadial_forMassEq_equals_Normal(true),
          m_bBndDataNeeded(false),
		  m_spCutElementHandler(cutElementHandler)
{
    m_vBF.clear();
	m_vRadialAtIP.clear();
	m_vRadialAtCo.clear();
}


template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
resize_local_indices(LocalVector& locU)
{
// 	1. get global indices
	m_ind = locU.get_indices();
	m_numFct = m_ind.num_fct();
	m_numCo = this->numCo();

// 2. get new local algebra with potentially more indices
	m_ind.resize_fct(m_numFct);
	// -> here, the original data, i.e. global index information of u.m_pIndex (-> 'u.m_pIndex.m_vvIndex[num_fct][num_dof]'),
	// is preserved, and for the new indices, global index = 0 is appended!
	for (size_t fct = 0; fct < m_numFct; ++fct)
		m_ind.resize_dof(fct, m_numCo);

	return;
}

template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::resize_local_indices(LocalVector& locU, size_t numCo)
{
// 	1. get global indices
	m_ind = locU.get_indices();
	m_numFct = m_ind.num_fct();
//	m_numCo = numCo;

// 2. get new local algebra with potentially more indices
	m_ind.resize_fct(m_numFct);
	// -> here, the original data, i.e. global index information of u.m_pIndex (-> 'u.m_pIndex.m_vvIndex[num_fct][num_dof]'),
	// is preserved, and for the new indices, global index = 0 is appended!
	for (size_t fct = 0; fct < m_numFct; ++fct)
		m_ind.resize_dof(fct, numCo);

	return;
}


template <int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
switch_order()
{
	UG_LOG("============= VOR switch_order ===============\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i  )
		UG_LOG("buffer_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

// A: switch order of 'm_vQuadriCorners_for2'
	std::vector<MathVector<dim> > buffer_vCornerCoords;
	buffer_vCornerCoords.push_back(m_vQuadriCorners_for2[1]);
	buffer_vCornerCoords.push_back(m_vQuadriCorners_for2[0]);

	m_vQuadriCorners_for2.clear();
	m_vQuadriCorners_for2.push_back(buffer_vCornerCoords[0]);
	m_vQuadriCorners_for2.push_back(buffer_vCornerCoords[1]);

// B: switch order of 'm_vQuadriOrigID':
	std::vector<size_t> buffer_vQuadriOrigID;
	buffer_vQuadriOrigID.push_back(this->m_vQuadriOrigID[1]);
	buffer_vQuadriOrigID.push_back(this->m_vQuadriOrigID[0]);

	this->m_vQuadriOrigID.clear();
	this->m_vQuadriOrigID.push_back(buffer_vQuadriOrigID[0]);
	this->m_vQuadriOrigID.push_back(buffer_vQuadriOrigID[1]);

	UG_LOG("============= NACH switch_order ===============\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i  )
		UG_LOG("m_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

}


// called by geo.update()!!
template <int TWorldDim>
bool InterfaceHandlerLocalParticle<TWorldDim>::
update(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords)
{
	bool do_update_local = false;

// needs to be cleared, since it can contain faces from former assembling
    m_vBF.clear();
    
// getting cut element modus
// REMARK: In base class, here: call of 'compute_element_modus()'
//         instead of 'get_element_modus()'.
	this->m_elemModus = this->get_element_modus(elem);  // element modus was computed via
                                                        // 'compute_element_modus()' during 'update_interface_data()'

	switch(this->m_elemModus)
	{
 		case INSIDE_DOM:	   if ( dim == 2 ) this->set_element_data(elem, vCornerCoords, ROID_TRIANGLE);
 							   if ( dim == 3 ) this->set_element_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case OUTSIDE_DOM: 	   if ( this->StdFV_assembling() ){
									this->set_element_data(elem, vCornerCoords, ROID_UNKNOWN);
  								}
								else
									this->set_element_data(elem, vCornerCoords, ROID_UNKNOWN); 
							   break;
		case CUT_BY_INTERFACE: this->compute_cut_element_data(elem);
								//if ( this->m_roid == ROID_PYRAMID )UG_THROW("PYRAMID\n");
							   do_update_local = true;
						  	   break;  // cut element assembling
/*		case CUT_BY_2_INTERFACE: if ( this->StdFV_assembling() )
									this->compute_cut_element_data(elem);
								else
									compute_cut_element_data_for2(elem);
								do_update_local = true;
								break;*/
		default:
			throw(UGError("Error in InterfaceHandlerLocalParticle::update(): switch(m_elemModus)!"));
	}

 	return do_update_local;

}

} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_PARTICLE_IMPL_H_ */
