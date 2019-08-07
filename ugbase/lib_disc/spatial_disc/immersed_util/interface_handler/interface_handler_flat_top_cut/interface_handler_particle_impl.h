/*
 * interface_handler_local_impl.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_FLAT_TOP_IMPL_H_
#define INTERFACE_HANDLER_FLAT_TOP_IMPL_H_



namespace ug{

///////////////////////////////////////////////////////////////
/// methods for class 'IInterfaceHandlerLocal'
///////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////
/// methods for class 'InterfaceHandlerLocalParticle'
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//	Constructor
///////////////////////////////////////////////////////////////
template <int TWorldDim>
InterfaceHandlerLocalParticle<TWorldDim>::
InterfaceHandlerLocalParticle(SmartPtr<CutElementHandlerFlatTop<dim> > cutElementHandler,
			   number fluidDensity, number fluidKinVisc) :
		  InterfaceHandlerLocalBase<dim>(cutElementHandler),
          m_prtIndex(-1),
          m_fluidDensity(fluidDensity),
		  m_fluidKinVisc(fluidKinVisc),
          m_numFct(0), m_numCo(0),
 		  m_bRadial_forMassEq_equals_Normal(true),
          m_bBndDataNeeded(false),
//          m_spInterfaceProvider(interfaceProvider),
		  m_spCutElementHandler(cutElementHandler)
{

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
set_prtIndex(GridObject* elem)
{
 	int isOutside = 0;

 //	collect all vertices of the element
	std::vector<Vertex*> vVertex;
	CollectVertices(vVertex, *this->m_spMG, elem);

	std::vector<int> vPrtIndex(vVertex.size(), -1);

    m_vPrtIndices.clear();

	int refIndex = -1;
	for(size_t i = 0; i < vVertex.size(); ++i)
	{
  		if ( is_FTVertex(vVertex[i]) )
  		{
  			if ( !m_spCutElementHandler->is_outsideFluid(vPrtIndex[i], vVertex[i]) )
				UG_THROW("in InterfaceHandlerLocalParticle::get_vPrtIndex(): inconsisent boolian for vertex: "
						"is_FTVertex = true, but is_outsideFLuid = false!\n");
 			refIndex = vPrtIndex[i];
            m_vPrtIndices[vVertex[i]] = refIndex;
 			isOutside++;
  		}
	}
    
	bool set_prtIndex = true;
	if ( isOutside > 0 )
	{
		for(size_t i = 0; i < vVertex.size(); ++i)
		{
 			if ( vPrtIndex[i] != refIndex && vPrtIndex[i] != -1)
			{
				set_prtIndex = false;
				UG_LOG("SPECIAL case: 1 element is cut by 2 DIFFERENT particles! EXIT...\n");
			}
		}

	}


	// set m_prtIndex to unique prtIndex found:
	if ( set_prtIndex ) m_prtIndex = refIndex;


}

template <int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
switch_order()
{
	UG_LOG("============= VOR switch_order ===============\n");
	for ( size_t i = 0; i < this->m_vQuadriOrigID.size(); ++i  )
		UG_LOG("buffer_vQuadriOrigID[" << i << "]: " << this->m_vQuadriOrigID[i] << "\n");

// A: switch 'buffer_vCornerCoords'
	std::vector<MathVector<dim> > buffer_vCornerCoords;
	buffer_vCornerCoords.push_back(m_vQuadriCorners_for2[1]);
	buffer_vCornerCoords.push_back(m_vQuadriCorners_for2[0]);

	m_vQuadriCorners_for2.clear();
	m_vQuadriCorners_for2.push_back(buffer_vCornerCoords[0]);
	m_vQuadriCorners_for2.push_back(buffer_vCornerCoords[1]);

// B: switch 'buffer_vQuadriOrigID':
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

/*
template <int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
compute_flat_top_data_for2(GridObject* elem)
{

	if ( dim == 2 ) { CollectCorners_FlatTop_2d_for2(elem); }
	if ( dim == 3 )
		UG_THROW("compute_flat_top_data_for2() not implemented for dim = 3!\n");


 // some checks:
	if ( this->m_vCornerCoords.size() == 0 )
		UG_THROW("this->m_vCornerCoords.size() = " << this->m_vCornerCoords.size() << "not possible!\n");

	if ( dim == 2 )
		if ( elem->reference_object_id() != ROID_TRIANGLE )
			UG_THROW("Discretisation only coded for triangular elements!\n");

	if ( dim == 3 )
		if ( elem->reference_object_id() != ROID_TETRAHEDRON )
			UG_THROW("Discretisation only coded for tetrahedral elements!\n");

 // output routines
	if ( 0 )
	{
		UG_LOG("_________________ compute_flat_top_data_for2()_________________\n");

		for ( size_t i = 0; i < this->m_vCornerCoords.size(); ++i )
			UG_LOG("this->m_vCornerCoords = " << this->m_vCornerCoords[i] << "\n");
		for ( size_t i = 0; i < this->m_vOriginalCornerID.size(); ++i )
			UG_LOG("Original: id = " << this->m_vOriginalCornerID[i] << "\n");
		UG_LOG("\n");
		for ( size_t i = 0; i < this->m_vInterfaceID.size(); ++i )
			UG_LOG("Interface: id = " << this->m_vInterfaceID[i] << "\n");
		UG_LOG("\n");
	}

}
*/

// see preprocess() of flat_top.h
// called by geo.update()!!
template <int TWorldDim>
bool InterfaceHandlerLocalParticle<TWorldDim>::
update(GridObject* elem, const MathVector<TWorldDim>* vCornerCoords)
{
	bool do_update_local = false;
// reset data
	m_vBF.clear();
   
// computing flat top modus
	this->m_elemModus = this->get_element_modus(elem); // computed via 'compute_element_modus()' during 'update_marker()'

	switch(this->m_elemModus)
	{
 		case INSIDE_DOM:	   if ( dim == 2 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TRIANGLE);
 							   if ( dim == 3 ) this->set_flat_top_data(elem, vCornerCoords, ROID_TETRAHEDRON);
							   break;	// usual assembling
		case OUTSIDE_DOM: 	   if ( this->StdFV_assembling() ){
			   	   	   	   	   	   	//compute_flat_top_data(elem);
									this->set_flat_top_data(elem, vCornerCoords, ROID_UNKNOWN); // => enforce update_local() for next element! (necessary, because hereafter num_SCVF == 0!!)
  								}
								else
									this->set_flat_top_data(elem, vCornerCoords, ROID_UNKNOWN); // => enforce update_local() for next element! (necessary, because hereafter num_SCVF == 0!!)
							   break;
		case CUT_BY_INTERFACE: this->compute_flat_top_data(elem);
								//if ( this->m_roid == ROID_PYRAMID )UG_THROW("PYRAMID\n");
							   do_update_local = true;
						  	   break;  // flat top assembling
/*		case CUT_BY_2_INTERFACE: if ( this->StdFV_assembling() )
									this->compute_flat_top_data(elem);
								else
									compute_flat_top_data_for2(elem);
								do_update_local = true;
								break;*/
		default:
			throw(UGError("Error in InterfaceHandlerLocalParticle::update(): switch(m_elemModus)!"));
	}

 	return do_update_local;

}

} // end namespace ug



#endif /* INTERFACE_HANDLER_FLAT_TOP_IMPL_H_ */
