/*
 * immersed_interface_base_impl.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef IMMERSED_INTERFACE_BASE_IMPL_
#define IMMERSED_INTERFACE_BASE_IMPL_


namespace ug{

///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'IInterfaceMapper'
///////////////////////////////////////////////////////////
template <typename TAlgebra>
IInterfaceMapper<TAlgebra>::
IInterfaceMapper(SmartPtr<IInterfaceHandlerLocal> localHandler)
	: m_spInterfaceHandlerLocal(localHandler)
{
}


template <typename TAlgebra>
void IInterfaceMapper<TAlgebra>::
set_identity_mat(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
 	const LocalIndices& rowInd = lmat.get_row_indices();

	for(size_t fct1=0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for(size_t dof1=0; dof1 < lmat.num_all_row_dof(fct1); ++dof1)
		{
			const size_t rowIndex = rowInd.index(fct1,dof1);
			const size_t rowComp = rowInd.comp(fct1,dof1);

			BlockRef(mat(rowIndex, rowIndex), rowComp, rowComp) = 1.0;

		}

}


template <typename TAlgebra>
void IInterfaceMapper<TAlgebra>::
add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM:
			set_identity_mat(mat, lmat, dd);
			break;
		case INSIDE_DOM:
			AddLocalMatrixToGlobal(mat, lmat);
			break;
		case CUT_BY_INTERFACE:
			add_local_mat_to_global_interface(mat, lmat, dd);
			break;
		case CUT_BY_2_INTERFACE:
			add_local_mat_to_global_interface_for2(mat, lmat, dd);
			break;
		default:
			throw(UGError("Error in IInterfaceMapper::add_local_mat_to_global()!"));

	}

}


template <typename TAlgebra>
void IInterfaceMapper<TAlgebra>::
add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM:
			// add nothing;
			break;
		case INSIDE_DOM:
			AddLocalVector(vec, lvec);
			break;
		case CUT_BY_INTERFACE:
			add_local_vec_to_global_interface(vec, lvec, dd);
			break;
		case CUT_BY_2_INTERFACE:
			add_local_vec_to_global_interface_for2(vec, lvec, dd);
			break;
		default:
			throw(UGError("Error in IInterfaceMapper::add_local_vec_to_global()!"));

	}
}


///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'IInterfaceBndCond'
///////////////////////////////////////////////////////////

template <typename TDomain>
IInterfaceBndCond<TDomain>::
IInterfaceBndCond(const char* functions, const char* subsets,
				  SmartPtr<IInterfaceHandlerLocal> localHandler)
		: IElemDisc<TDomain>(functions, subsets),
		  m_spInterfaceHandlerLocal(localHandler)
{
}

template <typename TDomain>
IInterfaceBndCond<TDomain>::
IInterfaceBndCond(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset,
				  SmartPtr<IInterfaceHandlerLocal> localHandler)
		: IElemDisc<TDomain>(vFct, vSubset),
		  m_spInterfaceHandlerLocal(localHandler)
{
}


///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'IImmersedInterface'
///////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
IImmersedInterface<TDomain, TAlgebra>::
IImmersedInterface(SmartPtr<IAssemble<TAlgebra> > ass,
				 const char* functions, const char* subsets,
				 SmartPtr<IInterfaceHandlerLocal> localHandler)
				 : m_spInterfaceHandlerLocal(localHandler),
				   m_spInterfaceMapper(new IInterfaceMapper<TDomain>(localHandler)),
				   m_spInterfaceBndCond(new IInterfaceBndCond<TDomain>(functions, subsets, localHandler))
{

	SmartPtr<AssemblingTuner<TAlgebra> > assAdapt = ass->ass_tuner();
	assAdapt->set_mapping(m_spInterfaceMapper.get());
 
}

template <typename TDomain, typename TAlgebra>
IImmersedInterface<TDomain, TAlgebra>::
IImmersedInterface(SmartPtr<IAssemble<TAlgebra> > ass,
				 const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset,
				 SmartPtr<IInterfaceHandlerLocal> localHandler)
				 : m_spInterfaceHandlerLocal(localHandler),
				   m_spInterfaceMapper(new IInterfaceMapper<TDomain>(localHandler)),
				   m_spInterfaceBndCond(new IInterfaceBndCond<TDomain>(vFct, vSubset, localHandler))
{

	SmartPtr<AssemblingTuner<TAlgebra> > assAdapt = ass->ass_tuner();
	assAdapt->set_mapping(m_spInterfaceMapper.get());
 
}


} // end ug namespace


#endif /* IMMERSED_INTERFACE_BASE_IMPL_ */
