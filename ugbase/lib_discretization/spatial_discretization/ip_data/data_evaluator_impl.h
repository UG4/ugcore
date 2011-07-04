/*
 * data_evaluator_impl.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR_IMPL__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR_IMPL__

namespace ug{

template <typename TElem>
bool
DataEvaluator::
prepare_elem_loop(local_index_type& ind, number time)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	const ReferenceObjectID id = reference_element_type::REFERENCE_OBJECT_ID;

//	remove ip series for all used IPData
	for(size_t i = 0; i < m_vConstData.size(); ++i)
		m_vConstData[i]->clear_ips();
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->clear_ips();
	for(size_t i = 0; i < m_vIDataExport.size(); ++i)
		m_vIDataExport[i]->clear_export_ips();
	for(size_t i = 0; i < m_vLinkerData.size(); ++i)
		m_vLinkerData[i]->clear_ips();

// 	set elem type in elem disc
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
		if(!(*m_pvElemDisc)[i]->set_geometric_object_type(id))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': "
					"Cannot set geometric object type for Disc " << i <<".\n");
			return false;
		}

// 	prepare loop (elem disc set local ip series here)
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
		if(!(*m_pvElemDisc)[i]->prepare_elem_loop())
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': "
					"Cannot prepare element loop.\n");
			return false;
		}

//	prepare data imports
	for(size_t i = 0; i < m_vIDataImport.size(); ++i)
	{
	//	set id for imports
		if(!m_vIDataImport[i]->set_geometric_object_type(id))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': "
					"Cannot set geometric object type for Import " << i <<".\n");
			return false;
		}
	}

//	prepare data exports
	for(size_t i = 0; i < m_vIDataExport.size(); ++i)
	{
	//	set id for imports
		if(!m_vIDataExport[i]->set_geometric_object_type(id))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': "
					"Cannot set geometric object type for Export " << i <<".\n");
			return false;
		}
	}

//	prepare data linker
	for(size_t i = 0; i < m_vLinkerDepend.size(); ++i)
	{
		if(!m_vLinkerDepend[i]->make_ready())
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_element': "
					"Linker not ready.\n");
			return false;
		}
	}

//	evaluate constant data
	for(size_t i = 0; i < m_vConstData.size(); ++i)
		m_vConstData[i]->compute();

//	we're done
	return true;
}


template <typename TElem>
bool
DataEvaluator::
prepare_elem(TElem* elem, local_vector_type& u, const local_index_type& ind)
{
//	prepare data imports
//	adjust lin defect array
	for(size_t i = 0; i < m_vIDataImport.size(); ++i)
		m_vIDataImport[i]->resize(ind, m_vMapImp[i]);

//	prepare data exports
//	adjust derivative array
	for(size_t i = 0; i < m_vIDataExport.size(); ++i)
		m_vIDataExport[i]->resize(ind, m_vMapExp[i]);

//	prepare data linker
//	adjust derivative array
	for(size_t i = 0; i < m_vLinkerDepend.size(); ++i)
		m_vLinkerDepend[i]->resize(ind, m_vMapLinker[i]);

// 	prepare element
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));

	//	prepare for elem disc
		if(!(*m_pvElemDisc)[i]->prepare_elem(elem, u, ind))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_element': "
					"Cannot prepare element for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}



} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR_IMPL__ */
