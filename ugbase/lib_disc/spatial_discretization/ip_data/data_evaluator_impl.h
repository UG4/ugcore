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
prepare_elem_loop(local_index_type& ind, number time, bool bMassPart)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	const ReferenceObjectID id = reference_element_type::REFERENCE_OBJECT_ID;

//	remove ip series for all used IPData
	for(size_t i = 0; i < m_vConstData.size(); ++i) m_vConstData[i]->clear();
	for(size_t i = 0; i < m_vPosData.size(); ++i)   m_vPosData[i]->clear();
	for(size_t i = 0; i < m_vDependentIPData.size(); ++i) m_vDependentIPData[i]->clear();

// 	set elem type in elem disc
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
		if(!(*m_pvElemDisc)[i]->set_roid(id))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': "
					"Cannot set geometric object type for Disc " << i <<".\n");
			return false;
		}

//	set geometric type at imports
	for(size_t i = 0; i < m_vStiffDataImport.size(); ++i)
		if(!m_vStiffDataImport[i]->set_roid(id))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': Cannot set "
					" geometric object type "<<id<<" for Import " << i <<
					" (Stiffness part).\n");
			return false;
		}
	if(bMassPart)
		for(size_t i = 0; i < m_vMassDataImport.size(); ++i)
			if(!m_vMassDataImport[i]->set_roid(id))
			{
				UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': Cannot set "
						" geometric object type "<<id<<" for Import " << i <<
						" (Mass part).\n");
				return false;
			}

//	set geometric type at exports
	for(size_t i = 0; i < m_vDataExport.size(); ++i)
		if(!m_vDataExport[i]->set_roid(id))
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_elem_loop': "
					"Cannot set geometric object type for Export " << i <<".\n");
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

//	check, that all dependent data is ready for evaluation
	for(size_t i = 0; i < m_vDependentIPData.size(); ++i)
		if(!m_vDependentIPData[i]->is_ready())
		{
			UG_LOG("ERROR in 'DataEvaluator::prepare_element': Dependent IPData "
					"(e.g. Linker or Export) is not ready for evaluation.\n");
			return false;
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
prepare_elem(TElem* elem, local_vector_type& u, const local_index_type& ind,
             bool bDeriv, bool bMassPart)
{
//	adjust lin defect array of imports
	if(bDeriv)
	{
		for(size_t i = 0; i < m_vStiffDataImport.size(); ++i)
			m_vStiffDataImport[i]->resize(ind, m_vStiffImpMap[i]);
		if(bMassPart)
			for(size_t i = 0; i < m_vMassDataImport.size(); ++i)
				m_vMassDataImport[i]->resize(ind, m_vMassImpMap[i]);
	}

//	adjust derivative array of exports
	for(size_t i = 0; i < m_vDependentIPData.size(); ++i)
		m_vDependentIPData[i]->resize(ind, m_vDependentMap[i]);

// 	prepare element
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));

		if(m_vbNeedLocTimeSeries[i])
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	prepare for elem disc
		if(!(*m_pvElemDisc)[i]->prepare_elem(elem, u))
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
