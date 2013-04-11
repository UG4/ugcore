/*
 * data_evaluator_impl.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__

namespace ug{

template <typename TDomain>
template <typename TElem>
void DataEvaluator<TDomain>::prepare_timestep_elem(TElem* elem, LocalVector& u)
{

// 	prepare timestep
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get map
		const FunctionIndexMapping& map = m_vElemDisc[PT_ALL][i]->map();

	//	access disc functions
		u.access_by_map(map);

		if(m_vElemDisc[PT_ALL][i]->local_time_series_needed())
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map);

	//	prepare timestep for elem disc
		try{
			m_vElemDisc[PT_ALL][i]->fast_prep_timestep_elem(elem, u);
		}
		UG_CATCH_THROW("DataEvaluator<TDomain>::prepare_timestep_element: "
						"Cannot prepare timestep on element for IElemDisc "<<i);
	}
}


template <typename TDomain>
template <typename TElem>
void DataEvaluator<TDomain>::
prepare_elem(TElem* elem, LocalVector& u, const LocalIndices& ind,
             bool bDeriv)
{
//	remember element
	m_pElem = elem;

	UG_ASSERT(m_vElemDisc[PT_ALL].size() > 0, "No elem discs, but assembling");

//	get corners
	ElemGlobCornerCoords<TDomain, TElem>& co_coord = Provider<ElemGlobCornerCoords<TDomain, TElem> >::get();
	co_coord.update(&m_vElemDisc[PT_ALL][0]->domain(), elem);
	m_vCornerCoords = co_coord.vGlobalCorner();

// 	prepare element
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get map
		const FunctionIndexMapping& map = m_vElemDisc[PT_ALL][i]->map();

	//	access disc functions
		u.access_by_map(map);

		if(m_vElemDisc[PT_ALL][i]->local_time_series_needed())
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map);

	//	prepare for elem disc
		try{
			m_vElemDisc[PT_ALL][i]->fast_prep_elem(elem, u);
		}
		UG_CATCH_THROW("DataEvaluator<TDomain>::prep_elem: "
						"Cannot prepare element for IElemDisc "<<i);
	}

//	adjust lin defect array of imports and derivative array of exports
//	INFO: This is place here, since the 'prepare_elem' method of an element
//			disc may change the number of integration points, even if the type
//			of the element (e.g. triangle, quad) stays the same. This is the
//			case for, e.g., the NeumannBoundary element disc.
	if(bDeriv)
	{
		for(size_t i = 0; i < m_vImport[PT_ALL][MASS].size(); ++i)
			m_vImport[PT_ALL][MASS][i]->update_dof_sizes(ind);
		for(size_t i = 0; i < m_vImport[PT_ALL][STIFF].size(); ++i)
			m_vImport[PT_ALL][STIFF][i]->update_dof_sizes(ind);
		for(size_t i = 0; i < m_vImport[PT_ALL][RHS].size(); ++i)
			m_vImport[PT_ALL][RHS][i]->update_dof_sizes(ind);

		for(size_t i = 0; i < m_vDependentData.size(); ++i)
			m_vDependentData[i]->update_dof_sizes(ind);
	}

	compute_elem_data(u, elem, m_vCornerCoords, bDeriv);
}

template <typename TDomain>
template <typename TElem>
void DataEvaluator<TDomain>::
finish_timestep_elem(TElem* elem, const number time, LocalVector& u)
{
// 	finish timestep
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get map
		const FunctionIndexMapping& map = m_vElemDisc[PT_ALL][i]->map();

	//	access disc functions
		u.access_by_map(map);

		if(m_vElemDisc[PT_ALL][i]->local_time_series_needed())
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map);

	//	finish timestep for elem disc
		try{
			m_vElemDisc[PT_ALL][i]->fast_fsh_timestep_elem(elem, time, u);
		}
		UG_CATCH_THROW("DataEvaluator<TDomain>::finish_timestep_element: "
						"Cannot finish timestep on element for IElemDisc "<<i);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__ */
