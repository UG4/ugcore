/*
 * data_evaluator_impl.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__

namespace ug{

template <typename TElem>
void DataEvaluator::
prepare_timestep_elem(TElem* elem, LocalVector& u)
{

// 	prepare timestep
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get map
		const FunctionIndexMapping& map = m_vElemDisc[PT_ALL][i].map;

	//	access disc functions
		u.access_by_map(map);

		if(m_vElemDisc[PT_ALL][i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map);

	//	prepare timestep for elem disc
		try{
			m_vElemDisc[PT_ALL][i].elemDisc->fast_prep_timestep_elem(elem, u);
		}
		UG_CATCH_THROW("DataEvaluator::prepare_timestep_element: "
						"Cannot prepare timestep on element for IElemDisc "<<i);
	}
}

template <typename TElem>
void DataEvaluator::
prepare_elem_loop(int si)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	const ReferenceObjectID id = reference_element_type::REFERENCE_OBJECT_ID;

// 	set elem type in elem disc
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
		try{m_vElemDisc[PT_ALL][i].elemDisc->set_roid(id, m_discPart);}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: "
						"Cannot set geometric object type for Disc " << i);
	}

	clear_positions_in_user_data();

// 	prepare loop (elem disc set local ip series here)
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
		try{m_vElemDisc[PT_ALL][i].elemDisc->fast_prep_elem_loop(id, si);}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: "
						"Cannot prepare element loop.");
	}

//	extract data imports and userdatas
	try{extract_imports_and_userdata(m_discPart);}
	UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: "
					"Cannot extract imports and userdata.");

//	set geometric type at imports
	for(size_t i = 0; i < m_vImport[PT_ALL][MASS].size(); ++i){
		try{m_vImport[PT_ALL][MASS][i].import->set_roid(id);}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: Cannot set  geometric "
						"object type "<<id<<" for Import "<<i<<" (Mass part).");
	}
	for(size_t i = 0; i < m_vImport[PT_ALL][STIFF].size(); ++i){
		try{m_vImport[PT_ALL][STIFF][i].import->set_roid(id);}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: Cannot set  geometric "
						"object type "<<id<<" for Import "<<i<<" (Stiff part).");
	}
	for(size_t i = 0; i < m_vImport[PT_ALL][RHS].size(); ++i){
		try{m_vImport[PT_ALL][RHS][i].import->set_roid(id);}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: Cannot set  geometric "
						"object type "<<id<<" for Import "<<i<<" (Rhs part).");
	}

//	set geometric type at dependent data
	for(size_t i = 0; i < m_vDependentData.size(); ++i){
		try{m_vDependentData[i]->set_roid(id);}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: "
						"Cannot set geometric object type for Export " << i);
		try{m_vDependentData[i]->check_setup();}
		UG_CATCH_THROW("DataEvaluator::prep_elem: Dependent UserData "<<i<<
		                " (e.g. Linker or Export) is not ready for evaluation.)");
	}

//	evaluate constant data
	for(size_t i = 0; i < m_vConstData.size(); ++i)
		m_vConstData[i]->compute(NULL, NULL, false);
}


template <typename TElem>
void DataEvaluator::
prepare_elem(TElem* elem, LocalVector& u, const LocalIndices& ind,
             bool bDeriv)
{
// 	prepare element
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get map
		const FunctionIndexMapping& map = m_vElemDisc[PT_ALL][i].map;

	//	access disc functions
		u.access_by_map(map);

		if(m_vElemDisc[PT_ALL][i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map);

	//	prepare for elem disc
		try{
			m_vElemDisc[PT_ALL][i].elemDisc->fast_prep_elem(elem, u);
		}
		UG_CATCH_THROW("DataEvaluator::prep_elem: "
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
			m_vImport[PT_ALL][MASS][i].import->set_dof_sizes(ind, m_vImport[PT_ALL][MASS][i].map);
		for(size_t i = 0; i < m_vImport[PT_ALL][STIFF].size(); ++i)
			m_vImport[PT_ALL][STIFF][i].import->set_dof_sizes(ind, m_vImport[PT_ALL][STIFF][i].map);
		for(size_t i = 0; i < m_vImport[PT_ALL][RHS].size(); ++i)
			m_vImport[PT_ALL][RHS][i].import->set_dof_sizes(ind, m_vImport[PT_ALL][RHS][i].map);

		for(size_t i = 0; i < m_vDependentData.size(); ++i)
			m_vDependentData[i]->set_dof_sizes(ind, m_vDependentMap[i]);
	}

	compute_elem_data(u, elem, bDeriv);
}

template <typename TElem>
void DataEvaluator::
finish_timestep_elem(TElem* elem, const number time, LocalVector& u)
{
// 	finish timestep
	for(size_t i = 0; i < m_vElemDisc[PT_ALL].size(); ++i)
	{
	//	get map
		const FunctionIndexMapping& map = m_vElemDisc[PT_ALL][i].map;

	//	access disc functions
		u.access_by_map(map);

		if(m_vElemDisc[PT_ALL][i].needLocTimeSeries)
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map);

	//	finish timestep for elem disc
		try{
			m_vElemDisc[PT_ALL][i].elemDisc->fast_fsh_timestep_elem(elem, time, u);
		}
		UG_CATCH_THROW("DataEvaluator::finish_timestep_element: "
						"Cannot finish timestep on element for IElemDisc "<<i);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__ */
