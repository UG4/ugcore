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
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));

		if(m_vbNeedLocTimeSeries[i])
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	prepare timestep for elem disc
		try{
			(*m_pvElemDisc)[i]->prepare_timestep_elem(elem, u);
		}
		UG_CATCH_THROW("DataEvaluator::prepare_timestep_element: "
						"Cannot prepare timestep on element for IElemDisc "<<i);
	}
}

template <typename TElem>
void DataEvaluator::
prepare_elem_loop(bool bMassPart)
{
//	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	reference object id
	const ReferenceObjectID id = reference_element_type::REFERENCE_OBJECT_ID;

//	copy function group in import/export of element discs
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
	{
		for(size_t imp = 0; imp < (*m_pvElemDisc)[i]->num_imports(); ++imp)
			(*m_pvElemDisc)[i]->get_import(imp).set_function_group(m_vElemDiscFctGrp[i]);

		for(size_t exp = 0; exp < (*m_pvElemDisc)[i]->num_exports(); ++exp)
			(*m_pvElemDisc)[i]->get_export(exp)->set_function_group(m_vElemDiscFctGrp[i]);
	}

//	extract data imports and ipdatas
	if(!extract_imports_and_ipdata(bMassPart))
		UG_THROW_FATAL("DataEvaluator::prepare_elem_loop: "
						"Cannot extract imports and ipdata.");

// 	set elem type in elem disc
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
		if(!(*m_pvElemDisc)[i]->set_roid(id))
			UG_THROW_FATAL("DataEvaluator::prepare_elem_loop: "
							"Cannot set geometric object type for Disc " << i);

// 	prepare loop (elem disc set local ip series here)
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
	{
		try{
			(*m_pvElemDisc)[i]->prepare_elem_loop();
		}
		UG_CATCH_THROW("DataEvaluator::prepare_elem_loop: "
						"Cannot prepare element loop.");
	}

//	copy function group in import/export of element discs
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
	{
		for(size_t imp = 0; imp < (*m_pvElemDisc)[i]->num_imports(); ++imp)
			(*m_pvElemDisc)[i]->get_import(imp).set_function_group(m_vElemDiscFctGrp[i]);

		for(size_t exp = 0; exp < (*m_pvElemDisc)[i]->num_exports(); ++exp)
			(*m_pvElemDisc)[i]->get_export(exp)->set_function_group(m_vElemDiscFctGrp[i]);
	}

//	extract data imports and ipdatas
	if(!extract_imports_and_ipdata(bMassPart))
		UG_THROW_FATAL("DataEvaluator::prepare_elem_loop: "
						"Cannot extract imports and ipdata.");

//	set geometric type at imports
	for(size_t i = 0; i < m_vStiffDataImport.size(); ++i)
		if(!m_vStiffDataImport[i]->set_roid(id))
			UG_THROW_FATAL("DataEvaluator::prepare_elem_loop: Cannot set "
					" geometric object type "<<id<<" for Import " << i <<
					" (Stiffness part).");

	if(bMassPart)
		for(size_t i = 0; i < m_vMassDataImport.size(); ++i)
			if(!m_vMassDataImport[i]->set_roid(id))
				UG_THROW_FATAL("DataEvaluator::prepare_elem_loop: Cannot set "
						" geometric object type "<<id<<" for Import " << i <<
						" (Mass part).");

//	set geometric type at exports
	for(size_t i = 0; i < m_vDataExport.size(); ++i)
		if(!m_vDataExport[i]->set_roid(id))
			UG_THROW_FATAL("DataEvaluator::prepare_elem_loop: "
							"Cannot set geometric object type for Export " << i);

//	check, that all dependent data is ready for evaluation
	for(size_t i = 0; i < m_vDependentIPData.size(); ++i)
		if(!m_vDependentIPData[i]->is_ready())
			UG_THROW_FATAL("DataEvaluator::prepare_element: Dependent IPData "
							"(e.g. Linker or Export) is not ready for evaluation.");

//	evaluate constant data
	for(size_t i = 0; i < m_vConstData.size(); ++i)
		m_vConstData[i]->compute();
}


template <typename TElem>
void DataEvaluator::
prepare_elem(TElem* elem, LocalVector& u, const LocalIndices& ind,
             bool bDeriv, bool bMassPart)
{
//	adjust lin defect array of imports and derivative array of exports
	if(bDeriv)
	{
		for(size_t i = 0; i < m_vStiffDataImport.size(); ++i)
			m_vStiffDataImport[i]->set_dof_sizes(ind, m_vStiffImpMap[i]);
		if(bMassPart)
			for(size_t i = 0; i < m_vMassDataImport.size(); ++i)
				m_vMassDataImport[i]->set_dof_sizes(ind, m_vMassImpMap[i]);

		for(size_t i = 0; i < m_vDependentIPData.size(); ++i)
			m_vDependentIPData[i]->set_dof_sizes(ind, m_vDependentMap[i]);
	}

// 	prepare element
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));

		if(m_vbNeedLocTimeSeries[i])
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	prepare for elem disc
		try{
			(*m_pvElemDisc)[i]->prepare_elem(elem, u);
		}
		UG_CATCH_THROW("DataEvaluator::prepare_element: "
						"Cannot prepare element for IElemDisc "<<i);
	}
}

template <typename TElem>
void DataEvaluator::
finish_timestep_elem(TElem* elem, const number time, LocalVector& u)
{
// 	finish timestep
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));

		if(m_vbNeedLocTimeSeries[i])
			for(size_t t=0; t < m_pLocTimeSeries->size(); ++t)
				m_pLocTimeSeries->solution(t).access_by_map(map(i));

	//	finish timestep for elem disc
		try{
			(*m_pvElemDisc)[i]->finish_timestep_elem(elem, time, u);
		}
		UG_CATCH_THROW("DataEvaluator::finish_timestep_element: "
						"Cannot finish timestep on element for IElemDisc "<<i);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EVALUATOR_IMPL__ */
