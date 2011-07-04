/*
 * data_evaluator_impl.h
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "data_evaluator.h"

namespace ug{


bool
DataEvaluator::
set_elem_discs(const std::vector<IElemDisc*>& vElemDisc,
               const FunctionPattern& fctPat,
               bool bNonRegularGrid)
{
//	remember current elem discs
	m_pvElemDisc = &vElemDisc;

//	set function pattern to common function group
	m_commonFctGroup.set_function_pattern(fctPat);
	m_commonFctGroup.add_all();

//	create FunctionIndexMapping for each Disc
	m_vMap.resize(m_pvElemDisc->size());
	for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
	{
	//	a function group for the elem disc
		FunctionGroup elemFctGrp;

	//	create function group of this elem disc
		if(!ConvertStringToFunctionGroup(elemFctGrp, fctPat,
		                                 (*m_pvElemDisc)[i]->symb_fcts()))
		{
			UG_LOG("ERROR in 'DataEvaluator::set_elem_discs': Cannot find "
					"some symbolic Function Name for disc "<<i<<".\n");
			return false;
		}

	//	create a mapping between all functions and the function group of this
	//	element disc.
		if(!CreateFunctionIndexMapping(m_vMap[i], elemFctGrp, m_commonFctGroup))
		{
			UG_LOG("ERROR in 'DataEvaluator::set_elem_discs': Cannot create "
					"Function Index Mapping for disc "<<i<<".\n");
			return false;
		}

	//	set function group in import/export
		for(size_t imp = 0; imp < (*m_pvElemDisc)[i]->num_imports(); ++imp)
			(*m_pvElemDisc)[i]->get_import(imp).set_function_group(elemFctGrp);
		for(size_t exp = 0; exp < (*m_pvElemDisc)[i]->num_exports(); ++exp)
			(*m_pvElemDisc)[i]->get_export(exp).set_function_group(elemFctGrp);
	}

//	setup for non-regular grid
	if(!set_non_regular_grid(bNonRegularGrid))
	{
		UG_LOG("ERROR in 'DataEvaluator::set_elem_discs':"
				"Cannot setup non-regular grid.\n");
		return false;
	}

	return extract_imports_and_ipdata();
}


void
DataEvaluator::
clear()
{
//	remove imports scheduled for evaluation of lin defect
	m_vMapImp.clear();
	m_vMapImpConn.clear();
	m_vIDataImport.clear();

// 	remove data scheduled for evaluation
	m_vConstData.clear();
	m_vPosData.clear();

	m_vDependData.clear();
	m_vMapDepend.clear();

	m_vMapExp.clear();
	m_vIDataExport.clear();

	m_vMapLinker.clear();
	m_vLinkerDepend.clear();
	m_vLinkerData.clear();
}


bool
DataEvaluator::
add_data_to_eval_data(std::vector<IIPData*>& vEvalData,
                      std::vector<IIPData*>& vTryingToAdd)
{
//	if empty, we're done
	if(vTryingToAdd.empty()) return true;

//	search for element in already scheduled data
	std::vector<IIPData*>::iterator it, itEnd;
	it = find(vEvalData.begin(), vEvalData.end(), vTryingToAdd.back());

//	if found, skip this data
	if(it != vEvalData.end())
	{
		vTryingToAdd.pop_back();
		return true;
	}

//	search if element already contained in list. Then, the element
//	did start the adding procedure before and a circle dependency
//	is found
	itEnd = vTryingToAdd.end(); itEnd--;
	it = find(vTryingToAdd.begin(), itEnd, *itEnd);

//	if found, return error of circle dependency
	if(it != itEnd)
	{
		UG_LOG("ERROR in 'DataEvaluator::add_data_to_eval_data':"
				" Circle dependency of data detected for IP Data.\n");
		return false;
	}

//	add all dependent datas
	IIPData* data = vTryingToAdd.back();
	for(size_t i = 0; i < data->num_needed_data(); ++i)
	{
	//	add each data separately
		vTryingToAdd.push_back(data->needed_data(i));
		if(!add_data_to_eval_data(vEvalData, vTryingToAdd))
			return false;
	}

//	add this data to the evaluation list
	vEvalData.push_back(data);

//	pop last one, since now added to eval list
	vTryingToAdd.pop_back();

//	we're done
	return true;
}


bool
DataEvaluator::
extract_imports_and_ipdata()
{
//	clear imports and ipdata
	clear();

//	queue for all ip data needed
	std::vector<IIPData*> vEvalData;
	std::vector<IIPData*> vTryingToAdd;

//	loop elem discs
	for(size_t d = 0; d < m_pvElemDisc->size(); ++d)
	{
	//	loop imports
		for(size_t i = 0; i < (*m_pvElemDisc)[d]->num_imports(); ++i)
		{
		//	get import
			IDataImport* iimp = &((*m_pvElemDisc)[d]->get_import(i));

		//	skip non-given data (no need for evaluation)
			if(!iimp->data_given()) continue;

		//	push export on stack of needed data
			vTryingToAdd.push_back(iimp->get_data());

		//	add data and all dependency to evaluation list
			if(!add_data_to_eval_data(vEvalData, vTryingToAdd))
			{
				UG_LOG("ERROR in DataEvaluator::extract_imports_and_ipdata:"
						" Circle dependency of data detected for IP Data.\n");
				return false;
			}
			UG_ASSERT(vTryingToAdd.empty(), "All must have been added.");

		//	done iff zero-derivative
			if(iimp->zero_derivative()) continue;

		//	schedule import for computation of lin defect:
			m_vIDataImport.push_back(iimp);

		//	cast to dependent data
			IDependentIPData* dependData =
					dynamic_cast<IDependentIPData*>(iimp->get_data());

		//	check success
			if(dependData == NULL)
			{
				UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
						" Data seems dependent, but cast failed.\n");
				return false;
			}

		//	create FuncMap
			FunctionIndexMapping map;
			if(!CreateFunctionIndexMapping(map,
										   dependData->get_function_group(),
										   m_commonFctGroup))
			{
				UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
						"Cannot create Function Index Mapping for disc.\n");
				return false;
			}


		//	Remember FuncMap for lin defect (the same map as for elem disc)
			m_vMapImp.push_back(m_vMap[d]);

		//  Remember FuncMap for connected data
			m_vMapImpConn.push_back(map);
		}
	}

//	loop all needed ip data and group it
	for(size_t i = 0; i < vEvalData.size(); ++i)
	{
		IIPData* ipData = vEvalData[i];

	//	detect constant import
		if(ipData->constant_data())
		{
		//	schedule for evaluation of constant data
			m_vConstData.push_back(ipData);
			continue;
		}

	//	detect position dependent import
		if(ipData->zero_derivative())
		{
		//	schedule for evaluation of position dependent data
			m_vPosData.push_back(ipData);
			continue;
		}

	//	cast to dependent data
		IDependentIPData* dependData =
				dynamic_cast<IDependentIPData*>(ipData);

	//	check success
		if(dependData == NULL)
		{
			UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
					" Data seems dependent, but cast failed.\n");
			return false;
		}

	//	create FuncMap
		FunctionIndexMapping map;
		if(!CreateFunctionIndexMapping(map,
		                               dependData->get_function_group(),
									   m_commonFctGroup))
		{
			UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
					"Cannot create Function Index Mapping for disc.\n");
			return false;
		}


	//	save as dependent data
		m_vDependData.push_back(ipData);
		m_vMapDepend.push_back(map);

	//	cast to data export
		IDataExport* exp = dynamic_cast<IDataExport*>(ipData);

	//	Data Export case
		if(exp != NULL)
		{
		//	schedule for evaluation of IDataExports
			m_vIDataExport.push_back(exp);

		// 	remember function map
			m_vMapExp.push_back(map);
		}
		else
	//	Linker case
		{
		//	schedule for evaluation of linker
			m_vLinkerData.push_back(ipData);
			m_vLinkerDepend.push_back(dependData);

		// 	remember function map
			m_vMapLinker.push_back(map);
		}
	}

//	we're done
	return true;
}



bool
DataEvaluator::
compute_elem_data(local_vector_type & u,
                  local_index_type& ind,
                  bool computeDeriv)
{
//	evaluate position data
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->compute();

// 	process dependent data
	for(size_t i = 0; i < m_vDependData.size(); ++i)
	{
	//	cast to data export
		IDataExport* exp = dynamic_cast<IDataExport*>(m_vDependData[i]);

	//	linker case
		if(exp == NULL)
		{
			m_vDependData[i]->compute(computeDeriv);
		}
	//	export case
		else
		{
			u.access_by_map(m_vMapDepend[i]);
			if(!exp->compute_export(u, computeDeriv))
			{
				UG_LOG("ERROR in 'DataEvaluator::compute_elem_data':"
						"Cannot compute data for Export " << i <<".\n");
				return false;
			}
		}
	}
//	we're done
	return true;
}


bool
DataEvaluator::
compute_lin_defect_JA(local_vector_type & u, local_index_type& ind)
{
//	compute linearized defect
	for(size_t i = 0; i < m_vIDataImport.size(); ++i)
	{
	//	skip imports that are not located in mass part
		if(m_vIDataImport[i]->in_mass_part()) continue;

	//	set correct access for import
		u.access_by_map(m_vMapImp[i]);

	//	compute linearization of defect
		if(!m_vIDataImport[i]->compute_lin_defect(u))
		{
			UG_LOG("ERROR in 'DataEvaluator::compute_elem_data':"
					"Cannot compute lin defect for Import " << i <<".\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool
DataEvaluator::
compute_lin_defect_JM(local_vector_type & u, local_index_type& ind)
{
//	compute linearized defect
	for(size_t i = 0; i < m_vIDataImport.size(); ++i)
	{
	//	skip imports that are not located in mass part
		if(!m_vIDataImport[i]->in_mass_part()) continue;

	//	set correct access for import
		u.access_by_map(m_vMapImp[i]);

	//	compute linearization of defect
		if(!m_vIDataImport[i]->compute_lin_defect(u))
		{
			UG_LOG("ERROR in 'DataEvaluator::compute_elem_data':"
					"Cannot compute lin defect for Import " << i <<".\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool
DataEvaluator::
add_coupl_JA(local_matrix_type& J, local_index_type& indRow)
{

	for(size_t i = 0; i < m_vIDataImport.size(); ++i)
	{
	//	skip imports that are not located in mass part
		if(m_vIDataImport[i]->in_mass_part()) continue;

	//	rows are given by import
	//	cols are given by export
		J.access_by_map(m_vMapImp[i], m_vMapImpConn[i]);

	//	add off diagonal coupling
		m_vIDataImport[i]->assemble_jacobian(J);
	}

//	we're done
	return true;
}


bool
DataEvaluator::
add_coupl_JM(local_matrix_type& J, local_index_type& indRow)
{
	for(size_t i = 0; i < m_vIDataImport.size(); ++i)
	{
	//	skip imports that are not located in mass part
		if(!m_vIDataImport[i]->in_mass_part()) continue;

	//	rows are given by import
	//	cols are given by export
		J.access_by_map(m_vMapImp[i], m_vMapImpConn[i]);

	//	add off diagonal coupling
		m_vIDataImport[i]->assemble_jacobian(J);
	}

//	we're done
	return true;
}



bool
DataEvaluator::
assemble_JA(local_matrix_type& A, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		A.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->assemble_JA(A, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_JA': "
					"Cannot assemble JA for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool
DataEvaluator::
assemble_JM(local_matrix_type& M, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		M.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->assemble_JM(M, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_JM': "
					"Cannot assemble JM for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool
DataEvaluator::
assemble_A(local_vector_type& d, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		d.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->assemble_A(d, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_A': "
					"Cannot assemble Defect (A) for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool
DataEvaluator::
assemble_M(local_vector_type& d, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		d.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->assemble_M(d, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_M': "
					"Cannot assemble Defect (M) for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}



bool
DataEvaluator::
assemble_rhs(local_vector_type& rhs)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		rhs.access_by_map(map(i));

	//	assemble rhs
		if(!(*m_pvElemDisc)[i]->assemble_f(rhs))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_rhs': "
					"Cannot assemble rhs for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool
DataEvaluator::
set_time_dependent(bool bTimeDep, number time, LocalVectorTimeSeries* locTimeSeries)
{
	bool bNeedLocTimeSeries = false;

	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
		bNeedLocTimeSeries
			|= (*m_pvElemDisc)[i]->set_time_dependent(bTimeDep, time, locTimeSeries);

	return bNeedLocTimeSeries;
}


bool
DataEvaluator::
set_non_regular_grid(bool bNonRegularGrid)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//  let disc use non-regular grid assemblings
		if(!(*m_pvElemDisc)[i]->treat_non_regular_grid(bNonRegularGrid))
		{
			UG_LOG("ERROR in 'DataEvaluator::set_non_regular_grid': "
					" Elem Disc " << i << " does not support non-regular"
				   " grids, but this is requested.\n");
			return false;
		}
	}

//	check if hanging dofs are really used
	m_bUseHanging = false;
	if(bNonRegularGrid)
		for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
			m_bUseHanging |= (*m_pvElemDisc)[i]->use_hanging();

//	we're done
	return true;
}


bool
DataEvaluator::
finish_element_loop()
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
		if(!(*m_pvElemDisc)[i]->finish_element_loop())
		{
			UG_LOG("ERROR in 'DataEvaluator::finish_element_loop': "
					"Cannot finish element loop for IElemDisc "<<i<<".\n");
			return false;
		}

//	we're done
	return true;
}


} // end namespace ug

