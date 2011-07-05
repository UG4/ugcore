/*
 * data_evaluator_impl.h
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "data_evaluator.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// DataEvaluator Setup
///////////////////////////////////////////////////////////////////////////////

bool DataEvaluator::set_elem_discs(const std::vector<IElemDisc*>& vElemDisc,
                                   const FunctionPattern& fctPat,
                                   bool bNonRegularGrid,
                                   bool bMassPart)
{
//	remember current elem discs
	m_pvElemDisc = &vElemDisc;

//	set function pattern to common function group
	m_commonFctGroup.set_function_pattern(fctPat);
	m_commonFctGroup.add_all();

//	create FunctionIndexMapping for each Disc
	m_vElemDiscMap.resize(m_pvElemDisc->size());
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

	//	copy function group in import/export of element discs
		for(size_t imp = 0; imp < (*m_pvElemDisc)[i]->num_imports(); ++imp)
			(*m_pvElemDisc)[i]->get_import(imp).set_function_group(elemFctGrp);

		for(size_t exp = 0; exp < (*m_pvElemDisc)[i]->num_exports(); ++exp)
			(*m_pvElemDisc)[i]->get_export(exp).set_function_group(elemFctGrp);

	//	create a mapping between all functions and the function group of this
	//	element disc.
		if(!CreateFunctionIndexMapping(m_vElemDiscMap[i], elemFctGrp,
		                               	   	   	   	   	   	   m_commonFctGroup))
		{
			UG_LOG("ERROR in 'DataEvaluator::set_elem_discs': Cannot create "
					"Function Index Mapping for disc "<<i<<".\n");
			return false;
		}
	}

//	setup for non-regular grid
	if(!set_non_regular_grid(bNonRegularGrid))
	{
		UG_LOG("ERROR in 'DataEvaluator::set_elem_discs':"
				"Cannot setup non-regular grid.\n");
		return false;
	}

	return extract_imports_and_ipdata(bMassPart);
}

void
DataEvaluator::
clear_extracted_data_and_mappings()
{
	m_vMassDataImport.clear();
	m_vStiffDataImport.clear();
	m_vMassImpMap.clear();
	m_vStiffImpMap.clear();
	m_vMassImpConnMap.clear();
	m_vStiffImpConnMap.clear();

	m_vConstData.clear();
	m_vPosData.clear();

	m_vDependentIPData.clear();
	m_vDependentMap.clear();

	m_vDataExport.clear();
	m_vExpMap.clear();

	m_vLinkerMap.clear();
	m_vDataLinker.clear();
}

bool DataEvaluator::add_data_to_eval_data(std::vector<IIPData*>& vEvalData,
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

bool DataEvaluator::extract_imports_and_ipdata(bool bMassPart)
{
//	check that elem disc given
	if(m_pvElemDisc == NULL)
	{
		UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata': No vector"
				" of IElemDisc* set. Cannot extract imports and exports.\n");
		return false;
	}

//	clear imports and ipdata
	clear_extracted_data_and_mappings();

//	queue for all ip data needed
	std::vector<IIPData*> vEvalData;
	std::vector<IIPData*> vTryingToAdd;

//	In the next loop we extract all need IPData:
//	We only process the DataImport if there has been set data to the import
//	since otherwise no evaluation is needed.
//	If there is data given, we get the connected IPData and add it to the vector
//	of EvaluationData. This simply adds the IPData to the queue for IPData, if
//	the data does not depend on other Data. But if the IPData itself has
//	dependencies to other IPData, this data is added first (in a recursive
//	process). Of coarse, no circle dependency between IPData is allowed.

//	loop elem discs
	for(size_t d = 0; d < m_pvElemDisc->size(); ++d)
	{
	//	loop imports
		for(size_t i = 0; i < (*m_pvElemDisc)[d]->num_imports(); ++i)
		{
		//	get import
			IDataImport* iimp = &((*m_pvElemDisc)[d]->get_import(i));

		//	skip, if in mass part but no mass part wanted
			if(!bMassPart && iimp->in_mass_part()) continue;

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
		//	check that queue is empty now, else some internal error occured
			if(!vTryingToAdd.empty())
			{
				UG_LOG("ERROR in DataEvaluator::extract_imports_and_ipdata:"
						" Internal Error, IPData queue not empty after adding.\n");
				return false;
			}
		}
	}

//	Now, we have processed all imports, that must be evaluated and have a long
//	vector of IPData that is connected to those imports. The IPData is already
//	sorted in this way: Data that depends on other data appears after the data
//	it depends on. This is important since we will schedule now the data for
//	evaluation and the data, that is needed by other data, will be computed
//	first. In addition, the data linker have to update their FunctionGroup and
//	must be sure that the data they depend on has already a correct FunctionGroup
//	set. This all is ensured by the (already produced) correct ordering.
//
//	In the next loop we process all IPData, that will be evaluated during
//	assembling (i.e. is connected to an Import). First, we check if the data
//	is constant. If so simply add it the the Constant Data vector; nothing more
//	has to be done here. Else we check if the data depends on the primary
//	unknowns. If this is not the case, the IPData must be a Position-dependent
//	data, but not constant. Thus, schedule it at the Position Data vector.
//	If the data depends on the primary unknowns we must proceed as follows:
//	First, we update the FunctionGroup of the Data, since it could be a linker
//	and having an incorrect FunctionGroup (iff the FunctionGroup of the data
//	the linker depends on has been changed). Then we create the function
//	mapping between the functions the linker depends on and the common Function
//	Group.

//	loop all needed ip data and group it
	for(size_t i = 0; i < vEvalData.size(); ++i)
	{
	//	get the ip data
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

	//	update function group of dependent data
		if(!dependData->update_function_group())
		{
			UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
					" Cannot update FunctinoGroup of IDependentData.\n");
			return false;
		}

	//	create FuncMap
		FunctionIndexMapping map;
		if(!CreateFunctionIndexMapping(map,
		                               dependData->get_function_group(),
									   m_commonFctGroup))
		{
			UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
					"Cannot create Function Index Mapping for IDependData.\n");
			return false;
		}

	//	now we have to remember the mapping and schedule the dependent ipdata for
	//	evaluation at the correct queue

	//	save as dependent data
		m_vDependentIPData.push_back(dependData);
		m_vDependentMap.push_back(map);

	//	cast to data export
		IDataExport* exp = dynamic_cast<IDataExport*>(ipData);

	//	Data Export case
		if(exp != NULL)
		{
		//	schedule for evaluation of IDataExports
			m_vDataExport.push_back(exp);

		// 	remember function map
			m_vExpMap.push_back(map);
		}
		else
	//	Linker case
		{
		//	schedule for evaluation of linker
			m_vDataLinker.push_back(dependData);

		// 	remember function map
			m_vLinkerMap.push_back(map);
		}
	}

//	In a second loop over the data imports, we schedule the DataImports for
//	evaluation and compute the correct FunctionMapping for the linearization
//	of the defect and the Data, the Import is connected to:
//	If the IPData does not depend on the primary unknowns, we're done. Else
//	we have to setup the Function mappings between the common function group
//	and the DataImport-FunctionGroup. This is simply the same function map as
//	for the element discretization, since the DataImport depends by definition
//	from and only from the primary variables of its associated IElemDisc.

//	loop elem discs
	for(size_t d = 0; d < m_pvElemDisc->size(); ++d)
	{
	//	loop imports
		for(size_t i = 0; i < (*m_pvElemDisc)[d]->num_imports(); ++i)
		{
		//	get import
			IDataImport* iimp = &((*m_pvElemDisc)[d]->get_import(i));

		//	skip, if in mass part but no mass part wanted
			if(!bMassPart && iimp->in_mass_part()) continue;

		//	skip non-given data (no need for evaluation)
			if(!iimp->data_given()) continue;

		//	done if and only if zero-derivative
			if(iimp->zero_derivative()) continue;

		//	get and cast dependent data
			IDependentIPData* dependData =
								dynamic_cast<IDependentIPData*>(iimp->get_data());

		//	check success
			if(dependData == NULL)
			{
				UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
						" Data seems dependent, but cast failed.\n");
				return false;
			}

		//	create FuncMap for data
		//	this is ok, since the function group has been updated in the
		//	previous loop over all needed data
			FunctionIndexMapping map;
			if(!CreateFunctionIndexMapping(map,
										   dependData->get_function_group(),
										   m_commonFctGroup))
			{
				UG_LOG("ERROR in 'DataEvaluator::extract_imports_and_ipdata':"
						"Cannot create Function Index Mapping for DependentData.\n");
				return false;
			}

		//	check if data is located in mass part or stiffness part
			if(iimp->in_mass_part())
			{
			//	schedule import for computation of lin defect:
				m_vMassDataImport.push_back(iimp);

			//	Remember FuncMap for lin defect (the same map as for elem disc)
				m_vMassImpMap.push_back(m_vElemDiscMap[d]);

			//  Remember FuncMap for connected data
				m_vMassImpConnMap.push_back(map);
			}
			else
			{
			//	schedule import for computation of lin defect:
				m_vStiffDataImport.push_back(iimp);

			//	Remember FuncMap for lin defect (the same map as for elem disc)
				m_vStiffImpMap.push_back(m_vElemDiscMap[d]);

			//  Remember FuncMap for connected data
				m_vMassImpConnMap.push_back(map);
			}
		}
	}

//	we're done
	return true;
}


bool DataEvaluator::set_time_dependent(bool bTimeDep, number time,
                                       LocalVectorTimeSeries* locTimeSeries)
{
	bool bNeedLocTimeSeries = false;

	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
		bNeedLocTimeSeries
			|= (*m_pvElemDisc)[i]->set_time_dependent(bTimeDep, time, locTimeSeries);

	return bNeedLocTimeSeries;
}


bool DataEvaluator::set_non_regular_grid(bool bNonRegularGrid)
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

///////////////////////////////////////////////////////////////////////////////
// Assemble routines
///////////////////////////////////////////////////////////////////////////////

bool DataEvaluator::compute_elem_data(local_vector_type & u, bool bDeriv)
{
//	evaluate position data
	for(size_t i = 0; i < m_vPosData.size(); ++i)
		m_vPosData[i]->compute();

// 	process dependent data:
//	We can not simply compute exports first, then Linker, because an export
//	itself could depend on other data if implemented somehow in the IElemDisc
//	(e.g. using data from some DataImport). Thus, we have to loop the soretd
//	vector of all dependent data (that is correctly sorted the way that always
//	needed data has previously computed). We look up, if a Export is given, if
//	so compute it, else compute the linker

//	loop all dependent data
	for(size_t i = 0; i < m_vDependentIPData.size(); ++i)
	{
	//	check if current solution is needed
		if(m_vDependentIPData[i]->comp_needs_sol())
		{
			if(!m_vDependentIPData[i]->compute(u, bDeriv))
			{
				UG_LOG("ERROR in 'DataEvaluator::compute_elem_data':"
						"Cannot compute data for Export " << i <<".\n");
				return false;
			}
		}
		else
			m_vDependentIPData[i]->compute(bDeriv);
	}

//	we're done
	return true;
}

bool DataEvaluator::ass_JA_elem(local_matrix_type& A, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		A.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->ass_JA_elem(A, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_JA': "
					"Cannot assemble JA for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::ass_JM_elem(local_matrix_type& M, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		M.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->ass_JM_elem(M, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_JM': "
					"Cannot assemble JM for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::ass_dA_elem(local_vector_type& d, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		d.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->ass_dA_elem(d, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_A': "
					"Cannot assemble Defect (A) for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::ass_dM_elem(local_vector_type& d, local_vector_type& u)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		u.access_by_map(map(i));
		d.access_by_map(map(i));

	//	assemble JA
		if(!(*m_pvElemDisc)[i]->ass_dM_elem(d, u))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_M': "
					"Cannot assemble Defect (M) for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::ass_rhs_elem(local_vector_type& rhs)
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
	{
	//	access disc functions
		rhs.access_by_map(map(i));

	//	assemble rhs
		if(!(*m_pvElemDisc)[i]->ass_rhs_elem(rhs))
		{
			UG_LOG("ERROR in 'DataEvaluator::assemble_rhs': "
					"Cannot assemble rhs for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::finish_elem_loop()
{
	for(size_t i = 0; i < (*m_pvElemDisc).size(); ++i)
		if(!(*m_pvElemDisc)[i]->finish_elem_loop())
		{
			UG_LOG("ERROR in 'DataEvaluator::finish_element_loop': "
					"Cannot finish element loop for IElemDisc "<<i<<".\n");
			return false;
		}

//	we're done
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Coupling
///////////////////////////////////////////////////////////////////////////////

bool DataEvaluator::compute_lin_defect_JA(local_vector_type & u)
{
//	compute linearized defect
	for(size_t i = 0; i < m_vStiffDataImport.size(); ++i)
	{
	//	set correct access for import
		u.access_by_map(m_vStiffImpMap[i]);

	//	compute linearization of defect
		if(!m_vStiffDataImport[i]->compute_lin_defect(u))
		{
			UG_LOG("ERROR in 'DataEvaluator::compute_elem_data': Cannot compute"
					" linearized defect for Import " << i <<" (Stiffness part).\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::compute_lin_defect_JM(local_vector_type & u)
{
//	compute linearized defect
	for(size_t i = 0; i < m_vMassDataImport.size(); ++i)
	{
	//	set correct access for import
		u.access_by_map(m_vMassImpMap[i]);

	//	compute linearization of defect
		if(!m_vMassDataImport[i]->compute_lin_defect(u))
		{
			UG_LOG("ERROR in 'DataEvaluator::compute_elem_data': Cannot compute"
					" linearized defect for Import " << i <<" (Mass part).\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool DataEvaluator::add_coupl_JA(local_matrix_type& J)
{
//	loop all imports located in the stiffness part
	for(size_t i = 0; i < m_vStiffDataImport.size(); ++i)
	{
	//	rows are given by import, cols are given by connected data
		J.access_by_map(m_vStiffImpMap[i], m_vStiffImpConnMap[i]);

	//	add off diagonal coupling
		m_vStiffDataImport[i]->assemble_jacobian(J);
	}

//	we're done
	return true;
}

bool DataEvaluator::add_coupl_JM(local_matrix_type& J)
{
//	loop all imports located in the mass part
	for(size_t i = 0; i < m_vMassDataImport.size(); ++i)
	{
	//	rows are given by import, cols are given by connected data
		J.access_by_map(m_vMassImpMap[i], m_vMassImpConnMap[i]);

	//	add off diagonal coupling
		m_vMassDataImport[i]->assemble_jacobian(J);
	}

//	we're done
	return true;
}

} // end namespace ug

