
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT_IMPL__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT_IMPL__

#include <iostream>
#include <cassert>

namespace ug{


//////////////////////////////
// Data Export Possibility
//////////////////////////////

template<typename TDataType, typename TAlgebra>
DataClassExportPossibility<TDataType, TAlgebra>::
~DataClassExportPossibility()
{
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		UG_ASSERT(delete_data_export(m_vCreatedDataExports[i]),
					"DataClassExportPossibility::~DataClassExportPossibility: Cannot delete Exports.");
	}
}

template<typename TDataType, typename TAlgebra>
DataExportItem*
DataClassExportPossibility<TDataType, TAlgebra>::
create_data_export()
{
	DataClassExport<TDataType, TAlgebra> * exp =
		new DataClassExport<TDataType, TAlgebra>(this->name(), this, m_pExportingClass, m_nrExport);
	m_vCreatedDataExports.push_back(dynamic_cast<DataExportItem*>(exp));

	// A Class export does only depend on one (its) system
	if(!exp->set_num_sys(1)) return NULL;

	// setting sys id and num_sh
	if(!exp->set_num_fct(m_numFct)) return NULL;
	for(size_t fct = 0; fct < m_numFct; ++fct)
	{
		if(!exp->set_num_dofs(fct, m_vNumDoFsPerFct[fct])) return NULL;
	}
	if(!exp->set_sys_id(m_sysId)) return NULL;

	return dynamic_cast<DataExportItem*>(exp);
}

template<typename TDataType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TAlgebra>::
set_num_fct(size_t num_fct)
{
	m_numFct = num_fct;
	m_vNumDoFsPerFct.resize(num_fct);
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TAlgebra> * exp =
			dynamic_cast<DataClassExport<TDataType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_num_fct(num_fct))
			{UG_LOG("DataClassExportPossibility::set_num_fct: Cannot set num_fcts in export.\n"); return false;}
	}
	return true;
}

template<typename TDataType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TAlgebra>::
set_num_dofs(size_t fct, size_t num_dofs)
{
	UG_ASSERT(fct < m_vNumDoFsPerFct.size(), "Invalid index. (fct = "<<fct<<", size = "<<m_vNumDoFsPerFct.size()<<")");
	m_vNumDoFsPerFct[fct] = num_dofs;
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TAlgebra> * exp =
			dynamic_cast<DataClassExport<TDataType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_num_dofs(fct, num_dofs))
			{UG_LOG("DataClassExportPossibility::set_num_dofs: Cannot set num_dofs in export.\n"); return false;}
	}
	return true;
}

template<typename TDataType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TAlgebra>::
set_sys_id(size_t sys_id)
{
	m_sysId = sys_id;
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TAlgebra> * exp =
			dynamic_cast<DataClassExport<TDataType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_sys_id(sys_id))
			{UG_LOG("DataClassExportPossibility::set_sys_id: Cannot set sys id in export.\n"); return false;}
	}
	return true;
}

template<typename TDataType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TAlgebra>::
set_local_solution(const local_vector_type& u)
{
	m_pSolution = &u;
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TAlgebra> * exp =
				dynamic_cast<DataClassExport<TDataType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_local_solution(u))
			{UG_LOG("DataClassExportPossibility::set_local_solution: Cannot set solution in export.\n"); return false;}
	}
	return true;
}

//////////////////////////////
// Data Class Export
//////////////////////////////


template<typename TDataType, typename TAlgebra>
bool
DataClassExport<TDataType, TAlgebra>::
set_num_fct(size_t num_fct)
{
	if(!DataExport<TDataType>::set_num_fct(num_fct, 0))
		{UG_LOG("DataClassExport::set_num_fct: Error while setting num_fct.\n"); return false;}
	return true;
}

template<typename TDataType, typename TAlgebra>
bool
DataClassExport<TDataType, TAlgebra>::
set_num_dofs(size_t fct, size_t num_dofs)
{
	if(!DataExport<TDataType>::set_num_dofs(fct, num_dofs, 0))
		{UG_LOG("DataClassExport::set_num_dofs: Error while setting num_dofs.\n"); return false;}
	return true;
}


template<typename TDataType, typename TAlgebra>
bool DataClassExport<TDataType, TAlgebra>::
set_sys_id(size_t sys_id)
{
	if(!DataExport<TDataType>::set_sys_id(sys_id, 0)) return false;
	this->m_vSysId[0] = sys_id;
	return true;
}

template<typename TDataType, typename TAlgebra>
bool DataClassExport<TDataType, TAlgebra>::
set_local_solution(const local_vector_type& u)
{
	UG_ASSERT(u.num_dofs() == this->num_dofs(),
			"Wrong number of unknowns in local vector. Must match the number set in this export.");
	// casting away constness, such that we can set the subfunctionmap as we like
	m_pSolution = const_cast<local_vector_type*>(&u);
	m_pSubFunctionMap = &(m_pSolution->get_indices().get_sub_function_map());
	return true;
}


template<typename TDataType, typename TAlgebra>
void
DataClassExport<TDataType, TAlgebra>::
compute(bool compute_derivatives)
{
	UG_ASSERT(this->m_vValue.size() == this->num_ip(),
				"Size of Data Array and Position Array must be equal. Internal error.");
	UG_ASSERT(this->num_sys() == 1,
			"An Export has exactly one system it depends on");

	m_pSolution->get_indices().access_sub_function_group(*m_pSubFunctionMap);

	switch(this->m_posDim)
	{
	case 1:
		m_pExportingClass->template data_export<TDataType>(m_nrExport, this->m_vValue, this->m_vvDerivatives[0],
															this->template get_positions<1>(), *m_pSolution, compute_derivatives);
		break;
	case 2:
		m_pExportingClass->template data_export<TDataType>(m_nrExport, this->m_vValue, this->m_vvDerivatives[0],
															this->template get_positions<2>(), *m_pSolution, compute_derivatives);
		break;
	case 3:
		m_pExportingClass->template data_export<TDataType>(m_nrExport, this->m_vValue, this->m_vvDerivatives[0],
															this->template get_positions<3>(), *m_pSolution, compute_derivatives);
		break;
	default: UG_LOG("Dimension " << this->m_posDim << " not supported.\n"); UG_ASSERT(0, "Dimension " << this->m_posDim << " not supported."); exit(1);
	}
}



} // end namespace ug

#endif
