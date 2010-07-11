
#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT_IMPL__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_CLASS_EXPORT_IMPL__

#include <iostream>
#include <cassert>

namespace ug{


//////////////////////////////
// Data Export Possibility
//////////////////////////////

template<typename TDataType, typename TPositionType, typename TAlgebra>
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
~DataClassExportPossibility()
{
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		UG_ASSERT(delete_data_export(m_vCreatedDataExports[i]),
					"DataClassExportPossibility::~DataClassExportPossibility: Cannot delete Exports.");
	}
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
DataExportItem*
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
create_data_export()
{
	DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
		new DataClassExport<TDataType, TPositionType, TAlgebra>(this->name(), this, m_pExportingClass, m_nrExport);
	m_vCreatedDataExports.push_back(dynamic_cast<DataExportItem*>(exp));

	// A Class export does only depend on one (its) system
	if(!exp->set_num_sys(1)) return NULL;

	// setting sys id and num_sh
	if(!exp->set_num_sh(m_numSh)) return NULL;
	if(!exp->set_sys_id(m_sysId)) return NULL;

	return dynamic_cast<DataExportItem*>(exp);
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
set_num_sh(size_t num_sh)
{
	m_numSh = num_sh;
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
			dynamic_cast<DataClassExport<TDataType, TPositionType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_num_sh(num_sh)) return false;
	}
	return true;
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
set_sys_id(size_t sys_id)
{
	m_sysId = sys_id;
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
			dynamic_cast<DataClassExport<TDataType, TPositionType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_sys_id(sys_id)) return false;
	}
	return true;
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
set_local_solution(const local_vector_type& u)
{
	m_pSolution = &u;
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
				dynamic_cast<DataClassExport<TDataType, TPositionType, TAlgebra> *>(m_vCreatedDataExports[i]);
		if(!exp->set_local_solution(u)) return false;
	}
	return true;
}

//////////////////////////////
// Data Class Export
//////////////////////////////

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool DataClassExport<TDataType, TPositionType, TAlgebra>::
set_num_sh(size_t num_sh)
{
	if(!DataExport<TDataType, TPositionType>::set_num_sh(num_sh, 0)) return false;
	this->m_vNumSh[0] = num_sh;
	return true;
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool DataClassExport<TDataType, TPositionType, TAlgebra>::
set_sys_id(size_t sys_id)
{
	if(!DataExport<TDataType, TPositionType>::set_sys_id(sys_id, 0)) return false;
	this->m_vSysId[0] = sys_id;
	return true;
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool DataClassExport<TDataType, TPositionType, TAlgebra>::
set_local_solution(const local_vector_type& u)
{
	UG_ASSERT(u.size() == this->m_vNumSh[0],
			"Wrong number of unknowns in local vector. Must match the number set in this export.");
	m_pSolution = &u;
	return true;
}


template<typename TDataType, typename TPositionType, typename TAlgebra>
void
DataClassExport<TDataType, TPositionType, TAlgebra>::
compute(bool compute_derivatives)
{
	UG_ASSERT(this->m_vValue.size() == this->m_vPosition.size(),
				"Size of Data Array and Position Array must be equal. Internal error.");
	UG_ASSERT(this->num_sys() == 1,
			"An Export has exactly one system it depends on");

	m_pExportingClass->template data_export<TDataType, TPositionType>(m_nrExport, this->m_vValue, this->m_vvvDerivatives[0],
																		this->m_vPosition, *m_pSolution, compute_derivatives);

}



} // end namespace ug

#endif
