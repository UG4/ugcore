
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
	UG_DLOG(LIB_DISC_LINKER, 2, "DataClassExportPossibility::~DataClassExportPossibility: Deleting Data Export Possibility " << this->name() << ".\n");
	for(size_t i = 0; i < m_createdDataExports.size(); ++i)
	{
		UG_ASSERT(delete_data_export(m_createdDataExports[i]), "DataClassExportPossibility::~DataClassExportPossibility: Cannot delete Exports.");
	}
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
DataExportItem*
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
create_data_export()
{
	DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
		new DataClassExport<TDataType, TPositionType, TAlgebra>(this->name(), this, m_pExportingClass, m_nrExport);
	m_createdDataExports.push_back(dynamic_cast<DataExportItem*>(exp));
	if(exp->set_num_sh(m_sys, m_num_sh) != true) return false;
	return dynamic_cast<DataExportItem*>(exp);
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
set_num_sh(size_t sys, size_t num_sh)
{
	m_sys = sys;
	m_num_sh = num_sh;
	for(size_t i = 0; i < m_createdDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
			dynamic_cast<DataClassExport<TDataType, TPositionType, TAlgebra> *>(m_createdDataExports[i]);
		if(exp->set_num_sh(sys, num_sh) != true) return false;
	}
	return true;
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool
DataClassExportPossibility<TDataType, TPositionType, TAlgebra>::
set_local_solution(const local_vector_type& u)
{
	m_u = &u;
	for(size_t i = 0; i < m_createdDataExports.size(); ++i)
	{
		DataClassExport<TDataType, TPositionType, TAlgebra> * exp =
				dynamic_cast<DataClassExport<TDataType, TPositionType, TAlgebra> *>(m_createdDataExports[i]);
		if(exp->set_local_solution(u) != true) return false;
	}
	return true;
}


template<typename TDataType, typename TPositionType, typename TAlgebra>
bool DataClassExport<TDataType, TPositionType, TAlgebra>::
set_num_sh(size_t sys, size_t num_sh)
{
	this->m_num_sys = 1;
	this->m_sys.resize(this->num_sys());
	this->m_num_sh.resize(this->num_sys());
	this->m_sys[0] = sys;
	this->m_num_sh[0] = num_sh;
	this->m_values.resize(this->num_ip());
	this->m_derivatives.resize(this->num_sys());
	for(size_t s = 0; s < this->num_sys(); ++s)
	{
		this->m_derivatives[s].resize(this->num_ip());
		for(size_t ip = 0; ip < this->m_derivatives[s].size(); ++ip)
		{
			this->m_derivatives[s][ip].resize(this->num_sh(s));
		}
	}

	// remove local values and indices, since they may be invalid
	this->m_u = NULL;
	return true;
}

template<typename TDataType, typename TPositionType, typename TAlgebra>
bool DataClassExport<TDataType, TPositionType, TAlgebra>::
set_local_solution(const local_vector_type& u)
{
	UG_ASSERT(u.size() == this->m_num_sh[0], "Wrong number of unknowns in local vector. Must match the number set in this export.");
	m_u = &u;
	return true;
}


template<typename TDataType, typename TPositionType, typename TAlgebra>
void
DataClassExport<TDataType, TPositionType, TAlgebra>::
compute(bool compute_derivatives)
{
	UG_ASSERT(this->m_values.size() == this->m_positions.size(), "Size of Data Array and Position Array must be equal. Internal error.");
//	UG_ASSERT(this->m_evalFunction != NULL, "No evaluation function set. Internal error.");
	UG_ASSERT(this->num_sys() == 1, "An Export has exactly one system it depends on");

	m_pExportingClass->template data_export<TDataType, TPositionType>(m_nrExport, this->m_values, this->m_derivatives[0], this->m_positions, *m_u, compute_derivatives);

// old
//	(m_ExportingClass->*m_evalFunction)(this->m_values, this->m_derivatives[0], this->m_positions, *m_u, compute_derivatives);
}



} // end namespace ug

#endif
