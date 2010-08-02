/*
 * elementdata_impl.h
 *
 *  Created on: 09.11.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT_IMPL__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT_IMPL__

#include <iostream>
#include <cassert>
#include <algorithm>

namespace ug{


//////////////////////////////
// Data Export
//////////////////////////////

template<typename TDataType>
DataExport<TDataType>::
~DataExport()
{
	// unregister all registered imports
	typename std::vector<DataImportItem*>::reverse_iterator rit;
	for(rit = m_vImportList.rbegin(); rit != m_vImportList.rend(); ++rit)
	{
		remove_data_import(*rit);
	}

	if(this->get_possibility_item() != NULL)
		this->get_possibility_item()->delete_data_export(dynamic_cast<DataExportItem*>(this));
};

template<typename TDataType>
bool DataExport<TDataType>::
add_data_import(DataImportItem* importItem)
{
	DataImport<TDataType>* Import = dynamic_cast<DataImport<TDataType>*>(importItem);

	if(Import == NULL)
	{
		UG_LOG("DataExport::add_data_import: Data type and/or position type of Export<->Import does not match. Cannot register.\n");
		return false;
	}

	// check, that import is not already linked to other export
	if(Import->is_linked())
	{
		UG_LOG("DataExport::add_data_import: DataImport already registered to an DataExport. Invalid operation.\n");
		return false;
	}

	// TODO: is this necessary, since already known, that Import is not linked?
	// Find requested Import in Import list
	typename std::vector<DataImportItem*>::iterator importIter;
	importIter = find(m_vImportList.begin(), m_vImportList.end(), importItem);
	if(importIter != m_vImportList.end())
	{
		UG_LOG("DataExport::add_data_import: DataImport already registered at this DataExport. Invalid operation.\n");
		return false;
	}

	// add import
	m_vImportList.push_back(importItem);

	// set positions if already positions set
	if(m_posDim > 0)
	{
		switch(m_posDim)
		{
		case 1: set_positions<1>(Import->get_positions<1>()); break;
		case 2: set_positions<2>(Import->get_positions<2>()); break;
		case 3: set_positions<3>(Import->get_positions<3>()); break;
		default: UG_LOG("Dimension " << m_posDim << " not supported.\n"); return false;
		}
	}

	// register export at import
	Import->m_pTypeExport = this;

	// register export base at import base
	Import->m_pExport = dynamic_cast<DataExportItem*>(this);

	return true;
}

template<typename TDataType>
bool DataExport<TDataType>::
remove_data_import(DataImportItem* importItem)
{
	DataImport<TDataType>* Import = dynamic_cast<DataImport<TDataType>*>(importItem);

	if(Import == NULL)
	{
		UG_LOG("DataExport::add_data_import: Data type and/or position type of Export<->Import does not match. Cannot register.\n");
		return false;
	}

	typename std::vector<DataImportItem*>::iterator importIter;
	importIter = find(m_vImportList.begin(), m_vImportList.end(), importItem);
	if(importIter == m_vImportList.end())
	{
		UG_LOG("DataExport::remove_data_import: DataImport not found. Cannot unregister.\n");
		return false;
	}

	Import->m_pExport = NULL;
	Import->m_pTypeExport = NULL;

	// erase import
	m_vImportList.erase(importIter);

	// if import list is void, delete positions and values
	if(m_vImportList.empty())
	{
		switch(m_posDim)
		{
		case 1: get_positions<1>().clear();
		case 2: get_positions<1>().clear();
		case 3: get_positions<1>().clear();
		default: UG_LOG("Dimension "<< this->m_posDim << "not supported.\n"); return false;
		}
		m_vValue.clear();
		for(size_t s = 0; s < m_vvvDerivatives.size(); ++s)
		{
			for(size_t ip = 0; ip < m_vvvDerivatives[s].size(); ++ip)
			{
				m_vvvDerivatives[s][ip].clear();
			}
			m_vvvDerivatives[s].clear();
		}
		m_vvvDerivatives.clear();
	}

	return true;
}


template <typename TDataType>
template <int dim>
bool DataExport<TDataType>::
set_positions(const std::vector<MathVector<dim> >& positions, bool overwrite)
{
	// get positions
	std::vector<MathVector<dim> >& vPosition = get_positions<dim>();

	// if no positions set
	if(vPosition.empty() || overwrite)
	{
		// copy positions and set same size for value list
		vPosition = positions;

		m_posDim = dim;
		m_numIp = vPosition.size();

		m_vValue.resize(num_ip());

		m_vvvDerivatives.resize(this->num_sys());
		for(size_t s = 0; s < this->num_sys(); ++s)
		{
			m_vvvDerivatives[s].resize(num_ip());
			for(size_t ip = 0; ip < m_vvvDerivatives[s].size(); ++ip)
			{
				m_vvvDerivatives[s][ip].resize(num_sh(s));
			}
		}

		//TODO : implement linker
/*		if(!set_linker_positions(vPosition))
			{UG_LOG("DataExport::set_positions: Error while setting positions in linker imports.\n"); return false;}
*/		return true;
	}

	// if already positions set, check, if positions are equal
	if(positions.size() != num_ip())
		{UG_LOG("DataExport::set_positions: Different number of positions as already set positions.\n");return false;}

	for(size_t ip = 0; ip < num_ip(); ++ip)
	{
		if(vPosition[ip] != positions[ip])
			{UG_LOG("DataExport::set_positions: Position " << ip << " is different from already set positions.\n"); return false;}
	}
	return true;
}

template<typename TDataType>
bool
DataExport<TDataType>::
equal(const DataExportItem& v) const
{
	const DataExport<data_type>* cast_v = dynamic_cast<const DataExport<data_type>*>(&v);

	// if type is not the same, they are nor equal
	if(cast_v == NULL) return false;

	// Compute if pointers to Base possibility are the same
	if(this->m_pPossibilityItem != cast_v->get_possibility_item()) return false;

	// number of positions have to match
	if(num_ip() != cast_v->num_ip()) return false;

	// check each position
	for(size_t ip = 0; ip < num_ip(); ++ip)
	{
		switch(m_posDim)
		{
		case 1:	if(position<1>(ip) != cast_v->position<1>(ip)) return false;
		case 2:	if(position<2>(ip) != cast_v->position<2>(ip)) return false;
		case 3:	if(position<3>(ip) != cast_v->position<3>(ip)) return false;
		default: UG_LOG("Dimension " << m_posDim << " not supported.\n"); UG_ASSERT(0, "Dimension " << m_posDim << " not supported."); return false;
		}
	}

	return true;
}

template<typename TDataType>
bool
DataExport<TDataType>::
print_positions() const
{
	if(m_posDim == 1)
	{
		std::vector<MathVector<1> >& vPosition1 = const_cast<DataExport<TDataType>*>(this)->get_positions<1>();
		for(size_t ip = 0; ip < num_ip(); ++ip)
			{UG_LOG(vPosition1[ip]); if(ip != num_ip()-1) UG_LOG(", ");}
		return true;
	}
	else if(m_posDim == 2)
	{
		std::vector<MathVector<2> >& vPosition2 = const_cast<DataExport<TDataType>*>(this)->get_positions<2>();
		for(size_t ip = 0; ip < num_ip(); ++ip)
			{UG_LOG(vPosition2[ip]); if(ip != num_ip()-1) UG_LOG(", ");}
		return true;
	}
	else if(m_posDim == 3)
	{
		std::vector<MathVector<3> >& vPosition3 = const_cast<DataExport<TDataType>*>(this)->get_positions<3>();
		for(size_t ip = 0; ip < num_ip(); ++ip)
			{UG_LOG(vPosition3[ip]); if(ip != num_ip()-1) UG_LOG(", ");}
		return true;
	}

	UG_LOG("Dimension  " << m_posDim << " not supported.\n");
	return false;
}

template<typename TDataType>
bool
DataExport<TDataType>::
print_values() const
{
	for(size_t ip = 0; ip < num_ip(); ++ip)
		{UG_LOG(operator[](ip));if(ip != num_ip()-1) UG_LOG(", ");}
	return true;
}

template<typename TDataType>
bool
DataExport<TDataType>::
print_derivatives(std::string offset) const
{
	UG_LOG("Depending on "<< this->num_sys()<< " systems.\n");
	for(size_t s = 0; s < this->num_sys(); ++s)
	{
		UG_LOG(offset << "w.r.t. sys = " << sys_id(s) << ": ");

		for(size_t k = 0; k < num_sh(s); ++k)
		{
			if(k==0) {UG_LOG(" k="<< k<< ": [ ");}
			else {UG_LOG(offset << "                 k="<< k<< ": [ ");};

			for(size_t ip = 0; ip < num_ip(); ++ip)
			{
				{
					UG_LOG(operator()(s,ip,k));
					if(ip != num_ip() - 1) UG_LOG(", ");
				}
			}
			UG_LOG(" ]\n");
		}
	}
	return true;
}

template<typename TDataType>
bool
DataExport<TDataType>::
print_info(std::string offset) const
{
	UG_LOG(offset << "IPs: [ "); if(!print_positions()) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Values: [ "); if(!print_values()) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Derivatives: "); if(!print_derivatives(offset)) return false;
	return true;
}


//////////////////////////////
// Data Import
//////////////////////////////

template <typename TDataType>
template <int dim>
bool
DataImport<TDataType>::
set_positions(const std::vector<MathVector<dim> >& pos, bool overwrite)
{
	// get positions
	std::vector<MathVector<dim> >& vPosition = get_positions<dim>();

	// copy new positions (elementwise)
	vPosition = pos;
	m_posDim = dim;
	m_numIp = vPosition.size();

	// adjust lin defect array for ip's
	for(size_t i = 0; i < m_vvLinearizedDefect.size(); ++i)
	{
		m_vvLinearizedDefect[i].resize(num_ip());
	}

	// if no export set, do nothing
	if(this->m_pExport == NULL) return true;

	// if export set is an DataExport
	typename ug::DataExport<TDataType>* Cast_Export =
			dynamic_cast< typename ug::DataExport<TDataType>* >(this->m_pExport);
	if(Cast_Export != NULL)
	{
		if(!Cast_Export->set_positions<dim>(vPosition, overwrite))
			{UG_LOG("DataImport::set_positions: Cannot set positions in export.\n"); return false;}
	}
	else
	{
		UG_LOG("Cast failed. By construction this should be impossible. Bug.\n"); return false;
	}

	// if export set is an DataLinker
	// TODO: implement

	return true;
}


template<typename TDataType>
bool
DataImport<TDataType>::
set_num_eq(size_t num_eq)
{
	m_numEq = num_eq;

	m_vvLinearizedDefect.resize(m_numEq);

	for(size_t i = 0; i < m_vvLinearizedDefect.size(); ++i)
	{
		m_vvLinearizedDefect[i].resize(num_ip());
	}

	return true;
}

template<typename TDataType>
bool
DataImport<TDataType>::
link_data_export(DataExportItem* exportItem)
{
	typename ug::DataExport<TDataType>* Export =
			dynamic_cast< typename ug::DataExport<TDataType>* >(exportItem);
	if(Export == NULL)
	{
		UG_LOG("DataImport::link_data_export: Data type and/or position type of Export<->Import does not match. Cannot register.\n");
		return false;
	}

	if(this->m_pExport != NULL)
	{
		UG_LOG("DataImport::link_data_export: Already an Export registered.\n");
		return false;
	}

	// register at export
	if(Export->add_data_import(this) != true)
	{
		UG_LOG("DataImport::link_data_export: Can not register at export.\n");
		this->m_pExport = NULL;
		return false;
	}

	// remember Export
	m_pTypeExport = Export;

	// remember export base
	this->m_pExport = exportItem;

	return true;
}


template<typename TDataType>
bool
DataImport<TDataType>::
clear_data_export()
{
	// if registered at export
	if(this->m_pExport != NULL)
	{
		if(!m_pTypeExport->remove_data_import(this))
		{UG_LOG("DataImport::clear_data_export: Can not unregister from export."); return false;}

		// reset export
		m_pTypeExport = NULL;
		this->m_pExport = NULL;
	}
	// if not registered: do nothing

	return true;
}


template<typename TDataType>
DataImport<TDataType>::
~DataImport()
{
	clear_data_export();
}

template<typename TDataType>
bool
DataImport<TDataType>::
print_positions() const
{
	if(m_posDim == 1)
	{
		std::vector<MathVector<1> >& vPosition1 = const_cast<DataImport<TDataType>*>(this)->get_positions<1>();
		for(size_t ip = 0; ip < num_ip(); ++ip)
			{UG_LOG(vPosition1[ip]); if(ip != num_ip()-1) UG_LOG(", ");}
		return true;
	}
	else if(m_posDim == 2)
	{
		std::vector<MathVector<2> >& vPosition2 = const_cast<DataImport<TDataType>*>(this)->get_positions<2>();
		for(size_t ip = 0; ip < num_ip(); ++ip)
			{UG_LOG(vPosition2[ip]); if(ip != num_ip()-1) UG_LOG(", ");}
		return true;
	}
	else if(m_posDim == 3)
	{
		std::vector<MathVector<3> >& vPosition3 = const_cast<DataImport<TDataType>*>(this)->get_positions<3>();
		for(size_t ip = 0; ip < num_ip(); ++ip)
			{UG_LOG(vPosition3[ip]); if(ip != num_ip()-1) UG_LOG(", ");}
		return true;
	}

	UG_LOG("Dimension " << m_posDim << " not supported.\n");
	return false;
}

template<typename TDataType>
bool
DataImport<TDataType>::
print_values() const
{
	if(!this->is_linked())
		{UG_LOG("Not linked."); return true;}

	for(size_t ip = 0; ip < num_ip(); ++ip)
		{UG_LOG(this->operator[](ip)); if(ip != num_ip()-1) UG_LOG(", ");}

	return true;
}

template<typename TDataType>
bool
DataImport<TDataType>::
print_derivatives(std::string offset) const
{
	if(!this->is_linked())
	{
		UG_LOG("Not linked.\n");
		return true;
	}

	UG_LOG("Depending on "<< this->num_sys()<< " systems.\n");
	for(size_t s = 0; s < this->num_sys(); ++s)
	{
		for(size_t k = 0; k < this->num_sh(s); ++k)
		{
			UG_LOG(offset << "w.r.t.: sys = " << this->sys_id(s) <<", k="<< k<< ": [ ");

			for(size_t ip = 0; ip < num_ip(); ++ip)
			{
				{
					UG_LOG(this->operator()(s, ip, k));
					if(ip != num_ip() - 1) UG_LOG(", ");
				}
			}
			UG_LOG(" ]\n");
		}
	}
	return true;
}

template<typename TDataType>
bool
DataImport<TDataType>::
print_info(std::string offset) const
{
	UG_LOG(offset << "IPs: [ "); if(!print_positions()) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Values: [ "); if(!print_values()) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Derivatives: "); if(!print_derivatives(offset)) return false;
	return true;
}

template <typename TDataType, typename TAlgebra>
std::ostream& operator<<(std::ostream& out, const DataImport<TDataType>& data)
{
	for(size_t ip = 0; ip < data.num_ip(); ++ip)
	{
		out << data[ip];
		if(ip != data.num_ip()-1) out << ", ";
	}
	return out;
}

} // End namespace ug

#endif /*__H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT_IMPL__*/
