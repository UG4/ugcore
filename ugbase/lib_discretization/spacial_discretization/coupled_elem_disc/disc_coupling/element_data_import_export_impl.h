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

namespace ug{


//////////////////////////////
// Data Export
//////////////////////////////

template<typename TDataType, typename TPositionType>
DataExport<TDataType, TPositionType>::
~DataExport()
{
	UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::~DataExport: Deleting Data Export " << this->name() << ".\n");
	// unregister all registered imports
	typename std::vector<DataImportItem*>::reverse_iterator rit;
	for(rit = m_importList.rbegin(); rit != m_importList.rend(); ++rit)
	{
		UG_DLOG(LIB_DISC_LINKER, 2, "Deleting Import....");
		remove_data_import(*rit);
		UG_DLOG(LIB_DISC_LINKER, 2, "done.\n");
	}
	UG_ASSERT(m_importList.empty(), "Import list is not empty, but DataExport will be destroyed. This results in undefined linked DataImports.");

	if(this->get_possibility_item() != NULL)
		this->get_possibility_item()->delete_data_export(dynamic_cast<DataExportItem*>(this));
};

template<typename TDataType, typename TPositionType>
bool DataExport<TDataType, TPositionType>::
add_data_import(DataImportItem* importItem)
{
	DataImport<TDataType, TPositionType>* Import = dynamic_cast<DataImport<TDataType, TPositionType>*>(importItem);

	if(Import == NULL)
	{
		UG_LOG("DataExport::add_data_import: Data type and/or position type of Export<->Import does not match. Cannot register.\n");
		return false;
	}

	// check, that import is not already linked to other export
	if(Import->is_linked() == true)
	{
		UG_LOG("DataExport::add_data_import: DataImport already registered to an DataExport. Invalid operation.\n");
		return false;
	}

	// TODO: is this necessary, since already known, that Import is not linked?
	// Find requested Import in Import list
	typename std::vector<DataImportItem*>::iterator importIter;
	importIter = find(m_importList.begin(), m_importList.end(), importItem);
	if(importIter != m_importList.end())
	{
		UG_LOG("DataExport::add_data_import: DataImport already registered at this DataExport. Invalid operation.\n");
		return false;
	}

	// add import
	m_importList.push_back(importItem);

	// set positions
	set_positions(Import->get_positions());

	// register export at import
	Import->m_Export = this;

	// register export base at import base
	Import->m_export = dynamic_cast<DataExportItem*>(this);

	UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::add_data_import: DataExport '" << this->name() <<"' registered at DataImport '" << importItem->name() << "'.\n");
	return true;
}

template<typename TDataType, typename TPositionType>
bool DataExport<TDataType, TPositionType>::
remove_data_import(DataImportItem* importItem)
{
	DataImport<TDataType, TPositionType>* Import = dynamic_cast<DataImport<TDataType, TPositionType>*>(importItem);

	if(Import == NULL)
	{
		UG_LOG("DataExport::add_data_import: Data type and/or position type of Export<->Import does not match. Cannot register.\n");
		return false;
	}

	typename std::vector<DataImportItem*>::iterator importIter;
	importIter = find(m_importList.begin(), m_importList.end(), importItem);
	if(importIter == m_importList.end())
	{
		UG_LOG("DataExport::remove_data_import: DataImport not found. Cannot unregister.\n");
		return false;
	}

	Import->m_export = NULL;
	Import->m_Export = NULL;

	// erase import
	m_importList.erase(importIter);

	// if import list is void, delete positions and values
	if(m_importList.empty())
	{
		m_positions.clear();
		m_values.clear();
		for(std::size_t s = 0; s < m_derivatives.size(); ++s)
		{
			for(std::size_t ip = 0; ip < m_derivatives[s].size(); ++ip)
			{
				m_derivatives[s][ip].clear();
			}
			m_derivatives[s].clear();
		}
		m_derivatives.clear();
	}

	UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::remove_data_import: DataExport '" << this->name() <<"' unregistered from DataImport '" << importItem->name() << "'.\n");
	return true;
}


template<typename TDataType, typename TPositionType>
bool DataExport<TDataType, TPositionType>::
set_positions(const std::vector<position_type>& positions, bool overwrite)
{
	// if no positions set
	if(m_positions.empty() || overwrite)
	{
		// copy positions and set same size for value list
		m_positions = positions;

		UG_DLOG(LIB_DISC_LINKER, 3, "DataExport::set_positions: Resize value array to size " << this->num_ip() << ".\n");
		m_values.resize(this->num_ip());

		UG_DLOG(LIB_DISC_LINKER, 3, "DataExport::set_positions: Resize derivatives array to size " << this->num_sys() <<" x " << this->num_ip() << " x num_sh(s) for each system s.\n");
		m_derivatives.resize(this->num_sys());
		for(std::size_t s = 0; s < this->num_sys(); ++s)
		{
			m_derivatives[s].resize(this->num_ip());
			for(std::size_t ip = 0; ip < m_derivatives[s].size(); ++ip)
			{
				m_derivatives[s][ip].resize(num_sh(s));
			}
		}

		UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::set_positions: Positions set.\n");

		if(set_linker_positions(m_positions) != true)
		{
			UG_LOG("DataExport::set_positions: Error while setting positions in linker imports.\n");
			return false;
		}
		UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::set_positions: Positions forwarded to Linker import.\n");

		return true;
	}

	// if already positions set, check, if positions are equal
	if(positions.size() != this->num_ip())
	{
		UG_LOG("DataExport::set_positions: Different number of positions as already set positions.\n");
		return false;
	}

	for(uint ip = 0; ip < this->num_ip(); ++ip)
	{
		if(m_positions[ip] != positions[ip])
		{
			UG_LOG("DataExport::set_positions: Position " << ip << " is different from already set positions.\n");
			return false;
		}
	}

	UG_DLOG(LIB_DISC_LINKER, 2, "DataExport::set_positions: Positions already set. New Import positions are equal.\n");
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataExport<TDataType, TPositionType>::
equal(const DataExportItem& v) const
{
	const DataExport<data_type, position_type>* cast_v = dynamic_cast<const DataExport<data_type, position_type>*>(&v);

	// if type is not the same, they are nor equal
	if(cast_v == NULL) return false;

	// Compute if pointers to Base possibility are the same
	if(this->m_possibilityItem != cast_v->get_possibility_item()) return false;

	// number of positions have to match
	if(this->num_ip() != cast_v->num_ip()) return false;

	// check each position
	for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
	{
		if(m_positions[ip] != cast_v->m_positions[ip]) return false;
	}

	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataExport<TDataType, TPositionType>::
print_positions() const
{
	for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
	{
		UG_LOG(m_positions[ip]);
		if(ip != this->num_ip() - 1) UG_LOG(", ");
	}
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataExport<TDataType, TPositionType>::
print_values() const
{
	for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
	{
		UG_LOG(operator[](ip));
		if(ip != this->num_ip() - 1) UG_LOG(", ");
	}
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataExport<TDataType, TPositionType>::
print_derivatives(std::string offset) const
{
	UG_LOG("Depending on "<< this->num_sys()<< " systems.\n");
	for(std::size_t s = 0; s < this->num_sys(); ++s)
	{
		UG_LOG(offset << "w.r.t. sys = " << sys(s) << ": ");

		for(std::size_t k = 0; k < num_sh(s); ++k)
		{
			if(k==0) {UG_LOG(" k="<< k<< ": [ ");}
			else {UG_LOG(offset << "                 k="<< k<< ": [ ");};

			for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
			{
				{
					UG_LOG(operator()(s,ip,k));
					if(ip != this->num_ip() - 1) UG_LOG(", ");
				}
			}
			UG_LOG(" ]\n");
		}
	}
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataExport<TDataType, TPositionType>::
print_info(std::string offset) const
{
	UG_LOG(offset << "IPs: [ "); if(print_positions() != true) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Values: [ "); if(print_values() != true) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Derivatives: "); if(print_derivatives(offset) != true) return false;
	return true;
}


//////////////////////////////
// Data Import Position
//////////////////////////////

template<typename TPositionType>
bool
DataImportPosition<TPositionType>::
set_positions(const std::vector<position_type>& pos, bool overwrite)
{
	// copy new positions (elementwise)
	m_positions = pos;
	return true;
}

//////////////////////////////
// Data Import
//////////////////////////////

template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
set_positions(const std::vector<position_type>& pos, bool overwrite)
{
	// copy new positions (elementwise)
	this->m_positions = pos;

	// adjust lin defect array for ip's
	for(std::size_t i = 0; i < m_linearized_defect.size(); ++i)
	{
		m_linearized_defect[i].resize(this->num_ip());
	}

	// if no export set, do nothing
	if(this->m_export == NULL) return true;

	// if export set is an DataExport
	typename ug::DataExport<TDataType, TPositionType>* Cast_Export = dynamic_cast< typename ug::DataExport<TDataType, TPositionType>* >(this->m_export);
	if(Cast_Export != NULL)
	{
		UG_DLOG(LIB_DISC_LINKER, 2, "DataImport::set_positions: DataExport set at this Import. Setting positions.\n");
		if(Cast_Export->set_positions(pos, overwrite) != true)
		{
			UG_LOG("DataImport::set_positions: Cannot set positions in export.\n");
			return false;
		}
	}

	// if export set is an DataLinker
	// TODO: implement

	UG_DLOG(LIB_DISC_LINKER, 2, "DataImport::set_positions: positions set.");
	return true;
}


template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
set_num_eq(std::size_t num_eq)
{
	m_num_eq = num_eq;

	m_linearized_defect.resize(m_num_eq);

	for(std::size_t i = 0; i < m_linearized_defect.size(); ++i)
	{
		m_linearized_defect[i].resize(this->num_ip());
	}

	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
link_data_export(DataExportItem* exportItem)
{
	typename ug::DataExport<TDataType, TPositionType>* Export = dynamic_cast< typename ug::DataExport<TDataType, TPositionType>* >(exportItem);
	if(Export == NULL)
	{
		UG_LOG("DataImport::link_data_export: Data type and/or position type of Export<->Import does not match. Cannot register.\n");
		return false;
	}

	if(this->m_export != NULL)
	{
		UG_LOG("DataImport::link_data_export: Already an Export registered.\n");
		return false;
	}

	// register at export
	if(Export->add_data_import(this) != true)
	{
		UG_LOG("DataImport::link_data_export: Can not register at export.\n");
		this->m_export = NULL;
		return false;
	}

	// remember Export
	m_Export = Export;

	// remember export base
	this->m_export = exportItem;

	UG_DLOG(LIB_DISC_LINKER, 2, "DataImport::link_data_export: DataExport '" << exportItem->name() << " to DataImport " << this->name() << ".\n");
	return true;
}


template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
clear_data_export()
{
	// if registered at export
	if(this->m_export != NULL)
	{
		if(m_Export->remove_data_import(this) != true)
		{
			UG_LOG("DataImport::clear_data_export: Can not unregister from export.");
			return false;
		}
		UG_DLOG(LIB_DISC_LINKER, 2, "DataImport::clear_data_export: DataExport removed from DataImport '" << this->name() << "'.\n");

		// reset export
		m_Export = NULL;
		this->m_export = NULL;
	}
	else
	{
		UG_DLOG(LIB_DISC_LINKER, 2, "DataImport::clear_data_export: No DataExport linked to DataImport '" << this->name() << "'.\n");
	}

	return true;
}


template<typename TDataType, typename TPositionType>
DataImport<TDataType, TPositionType>::
~DataImport()
{
	UG_DLOG(LIB_DISC_LINKER, 2, "DataImport::~DataImport: Deleting Data Import " << this->name() << ".\n");
	clear_data_export();
}

template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
print_positions() const
{
	for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
	{
		UG_LOG(this->m_positions[ip]);
		if(ip != this->num_ip() - 1) UG_LOG(", ");
	}
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
print_values() const
{
	if(!this->is_linked())
	{
		UG_LOG("Not linked.");
		return true;
	}

	for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
	{
		UG_LOG(this->operator[](ip));
		if(ip != this->num_ip() - 1) UG_LOG(", ");
	}
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
print_derivatives(std::string offset) const
{
	if(!this->is_linked())
	{
		UG_LOG("Not linked.\n");
		return true;
	}

	UG_LOG("Depending on "<< this->num_sys()<< " systems.\n");
	for(std::size_t s = 0; s < this->num_sys(); ++s)
	{
		for(std::size_t k = 0; k < this->num_sh(s); ++k)
		{
			UG_LOG(offset << "w.r.t.: sys = " << this->sys(s) <<", k="<< k<< ": [ ");

			for(std::size_t ip = 0; ip < this->num_ip(); ++ip)
			{
				{
					UG_LOG(this->operator()(s, ip, k));
					if(ip != this->num_ip() - 1) UG_LOG(", ");
				}
			}
			UG_LOG(" ]\n");
		}
	}
	return true;
}

template<typename TDataType, typename TPositionType>
bool
DataImport<TDataType, TPositionType>::
print_info(std::string offset) const
{
	UG_LOG(offset << "IPs: [ "); if(print_positions() != true) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Values: [ "); if(print_values() != true) return false; UG_LOG(" ]\n");
	UG_LOG(offset << "Derivatives: "); if(print_derivatives(offset) != true) return false;
	return true;
}

template <typename TDataType, typename TPositionType, typename TAlgebra>
std::ostream& operator<<(std::ostream& out, const DataImport<TDataType, TPositionType>& data)
{
	for(std::size_t ip = 0; ip < data.num_ip(); ++ip)
	{
		out << data[ip];
		if(ip != data.num_ip() - 1) out << ", ";
	}
	return out;
}

} // End namespace ug

#endif /*__H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT_IMPL__*/
