/*
 * elementdata.h
 *
 *  Created on: 30.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT__

#include <typeinfo>
#include <string>
#include <vector>

#include "common/common.h"

#include "lib_discretization/spacial_discretization/disc_coupling/element_data_items.h"

namespace ug{

// predeclaration
template <typename TDataType, typename TPositionType> class DataImport;
template <typename TDataType, typename TPositionType> class DataExport;


template <typename TDataType, typename TPositionType>
class DataExport : public DataExportItem{
	friend class DataImport<TDataType,TPositionType>;
	friend class DataContainer;

	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;

	protected:
		// Only Data Container can create an instance
		DataExport(std::string name, DataPossibilityItem* possibility) 	:
			DataExportItem(name,&typeid(TDataType),&typeid(TPositionType), possibility)
			{m_positions.clear(); m_values.clear(); m_derivatives.clear();};

	public:
		// number of values / positions / derivatives == number of integration points
		inline std::size_t num_ip() const {return m_positions.size();};

		// return value at ip (read only, are updated by compute())
		inline const data_type& operator[](std::size_t ip) const {
			UG_ASSERT(ip < m_values.size(), "Accessing Value array at not allocated position.");
			return m_values[ip];};

		// return vector of values
		inline const std::vector<data_type>& values() const	{return m_values;}

		// return value of derivative with respect to unknown k at ip (read only, are updated by compute())
		inline const data_type& operator()(std::size_t loc_sys, std::size_t ip, std::size_t k) const {
			UG_ASSERT(loc_sys < m_derivatives.size(), "Accessing Derivative array at not allocated position: sys = " << loc_sys);
			UG_ASSERT(ip < m_derivatives[loc_sys].size(), "Accessing Derivative array at not allocated position ip = " << ip);
			UG_ASSERT(k < m_derivatives[loc_sys][ip].size(), "Accessing Derivative array at not allocated position k = " << k);
			return m_derivatives[loc_sys][ip][k];};

		inline const std::vector<std::vector<data_type> >& derivatives(std::size_t loc_sys) const {
			UG_ASSERT(loc_sys < m_derivatives.size(), "Accessing Derivative array at not allocated position: sys = " << loc_sys);
			return m_derivatives[loc_sys];}

		// return position number i (read only, are induced by import)
		inline const position_type& position(std::size_t ip) const {
			UG_ASSERT(ip < m_positions.size(), "Accessing Position array at not allocated position.");
			return m_positions[ip];};

		// compute values at positions using evaluation function
		virtual void compute(bool compute_derivatives) = 0;

		// destructor
		~DataExport();


	public:
		// check if exports do
		// 1. use same evaluation function
		// 2. evaluate at exactly the same positions
		// (A necessary condition are same c++ template types, of course)
		virtual bool equal(const DataExportItem& v) const;

		// add an import
		virtual bool add_data_import(DataImportItem* importItem);

		// remove the import with this pointer
		virtual bool remove_data_import(DataImportItem* importItem);

	public:
		// set the positions
		// a) if first Import set values
		// b) if already set: check if positions are equal for other Imports
		bool set_positions(const std::vector<position_type>& positions, bool overwrite = false);

		// help function: is called to set values in linker imports as well
		virtual bool set_linker_positions(const std::vector<position_type>& positions) {return true;};

	public:
		// print infos to UG_LOG stream
		virtual bool print_positions() const;
		virtual bool print_values() const;
		virtual bool print_derivatives(std::string offset) const;
		virtual bool print_info(std::string offset) const;

	// Import (Linker) side
	public:
		// number of data exports linked by this linker
		virtual std::size_t num_data_exports() const {return 0;};

		// name of import i
		virtual std::string import_name(std::size_t) const {return "";};

		// add a Data Export number i
		virtual bool link_data_export(std::size_t i, DataExportItem* exportItem) {return false;};

		// remove Data Export number i
		virtual bool clear_data_export(std::size_t i) {return false;};

		// get registered export of slot i
		virtual const DataExportItem* get_data_export(std::size_t i) const {return NULL;};

		// return if an export is set at slot i
		virtual bool is_linked(std::size_t i) const {return false;};

		// return if all exports are set
		virtual bool is_linked() const {return false;};


	protected:
		// points (gives implicitly the number of ip's: num_ip = m_positions.size())
		std::vector<position_type> m_positions;

		// Data (size: (0, ... , num_ip-1))
		std::vector<data_type> m_values;

		// Data (size: (0, ..., num_sys-1) x (0, ... , num_ip-1) x (0, ... , num_sh-1)
		std::vector<std::vector<std::vector<data_type> > > m_derivatives;
};


template <typename TPositionType>
class DataImportPosition : public DataImportItem {
	public:
		typedef TPositionType position_type;

	public:
		DataImportPosition(std::string name, const std::type_info* dataType) :
			DataImportItem(name, dataType, &typeid(TPositionType))
			{m_positions.clear();};

		// number of values / positions / derivatives == number of integration points (num_ip)
		inline std::size_t num_ip() const {return m_positions.size();};

		// return position number i (read only, are induced by import)
		inline const position_type& position(std::size_t ip) const {
			UG_ASSERT(ip < m_positions.size(), "Accessing Position array at not allocated position.");
			return m_positions[ip];};

		// set new positions
		virtual bool set_positions(const std::vector<position_type>& pos, bool overwrite = false);

		// get positions
		const std::vector<position_type>& get_positions() const {return m_positions;}

	public:
		virtual bool print_positions() const {return false;}
		virtual bool print_values() const {return false;}
		virtual bool print_derivatives(std::string offset) const {return false;}
		virtual bool print_info(std::string offset) const {return false;}


	protected:
		// positions of this data import (saved here), gives implicitly num_ip = m_positions.size()
		// size: (0, ... , num_ip-1)
		std::vector<position_type> m_positions;


};


template <typename TDataType, typename TPositionType>
class DataImport : public DataImportPosition<TPositionType> {
	friend class DataExport<TDataType,TPositionType>;
	friend class DataContainer;

	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;

	public:
		DataImport(std::string name) :
			DataImportPosition<TPositionType>(name,&typeid(TDataType)),
			m_Export(NULL)
			{m_linearized_defect.clear();};

		// number of equations the data is provided for (num_eq)
		inline std::size_t num_eq() const {return m_linearized_defect.size();};

		// return value at ip (read only, are updated by compute())
		inline const data_type& operator[](std::size_t ip) const
			{UG_ASSERT(m_Export != NULL, "No Export set.");
			 return m_Export->operator[](ip);};

		// return vector of values
		inline const std::vector<data_type>& values() const
		{UG_ASSERT(m_Export != NULL, "No Export set.");
			return m_Export->values();
		}

		// return vector of values
		inline const std::vector<std::vector<data_type> >& derivatives(std::size_t loc_sys) const
		{UG_ASSERT(m_Export != NULL, "No Export set.");
			return m_Export->derivatives(loc_sys);
		}

		// return value of derivative with respect to unknown k at ip (read only, are updated by compute())
		inline const data_type& operator()(std::size_t sys, std::size_t ip, std::size_t k) const
			{UG_ASSERT(m_Export != NULL, "No Export set.");
			return m_Export->operator()(sys, ip, k);};

		// linearized defect (with respect to import) for the j'th local defect and position ip
		inline data_type& lin_defect(std::size_t j, std::size_t ip){
			UG_ASSERT(ip < m_linearized_defect[j].size(), "Accessing Global Index Array at position, that is not allocated.");
			UG_ASSERT(j < m_linearized_defect.size(), "Accessing Global Index Array at position, that is not allocated.");
			return m_linearized_defect[j][ip];
		}

		// linearized defect (with respect to import) for the j'th local defect and position ip
		inline const data_type& lin_defect(std::size_t j, std::size_t ip) const {
			UG_ASSERT(ip < m_linearized_defect[j].size(), "Accessing Global Index Array at position, that is not allocated.");
			UG_ASSERT(j < m_linearized_defect.size(), "Accessing Global Index Array at position, that is not allocated.");
			return m_linearized_defect[j][ip];
		}


	public:
		// set new positions
		virtual bool set_positions(const std::vector<position_type>& pos, bool overwrite = false);

		// set number of equations the defect is provided for
		bool set_num_eq(std::size_t num_eq);

		~DataImport();

	public:
		// register a data export (internally the import is registered at export as well)
		virtual bool link_data_export(DataExportItem* exportItem);

		// unregister the data export
		virtual bool clear_data_export();

	public:
		// print positions to UG_LOG stream
		virtual bool print_positions() const;
		virtual bool print_values() const;
		virtual bool print_derivatives(std::string offset) const;
		virtual bool print_info(std::string offset) const;

	protected:
		// number of equations the defect is provided for
		std::size_t m_num_eq;

		// linearization of defect of equation s for ip (defect is considered as function of this import, and linearized at given current solution)
		// (size: (0, ... , num_ip-1) x (0, ... , num_eq) )
		std::vector<std::vector<data_type> > m_linearized_defect;

		// Export, this Import is linked to
		DataExport<data_type, position_type>* m_Export;
};

// print current values of import to ostream
template <typename TDataType, typename TPositionType>
std::ostream& operator<<(std::ostream& out, const DataImport<TDataType,TPositionType>& data);

};
#include "element_data_import_export_impl.h"

#endif /* __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT__ */
