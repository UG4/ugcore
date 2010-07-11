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

#include "element_data_items.h"

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
			DataExportItem(name,&typeid(TDataType),&typeid(TPositionType), possibility),
			 m_bGlobalIPs(false)
			{m_vPosition.clear(); m_vValue.clear(); m_vvvDerivatives.clear();};

	public:
		// number of values / positions / derivatives == number of integration points
		inline size_t num_ip() const {return m_vPosition.size();};

		// return value at ip (read only, are updated by compute())
		inline const data_type& operator[](size_t ip) const {
			UG_ASSERT(ip < m_vValue.size(), "Accessing Value array at not allocated position.");
			return m_vValue[ip];};

		// return vector of values
		inline const std::vector<data_type>& values() const	{return m_vValue;}

		// return value of derivative with respect to unknown k at ip (read only, are updated by compute())
		inline const data_type& operator()(size_t loc_sys, size_t ip, size_t k) const {
			UG_ASSERT(loc_sys < m_vvvDerivatives.size(), "Accessing Derivative array at not allocated position: sys = " << loc_sys);
			UG_ASSERT(ip < m_vvvDerivatives[loc_sys].size(), "Accessing Derivative array at not allocated position ip = " << ip);
			UG_ASSERT(k < m_vvvDerivatives[loc_sys][ip].size(), "Accessing Derivative array at not allocated position k = " << k);
			return m_vvvDerivatives[loc_sys][ip][k];};

		// returns array of derivatives for system loc_sys
		inline const std::vector<std::vector<data_type> >& derivatives(size_t loc_sys) const {
			UG_ASSERT(loc_sys < m_vvvDerivatives.size(), "Accessing Derivative array at not allocated position: sys = " << loc_sys);
			return m_vvvDerivatives[loc_sys];}

		// return position number i (read only, are induced by import)
		inline const position_type& position(size_t ip) const {
			UG_ASSERT(ip < m_vPosition.size(), "Accessing Position array at not allocated position.");
			return m_vPosition[ip];};

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

		// specifies if local or global ips are needed to compute data
		// This is only needed for export identification, which will be performed on local ips
		// Then in computations, global ids must be given, if needed
		// As default, local ips are assumed
		void need_global_ids(bool bGlobalNeeded) {m_bGlobalIPs = bGlobalNeeded;}

		// help function: is called to set values in linker imports as well
		virtual bool set_linker_positions(const std::vector<position_type>& positions) {return true;};

	public:
		// print infos to UG_LOG stream
		virtual bool print_positions() const;
		virtual bool print_values() const;
		virtual bool print_derivatives(std::string offset) const;
		virtual bool print_info(std::string offset) const;

	protected:
		// sets the number of systems the export depends on
		virtual bool set_num_sys(size_t num_sys)
		{
			if(!DataExportItem::set_num_sys(num_sys)) return false;
			m_vvvDerivatives.resize(num_sys);
			for(size_t s = 0; s < this->num_sys(); ++s)
			{
				m_vvvDerivatives[s].resize(num_ip());
			}
			return true;
		}

		// sets the number of shape functions (num_sh) for the local system
		virtual bool set_num_sh(size_t num_sh, size_t loc_sys = 0)
		{
			if(!DataExportItem::set_num_sh(num_sh, loc_sys)) return false;
			for(size_t s = 0; s < num_sys(); ++s)
			{
				for(size_t ip = 0; ip < num_ip(); ++ip)
				{
					m_vvvDerivatives[s][ip].resize(num_sh);
				}
			}
			return true;
		}

	protected:
		// specifies if local or global ips are needed for this export to compute data
		bool m_bGlobalIPs;

		// points (gives implicitly the number of ip's: num_ip = m_vPosition.size())
		std::vector<position_type> m_vPosition;

		// Data (size: (0, ... , num_ip-1))
		std::vector<data_type> m_vValue;

		// Data (size: (0, ..., num_sys-1) x (0, ... , num_ip-1) x (0, ... , num_sh-1)
		std::vector<std::vector<std::vector<data_type> > > m_vvvDerivatives;
};


template <typename TPositionType>
class DataImportPosition : public DataImportItem {
	public:
		typedef TPositionType position_type;

	public:
		DataImportPosition(std::string name, const std::type_info* dataType) :
			DataImportItem(name, dataType, &typeid(TPositionType))
			{m_vPosition.clear();};

		// number of values / positions / derivatives == number of integration points (num_ip)
		inline size_t num_ip() const {return m_vPosition.size();};

		// return position number i (read only, are induced by import)
		inline const position_type& position(size_t ip) const {
			UG_ASSERT(ip < m_vPosition.size(), "Accessing Position array at not allocated position.");
			return m_vPosition[ip];};

		// set new positions
		virtual bool set_positions(const std::vector<position_type>& pos, bool overwrite = false);

		// get positions
		const std::vector<position_type>& get_positions() const {return m_vPosition;}

	public:
		virtual bool print_positions() const {return false;}
		virtual bool print_values() const {return false;}
		virtual bool print_derivatives(std::string offset) const {return false;}
		virtual bool print_info(std::string offset) const {return false;}

	protected:
		// positions of this data import (saved here), gives implicitly num_ip = m_vPosition.size()
		// size: (0, ... , num_ip-1)
		std::vector<position_type> m_vPosition;
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
			m_pTypeExport(NULL)
			{m_vvLinearizedDefect.clear();};

		// number of equations the data is provided for (num_eq)
		inline size_t num_eq() const {return m_vvLinearizedDefect.size();};

		// return value at ip (read only, are updated by compute())
		inline const data_type& operator[](size_t ip) const
			{UG_ASSERT(m_pTypeExport != NULL, "No Export set.");
			 return m_pTypeExport->operator[](ip);};

		// return vector of values
		inline const std::vector<data_type>& values() const
		{UG_ASSERT(m_pTypeExport != NULL, "No Export set.");
			return m_pTypeExport->values();
		}

		// return vector of values
		inline const std::vector<std::vector<data_type> >& derivatives(size_t loc_sys) const
		{UG_ASSERT(m_pTypeExport != NULL, "No Export set.");
			return m_pTypeExport->derivatives(loc_sys);
		}

		// return value of derivative with respect to unknown k at ip (read only, are updated by compute())
		inline const data_type& operator()(size_t loc_sys, size_t ip, size_t k) const
			{UG_ASSERT(m_pTypeExport != NULL, "No Export set.");
			return m_pTypeExport->operator()(loc_sys, ip, k);};

		// linearized defect (with respect to import) for the j'th local defect and position ip
		inline data_type& lin_defect(size_t j, size_t ip){
			UG_ASSERT(ip < m_vvLinearizedDefect[j].size(), "Accessing Global Index Array at position, that is not allocated.");
			UG_ASSERT(j < m_vvLinearizedDefect.size(), "Accessing Global Index Array at position, that is not allocated.");
			return m_vvLinearizedDefect[j][ip];
		}

		// linearized defect (with respect to import) for the j'th local defect and position ip
		inline const data_type& lin_defect(size_t j, size_t ip) const {
			UG_ASSERT(ip < m_vvLinearizedDefect[j].size(), "Accessing Global Index Array at position, that is not allocated.");
			UG_ASSERT(j < m_vvLinearizedDefect.size(), "Accessing Global Index Array at position, that is not allocated.");
			return m_vvLinearizedDefect[j][ip];
		}

		// add offdiagonal coupling to other system 's'
		virtual bool add_offdiagonal(FlexLocalMatrix& J, size_t loc_sys, number s_a)
		{
			for(size_t k = 0; k < this->num_sh(loc_sys); ++k)
			{
				for(size_t j = 0; j < this->num_eq(); ++j)
				{
					for(size_t ip = 0; ip < this->num_ip(); ++ip)
					{
						J(j, k) += s_a * (this->lin_defect(j, ip) * this->operator()(loc_sys, ip, k));
					}
				}
			}
			return true;
		}

	public:
		// set new positions
		virtual bool set_positions(const std::vector<position_type>& pos, bool overwrite = false);

		// set number of equations the defect is provided for
		bool set_num_eq(size_t num_eq);

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
		size_t m_numEq;

		// linearization of defect of equation s for ip (defect is considered as function of this import, and linearized at given current solution)
		// (size: (0, ... , num_ip-1) x (0, ... , num_eq) )
		std::vector<std::vector<data_type> > m_vvLinearizedDefect;

		// Export, this Import is linked to
		DataExport<data_type, position_type>* m_pTypeExport;
};

// print current values of import to ostream
template <typename TDataType, typename TPositionType>
std::ostream& operator<<(std::ostream& out, const DataImport<TDataType,TPositionType>& data);

};
#include "element_data_import_export_impl.h"

#endif /* __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT__ */
