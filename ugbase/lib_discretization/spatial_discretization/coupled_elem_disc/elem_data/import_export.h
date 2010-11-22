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

#include "common/math/ugmath.h"

#include "data_items.h"
#include "ip_derivative.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{

// predeclaration
template <typename TDataType> class DataImport;
template <typename TDataType> class DataExport;


template <typename TDataType>
class DataExport : public DataExportItem{
	friend class DataImport<TDataType>;
	friend class DataContainer;

	public:
		typedef TDataType data_type;

	protected:
		// Only Data Container can create an instance
		DataExport(std::string name, DataPossibilityItem* possibility) 	:
			DataExportItem(name,&typeid(TDataType), possibility), m_bGlobalIPs(false),
			m_posDim(-1), m_numIp(0)
			{m_vValue.clear(); m_vvDerivatives.clear();};

	public:
		// number of values / positions / derivatives == number of integration points
		inline size_t num_ip() const {return m_numIp;};

		// return value at ip (read only, are updated by compute())
		inline const data_type& operator[](size_t ip) const {
			UG_ASSERT(ip < m_vValue.size(), "Invalid index. (ip = "<<ip<<", size = "<<m_vValue.size()<<")");
			return m_vValue[ip];};

		// return vector of values
		inline const std::vector<data_type>& values() const	{return m_vValue;}

		// return value of derivative with respect to unknown k at ip (read only, are updated by compute())
		inline const data_type& operator()(size_t loc_sys, size_t ip, size_t fct, size_t dof) const {
			UG_ASSERT(loc_sys < m_vvDerivatives.size(), "Invalid index. (loc_sys = "<<loc_sys<<", size = "<<m_vvDerivatives.size()<<")");
			UG_ASSERT(ip < m_vvDerivatives[loc_sys].size(), "Invalid index. (ip = "<<ip<<", size = "<< m_vvDerivatives[loc_sys].size()<<")");
			return m_vvDerivatives[loc_sys][ip](fct, dof);};

		// returns array of derivatives for system loc_sys
		inline const std::vector<IPDerivative<data_type> >& derivatives(size_t loc_sys) const {
			UG_ASSERT(loc_sys < m_vvDerivatives.size(), "Invalid index. (loc_sys = "<<loc_sys<<", size = "<<m_vvDerivatives.size()<<")");
			return m_vvDerivatives[loc_sys];}

		// return position number i (read only, are induced by import)
		template <int dim>
		inline const MathVector<dim>& position(size_t ip) const {
			UG_ASSERT(ip < num_ip(), "Invalid index. (ip = "<<ip<<", size = "<<num_ip()<<")");
			return const_cast<DataExport<TDataType>*>(this)->get_positions<dim>()[ip];};

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
		template <int dim>
		bool set_positions(const std::vector<MathVector<dim> >& positions, bool overwrite = false);

		// specifies if local or global ips are needed to compute data
		// This is only needed for export identification, which will be performed on local ips
		// Then in computations, global ids must be given, if needed
		// As default, local ips are assumed
		void need_global_ids(bool bGlobalNeeded) {m_bGlobalIPs = bGlobalNeeded;}

		// help function: is called to set values in linker imports as well
		// TODO: What to do here ?!?!??
		virtual bool set_linker_positions(const std::vector<MathVector<2> >& positions) {return true;};

	public:
		// print infos to UG_LOG stream
		virtual bool print_positions() const;
		virtual bool print_values() const;
		virtual bool print_derivatives(std::string offset) const;
		virtual bool print_info(std::string offset) const;

		// sets the number of systems the export depends on
		virtual bool set_num_sys(size_t num_sys)
		{
			if(!DataExportItem::set_num_sys(num_sys)) return false;
			m_vvDerivatives.resize(num_sys);
			for(size_t s = 0; s < this->num_sys(); ++s)
			{
				m_vvDerivatives[s].resize(num_ip());
			}
			return true;
		}

		// sets the number of functions for the local system
		virtual bool set_num_fct(size_t num_fct, size_t loc_sys = 0)
		{
			if(!DataExportItem::set_num_fct(num_fct, loc_sys))
				{UG_LOG("DataExport::set_num_fct: Error while setting num_fct.\n"); return false;}
			for(size_t s = 0; s < num_sys(); ++s)
			{
				for(size_t ip = 0; ip < num_ip(); ++ip)
				{
					m_vvDerivatives[s][ip].set_num_fct(num_fct);
				}
			}
			return true;
		}

		// sets the number of dofs for a functions for the local system
		virtual bool set_num_dofs(size_t fct, size_t num_dofs, size_t loc_sys = 0)
		{
			if(!DataExportItem::set_num_dofs(fct, num_dofs, loc_sys))
				{UG_LOG("DataExport::set_num_dofs: Error while setting num_dofs.\n"); return false;}
			for(size_t s = 0; s < num_sys(); ++s)
			{
				for(size_t ip = 0; ip < num_ip(); ++ip)
				{
					m_vvDerivatives[s][ip].set_num_dofs(fct, num_dofs);
				}
			}
			return true;
		}

	protected:
		template <int dim>
		std::vector<MathVector<dim> >& get_positions()
		{
			static std::vector<MathVector<dim> > vPositions;
			return vPositions;
		}

	protected:
		// specifies if local or global ips are needed for this export to compute data
		bool m_bGlobalIPs;

		// current dimension for positions
		int m_posDim;

		// current number of ips
		size_t m_numIp;

		// Data (size: (0, ... , num_ip-1))
		std::vector<data_type> m_vValue;

		// Data (size: (0, ..., num_sys-1) x (0, ... , num_ip-1))
		std::vector<std::vector<IPDerivative<data_type> > > m_vvDerivatives;
};

template <typename TDataType>
class DataImport : public DataImportItem {
	friend class DataExport<TDataType>;
	friend class DataContainer;

	public:
		typedef TDataType data_type;

	public:
		DataImport(std::string name) :
			DataImportItem(name,&typeid(TDataType)),
			m_posDim(-1), m_numIp(0), m_pTypeExport(NULL)
			{m_vLinearizedDefect.clear();};

		// number of values / positions / derivatives == number of integration points (num_ip)
		inline size_t num_ip() const {return m_numIp;};

		// return position number i (read only, are induced by import)
		template <int dim>
		inline const MathVector<dim>& position(size_t ip) const {
			UG_ASSERT(ip < num_ip(), "Invalid index. (ip = "<<ip<<", size = "<<num_ip()<<")");
			return get_positions<dim>()[ip];};

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
		inline const std::vector<IPDerivative<data_type> >& derivatives(size_t loc_sys) const
		{UG_ASSERT(m_pTypeExport != NULL, "No Export set.");
			return m_pTypeExport->derivatives(loc_sys);
		}

		// return value of derivative with respect to unknown k at ip (read only, are updated by compute())
		inline const data_type& operator()(size_t loc_sys, size_t ip, size_t fct, size_t dof) const
			{UG_ASSERT(m_pTypeExport != NULL, "No Export set.");
			return m_pTypeExport->operator()(loc_sys, ip, fct, dof);};

		// linearized defect (with respect to import) for the j'th local defect and position ip
		inline data_type& lin_defect(size_t ip, size_t fct, size_t dof){
			UG_ASSERT(ip < m_vLinearizedDefect.size(), "Invalid index. (ip = "<<ip<<", size = "<<m_vLinearizedDefect.size()<<")");
			return m_vLinearizedDefect[ip](fct, dof);
		}

		// linearized defect (with respect to import) at position ip
		inline const data_type& lin_defect(size_t ip, size_t fct, size_t dof) const {
			UG_ASSERT(ip < m_vLinearizedDefect.size(), "Invalid index. (ip = "<<ip<<", size = "<<m_vLinearizedDefect.size()<<")");
			return m_vLinearizedDefect[ip](fct, dof);
		}

		// add offdiagonal coupling to other system 's'
		virtual bool add_offdiagonal(LocalMatrixBase& J, size_t loc_sys, number s_a)
		{
			{UG_LOG("DataImport::add_offdiagonal: Should not be called here.\n"); return false;}
			return false;
		}

	public:
		// set new positions
		// TODO: What to do here ?!?!??!? this was virtual
		template <int dim>
		bool set_positions(const std::vector<MathVector<dim> >& pos, bool overwrite = false);

		// set number of equations the defect is provided for
		bool set_num_eq_fct(size_t num_eq_fct)
		{
			if(!DataImportItem::set_num_eq_fct(num_eq_fct)) return false;
			for(size_t ip = 0; ip < m_vLinearizedDefect.size(); ++ip)
			{
				m_vLinearizedDefect[ip].set_num_fct(num_eq_fct);
			}
			return true;
		}

		// set number of equations the defect is provided for
		bool set_num_eq_dofs(size_t fct, size_t num_eq_dofs)
		{
			if(!DataImportItem::set_num_eq_dofs(fct, num_eq_dofs)) return false;
			for(size_t ip = 0; ip < m_vLinearizedDefect.size(); ++ip)
			{
				m_vLinearizedDefect[ip].set_num_dofs(fct, num_eq_dofs);

			}
			return true;
		}

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
		// get positions
		template <int dim>
		std::vector<MathVector<dim> >& get_positions() const
		{
			static std::vector<MathVector<dim> > vPositions;
			return vPositions;
		}

	protected:
		// current dimension for positions
		int m_posDim;

		// current number of ips
		size_t m_numIp;

		// linearization of defect of equation s for ip (defect is considered as function of this import, and linearized at given current solution)
		// (size: (0, ... , num_ip-1) x (0, ... , num_eq) )
		std::vector<IPDerivative<data_type> > m_vLinearizedDefect;

		// Export, this Import is linked to
		DataExport<data_type>* m_pTypeExport;
};

template <typename TDataType, typename TAlgebra>
class AlgebraDataImport : public DataImport<TDataType> {
	public:
		typedef TDataType data_type;

		typedef typename TAlgebra::matrix_type::value_type value_type;

	public:
		AlgebraDataImport(std::string name) : DataImport<TDataType>(name)
		{};

		// add offdiagonal coupling to other system 's'
		virtual bool add_offdiagonal(LocalMatrixBase& baseJ, size_t loc_sys, number s_a)
		{
			LocalMatrix<value_type>* J = dynamic_cast<LocalMatrix<value_type>*>(&baseJ);
			if(J == NULL)
				{UG_LOG("DataImport::add_offdiagonal: Cannot convert Matrix to type based matrix.\n"); return false;}

			for(size_t fct1 = 0; fct1 < this->num_eq_fct(); ++fct1)
			{
				for(size_t dof1 = 0; dof1 < this->num_eq_dofs(fct1); ++dof1)
				{
					for(size_t fct2 = 0; fct2 < this->num_fct(loc_sys); ++fct2)
					{
						for(size_t dof2 = 0; dof2 < this->num_dofs(fct2, loc_sys); ++dof2)
						{
							for(size_t ip = 0; ip < this->num_ip(); ++ip)
							{
								(*J)(fct1, dof1, fct2, dof2) += s_a * (this->lin_defect(ip, fct1, dof1) * this->operator()(loc_sys, ip, fct2, dof2));
							}
						}
					}
				}
			}
			return true;
		}
};

// print current values of import to ostream
template <typename TDataType>
std::ostream& operator<<(std::ostream& out, const DataImport<TDataType>& data);

} // end namespace ug

#include "import_export_impl.h"

#endif /* __H__LIB_DISCRETIZATION__ELEMENT_DATA_IMPORT_EXPORT__ */
