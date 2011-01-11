/*
 * data_export.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EXPORT__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EXPORT__

#include "user_data.h"
#include "../elem_disc/elem_disc_interface.h"

namespace ug{

/// Base class for Data Export
/**
 * An base class for all data exports
 */
template <typename TAlgebra>
class IDataExport
{
	public:
	///	Constructor
		IDataExport() : m_id(-1), m_pObj(NULL) {}
		~IDataExport()	{}

	/// clear ips
		virtual void clear_ips_virt() = 0;

	/// set	Function Group of functions
		void set_function_group(const FunctionGroup& fctGrp)
			{m_pFctGrp = &fctGrp;}

	///	Function Group of functions
		const FunctionGroup& get_function_group() const
		{
			UG_ASSERT(m_pFctGrp != NULL, "No func group set.");
			return *m_pFctGrp;
		}

	///	number of fuctions this export depends on
		size_t num_fct() const
		{
			UG_ASSERT(m_pFctGrp != NULL, "No func group set.");
			return m_pFctGrp->num_fct();
		}

	///	resize arrays
		virtual void resize(const LocalIndices& ind) = 0;

	protected:
	/// number of functions the data depends on
		const FunctionGroup* m_pFctGrp;

	protected:
	///	type of local vector
		typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;

	///	type of evaluation function
		typedef bool (IElemDisc<TAlgebra>::*ExportFunc)(const local_vector_type& u,
														bool compDeriv);

	///	function pointers for all elem types
		std::vector<ExportFunc>	m_vExportFunc;

	/// current Geom Object
		int m_id;

	///	elem disc
		IElemDisc<TAlgebra>* m_pObj;

	public:
	///	sets the geometric object type
		bool set_geometric_object_type(int id)
		{
			if(id < 0 || (size_t)id >= m_vExportFunc.size() || m_vExportFunc[id] == NULL)
			{
				UG_LOG("No or not all functions registered "
						"for object with reference object id " << id << ".\n");
				m_id = -1; return false;
			}
			else{m_id = id;	return true;}
		}

	///	register evaluation of linear defect for a element
		template <typename TFunc>
		void register_export_func(int id, IElemDisc<TAlgebra>* obj, TFunc func)
		{
		//	make sure that there is enough space
			if((size_t)id >= m_vExportFunc.size())
				m_vExportFunc.resize(id+1, NULL);

			m_vExportFunc[id] = (ExportFunc)func;

			if(m_pObj == NULL) m_pObj = obj;
			else if(m_pObj != obj)
				throw(UGFatalError("Exports assume to be used by on object for all functions."));
		}

	///	compute lin defect
		bool compute_export(const local_vector_type& u, bool compDeriv)
		{
			UG_ASSERT(m_id >=0, "ElemType id is not set correctly.");
			UG_ASSERT((size_t)m_id < m_vExportFunc.size(), "id "<<m_id<<" not registered");
			UG_ASSERT(m_vExportFunc[m_id] != NULL, "Func pointer is NULL");
			return (m_pObj->*(m_vExportFunc[m_id]))(u, compDeriv);
		}
};

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim, typename TAlgebra>
class DataExport : 	public IPData<TData, dim>,
					public IDataExport<TAlgebra>
{
	public:
		using IPData<TData, dim>::num_series;
		using IPData<TData, dim>::num_ip;
		using IPData<TData, dim>::local_ips;

	public:
	/// Constructor
		DataExport()
		{}

	///	clear ips
		virtual void clear_ips_virt(){this->clear_ips();}

	/// compute values (and derivatives iff compDeriv == true)
		virtual void compute(bool compDeriv = false)
			{throw(UGFatalError("Not Implemented."));}

	////////////////////////////////
	// Derivatives
	////////////////////////////////

	/// number of shapes for local function
		size_t num_sh(size_t s, size_t fct) const
		{
			const size_t ip = 0;
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			return m_vvvDeriv[s][ip][fct].size();
		}

	///	returns the derivative of the local function, at ip and for a dof
		const TData& deriv(size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			UG_ASSERT(dof < m_vvvDeriv[s][ip][fct].size(), "Invalid index.");
			return m_vvvDeriv[s][ip][fct][dof];
		}

	///	returns the derivatives of the local function, at ip
		TData* deriv(size_t s, size_t ip, size_t fct)
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			return &(m_vvvDeriv[s][ip][fct][0]);
		}

	///	resize lin defect arrays
		virtual void resize(const LocalIndices& ind)
		{
		//	resize num fct
			for(size_t s = 0; s < m_vvvDeriv.size(); ++s)
				for(size_t ip = 0; ip < m_vvvDeriv[s].size(); ++ip)
				{
				//	resize num fct
					m_vvvDeriv[s][ip].resize(ind.num_fct());

				//	resize dofs
					for(size_t fct = 0; fct < ind.num_fct(); ++fct)
						m_vvvDeriv[s][ip][fct].resize(ind.num_dofs(fct));
				}
		}

	protected:
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
		{
		//	adjust data arrays
			m_vvvDeriv.resize(vNumIP.size());
			for(size_t s = 0; s < vNumIP.size(); ++s)
				m_vvvDeriv[s].resize(vNumIP[s]);

		//	resize values
			IPData<TData, dim>::adjust_global_ips_and_data(vNumIP);
		}

	protected:

	///	Derivatives
	// Data (size: (0,...,num_series-1) x (0,...,num_ip-1) x (0,...,num_fct-1) x (0,...,num_sh(fct) )
		std::vector<std::vector<std::vector<std::vector<TData> > > > m_vvvDeriv;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EXPORT__ */
